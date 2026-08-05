# ---------------------------------------------------
# Radial velocity with the instrument zero point marginalized out
#
# This is the one RV type that did *not* fold into core `RadialVelocityObs`.
# Not because of the references — it takes `target`/`ref` like everything
# else, and v1's "StarAbsolute" name only ever meant "the star against the
# barycentre" — but because the marginalization couples every row of the
# instrument. There is no pointwise likelihood to hold out, which is a
# structural property of the estimator and survives any refactor.
# ---------------------------------------------------

using Bumper
using Octofitter: ObsContext, Priors, Derived, solutionat, refspec,
    likelihoodname, equal_length_cols, rv_cols, _system_number_type

"""
    MarginalizedRVObs(data; target, ref=Barycentre, name, variables=…,
                      trend_function=…)

Radial velocities [m/s] of `target` measured against `ref`, with this
instrument's zero point analytically marginalized out rather than sampled.

    rvs = MarginalizedRVObs(tab; target=A, ref=Barycentre, name="HARPS",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100)   # m/s, added in quadrature
        end)

`data` needs `:epoch` [MJD], `:rv` [m/s] and `:σ_rv` [m/s]. Declare no
`offset` — that is the parameter this type integrates over. Everything else
(`target`/`ref`, `jitter`, `trend_function`) means what it does on
[`RadialVelocityObs`](@ref).

Marginalizing removes one dimension per instrument, which usually samples
better. Reach for `RadialVelocityObs` instead when you want a specific prior
on the zero point, correlations between instruments' zero points, a
hierarchical model over them — or a Gaussian process, which this type does
not support.

!!! note "No cross-validation"
    The marginalization couples every point in the instrument, so this type
    has no pointwise likelihoods and `likeobj_from_epoch_subset` errors. See
    that method's message.

# Renamed from v8
`MarginalizedStarAbsoluteRVObs` — "Star" and "Absolute" were both just ref
choices, and are now spelled as ones.
"""
struct MarginalizedRVObs{TTable<:Table,TT,TR,TF} <: Octofitter.AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    target::TT
    ref::TR
    trend_function::TF
    name::String
end

function MarginalizedRVObs(observations;
                           target, ref=Barycentre,
                           name,
                           trend_function=Octofitter._rv_no_trend,
                           variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    table = Table(observations)
    equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    issubset(rv_cols, Tables.columnnames(table)) ||
        error("Expected columns $rv_cols")
    if any(>=(mjd("2050")), table.epoch) || any(<=(mjd("1950")), table.epoch)
        @warn "Epochs fell outside 1950–2050; the expected format is MJD. Double-check your input."
    end
    table = table[sortperm(vec(table.epoch))]
    (haskey(priors.priors, :offset) || haskey(derived.variables, :offset)) &&
        error("`MarginalizedRVObs` marginalizes the zero point analytically, so an " *
              "`offset` variable would be integrated over *and* added — remove it, or " *
              "use `RadialVelocityObs` to sample it explicitly.")

    t, r = refspec(target), refspec(ref)
    return MarginalizedRVObs{typeof(table),typeof(t),typeof(r),typeof(trend_function)}(
        table, priors, derived, t, r, trend_function, String(name))
end

export MarginalizedRVObs

Octofitter.refspecs(obs::MarginalizedRVObs) = (obs.target, obs.ref)

function Octofitter.likeobj_from_epoch_subset(obs::MarginalizedRVObs, obs_inds)
    error("""
    Data subsetting is not supported for MarginalizedRVObs, and cannot be.

    This type integrates the instrument zero point out of the likelihood
    analytically. That integral runs over every row at once — the marginal
    likelihood of a subset is not any product of per-row terms, so there is no
    pointwise likelihood to hold out or to weight. Cross-validation and PSIS-LOO
    both need one.

    This is a property of the estimator, not a missing feature. Use
    `RadialVelocityObs` with an explicit `offset` variable if you need to
    cross-validate; the two agree on the posterior up to that parameter.
    """)
end

# ---------------------------------------------------
# Forward model
#
# `radvel(sol, target, ref)` — the same one call `RadialVelocityObs` makes.
# v1 summed `radvel(sol_j, mass_j)` over every planet by hand, which is the
# reflex superposition this refactor deletes.
# ---------------------------------------------------

function Octofitter.simulate!(rv_model, obs::MarginalizedRVObs, ctx::ObsContext)
    target = Octofitter.ref(ctx, obs.target)
    reference = Octofitter.ref(ctx, obs.ref)
    # No offset term: that is the whole point of the type.
    @inbounds for i in eachindex(obs.table.epoch)
        rv_model[i] = radvel(solutionat(ctx, i), target, reference) +
                      obs.trend_function(ctx.θ_obs, obs.table.epoch[i])
    end
    return (; rv_model, epochs=obs.table.epoch)
end
function Octofitter.simulate(obs::MarginalizedRVObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    return Octofitter.simulate!(Vector{T}(undef, length(obs.table.epoch)), obs, ctx)
end

"""
    marginalized_offset(obs::MarginalizedRVObs, ctx) -> Real

The zero point this instrument's data most prefer under the current sample:
the inverse-variance weighted mean residual, which is the mode (and mean) of
the Gaussian being integrated over in `ln_like`.

Nothing in the likelihood needs it — the integral is closed form — but
plotting does, since the data have to be put on the model's scale somehow.
"""
function marginalized_offset(obs::MarginalizedRVObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    sim = Octofitter.simulate(obs, ctx)
    num = zero(T)
    den = zero(T)
    for i in eachindex(obs.table.epoch)
        var = obs.table.σ_rv[i]^2 + jitter^2
        num += (obs.table.rv[i] - sim.rv_model[i]) / var
        den += 1 / var
    end
    return num / den
end
export marginalized_offset

function Octofitter.ln_like(obs::MarginalizedRVObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    L = length(obs.table.epoch)
    ll = zero(T)

    @no_escape ctx.buf begin
        rv_model = @alloc(T, L)
        Octofitter.simulate!(rv_model, obs, ctx)

        # Marginalize the instrument zero point γ out under a flat prior,
        # following the Orvara paper (Brandt et al. 2021). With residuals rᵢ
        # (data − model, before any offset) and variances varᵢ,
        #
        #   A = Σ 1/varᵢ,   B = −2 Σ rᵢ/varᵢ,   C = Σ rᵢ²/varᵢ
        #
        # so that Σ (rᵢ − γ)²/varᵢ = A γ² + B γ + C. Completing the square and
        # integrating over γ ∈ ℝ gives ∫ exp(−½A(γ + B/2A)²) dγ = √(2π/A), and
        # hence the normalized marginal log-likelihood
        #
        #   log ∫ L(γ) dγ = −½ Σ log(2π varᵢ) − C/2 + B²/(8A) + ½log(2π) − ½log A
        #
        # **This differs from v1 by a factor of two.** v1 (and this port, up to
        # the fix) computed `−Σ log(2π varᵢ) + B²/(4A) − C − log A`, which is
        # exactly 2·(the above) − log(2π). The additive constant is harmless,
        # but the factor of two is not: it is not a constant, so it sharpened
        # this instrument's contribution against every other term in a joint
        # model and against the priors. Will's call (2026-08-05) is to fix it
        # going forward rather than preserve the v1 behaviour. Posteriors from
        # models using this likelihood will shift accordingly; models with a
        # single marginalized-RV instrument and no other data are affected the
        # least, joint models the most.
        A = zero(T)
        B = zero(T)
        C = zero(T)
        @inbounds for i in 1:L
            var = obs.table.σ_rv[i]^2 + jitter^2
            resid = obs.table.rv[i] - rv_model[i]
            A += 1 / var
            B -= 2resid / var
            C += resid^2 / var
            ll -= log(2π * var) / 2
        end
        ll += B^2 / (8A) - C / 2 + log(2π) / 2 - log(A) / 2
    end

    return ll
end

function Octofitter.generate_from_params(obs::MarginalizedRVObs, ctx::ObsContext; add_noise)
    sim = Octofitter.simulate(obs, ctx)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    rv = add_noise ? sim.rv_model .+ randn.() .* hypot.(obs.table.σ_rv, jitter) : sim.rv_model
    return MarginalizedRVObs(Table(; epoch=obs.table.epoch, rv, σ_rv=obs.table.σ_rv);
        target=obs.target, ref=obs.ref, obs.name,
        trend_function=obs.trend_function, variables=(obs.priors, obs.derived))
end

# ---------------------------------------------------
# Plotting protocol
#
# Same single channel as `RadialVelocityObs`; the only difference is where
# the zero point comes from, since there is no `offset` variable to read.
# ---------------------------------------------------

Octofitter.plotchannels(obs::MarginalizedRVObs) = (
    Octofitter.PlotChannel(:rv, "radial velocity", "m/s";
        query=Octofitter.ObservableQuery(PlanetOrbits.radvel, obs.target, obs.ref)),
)

function Octofitter.residuals(obs::MarginalizedRVObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    sim = Octofitter.simulate(obs, ctx)          # no offset in the model
    offset = marginalized_offset(obs, ctx)
    ep = collect(Float64, obs.table.epoch)
    # Remove the trend as well as the marginalized zero point, on both sides:
    # the `radvel` query the panel draws its curve from knows about neither, so
    # points that kept the trend would not lie on their own curve. (Same rule
    # as `RadialVelocityObs`; `simulate` puts the trend into `rv_model`.)
    trend = obs.trend_function.(Ref(ctx.θ_obs), obs.table.epoch)
    model_pure = sim.rv_model .- trend
    data_cal = obs.table.rv .- offset .- trend
    return (;
        rv=(; epoch=ep, data=collect(data_cal), model=collect(model_pure),
            resid=collect(obs.table.rv .- offset .- sim.rv_model),
            σ=collect(float.(obs.table.σ_rv)),
            σ_eff=collect(hypot.(obs.table.σ_rv, jitter)),
            use=trues(length(ep))),
    )
end
