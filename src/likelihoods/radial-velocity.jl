# ---------------------------------------------------
# Radial velocity
#
# One likelihood covers both kinds of radial velocity. Both are
# `radvel(sol, target, ref)`, differing only in which references they name:
#   star vs the system barycentre → the classic reflex signal
#   companion vs the host star    → relative RV
#
# The nuisance machinery — the correlated-noise model, the trend function, the
# held-out-data path — is shared, because none of it depends on which pair of
# bodies is being differenced.
# ---------------------------------------------------

"""
    RadialVelocityObs(data; target, ref=Barycentre, name, variables=…,
                      trend_function=…, gaussian_process=…)

Radial velocities [m/s] of `target` measured against `ref`, positive
receding.

    rvs = RadialVelocityObs(tab; target=A, ref=Barycentre, name="HARPS",
        variables=@variables begin
            offset ~ Normal(0, 100)      # m/s, instrument zero point
            jitter ~ LogUniform(0.1, 100) # m/s, added in quadrature
        end)

`data` needs `:epoch` [MJD], `:rv` [m/s] and `:σ_rv` [m/s].

One `RadialVelocityObs` per instrument: the offset and jitter are per
instrument, and that is exactly the granularity this object carries.

# Variables
  - `offset` [m/s] — instrument zero point, added to the model. **Optional**,
    for both absolute and relative RV. v8's `StarAbsoluteRVObs` injected
    `offset ~ Uniform(-1000, 1000)` when no variables block was given; v9
    never invents a prior, so a model that relied on that now has *no*
    offset unless it declares one.
  - `jitter` [m/s] — added in quadrature to `σ_rv`. Optional, same reasoning.

# Trends

    trend_function = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000)

Evaluated per epoch and added to the model alongside `offset`. Its
parameters are ordinary variables of this observation.

# Corrections and data provenance

Two keywords declare what your reduction already did, and one is an
experimentation hook. See the "How Octofitter Computes Orbits" manual page.

    secular_acceleration = :model | :data_corrected     # default :model

The perspective-acceleration drift — the barycentre's line-of-sight velocity
changing as its 3D motion swings the line of sight. `:model` adds it to the
prediction (its constant part washes into `offset`; the drift and curvature
are the payload: ~4.5 m/s/yr for a Barnard-class star, ~3 cm/s/yr for a
generic 30 km/s star at 30 pc). `:data_corrected` says your pipeline already
removed it, using frozen catalog values, so the model adds nothing.

Why this one is model-side by default while the barycentric correction is
not: the Earth-side terms depend only on the Earth ephemeris and the catalog
direction — zero fit parameters — so a pipeline removes them once, perfectly.
Perspective acceleration depends only on **μ, ϖ and rv, which are sampled
parameters here**, and its signature (a smooth trend with slight curvature) is
degenerate with a long-period companion. `:model` is also what makes a joint
RV + absolute-astrometry fit self-consistent: the same sampled (μ, ϖ) predict
both data types.

It requires an `AbsoluteFrame` to be non-zero; with `plx` alone the term is
definitionally zero and no error is raised. It is rejected at construction for
a **relative** series, with an explanation — a blanket "it cancels" is false,
and the half that does not cancel is large.

`secular_acceleration` is this observation's **only** provenance dimension.
Whether the Einstein term is modelled is a property of the whole fit, not of
one instrument — no pipeline can have removed the part of it that varies with
the orbit — so it lives on the system as `System(…; einstein_rv=:on|:off)`.

# Correlated noise

    gaussian_process = θ_obs -> GP(θ_obs.gp_η₁^2 * SqExponentialKernel() ∘
                                   ScaleTransform(1/θ_obs.gp_η₂))

Fits the residuals with a GP rather than assuming they are independent. The
callable is handed this observation's variables and must return a GP object;
the residuals it sees are at *this observation's own epochs, in table order*,
which the constructor has sorted.

Two backends are supported and neither is a dependency of Octofitter: an
AbstractGPs `GP` works through the duck-typed default path
([`gp_condition`](@ref)), and `OctofitterRadialVelocity` adds methods for its
vendored Celerite. Load `OctofitterRadialVelocity` for the Celerite kernels.

# Cross-validation

`likeobj_from_epoch_subset(obs, rows)` returns an observation whose `ln_like`
is the log-likelihood of `rows` alone. Without a GP that is just the subset;
with one the remaining rows are kept as the *conditioning* set and `rows`
become `held_out_table`, since a correlated model cannot score a point
without the points it is correlated with. Prediction is implemented for
Celerite only — the AbstractGPs case throws, exactly as in v1.
"""
struct RadialVelocityObs{TTable<:Table,THeld<:Table,TT,TR,GP,TF} <: AbstractObs
    table::TTable
    # Rows this observation *scores* while conditioning on `table`. Empty in
    # every ordinary model; populated only by `likeobj_from_epoch_subset`.
    held_out_table::THeld
    priors::Priors
    derived::Derived
    target::TT
    ref::TR
    gaussian_process::GP
    trend_function::TF
    name::String
    # `:model`, `:data_corrected`, or `:not_applicable` for a relative-RV
    # series, where the term this governs cancels. Resolved at construction.
    # The *only* provenance dimension an RV series has: whether the Einstein
    # term is modelled is a property of the whole fit, not of one instrument,
    # and lives on `System` as `einstein_rv`.
    secular_acceleration::Symbol
end

const rv_cols = (:epoch, :rv, :σ_rv)

# Default trend: nothing, in whatever number type the sample is being
# evaluated at (so it composes with ForwardDiff without widening).
_rv_no_trend(θ_obs, epoch) = zero(_system_number_type(θ_obs))

function RadialVelocityObs(observations;
                           target, ref=Barycentre,
                           name,
                           gaussian_process=nothing,
                           trend_function=_rv_no_trend,
                           held_out=nothing,
                           secular_acceleration=nothing,
                           variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    # Collect the columns: a multithreaded `CSV.read` returns `ChainedVector`s,
    # which the per-epoch indexing in the likelihood cannot take a row view of.
    table = materialize_cols(Table(observations))
    equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    issubset(rv_cols, Tables.columnnames(table)) ||
        error("Expected columns $rv_cols")
    if any(>=(mjd("2050")), table.epoch) || any(<=(mjd("1950")), table.epoch)
        @warn "Epochs fell outside 1950–2050; the expected format is MJD. Double-check your input."
    end
    table = table[sortperm(vec(table.epoch))]

    held_out_table = isnothing(held_out) ? empty(table) : Table(held_out)
    isempty(held_out_table) || issubset(rv_cols, Tables.columnnames(held_out_table)) ||
        error("Expected columns $rv_cols in the held-out table")

    # The GP is fit to residuals in table order, and the constructor has just
    # sorted them — but "sorted" is not "strictly increasing", and Celerite's
    # O(N) Cholesky assumes the latter (two rows at the same epoch also make
    # the covariance singular for any stationary kernel). Assert it here
    # rather than let it surface as a `PosDefException` that the sampler
    # guards below would quietly turn into `-Inf`.
    if !isnothing(gaussian_process) && !_rv_strictly_increasing(table.epoch)
        i = findfirst(k -> table.epoch[k+1] <= table.epoch[k], 1:length(table.epoch)-1)
        error("""
        `gaussian_process` requires strictly increasing epochs, but rows $i and $(i+1)
        of "$name" are both at epoch $(table.epoch[i]).

        Two measurements at the identical epoch carry no information a stationary
        kernel can separate. Merge them, or nudge one epoch by the exposure time.
        """)
    end

    t, r = refspec(target), refspec(ref)
    sec = _resolve_secular_acceleration(secular_acceleration, r, String(name))
    return RadialVelocityObs{typeof(table),typeof(held_out_table),typeof(t),typeof(r),
                             typeof(gaussian_process),typeof(trend_function)}(
        table, held_out_table, priors, derived, t, r,
        gaussian_process, trend_function, String(name), sec)
end

# Whole-system barycentre only: that is what makes a series *absolute*, and
# the perspective-acceleration drift is a property of the frame, not of any
# pair inside it.
_is_absolute_rv(::BarycentreSpec{()}) = true
_is_absolute_rv(@nospecialize(_)) = false

function _resolve_secular_acceleration(spec, r, name)
    if isnothing(spec)
        return _is_absolute_rv(r) ? :model : :not_applicable
    end
    spec in (:model, :data_corrected) || error(
        "`secular_acceleration` takes `:model` (the default for an absolute " *
        "series — Octofitter adds the propagated perspective-acceleration " *
        "drift to the prediction) or `:data_corrected` (your pipeline already " *
        "removed it, so the model adds nothing); got $(repr(spec)).")
    _is_absolute_rv(r) || _err_relative_secular(name)
    return spec
end

@noinline _err_relative_secular(name) = error("""
`secular_acceleration` is only meaningful for an absolute series (`ref=Barycentre`),
but "$name" measures one body against another.

Not because "it cancels" — half of it does and half of it does not, and the half
that does not is large:

  * What cancels is the **common-mode** barycentre drift: both bodies share the
    barycentre's line-of-sight velocity, so it drops out of the difference. That
    is the only part this keyword governs, which is why there is nothing here for
    it to declare.

  * What does **not** cancel is the **differential** term. Two bodies sit at
    different sky positions, so the barycentre's transverse space velocity
    projects differently onto each: ≈ 0.023 · μ["/yr] · a[AU] m/s — 24 cm/s for a
    Barnard-like host with a 1 AU companion, and ~239 m/s at 1000 AU. That term is
    already modelled, as correction #4 of PlanetOrbits' observing-geometry pass,
    computed from the sampled frame velocity and each body's fitted direction —
    which is also why no pipeline can have removed it for you.

So: drop the keyword. If this is a high-proper-motion system, make sure
`observing_geometry` is on for the differential term (`:auto` will keep it on for
data like this, since Δprediction ≫ σ).
""")

_rv_strictly_increasing(x) = all(k -> x[k] < x[k+1], 1:length(x)-1)

export RadialVelocityObs

refspecs(obs::RadialVelocityObs) = (obs.target, obs.ref)

# Held-out rows need solved states too, and they are not in `table`. Putting
# them at the end of `epochs(obs)` means `solutionat(ctx, L + k)` reaches
# them: v1 re-ran `orbitsolve` per held-out point because per-planet solution
# caches only covered the likelihood's own table, and there is no reason to
# pay that now that one trajectory covers the whole model.
epochs(obs::RadialVelocityObs) =
    isempty(obs.held_out_table) ? collect(Float64, obs.table.epoch) :
    vcat(collect(Float64, obs.table.epoch), collect(Float64, obs.held_out_table.epoch))

"""
    likeobj_from_epoch_subset(obs::RadialVelocityObs, rows)

An observation whose `ln_like` is the log-likelihood of `rows` of `obs`.

Without a GP the rows are independent, so that is literally the subset —
identical to what `RelAstromObs` does. With a GP they are not: the score of a
held-out point depends on the points it is correlated with, so the *other*
rows are kept as the conditioning set and `rows` move to `held_out_table`.
Both spellings answer the same question, which is what cross-validation
asks.

`rows = :` means "all of them", i.e. the ordinary full-data likelihood, with
nothing left to condition on.
"""
function likeobj_from_epoch_subset(obs::RadialVelocityObs, rows)
    inds = _rv_rowinds(obs.table, rows)
    keep = setdiff(1:length(obs.table.epoch), inds)
    if isnothing(obs.gaussian_process) || isempty(keep)
        return _rv_rebuild(obs, obs.table[inds], nothing)
    end
    return _rv_rebuild(obs, obs.table[keep], obs.table[inds])
end

_rv_rowinds(table, ::Colon) = collect(1:length(table.epoch))
_rv_rowinds(table, i::Integer) = [i]
_rv_rowinds(table, rows) = collect(rows)

_rv_rebuild(obs::RadialVelocityObs, table, held_out) = RadialVelocityObs(
    table; target=obs.target, ref=obs.ref, obs.name, held_out,
    gaussian_process=obs.gaussian_process, trend_function=obs.trend_function,
    secular_acceleration=(obs.secular_acceleration === :not_applicable ? nothing :
                          obs.secular_acceleration),
    variables=(obs.priors, obs.derived))

# ---------------------------------------------------
# Forward model
# ---------------------------------------------------

# Fill `out` with the modelled RV at `epochvec`, whose rows start at column
# `i0 + 1` of this observation's epoch map. Refs are resolved by the caller,
# outside any epoch loop.
@inline function _rv_model!(out, obs::RadialVelocityObs, ctx::ObsContext,
                            target, reference, offset, epochvec, i0)
    T = _system_number_type(ctx.θ_system)
    # Both loop-invariant, so the branches hoist out of the loop. One is a
    # whole-model setting (`System(…; einstein_rv=…)`, reaching us through the
    # context) and one is this series' own provenance declaration.
    spectroscopic = ctx.einstein_rv
    secular = obs.secular_acceleration === :model
    @inbounds for i in eachindex(epochvec)
        sol = solutionat(ctx, i0 + i)
        v = spectroscopic ? radvel(sol, target, reference) :
            velz(sol, target, reference) * _AU_PER_YEAR_TO_M_PER_S
        if secular
            v += _frame_drift_rv(sol, T)
        end
        out[i] = v + offset + obs.trend_function(ctx.θ_obs, epochvec[i])
    end
    return out
end

const _AU_PER_YEAR_TO_M_PER_S = PlanetOrbits.au2m / PlanetOrbits.year2sec_julian

# The frame's own propagated radial velocity relative to `ref_epoch` — the
# perspective-acceleration drift. Its constant part washes into `offset`; the
# drift and curvature are the payload (~4.5 m/s/yr for a Barnard-class star,
# ~3 cm/s/yr for a generic 30 km/s star at 30 pc).
#
# `fr.rv` is the right anchor under either setting of `barycentric_lighttime`:
# PlanetOrbits leaves the radial readout as the coordinate rate ḋ at the
# emission event (the spectroscopic convention) rather than converting it to
# an apparent rate the way it does the proper motions, so
# `frame_rv(ref_epoch) == rv` exactly on both paths and this difference is
# zero there by construction.
#
# Definitionally zero without an absolute frame, and that is not an error: the
# frame level chooses the physics, exactly as it does everywhere else in
# PlanetOrbits. The `:auto` correction report is the safety net that tells a
# user their high-pm star wants one.
@inline _frame_drift_rv(sol, ::Type{T}) where {T} =
    _frame_drift_rv(sol, PlanetOrbits.frame(sol), T)
@inline _frame_drift_rv(sol, fr::PlanetOrbits.AbsoluteFrame, ::Type{T}) where {T} =
    convert(T, PlanetOrbits.frame_rv(sol) - fr.rv)
@inline _frame_drift_rv(_, ::PlanetOrbits.AbstractFrame, ::Type{T}) where {T} = zero(T)

function simulate!(rv_model, obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    offset = hasproperty(ctx.θ_obs, :offset) ? ctx.θ_obs.offset : zero(T)
    target = ref(ctx, obs.target)
    reference = ref(ctx, obs.ref)
    _rv_model!(rv_model, obs, ctx, target, reference, offset, obs.table.epoch, 0)
    return (; rv_model, epochs=obs.table.epoch)
end
function simulate(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    return simulate!(Vector{T}(undef, length(obs.table.epoch)), obs, ctx)
end

function ln_like(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    ll = zero(T)
    L = length(obs.table.epoch)
    @no_escape ctx.buf begin
        rv_model = @alloc(T, L)
        resid = @alloc(T, L)
        rv_var = @alloc(T, L)
        simulate!(rv_model, obs, ctx)
        @inbounds for i in 1:L
            resid[i] = obs.table.rv[i] - rv_model[i]
            rv_var[i] = obs.table.σ_rv[i]^2 + jitter^2
        end
        if isnothing(obs.gaussian_process)
            @inbounds for i in 1:L
                ll -= (resid[i]^2 / rv_var[i] + log(2π * rv_var[i])) / 2
            end
        else
            ll = _rv_gp_ln_like(obs, ctx, resid, rv_var, T)
        end
    end
    return ll
end

# The GP branch lives out of line for two reasons: it needs its own
# `@no_escape` for the held-out buffers, and `@no_escape` forbids the early
# returns the `-Inf` guards want.
function _rv_gp_ln_like(obs::RadialVelocityObs, ctx::ObsContext, resid, rv_var, ::Type{T}) where {T}
    θ_obs = ctx.θ_obs
    # These guards are ugly and load-bearing. A sampler proposing a kernel
    # with a non-positive length scale, or hyperparameters that make the
    # covariance numerically indefinite, must see a rejected sample rather
    # than a crashed chain.
    local fx
    conditioned = true
    try
        gp = @inline obs.gaussian_process(θ_obs)
        fx = gp_condition(gp, obs.table.epoch, rv_var)
    catch err
        if err isa DomainError || err isa PosDefException || err isa ArgumentError
            conditioned = false
        else
            rethrow(err)
        end
    end
    conditioned || return convert(T, -Inf)

    ll = zero(T)
    try
        if isempty(obs.held_out_table)
            ll += gp_ln_like(fx, resid)
        else
            ll += _rv_gp_heldout_ln_like(obs, ctx, fx, resid, T)
        end
    catch err
        if err isa PosDefException || err isa DomainError
            ll = convert(T, -Inf)
        else
            rethrow(err)
        end
    end
    return ll
end

function _rv_gp_heldout_ln_like(obs::RadialVelocityObs, ctx::ObsContext, fx, resid, ::Type{T}) where {T}
    θ_obs = ctx.θ_obs
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
    offset = hasproperty(θ_obs, :offset) ? θ_obs.offset : zero(T)
    L = length(obs.table.epoch)
    Lh = length(obs.held_out_table.epoch)
    target = ref(ctx, obs.target)
    reference = ref(ctx, obs.ref)
    ll = zero(T)
    @no_escape ctx.buf begin
        model_h = @alloc(T, Lh)
        _rv_model!(model_h, obs, ctx, target, reference, offset,
            obs.held_out_table.epoch, L)
        # The GP's own predictive variance at the held-out epochs, plus that
        # point's white noise. Conditioning is on `resid`, the training rows.
        pred, pred_var = gp_predict(fx, resid, obs.held_out_table.epoch)
        @inbounds for k in 1:Lh
            r = obs.held_out_table.rv[k] - model_h[k]
            σ² = obs.held_out_table.σ_rv[k]^2 + jitter^2
            ll += logpdf(Normal(pred[k], sqrt(pred_var[k] + σ²)), r)
        end
    end
    return ll
end

function generate_from_params(obs::RadialVelocityObs, ctx::ObsContext; add_noise)
    sim = simulate(obs, ctx)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    rv = add_noise ? sim.rv_model .+ randn.() .* hypot.(obs.table.σ_rv, jitter) : sim.rv_model
    return _rv_rebuild(obs, Table(; epoch=obs.table.epoch, rv, σ_rv=obs.table.σ_rv), nothing)
end

# ---------------------------------------------------
# Gaussian-process backend hooks
#
# `RadialVelocityObs` lives in core Octofitter, which depends on neither
# AbstractGPs nor the vendored Celerite. The three functions below are the
# entire interface a backend has to satisfy; the defaults implement the
# AbstractGPs shape by duck typing (`gp(x, σ²)` then `logpdf`), which needs no
# package at all, and `OctofitterRadialVelocity` adds Celerite methods.
#
# v1 achieved the same thing with `gp isa Celerite.CeleriteGP` branches inside
# the likelihood, which only worked because the likelihood itself lived in the
# subpackage.
# ---------------------------------------------------

"""
    gp_condition(gp, epochs, σ²) -> fx

Attach the observation's epochs and per-point white-noise variances to the
user's GP, giving whatever object [`gp_ln_like`](@ref) and
[`gp_predict`](@ref) consume.

The default is AbstractGPs' finite-GP call, `gp(epochs, σ²)`. A backend that
computes a factorization in place (Celerite) returns the GP itself.
"""
gp_condition(gp, epochs, σ²) = gp(epochs, σ²)

"""
    gp_ln_like(fx, residuals) -> Real

Log-likelihood of `residuals` under the conditioned GP from
[`gp_condition`](@ref). Defaults to `logpdf`.
"""
gp_ln_like(fx, residuals) = logpdf(fx, residuals)

"""
    gp_predict(fx, residuals, epochs) -> (mean, var)

Posterior predictive mean and variance at `epochs`, given `residuals` at the
conditioning epochs. Two callers: cross-validation scores held-out points
against it, and [`noisemodel`](@ref) turns it into the correlated-noise band a
plot draws.

`OctofitterRadialVelocity` implements both shipped backends — Celerite, and
AbstractGPs via `posterior`/`mean_and_var`. (v8 had the Celerite case only, so
cross-validating an AbstractGPs-correlated fit threw; that hole is closed.)
There is no default, because there is no way to guess it: a backend that
implements [`gp_condition`](@ref) and [`gp_ln_like`](@ref) must implement this
too if its fits are to be cross-validated or plotted.
"""
gp_predict(fx, residuals, epochs) = error("""
Posterior prediction is not implemented for GP backend $(typeof(fx)).

Cross-validating a GP-correlated `RadialVelocityObs` needs the GP's predictive
distribution at the held-out epochs, and plotting its noise band needs the same
thing on the plot grid. Add a method:

    Octofitter.gp_predict(fx::$(typeof(fx)), residuals, epochs) = (mean, variance)

`OctofitterRadialVelocity` does this for Celerite and for AbstractGPs.
""")

# --- the `:auto` correction report -------------------------------------------

has_correction_impact(::Type{<:RadialVelocityObs}) = true
correction_impact(obs::RadialVelocityObs, a::ObsContext, b::ObsContext) =
    _simulate_impact(simulate(obs, a), simulate(obs, b), (:rv_model,),
                     _tightest(obs.table.σ_rv))

# Two magnitudes worth reporting for every RV series, neither of which any
# flag governs: the Einstein term (always modelled — no pipeline can have
# removed the part of it that varies) and the perspective-acceleration drift
# (declared per series, and zero without an absolute frame). Both as
# peak-to-peak over the series, which is the part an offset cannot absorb.
function correction_advisories(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    target = ref(ctx, obs.target)
    reference = ref(ctx, obs.ref)
    σ = _tightest(obs.table.σ_rv)
    lo_e = Inf; hi_e = -Inf
    lo_s = Inf; hi_s = -Inf
    @inbounds for i in eachindex(obs.table.epoch)
        sol = solutionat(ctx, i)
        e = float(radvel(sol, target, reference) -
                  velz(sol, target, reference) * _AU_PER_YEAR_TO_M_PER_S)
        lo_e = min(lo_e, e); hi_e = max(hi_e, e)
        s = float(_frame_drift_rv(sol, T))
        lo_s = min(lo_s, s); hi_s = max(hi_s, s)
    end
    adv = (("Einstein term, peak-to-peak [m/s]", hi_e - lo_e, σ),)
    obs.secular_acceleration === :model || return adv
    return (adv..., ("secular-acceleration drift, peak-to-peak [m/s]", hi_s - lo_s, σ))
end
