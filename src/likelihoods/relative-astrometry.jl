# ---------------------------------------------------
# Relative astrometry
#
# What this file no longer contains is the point: v1 rebuilt the companion's
# apparent position by superposing the reflex motions of every *inner*
# companion, selecting them at runtime by comparing semi-major axes. That
# loop appeared here and in five other likelihoods, it broke for crossing or
# equal-`a` orbits, and it had no way to express a moon.
#
# `raoff(sol, target, ref)` is now the whole model.
# ---------------------------------------------------

const astrom_cols_radec = (:epoch, :ra, :dec, :σ_ra, :σ_dec)
const astrom_cols_seppa = (:epoch, :pa, :sep, :σ_pa, :σ_sep)

"""
    RelAstromObs(data; target, ref, name, variables=@variables begin end)

Relative astrometry: the sky-plane offset of `target` from `ref` [mas].

    astrom = RelAstromObs(tab; target=b, ref=A, name="GPI")

`data` needs an `:epoch` column [MJD] plus either `:ra`, `:dec`, `:σ_ra`,
`:σ_dec` (optionally `:cor`) or `:pa`, `:sep`, `:σ_pa`, `:σ_sep`; angles are
mas, position angles radians.

Both references take the full grammar — a body, `Barycentre(A, b)`,
`Photocentre` — so an inner binary's photocentre versus an outer companion is
spelled directly rather than assembled from per-planet terms.

# Variables
  - `jitter` [mas] — added in quadrature to both components.
  - `platescale` — multiplicative scale on the measured separation.
  - `northangle` [rad] — rotation of the measured position angle.
"""
struct RelAstromObs{TTable<:Table,TDist,TT,TR} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    pointwise::TDist
    target::TT
    ref::TR
    name::String
end

function RelAstromObs(observations;
                      target, ref,
                      name,
                      variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    table = Table(observations)
    equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    cols = Tables.columnnames(table)
    (issubset(astrom_cols_radec, cols) || issubset(astrom_cols_seppa, cols)) ||
        error("Expected columns $astrom_cols_radec or $astrom_cols_seppa")
    if any(>=(mjd("2050")), table.epoch) || any(<=(mjd("1950")), table.epoch)
        @warn "Epochs fell outside 1950–2050; the expected format is MJD. Double-check your input."
    end
    table = table[sortperm(vec(table.epoch))]

    if hasproperty(table, :pa) && hasproperty(table, :sep)
        σ₁, σ₂ = table.σ_pa, table.σ_sep
        (any(>=(2pi), table.pa) || any(<=(-2pi), table.pa)) &&
            @warn "Position angles fell outside [-2π, 2π]; the expected format is radians."
    else
        σ₁, σ₂ = table.σ_ra, table.σ_dec
    end

    # Pre-factorize the per-point 2×2 covariance; Distributions.jl caches the
    # factorization inside MvNormal.
    pointwise = if hasproperty(table, :cor)
        any(abs.(table.cor) .> 1 - 1e-5) && error("Correlation values are not well-specified")
        broadcast(σ₁, σ₂, table.cor) do a, b, c
            MvNormal(@SArray[a^2 c*a*b; c*a*b b^2])
        end
    else
        broadcast((a, b) -> MvNormal(Diagonal(@SArray[a^2, b^2])), σ₁, σ₂)
    end
    pointwise = (pointwise...,)

    t, r = refspec(target), refspec(ref)
    return RelAstromObs{typeof(table),typeof(pointwise),typeof(t),typeof(r)}(
        table, priors, derived, pointwise, t, r, String(name))
end

export RelAstromObs

refspecs(obs::RelAstromObs) = (obs.target, obs.ref)

likeobj_from_epoch_subset(obs::RelAstromObs, inds) = RelAstromObs(
    obs.table[inds, :, 1]; target=obs.target, ref=obs.ref, obs.name,
    variables=(obs.priors, obs.derived))

"""
    simulate(obs::RelAstromObs, ctx) -> (; ra_model, dec_model, epochs)

Model astrometry at this observation's epochs, allocating. `simulate!` fills
caller storage instead; both share one implementation with `ln_like`.
"""
function simulate!(ra_model, dec_model, obs::RelAstromObs, ctx::ObsContext)
    target = ref(ctx, obs.target)
    reference = ref(ctx, obs.ref)
    @inbounds for i in eachindex(obs.table.epoch)
        sol = solutionat(ctx, i)
        ra_model[i] = raoff(sol, target, reference)
        dec_model[i] = decoff(sol, target, reference)
    end
    return (; ra_model, dec_model, epochs=obs.table.epoch)
end
function simulate(obs::RelAstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    L = length(obs.table.epoch)
    return simulate!(Vector{T}(undef, L), Vector{T}(undef, L), obs, ctx)
end

function ln_like(obs::RelAstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
    platescale = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : one(T)
    northangle = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : zero(T)
    seppa = hasproperty(obs.table, :pa) && hasproperty(obs.table, :sep)

    ll = zero(T)
    L = length(obs.table.epoch)
    @no_escape begin
        ra_model = @alloc(T, L)
        dec_model = @alloc(T, L)
        simulate!(ra_model, dec_model, obs, ctx)

        for i in eachindex(obs.table.epoch)
            ram, decm = ra_model[i], dec_model[i]
            if seppa
                ρ = hypot(ram, decm)
                pa = atan(ram, decm)
                pa_diff = (obs.table.pa[i] + northangle - pa + π) % 2π - π
                pa_diff = pa_diff < -π ? pa_diff + 2π : pa_diff
                resid1 = pa_diff
                resid2 = obs.table.sep[i] * platescale - ρ
                σ₁, σ₂ = obs.table.σ_pa[i], obs.table.σ_sep[i]
            else
                pa_dat = atan(obs.table.dec[i], obs.table.ra[i]) + northangle
                sep_dat = hypot(obs.table.dec[i], obs.table.ra[i]) * platescale
                resid1 = sep_dat * cos(pa_dat) - ram
                resid2 = sep_dat * sin(pa_dat) - decm
                σ₁, σ₂ = obs.table.σ_ra[i], obs.table.σ_dec[i]
            end

            if iszero(jitter)
                ll += logpdf(obs.pointwise[i], @SVector[resid1, resid2])
            else
                s₁ = hypot(σ₁, jitter)
                s₂ = hypot(σ₂, jitter)
                Σ = if hasproperty(obs.table, :cor)
                    c = obs.table.cor[i]
                    @SArray[s₁^2 c*s₁*s₂; c*s₁*s₂ s₂^2]
                else
                    Diagonal(@SArray[s₁^2, s₂^2])
                end
                ll += logpdf(MvNormal(Σ), @SVector[resid1, resid2])
            end
        end
    end
    return ll
end

function generate_from_params(obs::RelAstromObs, ctx::ObsContext; add_noise)
    sim = simulate(obs, ctx)
    epoch = obs.table.epoch
    θ_obs = ctx.θ_obs
    platescale = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : 1.0
    northangle = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : 0.0
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : 0.0

    sep_model = hypot.(sim.ra_model, sim.dec_model)
    pa_model = atan.(sim.ra_model, sim.dec_model)
    if hasproperty(obs.table, :pa) && hasproperty(obs.table, :sep)
        sep = sep_model ./ platescale
        pa = pa_model .- northangle
        if add_noise
            sep = sep .+ randn.() .* obs.table.σ_sep
            pa = pa .+ randn.() .* obs.table.σ_pa
        end
        tab = Table(; epoch, sep, pa, σ_sep=obs.table.σ_sep, σ_pa=obs.table.σ_pa)
    else
        ra = (sep_model ./ platescale) .* sin.(pa_model .- northangle)
        dec = (sep_model ./ platescale) .* cos.(pa_model .- northangle)
        if add_noise
            ra = ra .+ randn.() .* hypot.(obs.table.σ_ra, jitter)
            dec = dec .+ randn.() .* hypot.(obs.table.σ_dec, jitter)
        end
        tab = Table(; epoch, ra, dec, σ_ra=obs.table.σ_ra, σ_dec=obs.table.σ_dec)
    end
    return RelAstromObs(tab; target=obs.target, ref=obs.ref, obs.name,
        variables=(obs.priors, obs.derived))
end
