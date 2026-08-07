# ---------------------------------------------------
# `HipparcosIADObs` — standalone Hipparcos IAD likelihood
#
# The reader and sky-path reconstruction already crossed to v2 in
# `hipparcos-iad.jl`; this is the likelihood that fits the per-transit
# abscissa residuals on their own. `G23HObs` carries the same data as one
# channel, but a Hipparcos-only fit is a real use case and stays available.
#
# The refs replace v1's `fluxratio_hip`, and the `iad_Δra`/`iad_Δdec`/
# `iad_Δplx`/`iad_Δpmra`/`iad_Δpmdec` block is to be recognized here as a
# *general* five-parameter frame-offset facility (Hipparcos today, HST FGS
# next) rather than a Hipparcos-specific one. See
# `design/observation-types-migration.md` §3.5, §6.6.
#
# What the v2 port changed, and why
# ---------------------------------
#
# 1. `fluxratio_hip` is no longer a required positional `Product` indexed by
#    planet declaration order. The observation names its `host` and its
#    `companions`, and the Hp-band ratios default to the bodies' own
#    `flux_Hp` variables (see `_g23h_fluxratios`). An explicit vector still
#    overrides, which is what per-draw resolved-flag gating needs.
#
# 2. The forward model is the tangent-plane one `G23HObs`'s `:iad_hip`
#    channel uses — a five-parameter frame solution plus the source's own
#    excursion, projected on scan and compared against
#    `proj_meas_alongscan`. v1 instead built the model position from an
#    `AbsoluteVisual` orbit's compensated (α, δ, ϖ(t)) and took the
#    perpendicular distance to the two-point abscissa line, while the data it
#    compared against had been reconstructed with the *tangent-plane* model —
#    the two halves were on different conventions. The scalar projection is
#    algebraically the same distance (the line's direction is (−sinϕ, cosϕ)),
#    so what actually changes is that both halves are now on one convention.
#    Results are therefore NOT bit-identical to v1.
#
# 3. There is no per-planet superposition loop. The Hipparcos grating
#    response is applied once over every companion (`_hippacentre!`), which
#    is the whole point of BINARYS: the modulated response of several
#    companions is not the sum of their individual responses.
# ---------------------------------------------------

"""
    HipparcosIADObs(; hip_id, host, companions=(), ref=Barycentre, …)

The Hipparcos per-transit abscissa residuals (van Leeuwen 2007 intermediate
astrometric data) as a standalone likelihood.

[`G23HObs`](@ref) carries the same data as its `:iad_hip` channel, alongside
the Gaia catalog channels; reach for this type when Hipparcos is all you
have, or when you want the abscissae without the catalog proper motions.

# Source membership

    host=A, companions=(b, c)

`host` is the star the Hipparcos entry is centred on and `companions` are the
bodies whose light may modulate its abscissa. `ref` is what the host's reflex
is measured against — `Barycentre` by default.

Companion Hp-band flux ratios come from the same three-tier lookup
[`G23HObs`](@ref) documents: this observation's own `fluxratio_hip` vector
(the per-draw override, which is what marginalized resolvedness needs), then
a system-level vector of that name, then the bodies' own `flux_Hp`
variables as `flux_Hp(companion) / flux_Hp(host)`. A model with no `flux_Hp`
anywhere and no vector leaves every companion dark, which is the right
answer when the companions contribute no light.

# The forward model, and what it constrains

The measured abscissa carries the Hipparcos catalog's *whole* sky path, so
the model has to reproduce it. That is the five-parameter frame offset of
`skypath.jl` — [`FrameOffset`](@ref) — plus the source's own excursion:

    b_i = (Δα + Δt·μα✱ + Δα_src)·cosϕ + (Δδ + Δt·μδ + Δδ_src)·sinϕ + ϖ·f_AL

with the residual `|proj_meas_alongscan − b|` compared against the
renormalized per-transit σ, inflated by the BINARYS first-harmonic factor.

By default the parallax `ϖ` is the **system's own `plx`** with no offset
(`iad_Δplx = 0`), while the position and proper motion are free nuisances
with wide priors. So a Hipparcos-only fit constrains the parallax and the
companion's reflex, and marginalizes the four linear terms — the standard
setup, and the reason the frame block is offsets rather than fixed values.
Override `variables=` to change that; `iad_Δplx ~ Uniform(-10, 10)` makes the
parallax a nuisance too, exactly as `G23HObs` treats it.

# Data

Give `hip_id` and the IAD is loaded through `hipparcos_iad` from the
`Hipparcos_IAD` data dependency, or pass an already-loaded `iad=(; table,
hip_sol)` (which is how the test suite runs offline). `renormalize`,
`attempt_correction` and `is_van_leeuwen` are forwarded to the loader.
`recalibrate=true` additionally applies the Brandt, Michalik & Brandt shift
(+0.140 mas on the residuals, 2.25 mas of extra dispersion) that `G23HObs`
applies unconditionally; it is off here so that this type reproduces the
catalog data as published unless asked otherwise.

# Variables

Defaults: `hip_iad_jitter` (excess per-transit dispersion, added in
quadrature), `iad_Δra`, `iad_Δdec`, `iad_Δpmra`, `iad_Δpmdec` and
`iad_Δplx = 0`, with `iad_pmra`/`iad_pmdec` derived by adding the offsets to
the Hipparcos catalog proper motion.
"""
struct HipparcosIADObs{TTable,THipSol,THost,TComp,TRef} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    hip_sol::THipSol
    name::String
    host::THost
    companions::TComp
    ref::TRef
end

export HipparcosIADObs

likelihoodname(obs::HipparcosIADObs) = obs.name
refspecs(obs::HipparcosIADObs) = (obs.host, obs.companions..., obs.ref)
_refdesc(obs::HipparcosIADObs) = _blend_refdesc(obs.host, obs.companions, obs.ref)

function HipparcosIADObs(;
        host,
        companions=(),
        ref=Barycentre,
        hip_id=nothing,
        iad=nothing,
        catalog=nothing,
        renormalize::Bool=true,
        attempt_correction::Bool=true,
        is_van_leeuwen::Bool=true,
        recalibrate::Bool=false,
        variables::Union{Nothing,Tuple{Priors,Derived}}=nothing,
        name::AbstractString="Hipparcos IAD",
    )
    isnothing(hip_id) && isnothing(iad) && error(
        "HipparcosIADObs needs either `hip_id` (to load the intermediate astrometric " *
        "data) or `iad=(; table, hip_sol)` (an already-loaded one).")
    loaded = isnothing(iad) ?
             hipparcos_iad(; hip_id, catalog, renormalize, attempt_correction,
                 is_van_leeuwen) : iad
    table = loaded.table
    hip_sol = loaded.hip_sol
    if recalibrate
        # Copy the columns first. `FlexTable(::Table)` wraps the *same* column
        # vectors, so recalibrating in place would shift the caller's own
        # table by +0.140 mas — and a second observation built from the same
        # preloaded `iad=` would then shift it again.
        ft = FlexTable(map(copy, Tables.columntable(table)))
        hipparcos_recalibrate!(ft)
        table = Table(ft)
    end

    hostspec = refspec(host)
    compspecs = map(refspec, Tuple(companions))
    refspec_ = refspec(ref)

    if isnothing(variables)
        variables = @variables begin
            hip_iad_jitter ~ LogUniform(0.001, 100)
            iad_Δra ~ Uniform(-1000, 1000)
            iad_Δdec ~ Uniform(-1000, 1000)
            # Fixed at zero by default: the parallax anchor is the system's own
            # `plx`, so leaving this free would turn the one thing a
            # Hipparcos-only fit is best at measuring into a nuisance.
            iad_Δplx = 0.0
            iad_Δpmra ~ Uniform(-1000, 1000)
            iad_Δpmdec ~ Uniform(-1000, 1000)
            iad_pmra = $(hip_sol.pm_ra) + iad_Δpmra
            iad_pmdec = $(hip_sol.pm_de) + iad_Δpmdec
        end
    end
    (priors, derived) = variables

    return HipparcosIADObs{typeof(table),typeof(hip_sol),typeof(hostspec),
        typeof(compspecs),typeof(refspec_)}(
        table, priors, derived, hip_sol, String(name), hostspec, compspecs, refspec_)
end

# Row subsets are transits, so this is the plain table restriction — no
# channel bookkeeping, unlike `G23HObs`.
# `table[inds]`, not `table[inds, :]`: the latter returns a 2-D table, which no
# longer matches the field's type parameter.
function likeobj_from_epoch_subset(obs::HipparcosIADObs, obs_inds)
    table = obs.table[obs_inds]
    return HipparcosIADObs{typeof(table),typeof(obs.hip_sol),typeof(obs.host),
        typeof(obs.companions),typeof(obs.ref)}(
        table, obs.priors, obs.derived, obs.hip_sol, obs.name,
        obs.host, obs.companions, obs.ref)
end

# ──────────────────────────────────────────────────────────────────────
# The forward model
# ──────────────────────────────────────────────────────────────────────

"""
    _hipiad_model!(resid, σ_infl, Δα, Δδ, obs, ctx, T)

Fill `resid` with the signed along-scan residual `measured − model` [mas] per
transit, and `σ_infl` with the BINARYS first-harmonic σ-inflation factor.
`Δα`/`Δδ` receive the source's own sky excursion, which is the Hipparcos
grating response of the host and its companions — *not* a photocentre, see
`_hippacentre!`.

Buffers must arrive zeroed (`Δα`, `Δδ`) and at one (`σ_infl`); they are
accumulated into, as `_hippacentre!` requires.
"""
function _hipiad_model!(resid, σ_infl, Δα, Δδ, obs::HipparcosIADObs,
                        ctx::ObsContext, ::Type{T}) where {T}
    θ_obs = ctx.θ_obs
    sys = ctx.system
    n = length(obs.table.epoch)

    # Resolve every reference once, outside the transit loop.
    hostref = ref(ctx, obs.host)
    reference = ref(ctx, obs.ref)
    comps = resolverefs(ctx, obs.companions)
    cidx = map(c -> c.idx, comps)
    masses = sys.masses
    # Test the primal: a differentiated zero mass is a Dual whose value is
    # zero but whose partials are not, so `iszero` on it is false.
    active = map(c -> !iszero(PlanetOrbits._primal(masses[c.idx])), comps)

    f_hip = _g23h_fluxratios(θ_obs, ctx.θ_system, :fluxratio_hip, Val(:Hp), sys,
        hostref.idx, cidx, active, T)
    _hippacentre!(Δα, Δδ, σ_infl, ctx, 1:n, obs.table.cosϕ, obs.table.sinϕ,
        hostref, reference, comps, f_hip, active, HIPPARCOS_GRID_STEP_ARCSEC)

    # The abscissae are on Hipparcos' frame, not the model's, so they get the
    # five-parameter frame offset of `skypath.jl`. The parallax anchors on the
    # system's own `plx` — with `iad_Δplx = 0` by default, that is what makes
    # these data constrain the parallax rather than a nuisance around it.
    plx_anchor = hasproperty(ctx.θ_system, :plx) ? ctx.θ_system.plx : obs.hip_sol.plx
    off = frame_offset(θ_obs, plx_anchor, T)

    inv_jy = inv(julian_year)
    ep = obs.table.epoch
    cϕ = obs.table.cosϕ
    sϕ = obs.table.sinϕ
    pf = obs.table.parallaxFactorAlongScan
    meas = obs.table.proj_meas_alongscan
    @inbounds @simd for i in eachindex(ep, resid)
        Δt = (ep[i] - hipparcos_catalog_epoch_mjd) * inv_jy
        resid[i] = meas[i] - frame_offset_alongscan(off, Δt, cϕ[i], sϕ[i], pf[i],
            Δα[i], Δδ[i])
    end
    return nothing
end

"""
    simulate(obs::HipparcosIADObs, ctx) -> NamedTuple

Per-transit model quantities at the sample in `ctx`: `resid` (signed
`measured − model` along-scan residual [mas]), `σ_inflation` (the BINARYS
first-harmonic factor multiplying each transit's formal σ), and the source's
own sky excursion `Δα_mas`/`Δδ_mas`.
"""
function simulate(obs::HipparcosIADObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    n = length(obs.table.epoch)
    resid = zeros(T, n)
    σ_inflation = ones(T, n)
    Δα_mas = zeros(T, n)
    Δδ_mas = zeros(T, n)
    _hipiad_model!(resid, σ_inflation, Δα_mas, Δδ_mas, obs, ctx, T)
    return (; resid, σ_inflation, Δα_mas, Δδ_mas)
end

function ln_like(obs::HipparcosIADObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    ll = zero(T)
    @no_escape ctx.buf begin
        n = length(obs.table.epoch)
        resid = @alloc(T, n); fill!(resid, zero(T))
        Δα = @alloc(T, n); fill!(Δα, zero(T))
        Δδ = @alloc(T, n); fill!(Δδ, zero(T))
        σ_infl = @alloc(T, n); fill!(σ_infl, one(T))
        _hipiad_model!(resid, σ_infl, Δα, Δδ, obs, ctx, T)

        jitter = hasproperty(ctx.θ_obs, :hip_iad_jitter) ?
                 T(ctx.θ_obs.hip_iad_jitter) : zero(T)
        jitter² = jitter * jitter
        half_log2π = T(0.5) * log(T(2π))
        rej = obs.table.reject
        sres = obs.table.sres_renorm
        @inbounds for i in 1:n
            rej[i] && continue
            # The per-transit σ is inflated by the BINARYS first-harmonic
            # factor: 1 in the unresolved limit, growing where the binary
            # modulation reduces the signal amplitude — and the residual noise
            # scales the same way.
            s = sres[i] * σ_infl[i]
            σ² = s * s + jitter²
            r = resid[i]
            ll += -T(0.5) * (r * r / σ² + log(σ²)) - half_log2π
        end
    end
    return isnan(ll) ? convert(T, -Inf) : ll
end

"""
    generate_from_params(obs::HipparcosIADObs, ctx; add_noise)

A new `HipparcosIADObs` whose abscissa residuals are the ones this sample
predicts. `res` is shifted by the model residual and `proj_meas_alongscan` is
rebuilt from it, so the regenerated observation is evaluated by exactly the
same code path — the reconstructed reference sky path and the scan geometry
are data and do not move.
"""
function generate_from_params(obs::HipparcosIADObs, ctx::ObsContext; add_noise)
    sim = simulate(obs, ctx)
    n = length(obs.table.epoch)
    jitter = hasproperty(ctx.θ_obs, :hip_iad_jitter) ?
             Float64(ctx.θ_obs.hip_iad_jitter) : 0.0
    # r = meas − model and meas = res + (reference path)·scan, so putting the
    # data on the model is res − r; the reference-path term is unchanged.
    res = collect(Float64, obs.table.res)
    @inbounds for i in 1:n
        res[i] -= Float64(sim.resid[i])
        if add_noise
            res[i] += randn() *
                      hypot(obs.table.sres_renorm[i] * Float64(sim.σ_inflation[i]), jitter)
        end
    end
    cols = Tables.columntable(obs.table)
    proj = @. res + cols.Δα✱ * cols.cosϕ + cols.Δδ * cols.sinϕ
    table = Table(merge(cols, (; res, proj_meas_alongscan=proj)))
    return HipparcosIADObs{typeof(table),typeof(obs.hip_sol),typeof(obs.host),
        typeof(obs.companions),typeof(obs.ref)}(
        table, obs.priors, obs.derived, obs.hip_sol, obs.name,
        obs.host, obs.companions, obs.ref)
end

# Along-scan residuals are `measured − model`, so a change in the model shows
# up in them one-for-one.
has_correction_impact(::Type{<:HipparcosIADObs}) = true
correction_impact(obs::HipparcosIADObs, a::ObsContext, b::ObsContext) =
    _simulate_impact(simulate(obs, a), simulate(obs, b), (:resid,),
                     _tightest(obs.table.sres_renorm))
