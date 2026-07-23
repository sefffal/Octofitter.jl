#=
GaiaDR3CompanionObs — a planet-level likelihood for a *resolved, bound stellar
companion* that has its own published Gaia DR3 five-parameter astrometric
solution (position, parallax, proper motion).

Physical picture
================
The user has a system consisting of

  * a (possibly dark) central object — modelled as the **host star** of the
    `System`, whose mass (`system.M`) and barycentric astrometry
    (`ra, dec, plx, pmra, pmdec`) are the quantities we want to constrain; and
  * a visible, gravitationally bound companion that Gaia has *resolved* and
    published a five-parameter solution for — modelled as a `Planet` orbiting
    the host.

Because the companion is bound, its measured Gaia astrometry is the sum of the
system barycentre's space motion and the companion's reflex orbit about the
barycentre. Folding in the published Gaia solution therefore constrains the
companion's orbit and, given an assumed companion mass, the mass (and location)
of the central object — even when the central object is completely dark.

This is the photometric/astrometric analogue of measuring an unseen mass from
the motion of a visible companion (cf. the Hipparcos-Gaia acceleration method,
but using a *single, fully resolved* secondary's own catalogue entry).

Interface notes (answering "how do the system-level astrometric variables reach
a planet-level likelihood?")
================
In Octofitter each planet's orbit is constructed as

    orbit_i = basis(; merge(θ_system, θ_system.planets[i])...)

so the system-level `ra, dec, plx, pmra, pmdec, ref_epoch, rv` are baked into
the per-planet `AbsoluteVisual` orbit. A planet-level likelihood receives that
fully-constructed orbit (and its pre-solved solutions) plus `θ_system`, so it
can predict the companion's *absolute* sky path with no extra plumbing. This is
why the observation can simply be attached to the `Planet` that represents the
companion.

The orbit basis **must** be `AbsoluteVisual{...}` so that the rigorous 3D
non-linear stellar-motion propagation (perspective acceleration, changing
parallax, light-travel time) is engaged.

The companion mass enters through the planet's standard `mass` variable (Jupiter
masses, as elsewhere in Octofitter); the assumed-mass prior is therefore placed
in the `Planet`'s `@variables` block (or threaded from a `System` variable). The
total orbital mass `system.M` should be defined as `M_central + M_companion`.
=#

const _GAIA_DR3_CORR_FIELDS = (
    :ra_dec_corr, :ra_parallax_corr, :ra_pmra_corr, :ra_pmdec_corr,
    :dec_parallax_corr, :dec_pmra_corr, :dec_pmdec_corr,
    :parallax_pmra_corr, :parallax_pmdec_corr, :pmra_pmdec_corr,
)

"""
    GaiaDR3CompanionObs(; gaia_id_dr3=..., kwargs...)

A planet-level likelihood that folds in the published Gaia DR3 five-parameter
astrometric solution of a **resolved, bound stellar companion**. Attach it to
the `Planet` that represents the companion. The planet's orbit basis must be
`AbsoluteVisual{...}`.

# Specifying the data
Provide *either*

  * `gaia_id_dr3` — the Gaia DR3 `source_id` of the companion. The full
    catalogue row is queried from the Gaia archive (cached locally), and the
    published position/parallax/proper-motion and their full 5×5 covariance are
    used automatically; **or**

  * the five parameters and their uncertainties directly:
    `ra`, `dec` (degrees), `parallax`, `pmra`, `pmdec` (mas, mas/yr), with
    `ra_error`, `dec_error` (mas; these are σ of α*=α·cosδ and δ as Gaia
    reports), `parallax_error`, `pmra_error`, `pmdec_error`. Correlation
    coefficients (`ra_dec_corr`, …) are optional and default to 0.

# Keyword arguments
- `inflation_factor=1.37`: multiplicative inflation applied to all published
  astrometric uncertainties before forming the likelihood (Brandt 2019). Set to
  `1.0` to use the catalogue uncertainties as-is.
- `include_rv=false`: if `true`, additionally constrain the companion's
  barycentric radial velocity against the Gaia DR3 `radial_velocity`
  (km/s). Requires `radial_velocity`/`radial_velocity_error` to be available
  (queried or supplied).
- `ref_epoch`: MJD of the catalogue solution. Defaults to the Gaia DR3 reference
  epoch (J2016.0).
- `rv_inflation_factor=1.0`: separate inflation for the RV uncertainty.
- `name="GaiaDR3"`.
- `variables`: optional observation-level `@variables` block. An
  `astrometric_jitter` variable (mas), if present, is added in quadrature to the
  position uncertainties.

# Example
```julia
companion = Planet(
    name = "B",
    basis = AbsoluteVisual{KepOrbit},
    observations = [
        GaiaDR3CompanionObs(gaia_id_dr3 = 1234567890123456789),
    ],
    variables = @variables begin
        mass ~ truncated(Normal(520, 50), lower=0)  # companion mass [Mjup] (assumed)
        a ~ LogUniform(1, 1000)
        e ~ Uniform(0, 0.99)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        M = system.M
        tp = θ_at_epoch_to_tperi(θ, 57388.5; M, e, a, i, ω, Ω)
    end
)

sys = System(
    name = "darkmass",
    companions = [companion],
    observations = [],
    variables = @variables begin
        M_central ~ LogUniform(0.1, 100)            # central (dark) mass [Msol] — constrained
        M = M_central + companions.B.mass*Octofitter.mjup2msol
        plx  ~ Uniform(1, 100)
        pmra  ~ Uniform(-500, 500)
        pmdec ~ Uniform(-500, 500)
        rv = 0.0
        ra  ~ Normal(<catalog_ra>,  0.1)
        dec ~ Normal(<catalog_dec>, 0.1)
        ref_epoch = 57388.5
    end
)
```
"""
struct GaiaDR3CompanionObs{TTable<:Table,TCat,TDist,TRV} <: AbstractObs
    table::TTable
    source_id::Union{Int,Nothing}
    gaia_sol::TCat              # NamedTuple of catalogue values (μ, σ, corr, rv)
    dist::TDist                 # 5-parameter MvNormal [Δα*, Δδ, ϖ, μα*, μδ] (inflated)
    σ_vec::SVector{5,Float64}   # inflated σ (for jitter path)
    C_mat::SMatrix{5,5,Float64,25}
    include_rv::Bool
    rv_dist::TRV                # Normal or Nothing
    ref_epoch::Float64          # MJD of catalogue solution
    fd_step_days::Float64       # finite-difference half-step for proper motion
    priors::Priors
    derived::Derived
    name::String
end
const GaiaDR3CompanionLikelihood = GaiaDR3CompanionObs
export GaiaDR3CompanionObs, GaiaDR3CompanionLikelihood

function GaiaDR3CompanionObs(;
    gaia_id_dr3=nothing,
    ra=nothing, dec=nothing, parallax=nothing, pmra=nothing, pmdec=nothing,
    ra_error=nothing, dec_error=nothing, parallax_error=nothing,
    pmra_error=nothing, pmdec_error=nothing,
    radial_velocity=nothing, radial_velocity_error=nothing,
    ref_epoch=nothing,
    inflation_factor=1.37,
    rv_inflation_factor=1.0,
    include_rv=false,
    fd_step_days=180.0,
    variables::Tuple{Priors,Derived}=(@variables begin end),
    name="GaiaDR3",
    corr=NamedTuple(),
)
    (priors, derived) = variables

    if !isnothing(gaia_id_dr3)
        source_id = gaia_id_dr3
        cat = _query_gaia_dr3(; gaia_id=gaia_id_dr3)
        ra            = cat.ra
        dec           = cat.dec
        parallax      = cat.parallax
        pmra          = cat.pmra
        pmdec         = cat.pmdec
        ra_error      = cat.ra_error
        dec_error     = cat.dec_error
        parallax_error= cat.parallax_error
        pmra_error    = cat.pmra_error
        pmdec_error   = cat.pmdec_error
        corr = NamedTuple{_GAIA_DR3_CORR_FIELDS}(
            map(f -> getproperty(cat, f), _GAIA_DR3_CORR_FIELDS)
        )
        if include_rv
            if isnothing(radial_velocity) && hasproperty(cat, :radial_velocity)
                radial_velocity = cat.radial_velocity
            end
            if isnothing(radial_velocity_error) && hasproperty(cat, :radial_velocity_error)
                radial_velocity_error = cat.radial_velocity_error
            end
        end
    else
        source_id = nothing
        # User-specified solution; require the five parameters and uncertainties.
        for (nm, v) in (
            (:ra, ra), (:dec, dec), (:parallax, parallax), (:pmra, pmra), (:pmdec, pmdec),
            (:ra_error, ra_error), (:dec_error, dec_error),
            (:parallax_error, parallax_error), (:pmra_error, pmra_error),
            (:pmdec_error, pmdec_error),
        )
            if isnothing(v)
                throw(ArgumentError("Provide `gaia_id_dr3`, or supply all of ra, dec, parallax, pmra, pmdec and their *_error. Missing: $nm"))
            end
        end
    end

    if isnothing(ref_epoch)
        ref_epoch = meta_gaia_DR3.ref_epoch_mjd
    end

    # Helper to read a (possibly missing) correlation coefficient.
    getcorr(s) = haskey(corr, s) ? Float64(getproperty(corr, s)) : 0.0

    # Data vector ordering: [Δα* (mas), Δδ (mas), ϖ (mas), μα* (mas/yr), μδ (mas/yr)].
    # Position residuals are taken relative to the catalogue position, so the
    # data values for Δα*/Δδ are 0 and the model carries the predicted offset.
    μ = Float64[0.0, 0.0, parallax, pmra, pmdec]
    σ = inflation_factor .* Float64[ra_error, dec_error, parallax_error, pmra_error, pmdec_error]

    C = @SMatrix [
        1.0                      getcorr(:ra_dec_corr)    getcorr(:ra_parallax_corr)  getcorr(:ra_pmra_corr)      getcorr(:ra_pmdec_corr)
        getcorr(:ra_dec_corr)    1.0                      getcorr(:dec_parallax_corr) getcorr(:dec_pmra_corr)     getcorr(:dec_pmdec_corr)
        getcorr(:ra_parallax_corr) getcorr(:dec_parallax_corr) 1.0                   getcorr(:parallax_pmra_corr) getcorr(:parallax_pmdec_corr)
        getcorr(:ra_pmra_corr)   getcorr(:dec_pmra_corr)  getcorr(:parallax_pmra_corr) 1.0                        getcorr(:pmra_pmdec_corr)
        getcorr(:ra_pmdec_corr)  getcorr(:dec_pmdec_corr) getcorr(:parallax_pmdec_corr) getcorr(:pmra_pmdec_corr) 1.0
    ]
    σvec = SVector{5,Float64}(σ)
    Σ = Diagonal(σvec) * C * Diagonal(σvec)
    dist = MvNormal(μ, Hermitian(Matrix(Σ)))

    # Optional radial-velocity constraint.
    if include_rv
        if isnothing(radial_velocity) || isnothing(radial_velocity_error)
            throw(ArgumentError("`include_rv=true` but no Gaia radial_velocity / radial_velocity_error is available. Supply them or use a source with a published Gaia RV."))
        end
        rv_dist = Normal(Float64(radial_velocity), rv_inflation_factor * Float64(radial_velocity_error))
    else
        rv_dist = nothing
    end

    gaia_sol = (;
        ra, dec, parallax, pmra, pmdec,
        ra_error, dec_error, parallax_error, pmra_error, pmdec_error,
        radial_velocity, radial_velocity_error,
        ref_epoch,
    )

    # Three epochs around the reference epoch: central finite difference for the
    # instantaneous proper motion, central epoch for position+parallax.
    table = Table(epoch=[ref_epoch - fd_step_days, ref_epoch, ref_epoch + fd_step_days])

    return GaiaDR3CompanionObs(
        table, source_id, gaia_sol, dist, σvec, C,
        include_rv, rv_dist, Float64(ref_epoch), Float64(fd_step_days),
        priors, derived, name,
    )
end

function likelihoodname(obs::GaiaDR3CompanionObs)
    return obs.name
end

# Not meaningfully subsettable (a single catalogue solution); return as-is so
# cross-validation machinery does not error.
function likeobj_from_epoch_subset(obs::GaiaDR3CompanionObs, obs_inds)
    return obs
end

# Predicted offset (Δα* mas, Δδ mas, parallax mas) of the *companion* at a given
# orbit solution, relative to the catalogue reference position.
@inline function _dr3_companion_offset(sol, M_companion_msol, ref_ra, ref_dec)
    α = sol.compensated.ra2     # barycentre RA at epoch [deg]
    δ = sol.compensated.dec2    # barycentre Dec at epoch [deg]
    # Companion position relative to the barycentre [mas].
    # raoff(sol) is the companion relative to the host; raoff(sol, M_companion)
    # is the host relative to the barycentre; their sum is the companion
    # relative to the barycentre (= (M_host/M_tot)·separation).
    Δα_comp = raoff(sol) + raoff(sol, M_companion_msol)
    Δδ_comp = decoff(sol) + decoff(sol, M_companion_msol)
    # Barycentre offset from the catalogue reference position [mas].
    Δα_bary = (α - ref_ra) * 3.6e6 * cosd(ref_dec)
    Δδ_bary = (δ - ref_dec) * 3.6e6
    return (Δα_bary + Δα_comp, Δδ_bary + Δδ_comp, sol.compensated.parallax2)
end

"""
    simulate(obs::GaiaDR3CompanionObs, ...) -> NamedTuple

Predict the companion's Gaia five-parameter solution (and barycentric RV) from
the current orbit. Returns `(; Δα_mas, Δδ_mas, parallax, pmra, pmdec, rv_kms)`
where the first two are offsets from the catalogue position.
"""
function simulate(
    obs::GaiaDR3CompanionObs,
    θ_system, θ_planet, θ_obs,
    orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start,
)
    orbit = orbits[i_planet]
    if !(orbit isa PlanetOrbits.AbsoluteVisualOrbit)
        error("GaiaDR3CompanionObs requires an AbsoluteVisual{...} orbit basis so that the full 3D stellar-motion propagation is engaged (got $(typeof(orbit))).")
    end
    if !hasproperty(θ_planet, :mass)
        error("GaiaDR3CompanionObs requires the companion Planet to define a `mass` variable (Jupiter masses).")
    end

    M_companion_msol = θ_planet.mass * mjup2msol
    ref_ra = obs.gaia_sol.ra
    ref_dec = obs.gaia_sol.dec

    i0 = orbit_solutions_i_epoch_start
    sol_m = orbit_solutions[i_planet][1 + i0]   # ref_epoch - fd_step
    sol_0 = orbit_solutions[i_planet][2 + i0]   # ref_epoch
    sol_p = orbit_solutions[i_planet][3 + i0]   # ref_epoch + fd_step

    (Δα_m, Δδ_m, _)    = _dr3_companion_offset(sol_m, M_companion_msol, ref_ra, ref_dec)
    (Δα_0, Δδ_0, plx)  = _dr3_companion_offset(sol_0, M_companion_msol, ref_ra, ref_dec)
    (Δα_p, Δδ_p, _)    = _dr3_companion_offset(sol_p, M_companion_msol, ref_ra, ref_dec)

    # Central finite difference for instantaneous proper motion [mas/yr].
    dt_yr = (obs.table.epoch[3] - obs.table.epoch[1]) / julian_year
    pmra_model = (Δα_p - Δα_m) / dt_yr
    pmdec_model = (Δδ_p - Δδ_m) / dt_yr

    # Companion barycentric radial velocity [m/s] -> [km/s].
    rv_kms = if obs.include_rv
        M_tot = totalmass(orbit)
        M_host = M_tot - M_companion_msol
        kep_rv = radvel(sol_0.sol)  # secondary relative to host [m/s]
        (sol_0.compensated.rv2 + (M_host / M_tot) * kep_rv) / 1000
    else
        zero(plx)
    end

    return (; Δα_mas=Δα_0, Δδ_mas=Δδ_0, parallax=plx, pmra=pmra_model, pmdec=pmdec_model, rv_kms)
end

function ln_like(obs::GaiaDR3CompanionObs, ctx::PlanetObservationContext)
    (; θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start) = ctx
    T = _system_number_type(θ_system)

    sim = simulate(obs, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)

    model = SVector{5,T}(sim.Δα_mas, sim.Δδ_mas, sim.parallax, sim.pmra, sim.pmdec)

    # Optional astrometric jitter (mas), added in quadrature to position errors.
    jitter = hasproperty(θ_obs, :astrometric_jitter) ? θ_obs.astrometric_jitter : zero(T)

    ll = zero(T)
    if jitter == 0
        ll += logpdf(obs.dist, model)
    else
        σ = SVector{5,T}(
            hypot(obs.σ_vec[1], jitter),
            hypot(obs.σ_vec[2], jitter),
            obs.σ_vec[3], obs.σ_vec[4], obs.σ_vec[5],
        )
        Σ = Diagonal(σ) * obs.C_mat * Diagonal(σ)
        ll += logpdf(MvNormal(SVector{5,T}(0, 0, obs.gaia_sol.parallax, obs.gaia_sol.pmra, obs.gaia_sol.pmdec), Hermitian(Σ)), model)
    end

    if obs.include_rv && !isnothing(obs.rv_dist)
        ll += logpdf(obs.rv_dist, sim.rv_kms)
    end

    return ll
end

# Generate a synthetic catalogue solution from the current parameters (used for
# simulation-based testing / injection-recovery).
function generate_from_params(obs::GaiaDR3CompanionObs, ctx::PlanetObservationContext; add_noise)
    (; θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start) = ctx
    sim = simulate(obs, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)

    # Convert the predicted Δα*/Δδ offsets back to absolute catalogue ra/dec.
    ref_ra = obs.gaia_sol.ra
    ref_dec = obs.gaia_sol.dec
    new_ra  = ref_ra  + sim.Δα_mas / 3.6e6 / cosd(ref_dec)
    new_dec = ref_dec + sim.Δδ_mas / 3.6e6
    new_plx = sim.parallax
    new_pmra = sim.pmra
    new_pmdec = sim.pmdec

    σ = obs.σ_vec
    if add_noise
        new_ra  += randn() * σ[1] / 3.6e6 / cosd(ref_dec)
        new_dec += randn() * σ[2] / 3.6e6
        new_plx += randn() * σ[3]
        new_pmra += randn() * σ[4]
        new_pmdec += randn() * σ[5]
    end

    rv = obs.gaia_sol.radial_velocity
    rv_err = obs.gaia_sol.radial_velocity_error
    if obs.include_rv
        rv = sim.rv_kms + (add_noise ? randn() * (rv_err === nothing ? 0.0 : rv_err) : 0.0)
    end

    return GaiaDR3CompanionObs(;
        ra=new_ra, dec=new_dec, parallax=new_plx, pmra=new_pmra, pmdec=new_pmdec,
        ra_error=obs.gaia_sol.ra_error, dec_error=obs.gaia_sol.dec_error,
        parallax_error=obs.gaia_sol.parallax_error,
        pmra_error=obs.gaia_sol.pmra_error, pmdec_error=obs.gaia_sol.pmdec_error,
        radial_velocity=rv, radial_velocity_error=rv_err,
        ref_epoch=obs.ref_epoch,
        inflation_factor=1.0,  # σ already stored inflated; avoid double-inflation
        include_rv=obs.include_rv,
        fd_step_days=obs.fd_step_days,
        variables=(obs.priors, obs.derived),
        name=obs.name,
    )
end
