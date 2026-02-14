# const DR4_REFERENCE_EPOCH = 2457936.875 # J2017.5 

"""
    GaiaDR4Astrom(
        observations_table,
        gaia_id=1234567890,
        variables=@variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)  # mas
        end,
        name="GaiaDR4",
        primary_star_perturbation=false
    )

A likelihood for fitting Gaia DR4 Individual Astronomical Data (IAD).
This includes along-scan astrometric measurements from Gaia.

The `astrometric_jitter` variable should be defined in the variables block
to account for systematic astrometric uncertainties.

If `primary_star_perturbation=true`, the observation variables (`ra_offset_mas`,
`pmra`, `pmdec`) are interpreted as the primary star's motion rather than the
barycentre's. The linear component of the companion perturbation is analytically
removed so that only the non-linear residual (acceleration/curvature) enters the
model. This breaks the degeneracy between proper motion and companion parameters
for wide-orbit companions.
"""
struct GaiaDR4AstromObs{TTable<:Table,TSol<:NamedTuple} <: Octofitter.AbstractObs
    table::TTable
    gaia_id::Int
    gaia_sol::TSol
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
    primary_star_perturbation::Bool
    detrend_Δt::Vector{Float64}
    detrend_inv_N::Float64
    detrend_inv_sum_Δt²::Float64
end
const GaiaDR4Astrom = GaiaDR4AstromObs
export GaiaDR4AstromObs, GaiaDR4Astrom

function GaiaDR4AstromObs(
    observations_table;
    gaia_id,
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin end),
    name="GaiaDR4",
    primary_star_perturbation::Bool=false
)
    (priors, derived) = variables
    table = Table(observations_table)
    if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
        table = Table(table; epoch=jd2mjd.(table.obs_time_tcb))
    end
    xyz = Table(Octofitter.geocentre_position_query.(table.epoch))
    table = Table(table; xyz)
    gaia_sol = Octofitter._query_gaia_dr3(; gaia_id=gaia_id)

    # Precompute detrending coefficients for primary_star_perturbation mode.
    # These allow O(N) linear detrending of the companion perturbation at each
    # MCMC step, removing the constant+slope components that are degenerate with
    # the fitted position and proper motion.
    epochs = table.epoch
    mean_epoch = sum(epochs) / length(epochs)
    Δt = collect((epochs .- mean_epoch) ./ 365.25)
    inv_N = 1.0 / length(epochs)
    inv_sum_Δt² = 1.0 / sum(Δt .^ 2)

    return GaiaDR4AstromObs{typeof(table),typeof(gaia_sol)}(
        table, gaia_id, gaia_sol, priors, derived, name,
        primary_star_perturbation, Δt, inv_N, inv_sum_Δt²
    )
end
function Octofitter.likeobj_from_epoch_subset(obs::GaiaDR4AstromObs, obs_inds)
    
    # Due to TypedTables bug, the line below creates a "matrix" table that isn't the same type as the input.
    # table = typeof(obs.table)(obs.table[setdiff(1:size(obs.table,1), obs_inds),:,1])
    # table = Table(collect(eachrow(obs.table))[setdiff(1:size(obs.table,1), obs_inds)]...)
    table = Table(first(eachcol(obs.table[setdiff(1:size(obs.table,1), obs_inds)])))
    return GaiaDR4AstromObs(
        table;
        gaia_id=obs.gaia_id,
        variables=(obs.priors, obs.derived),
        name=obs.name,
        primary_star_perturbation=obs.primary_star_perturbation
    )

end
#  likelihood function
function Octofitter.ln_like(
    likeobj::GaiaDR4AstromObs,
    ctx::SystemObservationContext
)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)


    # We separate out the the likelihood function into a `simulate` phase, 
    # and then use the results here to compute the lieklihood.
    # This way we can re-use the simulator for plots, generating fake data,
    # etc.

    # Get astrometric jitter from observation variables
    astrometric_jitter = hasproperty(θ_obs, :astrometric_jitter) ? θ_obs.astrometric_jitter : zero(T)
    astrometric_var = astrometric_jitter^2

    # Use Bumper.jl to avoid memory allocation overhead
    @no_escape begin
        N = size(likeobj.table,1)
        # Allocate memory to simulate the along scan residuals
        ra_offset_buffer = @alloc(T, N)
        dec_offset_buffer = @alloc(T, N)
        centroid_pos_al_model_buffer = @alloc(T, N)
        # Extra buffers for primary_star_perturbation detrending
        pert_ra_buffer = likeobj.primary_star_perturbation ? @alloc(T, N) : nothing
        pert_dec_buffer = likeobj.primary_star_perturbation ? @alloc(T, N) : nothing
        centroid_pos_al_model = Octofitter.simulate(
            likeobj,
            θ_system,
            θ_obs,
            orbits,
            orbit_solutions,
            orbit_solutions_i_epoch_start,
            ra_offset_buffer,
            dec_offset_buffer,
            centroid_pos_al_model_buffer,
            pert_ra_buffer,
            pert_dec_buffer,
        )
        # TODO: do we want to fit for a correlation parameter within each visibility window?
        # Extract the along-scan residuals from the NamedTuple returned by simulate
        along_scan_model = centroid_pos_al_model.along_scan_residuals_buffer
        for i in eachindex(likeobj.table.centroid_pos_al, along_scan_model)
            # There's an "outlier flag" -- I assume we should ignore ones that are flagged?
            if hasproperty(likeobj.table, :outlier_flag) && likeobj.table.outlier_flag[i] > 0
                continue
            end
            σ = sqrt(astrometric_var + likeobj.table.centroid_pos_error_al[i]^2)
            ll += logpdf(Normal(likeobj.table.centroid_pos_al[i], σ), along_scan_model[i])
        end
    end

    return ll
end

function Octofitter.simulate(
    likeobj::GaiaDR4AstromObs,
    θ_system,
    θ_obs,
    orbits,
    orbit_solutions,
    orbit_solutions_i_epoch_start,
    ra_offset_buffer=zeros(size(likeobj.table.epoch)),
    dec_offset_buffer=zeros(size(likeobj.table.epoch)),
    along_scan_residuals_buffer=zeros(size(likeobj.table.epoch)),
    pert_ra_buffer=likeobj.primary_star_perturbation ? zeros(size(likeobj.table.epoch)) : nothing,
    pert_dec_buffer=likeobj.primary_star_perturbation ? zeros(size(likeobj.table.epoch)) : nothing,
)

    T = Octofitter._system_number_type(θ_system)

    # All planets in the system have orbits defined with the same ra, dec, and proper motion,
    # since these are properties of the system.
    orbit = first(orbits)
    if length(orbits) > 1
        for i in eachindex(orbits)[2:end]
            if orbits[i].ra != orbit.ra ||
               orbits[i].dec != orbit.dec ||
               orbits[i].pmra != orbit.pmra ||
               orbits[i].pmdec != orbit.pmdec
                error("Planet orbits do not have matching ra, dec, pmra, and pmdec.")
            end
        end
    end

    # Compute position + proper motion at each epoch
    for i in eachindex(likeobj.table.epoch)

        orbitsol = first(orbit_solutions)[i+orbit_solutions_i_epoch_start]

        # fast path, not accounting for higher order effects
        if !(orbitsol isa PlanetOrbits.OrbitSolutionAbsoluteVisual)
            # Careful of units: proper motion is defined in mas per *Julian year*
            ra_offset_buffer[i] = (θ_obs.ra_offset_mas + θ_obs.pmra*(likeobj.table.epoch[i]-θ_obs.ref_epoch)/365.25)
            dec_offset_buffer[i] = (θ_obs.dec_offset_mas + θ_obs.pmdec*(likeobj.table.epoch[i]-θ_obs.ref_epoch)/365.25)

        else
            # Our orbit solutions calculate an updated Ra and Dec for the system's barycentre
            # with various non-linear corrections (see PlanetOrbits.jl for more details)
            α = orbitsol.compensated.ra2
            δ = orbitsol.compensated.dec2
            plx_at_epoch = orbitsol.compensated.parallax2 # Parallax distance may be changing

            ra_offset_buffer[i] = (α - likeobj.gaia_sol.ra)*60*60*1000*cosd(δ)
            dec_offset_buffer[i] = (δ - likeobj.gaia_sol.dec)*60*60*1000

        end
    end

    # Add perturbations from all planets.
    # In primary_star_perturbation mode, we remove the linear trend (constant + slope)
    # of the perturbation so that only the non-linear residual (curvature/acceleration)
    # enters the model. This breaks the degeneracy between proper motion and companion
    # parameters for wide orbits.
    if likeobj.primary_star_perturbation
        # Pass 1: accumulate raw perturbations and running sums for linear fit
        sum_pert_ra = zero(T);  dot_pert_ra = zero(T)
        sum_pert_dec = zero(T); dot_pert_dec = zero(T)
        for i in eachindex(likeobj.table.epoch)
            pert_ra = zero(T)
            pert_dec = zero(T)
            for planet_i in eachindex(orbits)
                sol = orbit_solutions[planet_i][i+orbit_solutions_i_epoch_start]
                pert_ra += raoff(sol, θ_system.planets[planet_i].mass * mjup2msol)
                pert_dec += decoff(sol, θ_system.planets[planet_i].mass * mjup2msol)
            end
            pert_ra_buffer[i] = pert_ra
            pert_dec_buffer[i] = pert_dec
            sum_pert_ra += pert_ra;   dot_pert_ra += likeobj.detrend_Δt[i] * pert_ra
            sum_pert_dec += pert_dec; dot_pert_dec += likeobj.detrend_Δt[i] * pert_dec
        end

        # Best-fit linear coefficients: pert ≈ mean + slope * Δt
        mean_pert_ra  = sum_pert_ra * likeobj.detrend_inv_N
        slope_pert_ra = dot_pert_ra * likeobj.detrend_inv_sum_Δt²
        mean_pert_dec  = sum_pert_dec * likeobj.detrend_inv_N
        slope_pert_dec = dot_pert_dec * likeobj.detrend_inv_sum_Δt²

        # Pass 2: add only the non-linear residual
        for i in eachindex(likeobj.table.epoch)
            ra_offset_buffer[i] += pert_ra_buffer[i] - mean_pert_ra - slope_pert_ra * likeobj.detrend_Δt[i]
            dec_offset_buffer[i] += pert_dec_buffer[i] - mean_pert_dec - slope_pert_dec * likeobj.detrend_Δt[i]
        end
    else
        # Original mode: add full perturbation (barycentric parameterization)
        for i in eachindex(likeobj.table.epoch)
            for planet_i in eachindex(orbits)
                sol = orbit_solutions[planet_i][i+orbit_solutions_i_epoch_start]
                ra_offset_buffer[i] += raoff(sol, θ_system.planets[planet_i].mass * mjup2msol)
                dec_offset_buffer[i] += decoff(sol, θ_system.planets[planet_i].mass * mjup2msol)
            end
        end
    end

    # Project onto along-scan direction and add parallax
    for i in eachindex(likeobj.table.epoch)
        s, c = sincos(likeobj.table.scan_pos_angle[i])
        along_scan_residuals_buffer[i] =
            ra_offset_buffer[i]*s +
            dec_offset_buffer[i]*c +
            θ_system.plx*likeobj.table.parallax_factor_al[i]
    end

    return (;
        along_scan_residuals_buffer,      # η = Δα* sin(ψ) + Δδ cos(ψ) — what we currently return
        ra_offset_buffer,       # Δα* in mas (total: PM + parallax + planets)
        dec_offset_buffer,      # Δδ in mas (total: PM + parallax + planets)
    )
end



# Generate new astrometry observations
function Octofitter.generate_from_params(obs::GaiaDR4AstromObs, ctx::SystemObservationContext; add_noise)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx
    sim_result = simulate(obs, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    # Extract the along-scan residuals from the NamedTuple returned by simulate
    along_scan_residuals = sim_result.along_scan_residuals_buffer

    new_table = deepcopy(obs.table)
    if add_noise
        for i in eachindex(new_table.centroid_pos_al)
            σ = new_table.centroid_pos_error_al[i]
            new_table.centroid_pos_al[i] = along_scan_residuals[i] + randn() * σ
        end
    else
        new_table.centroid_pos_al .= along_scan_residuals
    end

    # Use inner constructor to preserve gaia_sol and avoid network query
    return GaiaDR4AstromObs{typeof(new_table), typeof(obs.gaia_sol)}(
        new_table,
        obs.gaia_id,
        obs.gaia_sol,
        obs.priors,
        obs.derived,
        obs.name,
        obs.primary_star_perturbation,
        obs.detrend_Δt,
        obs.detrend_inv_N,
        obs.detrend_inv_sum_Δt²,
    )
end