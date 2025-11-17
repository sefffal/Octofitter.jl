# const DR4_REFERENCE_EPOCH = 2457936.875 # J2017.5 

"""
    GaiaDR4Astrom(
        observations_table,
        gaia_id=1234567890,
        variables=@variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)  # mas
        end,
        name="GaiaDR4"
    )

A likelihood for fitting Gaia DR4 Individual Astronomical Data (IAD). 
This includes along-scan astrometric measurements from Gaia.

The `astrometric_jitter` variable should be defined in the variables block
to account for systematic astrometric uncertainties.
"""
struct GaiaDR4AstromObs{TTable<:Table,TSol<:NamedTuple} <: Octofitter.AbstractObs
    table::TTable
    gaia_sol::TSol
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
end
const GaiaDR4Astrom = GaiaDR4AstromObs
export GaiaDR4AstromObs, GaiaDR4Astrom

function GaiaDR4AstromObs(
    observations_table;
    gaia_id,
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin end),
    name="GaiaDR4"
)
    (priors, derived) = variables
    table = Table(observations_table)
    gaia_sol = Octofitter._query_gaia_dr3(; gaia_id=gaia_id)
    return GaiaDR4AstromObs{typeof(table),typeof(gaia_sol)}(table, gaia_sol, priors, derived, name)
end
function Octofitter.likeobj_from_epoch_subset(obs::GaiaDR4AstromObs, obs_inds)
    return GaiaDR4AstromObs(
        obs.table[obs_inds,:,1]...,
        gaia_id=obs.gaia_sol.source_id,
        variables=(obs.priors, obs.derived),
        name=obs.name
    )

end
#  likelihood function
function Octofitter.ln_like(
    likeobj::GaiaDR4AstromObs,
    (;θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)::SystemObservationContext
)
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
        # Allocate memory to simulate the along scan residuals
        centroid_pos_al_model_buffer = @alloc(T, size(likeobj.table,1))
        centroid_pos_al_model = Octofitter.simulate(
            likeobj,
            θ_system,
            θ_obs,
            orbits,
            orbit_solutions,
            orbit_solutions_i_epoch_start,
            centroid_pos_al_model_buffer
        )
        # TODO: do we want to fit for a correlation parameter within each visibility window?
        for i in eachindex(likeobj.table.centroid_pos_al, centroid_pos_al_model)
            # There's an "outlier flag" -- I assume we should ignore ones that are flagged?
            if likeobj.table.outlier_flag[i] > 0
                continue
            end
            σ = sqrt(astrometric_var + likeobj.table.centroid_pos_error_al[i]^2)
            ll += logpdf(Normal(likeobj.table.centroid_pos_al[i], σ), centroid_pos_al_model[i])
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
    along_scan_residuals_buffer=zeros(size(likeobj.table.epoch))
)

    T = Octofitter._system_number_type(θ_system)

    # All planets inthe system have orbits defined with the same ra, dec, and proper motion,
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


    for i in eachindex(likeobj.table.epoch)

        orbitsol = first(orbit_solutions)[i+orbit_solutions_i_epoch_start]

        # fast path, not accounting for higher order effects
        if !(orbitsol isa PlanetOrbits.OrbitSolutionAbsoluteVisual)
            # Our orbit solutions calculate an updated Ra and Dec for the system's barycentre
            # with various non-linear corrections (see PlanetOrbits.jl for more details)
            # Careful of units: proper motion is defined in mas per *Julian year*
            
            # Faster to compute both sine and cos at the same time
            s, c = sincosd(likeobj.table.scan_pos_angle[i])
            
            astrometric_measurement_mas = 
                θ_system.ra_offset_mas + θ_system.pmra*(likeobj.table.epoch[i]-θ_system.ref_epoch)/365.25*s + 
                θ_system.dec_offset_mas + θ_system.pmdec*(likeobj.table.epoch[i]-θ_system.ref_epoch)/365.25*c + 
                θ_system.plx*likeobj.table.parallax_factor_al[i]
        else
            # Our orbit solutions calculate an updated Ra and Dec for the system's barycentre
            # with various non-linear corrections (see PlanetOrbits.jl for more details)
            α = orbitsol.compensated.ra2
            δ = orbitsol.compensated.dec2
            plx_at_epoch = orbitsol.compensated.parallax2 # Parallax distance may be changing

            # Faster to compute both sine and cos at the same time
            s, c = sincosd(likeobj.table.scan_pos_angle[i])
            astrometric_measurement_mas = 
                (α - likeobj.gaia_sol.ra)*60*60*1000*cosd(δ)*s + 
                (δ - likeobj.gaia_sol.dec)*60*60*1000*c + 
                plx_at_epoch*likeobj.table.parallax_factor_al[i]
        end

        # Add perturbations from all planets
        # TODO: these perturbations are added with a linear approximation. 
        # For truly huge companions or extreme accuracy, an angle formula should be used instead
        for planet_i in eachindex(orbits)
            sol = orbit_solutions[planet_i][i+orbit_solutions_i_epoch_start]
            # Add perturbation from planet
            x = raoff(sol, θ_system.planets[planet_i].mass * Octofitter.mjup2msol)
            y = decoff(sol, θ_system.planets[planet_i].mass * Octofitter.mjup2msol)
            astrometric_measurement_mas += x*sind(likeobj.table.scan_pos_angle[i]) + y*cosd(likeobj.table.scan_pos_angle[i]) # yes, this convention is correct
        end
        along_scan_residuals_buffer[i] = astrometric_measurement_mas
    end

    return along_scan_residuals_buffer
end



# Generate new astrometry observations
function Octofitter.generate_from_params(like::GaiaDR4AstromObs, θ_system, θ_obs, orbits, solutions, sol_start_i)
    along_scan_residuals_buffer_sim = simulate(like, θ_system, θ_obs, orbits, solutions, sol_start_i)

    new_table = deepcopy(like.table)
    new_table.centroid_pos_al .= along_scan_residuals_buffer_sim

    return GaiaDR4AstromObs(
        new_table,
        gaia_id=like.gaia_sol.source_id,
        variables=(like.priors, like.derived),
        name=like.name
    )
end
