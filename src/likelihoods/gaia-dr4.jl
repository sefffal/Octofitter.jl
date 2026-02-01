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
    gaia_id::Int
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
    if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
        table = Table(table; epoch=jd2mjd.(table.obs_time_tcb))
    end
    xyz = Table(Octofitter.geocentre_position_query.(table.epoch))
    table = Table(table; xyz)
    gaia_sol = Octofitter._query_gaia_dr3(; gaia_id=gaia_id)
    return GaiaDR4AstromObs{typeof(table),typeof(gaia_sol)}(table,gaia_id, gaia_sol, priors, derived, name)
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
        name=obs.name
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
        # Allocate memory to simulate the along scan residuals
        ra_offset_buffer = @alloc(T, size(likeobj.table,1))
        dec_offset_buffer = @alloc(T, size(likeobj.table,1))
        centroid_pos_al_model_buffer = @alloc(T, size(likeobj.table,1))
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

        # Faster to compute both sine and cos at the same time
        s, c = sincos(likeobj.table.scan_pos_angle[i])
        
        # fast path, not accounting for higher order effects
        if !(orbitsol isa PlanetOrbits.OrbitSolutionAbsoluteVisual)
            # Our orbit solutions calculate an updated Ra and Dec for the system's barycentre
            # with various non-linear corrections (see PlanetOrbits.jl for more details)
            # Careful of units: proper motion is defined in mas per *Julian year*
            ra_offset_buffer[i] = (θ_obs.ra_offset_mas + θ_obs.pmra*(likeobj.table.epoch[i]-θ_obs.ref_epoch)/365.25)
            dec_offset_buffer[i] = (θ_obs.dec_offset_mas + θ_obs.pmdec*(likeobj.table.epoch[i]-θ_obs.ref_epoch)/365.25)

        else
            # Our orbit solutions calculate an updated Ra and Dec for the system's barycentre
            # with various non-linear corrections (see PlanetOrbits.jl for more details)
            α = orbitsol.compensated.ra2
            δ = orbitsol.compensated.dec2
            plx_at_epoch = orbitsol.compensated.parallax2 # Parallax distance may be changing

            # Faster to compute both sine and cos at the same time
            ra_offset_buffer[i] = (α - likeobj.gaia_sol.ra)*60*60*1000*cosd(δ)
            dec_offset_buffer[i] = (δ - likeobj.gaia_sol.dec)*60*60*1000

        end

        # Add perturbations from all planets
        # TODO: these perturbations are added with a linear approximation. 
        # For truly huge companions or extreme accuracy, an angle formula should be used instead
        for planet_i in eachindex(orbits)
            sol = orbit_solutions[planet_i][i+orbit_solutions_i_epoch_start]
            # Add perturbation from planet
            ra_offset_buffer[i] += raoff(sol, θ_system.planets[planet_i].mass * mjup2msol)
            dec_offset_buffer[i] += decoff(sol, θ_system.planets[planet_i].mass * mjup2msol)

        end
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
        obs.name
    )
end