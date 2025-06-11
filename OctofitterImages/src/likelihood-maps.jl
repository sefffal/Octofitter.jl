
const likemap_cols = (:map, :epoch, :platescale,)

"""
    LogLikelihoodMap(
        table,
        name="likemap",
        variables=@variables begin
        end
    )

One or more maps of likelihood vs Δright-ascension and Δdeclination computed with some other tool.
The platescale column maps pixels in the matrix to separation in milliarcseconds.
Pass a Table with the following columns:
$likemap_cols

For example:
```julia
likemap_dat = Table(;
    epoch = [1234.0, 1584.7],
    map = [readfits("map1.fits"), readfits("map2.fits")],
    platescale = [19.4, 19.4]
)

LogLikelihoodMap(
    likemap_dat,
    name="GRAVITY",
    variables=@variables begin
        platescale = 1.0               # Platescale multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        northangle = 0.0               # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
    end
)
```
Epoch is in MJD.
Platescale is in mas/px.
"""
struct LogLikelihoodMap{TTable<:Table} <: Octofitter.AbstractLikelihood
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
    function LogLikelihoodMap(
        table;
        name::String="likemap",
        variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin end)
    )
        (priors, derived) = variables
        if !issubset(likemap_cols, columnnames(table))
            error("Expected columns $likemap_cols")
        end
        if !in(:fillvalue, columnnames(table))
            @info ":fillvalue not provided. Filling Inf/NaN with minimum log-like value"
            fillvalue = map(eachrow(table)) do row
                img = row[].map
                m = minimum(filter(isfinite, img))
                @info "row:" epoch=row[].epoch fillvalue=m
                return m
            end
            table = Table(table; fillvalue)
        end
        # Create linear interpolators over the input likemaps
        mapinterp = map(eachrow(table)) do row
            img = row[].map
            data = collect(img)
            data[.!isfinite.(data)] .= row[].fillvalue
            LinearInterpolation(parent.(dims(img)), data, extrapolation_bc=convert(eltype(img), row[].fillvalue))
        end
        table = Table(table; mapinterp)
        return new{typeof(table)}(table, priors, derived, name)
    end
end
# Legacy constructor for backward compatibility
LogLikelihoodMap(observations::NamedTuple...; kwargs...) = LogLikelihoodMap(Table(observations...); kwargs...)
export LogLikelihoodMap

function Octofitter.likeobj_from_epoch_subset(obs::LogLikelihoodMap, obs_inds)
    return LogLikelihoodMap(obs.table[obs_inds,:,1]...)
end

"""
Likelihood of there being planets in a sequence of likemaps.
"""
function Octofitter.ln_like(likemaps::LogLikelihoodMap, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)

    # Resolve the combination of system and planet parameters
    # as a Visual{KepOrbit} object. This pre-computes
    # some factors used in various calculations.
    # elements = construct_elements(θ_system, θ_planet)
    this_orbit = orbits[i_planet]


    likemaps_table = likemaps.table
    T = Octofitter._system_number_type(θ_planet)
    ll = zero(T)
    
    # Get calibration parameters from observation variables
    platescale_multiplier = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : one(T)
    northangle_offset = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : zero(T)
    
    for i_epoch in eachindex(likemaps_table.epoch)

        sol = orbit_solutions[i_planet][i_epoch+orbit_solutions_i_epoch_start]
        @assert isapprox(likemaps.table.epoch[i_epoch], PlanetOrbits.soltime(sol), rtol=1e-2)
        ra_host_perturbation = zero(T)
        dec_host_perturbation = zero(T)
        for (i_other_planet, key) in enumerate(keys(θ_system.planets))
            orbit_other = orbits[i_other_planet]
            # Only account for inner planets with non-zero mass
            if semimajoraxis(orbit_other) < semimajoraxis(this_orbit)
                θ_planet′ = θ_system.planets[key]
                if !hasproperty(θ_planet′, :mass)
                    continue
                end
                mass_other = θ_planet′.mass*Octofitter.mjup2msol
                sol′ = orbit_solutions[i_other_planet][i_epoch + orbit_solutions_i_epoch_start]
                # Note about `total mass`: for this to be correct, user will have to specify
                # `M` at the planet level such that it doesn't include the outer planets.

                ra_host_perturbation += raoff(sol′, mass_other)
                dec_host_perturbation += decoff(sol′, mass_other)

                @assert isapprox(likemaps.table.epoch[i_epoch], PlanetOrbits.soltime(sol), rtol=1e-2)
                @assert isapprox(likemaps.table.epoch[i_epoch], PlanetOrbits.soltime(sol′), rtol=1e-2)
            end
        end


        # Apply north angle rotation and platescale correction
        ra_raw = raoff(sol) - ra_host_perturbation
        dec_raw = decoff(sol) - dec_host_perturbation
        
        # Apply north angle rotation
        cos_θ = cos(northangle_offset)
        sin_θ = sin(northangle_offset)
        ra_rotated = ra_raw * cos_θ - dec_raw * sin_θ
        dec_rotated = ra_raw * sin_θ + dec_raw * cos_θ

        # Note the x reversal between RA and image coordinates
        x = -ra_rotated
        y = +dec_rotated

        # Get the log-likelihood at this position
        # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
        platescale = likemaps_table.platescale[i_epoch] * platescale_multiplier
        χ² = likemaps_table.mapinterp[i_epoch](x/platescale, y/platescale)

        # When we get a position that falls outside of our available
        # data (e.g. under the coronagraph) we cannot say anything
        # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
        # of zero.
        if !isfinite(χ²)
            ll += likemaps_table.fillvalue[i_epoch]
        end
        ll += χ²
    end

    return ll
end



# Generate new likemaps
# This isn't currently possible unless we can eg. find a way to plump the system log-likelihood into here.
# function Octofitter.generate_from_params(like::LogLikelihoodMap, θ_planet,  orbit::PlanetOrbits.AbstractOrbit)

#     newrows = map(like.table) do row
#         (;map, epoch, platescale) = row

#         injected = copy(map)
    
#         # Generate new astrometry point
#         os = orbitsolve(orbit, epoch)

#         ra = raoff(os)
#         dec = decoff(os)

#         phot = θ_planet[band]

#         # TODO: verify signs
#         dx = ra/platescale
#         dy = -dec/platescale
#         translation_tform = Translation(
#             mean(axes(psf,1))-mean(axes(image,1))+mean(dims(image,1))+dx,
#             mean(axes(psf,2))-mean(axes(image,2))+mean(dims(image,2))+dy
#         )
#         # TBD if we want to support rotations for handling negative sidelobes.

#         injected .+= psf_scaled

#         return merge(row, (;image=injected))
#     end

#     return LogLikelihoodMap(newrows)
# end
