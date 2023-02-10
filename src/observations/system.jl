



# # Overall log likelihood of the system given the parameters θ_system
# function ln_like(system::System, θ_system)
#     # Take some care to ensure type stability when using e.g. ForwardDiff
#     ll = zero(typeof(first(θ_system)))

#     # Fail fast if we have a negative stellar mass.
#     # Users should endeavour to use priors on e.g. stellar mass
#     # that are strictly positive, otherwise we are reduced to rejection sampling!
#     if hasproperty(θ_system, :M) && θ_system.M <= 0
#         return oftype(ll, -Inf)
#     end

#     # Go through each planet in the model and add its contribution
#     # to the ln-likelihood.
#     # out_of_bounds = Base.RefValue{Bool}(false)
#     elements = map(eachindex(system.planets)) do i
#         planet = system.planets[i]
#         θ_planet = θ_system.planets[i]

#         # Like negative stellar mass, users should use priors with supports
#         # that do not include these invalid values. But if we see them,
#         # give zero likelihood right away instead of an inscrutable error
#         # from some code expecting these invariants to hold.
#         if (hasproperty(θ_planet, :a) && θ_planet.a <= 0) ||
#             (hasproperty(θ_planet, :e) && !(0 <= θ_planet.e < 1))
#             out_of_bounds[] = true
#         end

#         # Resolve the combination of system and planet parameters
#         # as a orbit object. The type of orbitobject is stored in the 
#         # Planet type. This pre-computes some factors used in various calculations.
#         kep_elements = construct_elements(orbittype(planet), θ_system, θ_planet)

#         return kep_elements
#     end
#     # # Fail fast if out of bounds for one of the planets
#     # if out_of_bounds[]
#     #     return oftype(ll, -Inf) # Ensure type stability
#     # end

#     # Loop through the planets from the outside in. 
#     # Try to do this sorting in a non-allocating way.
#     # This way we have the option to account for each planets influence on the outer planets
#     # sma = map(elements) do elem
#     #     return elem.a
#     # end
#     # planet_sma_asc_ii = sortperm(SVector(sma))

#     # The above sorting is not currently used, so need to perform it.
#     planet_sma_asc_ii = 1:length(elements)


#     # Handle all observations attached to planets in order of semi-major axis
#     for j in eachindex(planet_sma_asc_ii)
#         i = planet_sma_asc_ii[j]
#         # Planet model
#         planet = system.planets[i]
#         # Parameters specific to this planet
#         θ_planet = θ_system.planets[i]
#         # Cached VisualOrbit with precomputed factors, etc.
#         planet_elements = elements[i]
#         # kep_elements, but for all planets interior to this one (given the current parameters)
#         # interior_planets = kep_elements[begin:min(end,i)]

#         # Loop through observations
#         for obs in planet.observations
#             ll += ln_like(obs, θ_planet, planet_elements, θ_system.planets, elements)
#         end
#     end

#     if !isfinite(ll)
#         return ll
#     end

#     # Loop through and add contribution of all observation types associated with this system as a whole
#     for obs in system.observations
#         ll += ln_like(obs, θ_system, elements)
#     end

#     return ll
# end


# # Overall log likelihood of the system given the parameters θ_system
# function ln_like(system::System, θ_system)
#     # Take some care to ensure type stability when using e.g. ForwardDiff
#     ll = zero(typeof(first(θ_system)))

#     # Should unroll this loop with RuntimeGeneratedFunction
#     # for i in eachindex(system.planets)
#     let i = 1
#         planet = system.planets[i]
#         θ_planet = θ_system.planets[i]

#         # Resolve the combination of system and planet parameters
#         # as a orbit object. The type of orbitobject is stored in the 
#         # Planet type. This pre-computes some factors used in various calculations.
#         planet_elements = construct_elements(orbittype(planet), θ_system, θ_planet)
#         # planet_elements = VisualOrbit((;
#         #     M=θ_system.M,
#         #     plx=θ_system.plx,
#         #     a=θ_planet.a,
#         #     i=θ_planet.i,
#         #     ω=θ_planet.ω,
#         #     Ω=θ_planet.Ω,
#         #     e=θ_planet.e,
#         #     τ=θ_planet.τ,
#         # ))
        
#         # Loop through observations ( could unroll this loop with RuntimeGeneratedFunction)
#         for obs in planet.observations
#             ll += ln_like(obs, θ_planet, planet_elements)
#         end
#     end
#     # Could unroll this loop

#     # # Loop through and add contribution of all observation types associated with this system as a whole
#     # for obs in system.observations
#     #     ll += ln_like(obs, θ_system)
#     # end

#     return ll
# end

function make_ln_like(system::System, θ_system)

    planet_exprs = Expr[]
    planet_keys = Symbol[]
    for i in eachindex(system.planets)
        planet = system.planets[i]
        θ_planet = θ_system.planets[i]
        key = Symbol("planet_$i")

        likelihood_exprs = map(eachindex(planet.observations)) do obs
            :(
                ll += ln_like(
                    system.planets[$(Meta.quot(i))].observations[$obs],
                    θ_system.planets.$i,
                    $(key)
                )
            )
        end

        planet_contruction = quote
            $key = $(construct_elements)($(orbittype(planet)), θ_system, θ_system.planets.$i)
            $(likelihood_exprs...)
        end
        push!(planet_exprs, planet_contruction)
        push!(planet_keys, key)
    end


    sys_exprs = map(system.observations) do obs
        :(ll += ln_like($obs, θ_system, elems))
    end

    return @RuntimeGeneratedFunction(:(function (system::System, θ_system)

        ll = zero(first(θ_system))

        # Construct all orbit elements and evaluate all their individual observation likelihoods
        $(planet_exprs...)

        # Construct a tuple of existing planet orbital elements
        elems = tuple($(planet_keys...))
        
        $(sys_exprs...)

        return ll
    end))


end




# Generate calibration data
"""
    generate(system, newparameters=drawfrompriors(system))

Generate a new system and observations from an existing model.
"""
function generate(system::System, θ_newsystem = drawfrompriors(system))

    # Generate new orbits for each planet in the system
    elements = map(eachindex(system.planets)) do i
        planet = system.planets[i]
        θ_newplanet = θ_newsystem.planets[i]

        if (hasproperty(θ_newplanet, :a) && θ_newplanet.a <= 0) ||
            (hasproperty(θ_newplanet, :e) && !(0 <= θ_newplanet.e < 1))
            out_of_bounds[] = true
        end

        neworbit = Octofitter.construct_elements(Octofitter.orbittype(planet), θ_newsystem, θ_newplanet)

        return neworbit
    end

    # Generate new observations for each planet in the system
    newplanets = map(1:length(system.planets)) do i
        planet = system.planets[i]
        elem = elements[i]

        newplanet_obs = map(planet.observations) do obs
            return genobs(obs, elem, planet)
        end
        newplanet = Planet{Octofitter.orbittype(planet)}(planet.priors, planet.derived, newplanet_obs..., name=planet.name)
    end

    # Generate new observations for the star
    newstar_obs = map(system.observations) do obs
        return genobs(obs, collect(elements), θ_newsystem)
    end

    # Generate new system
    newsystem = System(system.priors, system.derived, newstar_obs..., newplanets..., name=system.name)

    return newsystem
end
export generate
