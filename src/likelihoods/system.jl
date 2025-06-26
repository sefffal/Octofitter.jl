using Bumper




function make_ln_like(system::System, θ_system)

    # We want to gather up all observation epochs in a standardized order.
    # That way we can loop through and solve each orbit at each epoch.

    # We assume that all planets will need to be solved at every general system observation epoch
    # We assume that each planet will need to be solved at every epoch it has attached observations
    # We'll solve them once in a hot (possibly multi-threaded) loop, and pass views into each likelihood
    # function.

    # I want a vector of solution types per object
    # The solutions will be ordered first, each epoch in the system observations,
    # then, each epoch in that planet's observations.
    # Note that in practice MA is always fixed when solving; it's only e that we might want the gradient for.
    all_epochs = Float64[]
    epoch_start_index_mapping = Dict{Any,Int}()
    j = 1
    for obs in system.observations
        if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
            # TODO: deal with HGCA
            epoch_start_index_mapping[obs] = j
            j += length(obs.table.epoch)
            append!(all_epochs, obs.table.epoch)
        end
    end
    for i in 1:length(system.planets)
        for like in system.planets[i].observations
            if hasproperty(like, :table) && hasproperty(like.table, :epoch)
                epoch_start_index_mapping[like] = j
                j += length(like.table.epoch)
                append!(all_epochs, like.table.epoch)
            end
        end
    end

    planet_sol_keys = Symbol[]
    for i in eachindex(system.planets)
        sols_key = Symbol("sols_planet_$i")
        push!(planet_sol_keys, sols_key)
    end
    # TODO: this seems way overcomplicated? Just need a list of symbols
    # interpolated in.
    solutions_list = :(tuple($((:($sols_key) for sols_key in planet_sol_keys)...)))
    
    planet_keys = Symbol[]
    planet_construction_exprs = Expr[]
    planet_like_exprs = Expr[]
    planet_orbit_solution_exprs = Expr[]
    # Add planet declarations here
    planet_declarations = Expr[]
    j = 0
    for i in 1:length(system.planets)
        planet = system.planets[i]
        OrbitType = _planet_orbit_type(planet)
        key = Symbol("planet_$i")
        sols_key = planet_sol_keys[i]
        
        # Add declaration for this planet
        push!(planet_declarations, :($key = nothing))
        
        likelihood_exprs = map(enumerate(planet.observations)) do (i_like, like)
            i_epoch_start = get(epoch_start_index_mapping, like, 0)
            # Get the normalized observation name to access θ_obs
            if hasproperty(like, :name)
                obs_name = normalizename(likelihoodname(like))
                expr = :(
                    $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + ln_like(
                        system.planets[$(Meta.quot(i))].observations[$i_like],
                        θ_system,
                        θ_system.planets[$i],
                        hasproperty(θ_system.planets[$i].observations, $(Meta.quot(obs_name))) ? 
                            θ_system.planets[$i].observations.$obs_name :
                            (;),
                        elems,
                        ($solutions_list), # all orbit solutions
                        $i, # This planet index into orbit solutions
                        $(i_epoch_start-1) # start epoch index
                    );
                )
            else
                expr = :(
                    $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + ln_like(
                        system.planets[$(Meta.quot(i))].observations[$i_like],
                        θ_system,
                        θ_system.planets[$i],
                        (;),  # θ_obs
                        elems,
                        ($solutions_list), # all orbit solutions
                        $i, # This planet index into orbit solutions
                        $(i_epoch_start-1) # start epoch index
                    );
                )
            end
            j+=1
            return expr
        end

        likelihood_expr = quote
            $(likelihood_exprs...)
        end

        planet_contruction = quote
            $key = $(OrbitType)(;merge(θ_system, θ_system.planets[$i])...)
        end
    
        if isempty(all_epochs)  
            orbit_sol_expr = quote
                $sols_key = ()
            end
        else
            orbit_sol_expr = quote
                # Pre-solve kepler's equation for all epochs
                # epochs = Vector{Float64}(undef, $(length(epochs_planet_i)))
                epochs = @alloc(Float64, $(length(all_epochs)))
                $((
                    :(epochs[$j] = $(all_epochs[j]))
                    for j in 1:length(all_epochs)
                )...)

                sol0 = orbitsolve($key, first(epochs))
                # $sols_key = Vector{typeof(sol0)}(undef, length(epochs))
                $sols_key = @alloc(typeof(sol0), length(epochs))
                $sols_key[begin] = sol0
                $_kepsolve_all!(view($sols_key, 2:length(epochs)), $key, view(epochs, 2:length(epochs)))
            end
        end
        push!(planet_keys, key)
        push!(planet_construction_exprs, planet_contruction)
        push!(planet_like_exprs, likelihood_expr)
        push!(planet_orbit_solution_exprs, orbit_sol_expr)
    end


    
     sys_exprs = map(eachindex(system.observations)) do i
        like = system.observations[i]
        i_epoch_start = get(epoch_start_index_mapping, like, 0)
        # Get the normalized observation name to access θ_obs
        obs_name = normalizename(likelihoodname(like))
        expr = :(
            $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + ln_like(
                system.observations[$i], 
                θ_system, 
                hasproperty(θ_system.observations, $(Meta.quot(obs_name))) ? 
                    θ_system.observations.$obs_name :
                    (;),
                elems,
                ($solutions_list),
                $(i_epoch_start-1)
            );
            # if !isfinite($(Symbol("ll$(j+1)")))
            #     println("invalid likelihood value encountered")
            # end
        )
        j+=1
        return expr
    end

    return @RuntimeGeneratedFunction(:(function (system::System, θ_system)
        T = _system_number_type(θ_system)
        ll0 = zero(T)
        
        # Declare all planet variables before the try block
        $(planet_declarations...)

        # Try-catch only for planet construction
        try
            # Construct all orbit elements
            $(planet_construction_exprs...)
        catch err
            # Return -Inf if planet construction fails
            return convert(T, -Inf)
        end
        
        ll_out = @no_escape begin

            # Construct a tuple of existing planet orbital elements
            elems = tuple($(planet_keys...))

            # Solve all orbits
            $(planet_orbit_solution_exprs...)

            # evaluate all their individual observation likelihoods
            $(planet_like_exprs...)
            
            # And evaluate the overall system likelihoods
            $(sys_exprs...)

            $(Symbol("ll$j"))
        end

        return ll_out
    end))
end

const _kepsolve_use_threads = Ref(false)
function _kepsolve_all!(solutions, orbit, epochs)
    if _kepsolve_use_threads[]
        return _kepsolve_all_multithread!(solutions, orbit, epochs)
    else
        return _kepsolve_all_singlethread!(solutions, orbit, epochs)
    end
end
function _kepsolve_all_singlethread!(solutions, orbit, epochs)
    for epoch_i in eachindex(epochs)
        solutions[epoch_i] = orbitsolve(orbit, epochs[epoch_i])
    end
    return solutions
end

function _kepsolve_all_multithread!(solutions, orbit, epochs)
    Threads.@threads for epoch_i in 1:length(epochs)
        solutions[epoch_i] = orbitsolve(orbit, epochs[epoch_i])
    end
    return solutions
end


# Generate calibration data
"""
    generate_from_params(system, θ=drawfrompriors(system))

Generate a new system and observations from an existing model.
"""
function generate_from_params(system::System, θ_newsystem = drawfrompriors(system))

    # Generate new orbits for each planet in the system
    orbits = map(eachindex(system.planets)) do i
        planet = system.planets[i]
        neworbit = Octofitter.construct_elements(θ_newsystem, planet.name)
        return neworbit
    end

    # Generate new observations for each planet in the system
    newplanets = map(1:length(system.planets)) do i
        planet = system.planets[i]
        orbit = orbits[i]
        θ_newplanet = θ_newsystem.planets[i]
        newplanet_obs = map(planet.observations) do obs
            return generate_from_params(obs, θ_newplanet, orbit)
        end
        newplanet = Planet(
            variables=(planet.priors, planet.derived),
            basis=Octofitter.orbittype(planet),
            likelihoods=newplanet_obs,
            name=planet.name
        )
        return newplanet
    end

    # Generate new observations for the star
    newstar_obs = map(system.observations) do obs
        return generate_from_params(obs, θ_newsystem, collect(orbits))
    end

    # Generate new system
    newsystem = System(
        variables=(system.priors, system.derived),
        likelihoods=newstar_obs,
        companions=newplanets,
        name=system.name
    )

    return newsystem
end
export generate_from_params

