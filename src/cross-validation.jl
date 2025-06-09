
function run_all_functions(func_tuple::Tuple, systems_tuple::Tuple, params)
    return ntuple(i -> func_tuple[i](systems_tuple[i], params), length(func_tuple))
end

function pointwise_like(model, chain)
    # Resolve the array back into the nested named tuple structure used internally
    @info "Resolving chain"
    sample_nts = mcmcchain2result(model, chain)

    # Create system objects for each observation
    @info "Creating systems for each data point"
    per_obs_systems, epochs = generate_system_per_epoch(model.system)

    # Convert these system definitions into ln_like functions
    @info "Compiling likelihood functions"
    ln_likes = map(per_obs_systems) do system
        Octofitter.make_ln_like(system, first(sample_nts))
    end

    # Convert to tuples for better type stability and compile-time specialization
    ln_likes_tuple = Tuple(ln_likes)
    systems_tuple = Tuple(per_obs_systems)
    
    # Preallocate output array
    LL_out = zeros(length(sample_nts), length(per_obs_systems))
    
    # Disable threading in the solver -- we will thread in the outer loop
    Octofitter._kepsolve_use_threads[] = false

    @info "Calculating..."
    
    # Compute likelihoods for each sample
    Threads.@threads for i_sample in 1:size(LL_out, 1)
        @printf("Sample %5d ", i_sample)
        # Use the specialized function to get all likelihoods at once -- this compiles into faster code
        # than just looping through them here (triggers a bunch of dynamic dispatches)
        likelihoods = run_all_functions(ln_likes_tuple, systems_tuple, sample_nts[i_sample])
        
        # Store results in output array
        for i_obs in 1:size(LL_out, 2)
            print('.')
            LL_out[i_sample, i_obs] = likelihoods[i_obs]
        end
        println()
    end
    
    return LL_out, epochs
end



# Generate calibration data
"""
    prior_only_model(system, θ=drawfrompriors(system))

Creates a copy of a system model that is stripped of observations. The result is 
a model that only samples from the priors. This can be used eg. for tempering.
"""
function prior_only_model(system::System;exclude_all=false)


    # Generate new observations for each planet in the system
    newplanets = map(1:length(system.planets)) do i
        planet = system.planets[i]
        # newplanet_obs = exclude_all ? AbstractLikelihood[] : filter(_isprior, planet.observations)
        newplanet_obs = map(planet.observations) do obs
            if exclude_all || !_isprior(obs)
                # We have to make sure we have the same variables even if we drop the data
                return BlankLikelihood((obs.priors,obs.derived),obs.instrument_name)
            end
        end
        newplanet_obs = filter(!isnothing, newplanet_obs)
        newplanet = Planet(
            basis=Octofitter.orbittype(planet),
            variables=(planet.priors, planet.derived),
            likelihoods=newplanet_obs,
            name=planet.name
        )
        return newplanet
    end

    newsys_obs = map(system.observations) do obs
        if exclude_all || !_isprior(obs)
            # We have to make sure we have the same variables even if we drop the data
            return BlankLikelihood((obs.priors,obs.derived),obs.instrument_name)
        end
    end
    # Generate new system with observations
    newsystem = System(
        variables=(system.priors, system.derived),
        likelihoods=newsys_obs,
        companions=newplanets,
        name=system.name
    )

    return newsystem
end
export prior_only_model


"""
Given an existing model with N likelihood objects, return N copies
that each drop one of the likehood objects.
"""
function generate_kfold_systems(system)
    # We will number observations by index in their table (if there is a table)
    # (if no table, just discard observation object). Start with observations of star,
    # then observations attached to each planet.
    # Yikes!
    likeobj_count = _count_likeobj(system)
    
    per_obs_systems = map(1:likeobj_count) do i_obs
        j_obs = 0
        to_include_system = []
        to_include_planets = []
        for obs in system.observations
            if _isprior(obs)
                push!(to_include_system, obs)
            else
                j_obs += 1
                if i_obs != j_obs
                    push!(to_include_system, likeobj_from_epoch_subset(obs, :))
                end
            end
            # TODO, if prior always include
        end
        for planet in system.planets
            to_include_planet = []
            push!(to_include_planets, to_include_planet)
            for obs in planet.observations
                if _isprior(obs)
                    push!(to_include_planet, obs)
                else
                    j_obs += 1
                    if i_obs != j_obs
                        obs_row = likeobj_from_epoch_subset(obs, :)
                        push!(to_include_planet, likeobj_from_epoch_subset(obs, :))
                    end
                end
            end
        end

        planets_new = map(system.planets, to_include_planets) do planet, like_objs
            planet = Planet{orbittype(planet)}(
                planet.priors,
                planet.derived,
                like_objs...;
                planet.name
            )
        end

        new_name = Symbol(string(system.name, "_kfold_$i_obs"))
        return System(
            system.priors,
            system.derived,
            to_include_system...,
            planets_new...;
            name=new_name
        )
    end

    return per_obs_systems
end


"""
Given an existing model with N likelihood objects, run a user-provided callback
to decide which likelihood objects to retain
"""
function generate_system_filtered_like(user_predicate, system)
    # We will number observations by index in their table (if there is a table)
    # (if no table, just discard observation object). Start with observations of star,
    # then observations attached to each planet.
    # Yikes!
    likeobj_count = _count_likeobj(system)
    
    j_obs = 0
    to_include_system = []
    to_include_planets = []
    included_j_obs = Int[]
    for obs in system.observations
        # if _isprior(obs)
            # push!(to_include_system, obs)
        # else
            j_obs += 1
            if user_predicate(obs)
                push!(to_include_system, likeobj_from_epoch_subset(obs, :))
                push!(included_j_obs, j_obs)
            end
        # end
    end
    for planet in system.planets
        to_include_planet = []
        push!(to_include_planets, to_include_planet)
        for obs in planet.observations
            # if _isprior(obs)
                # push!(to_include_planet, obs)
            # else
                j_obs += 1
                if user_predicate(obs)
                    obs_row = likeobj_from_epoch_subset(obs, :)
                    push!(to_include_planet, likeobj_from_epoch_subset(obs, :))
                    push!(included_j_obs, j_obs)
                end
            # end
        end
    end


    planets_new = map(system.planets, to_include_planets) do planet, like_objs
        # println.(typeof.(like_objs))
        planet = Planet{orbittype(planet)}(
            planet.priors,
            planet.derived,
            like_objs...;
            planet.name
        )
    end

    new_name = Symbol(string(system.name, "_filt_obs_", join(included_j_obs, '_')))
    return System(
        system.priors,
        system.derived,
        to_include_system...,
        planets_new...;
        name=new_name
    )

end

"""
Given an existing model with N likelihood objects, return N copies
that each have only one likelihood attached
"""
function generate_systems_per_like(system)
    # We will number observations by index in their table (if there is a table)
    # (if no table, just discard observation object). Start with observations of star,
    # then observations attached to each planet.
    # Yikes!
    likeobj_count = _count_likeobj(system)
    
    per_obs_systems = map(1:likeobj_count) do i_obs
        j_obs = 0
        to_include_system = []
        to_include_planets = []
        for obs in system.observations
            if _isprior(obs)
                push!(to_include_system, obs)
            else
                j_obs += 1
                if i_obs == j_obs
                    push!(to_include_system, likeobj_from_epoch_subset(obs, :))
                end
            end
            # TODO, if prior always include
        end
        for planet in system.planets
            to_include_planet = []
            push!(to_include_planets, to_include_planet)
            for obs in planet.observations
                if _isprior(obs)
                    push!(to_include_planet, obs)
                else
                    j_obs += 1
                    if i_obs == j_obs
                        obs_row = likeobj_from_epoch_subset(obs, :)
                        push!(to_include_planet, likeobj_from_epoch_subset(obs, :))
                    end
                end
            end
        end

        planets_new = map(system.planets, to_include_planets) do planet, like_objs
            planet = Planet{orbittype(planet)}(
                planet.priors,
                planet.derived,
                like_objs...;
                planet.name
            )
        end

        new_name = Symbol(string(system.name, "_like_$i_obs"))
        return System(
            system.priors,
            system.derived,
            to_include_system...,
            planets_new...;
            name=new_name
        )
    end

    return per_obs_systems
end

"""
Given an existing model with N epochs of data spread across any number of likelihood objects, return N copies
that only contain data from a given epoch. Non-tabular/_isprior observations are included in all systems.
"""
function generate_system_per_epoch(system)
    
    # First, collect all tabular observations and count how many total epochs we have
    tabular_observations = []
    
    # Collect system-level tabular observations
    for obs in system.observations
        if hasproperty(obs, :table)
            push!(tabular_observations, (obs=obs, planet_idx=nothing, row_count=Tables.rowcount(obs.table)))
        end
    end
    
    # Collect planet-level tabular observations
    for (planet_idx, planet) in enumerate(system.planets)
        for obs in planet.observations
            if hasproperty(obs, :table)
                push!(tabular_observations, (obs=obs, planet_idx=planet_idx, row_count=Tables.rowcount(obs.table)))
            end
        end
    end
    
    # Calculate total epoch count from tabular observations
    total_epochs = sum(item -> item.row_count, tabular_observations)
    
    if total_epochs == 0
        return System[], Float64[]
    end
    
    epochs = Float64[]

    # Create a system for each epoch
    per_obs_systems = map(1:total_epochs) do i_obs
        # Initialize collections for observations
        to_include_system = []
        
        # Always include non-tabular system observations in every system
        for obs in system.observations
            if !hasproperty(obs, :table)
                push!(to_include_system, obs)
            end
        end
        
        # Track current epoch index
        current_epoch = 0
        
        # Add the matching tabular observation for this epoch
        for tab_obs in tabular_observations
            # Skip planet observations for now
            if tab_obs.planet_idx !== nothing
                continue
            end
            
            # For each row in the table
            for k_row in 1:tab_obs.row_count
                current_epoch += 1
                if current_epoch == i_obs
                    # This is the matching epoch, add just this row
                    o = likeobj_from_epoch_subset(tab_obs.obs, k_row)
                    push!(epochs, tab_obs.obs.table.epoch[k_row])
                    push!(to_include_system, o)
                    break
                end
            end
        end
        
        # Process all planets
        planets_new = map(1:length(system.planets)) do p_idx
            planet = system.planets[p_idx]
            to_include_planet = []
            
            # Always include non-tabular observations for this planet
            for obs in planet.observations
                if !hasproperty(obs, :table)
                    push!(to_include_planet, obs)
                end
            end
            
            # Add matching tabular observation for this planet if it exists
            for tab_obs in tabular_observations
                if tab_obs.planet_idx == p_idx
                    for k_row in 1:tab_obs.row_count
                        current_epoch += 1
                        if current_epoch == i_obs
                            # This is the matching epoch, add just this row
                            o = likeobj_from_epoch_subset(tab_obs.obs, k_row)
                            push!(epochs, tab_obs.obs.table.epoch[k_row])
                            push!(to_include_planet, o)
                            break
                        end
                    end
                end
            end
            
            # Create new planet with included observations
            return Planet{orbittype(planet)}(
                planet.priors,
                planet.derived,
                to_include_planet...;
                name=planet.name
            )
        end
        
        # Create the new system with all planets
        new_name = Symbol(string(system.name, "_epoch_$i_obs"))
        return System(
            system.priors,
            system.derived,
            to_include_system...,
            planets_new...;
            name=new_name
        )
    end

    return per_obs_systems, epochs
end