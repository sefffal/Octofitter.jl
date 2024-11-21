

"""
    pointwise_like(model, chain)

Given a model and posterior sample chain split out the posterior probability into prior and columns for each datapoint.
"""
function pointwise_like(model, chain)

    # Resolve the array back into the nested named tuple structure used internally.
    sample_nts = mcmcchain2result(model, chain,)

    # Got this part done!
    # Now need to construct a list of System objects from model, where each one only
    # has one datapoint.
    # Then loop through sample_nts * System_objects' and compute ln_like.

    per_obs_systems = generate_system_per_epoch(model.system)

    # Convert these system definitions into ln_like functions
    ln_likes = map(per_obs_systems) do system
        Octofitter.make_ln_like(system, first(sample_nts))
    end

    return broadcast(sample_nts, reshape(per_obs_systems,1,:), reshape(ln_likes,1,:)) do sample_nt, system, ln_like
        loglike = ln_like(system, sample_nt)
        return loglike
    end

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
        newplanet_obs = exclude_all ? AbstractLikelihood[] : filter(_isprior, planet.observations)
        newplanet = Planet{Octofitter.orbittype(planet)}(planet.priors, planet.derived, newplanet_obs..., name=planet.name)
        return newplanet
    end

    newsys_obs = exclude_all ? AbstractLikelihood[] : filter(_isprior, system.observations)

    # Generate new system with observations
    newsystem = System(system.priors, system.derived, newsys_obs..., newplanets..., name=system.name)

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
Given an existing model with N epochs of data spread across any number of likelihood objects, return N copies
that only contain data from a given epoch.
"""
function generate_system_per_epoch(system)
    
    # We will number observations by index in their table (if there is a table)
    # (if no table, just discard observation object). Start with observations of star,
    # then observations attached to each planet.
    observation_count = _count_epochs(system)
    
    per_obs_systems = map(1:observation_count) do i_obs
        j_obs = 0
        for obs in system.observations
            if hasproperty(obs, :table)
                for k_obs in 1:Tables.rowcount(obs.table)
                    j_obs += 1
                    if i_obs == j_obs
                        obs_row = likeobj_from_epoch_subset(obs, k_obs)
                        return System(
                            system.priors,
                            system.derived,
                            obs_row;
                            system.name
                        )
                    end
                end
            else
                j_obs += 1
                if i_obs == j_obs
                    return System(
                        system.priors,
                        system.derived,
                        obs;
                        system.name
                    )
                end
            end
        end
        for planet in system.planets
            for obs in planet.observations
                if hasproperty(obs, :table)
                    for k_obs in 1:Tables.rowcount(obs.table)
                        j_obs += 1
                        if i_obs == j_obs
                            obs_row = likeobj_from_epoch_subset(obs, k_obs)
                            planet = Planet{orbittype(planet)}(
                                planet.priors,
                                planet.derived,
                                (obs_row,);
                                planet.name
                            )
                            return System(
                                system.priors,
                                system.derived,
                                planet;
                                system.name
                            )
                        end
                    end
                else
                    j_obs += 1
                    if i_obs == j_obs
                        planet = Planet{orbittype(planet)}(
                            planet.priors,
                            planet.derived,
                            (obs,);
                            planet.name
                        )
                        return System(
                            system.priors,
                            system.derived,
                            planet;
                            system.name
                        )
                    end
                end
            end
        end
    end

    return per_obs_systems
end