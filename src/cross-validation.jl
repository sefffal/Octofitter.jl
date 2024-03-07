

"""
Given a model and original samples? split out the posterior probability into prior and columns for each datapoint.
"""
function pointwise_like(model, chain)

    # Resolve the array back into the nested named tuple structure used internally.
    sample_nts = mcmcchain2result(model, chain,)

    # Got this part done!
    # Now need to construct a list of System objects from model, where each one only
    # has one datapoint.
    # Then loop through sample_nts * System_objects' and compute ln_like.

    # We will number observations by index in their table (if there is a table)
    # (if no table, just discard observation object). Start with observations of star,
    # then observations attached to each planet.
    # Yikes!

    observation_count = 0
    for obj in [model.system; model.system.planets...]
        for obs in obj.observations
            if hasproperty(obs, :table)
                observation_count += Tables.rowcount(obs.table)
            else
                observation_count += 1
            end
        end
    end
    
    per_obs_systems = map(1:observation_count) do i_obs
        j_obs = 0
        for obs in model.system.observations
            # TODO: this has not been tested
            if hasproperty(obs, :table)
                for k_obs in 1:Tables.rowcount(obs.table)
                    j_obs += 1
                    if i_obs == j_obs
                        # TODO WARNING! internals access
                        T = Main.eval(typeof(obs).name.name)
                        obs_row = T(obs.table[k_obs,:])
                        return System(
                            model.system.priors,
                            model.system.derived,
                            obs_row;
                            model.system.name
                        )
                    end
                end
            else
                j_obs += 1
                if i_obs == j_obs
                    return System(
                        model.system.priors,
                        model.system.derived,
                        obs;
                        model.system.name
                    )
                end
            end
        end
        # TODO: planets
        for planet in model.system.planets
            for obs in planet.observations
                if hasproperty(obs, :table)
                    for k_obs in 1:Tables.rowcount(obs.table)
                        j_obs += 1
                        if i_obs == j_obs
                            # TODO WARNING! internals access
                            T = Main.eval(typeof(obs).name.name)
                            obs_row = T(obs.table[k_obs,:])
                            planet = Planet{orbittype(planet)}(
                                planet.priors,
                                planet.derived,
                                (obs_row,);
                                planet.name
                            )
                            return System(
                                model.system.priors,
                                model.system.derived,
                                planet;
                                model.system.name
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
                            model.system.priors,
                            model.system.derived,
                            planet;
                            model.system.name
                        )
                    end
                end
            end
        end
    end


    # Convert these system definitions into ln_like functions
    ln_likes = map(per_obs_systems) do system
        Octofitter.make_ln_like(system, first(sample_nts))
    end


    return broadcast(sample_nts, reshape(per_obs_systems,1,:), reshape(ln_likes,1,:)) do sample_nt, system, ln_like
        loglike = ln_like(system, sample_nt)
        return loglike
    end

end