#=
This file contains functions for creating prior- and posterior-predictive distributions.
That is, given uncertainty in model parameters θ⃗, generate a distribution of data values
y⃗.


How do we want to represent uncertainty in the observations? The same was as when we generate
a system! Just concatenate them together such that the epochs match.
=#
using Random

"""
    obstable1, obstable2, = Octofitter.prior_predictive_distribution(model.system, 1000)

This function creates a prior-predictive distribution from a system.
That is, it draws parameter values from the model priors and uses 
those to generate new observation tables.

"""
function prior_predictive_distribution(system::System, N::Integer; rng=Random.default_rng())
    # Sample parameters from priors
    θ_newsystem_N = sample_priors(rng, system, N)
    arr2nt = make_arr2nt(system)
    return predictive_distribution(system, arr2nt.(θ_newsystem_N))
end

function posterior_predictive_distribution(system::System, chains)
    # use parameters already sampled from posterior
    # TODO: need to reconstruct from chains... can use planet keys from system yuck
end

function predictive_distribution(system::System, θ_system_N::AbstractArray)

    # Use these to generate new systems complete with observations
    new_systems = map(θ_system_N) do θ_system
        generate_from_params(system, θ_system)
    end

    # Concatenate the results together
    observations = AbstractLikelihood[]

    for obs_i in eachindex(first(new_systems).observations)
        concat_obs = vcat(filter(!isnothing, map(sys->Table(sys.observations[obs_i]), new_systems)))
        push!(observations, concat_obs)
    end
    for planet_i in eachindex(first(new_systems).planets)
        planet = system.planets[planet_i]
        for (obs_i, obs) in enumerate(planet.observations)
            concat_obs_tables = vcat(filter(!isnothing, map(sys->Table(sys.planets[planet_i].observations[obs_i]), new_systems)))
            if length(concat_obs_tables) > 0
                # THis is hacky since it relies on undocumented Base internals to get the type name.
                # This should be possible via overloading similar?
                concatenated_observations = Base.typename(typeof(obs)).wrapper(vcat(concat_obs_tables...))
                push!(observations, concatenated_observations)
            end
    
        end
    end

    return observations


end


