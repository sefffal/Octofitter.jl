struct PlanetOrderPrior{T} <: Octofitter.AbstractLikelihood
    planets::T
    function PlanetOrderPrior(planets::Planet...) 
        tup = tuple(planets...)
        return new{typeof(tup)}(tup)
    end
end
export PlanetOrderPrior
Octofitter._isprior(::PlanetOrderPrior) = true
likelihoodname(likeobj::PlanetOrderPrior) = "PlanetOrderPrior_"*join(string.(p.name for p in likeobj.planets),"_")
function Octofitter.ln_like(prior::PlanetOrderPrior, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)
  
    if length(orbits) == 1
        return ll
    end

    for (planet_a, planet_b) in zip(prior.planets[1:end-1], prior.planets[2:end])
        # Using name, find index of planet in θ_system. 
        # Using index, find orbit.
        # Ensure orbits are ordered correctly.
        i_a = findfirst(==(planet_a.name), keys(θ_system.planets))
        i_b = findfirst(==(planet_b.name), keys(θ_system.planets))
        orbit_a = orbits[i_a]
        orbit_b = orbits[i_b]
        if semimajoraxis(orbit_a) >= semimajoraxis(orbit_b)
            ll += convert(T, -Inf)
        end
    end

    return ll
end