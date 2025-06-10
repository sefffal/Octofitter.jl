struct LimitClosestApproachAUPrior <: Octofitter.AbstractLikelihood
    hard_closest_approach_au::Float64
    soft_closest_approach_au::Float64
end
LimitClosestApproachAUPrior(soft_closest_approach_au) = LimitClosestApproachAUPrior(0, soft_closest_approach_au)
NonCrossingPrior() = LimitClosestApproachAUPrior(0.0, 0.0)
export NonCrossingPrior, LimitClosestApproachAUPrior
Octofitter._isprior(::LimitClosestApproachAUPrior) = true
function Octofitter.ln_like(prior::LimitClosestApproachAUPrior, θ_system,θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)
  
    if length(orbits) == 1
        return ll
    end

    # TODO: would be nice to find a way to make this non-allocating
    orbits_sorted = sort(collect(orbits), by=orbit->semimajoraxis(orbit))
    for (orbit_a, orbit_b) in zip(orbits_sorted[1:end-1], orbits_sorted[2:end])
        # orbit_b has sma > orbit_a by construction
        sep_farthest_inner = apoapsis(orbit_a)
        sep_closest_outer = periapsis(orbit_b)
        closest_approach = sep_closest_outer - sep_farthest_inner
        # hard cutoff
        if closest_approach <= prior.hard_closest_approach_au
            ll = convert(T, -Inf)
            break
        end
        # Apply soft constraint using repulsive potential if soft_closest_approach_au is > than hard_closest_approach_au
        if closest_approach < prior.soft_closest_approach_au
            ll -= 1 / (closest_approach - prior.soft_closest_approach_au)^2
        end
    end

    return ll::T
end