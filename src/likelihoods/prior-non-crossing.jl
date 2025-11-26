struct LimitClosestApproachAUPrior <: Octofitter.AbstractObs
    hard_closest_approach_au::Float64
    soft_closest_approach_au::Float64
end
LimitClosestApproachAUPrior(soft_closest_approach_au) = LimitClosestApproachAUPrior(0, soft_closest_approach_au)
NonCrossingPrior() = LimitClosestApproachAUPrior(0.0, 0.0)
likelihoodname(likeobj::LimitClosestApproachAUPrior) = "LimitClosestApproachAUPrior"
export NonCrossingPrior, LimitClosestApproachAUPrior
Octofitter._isprior(::LimitClosestApproachAUPrior) = true
function Octofitter.ln_like(prior::LimitClosestApproachAUPrior, ctx::SystemObservationContext)
    (; θ_system, orbits) = ctx
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




struct HillStabilityPrior <: Octofitter.AbstractObs
end
likelihoodname(likeobj::HillStabilityPrior) = "HillStabilityPrior"
export HillStabilityPrior
Octofitter._isprior(::HillStabilityPrior) = true
function Octofitter.ln_like(prior::HillStabilityPrior, ctx::SystemObservationContext)
    (; θ_system, orbits) = ctx
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)
  
    if length(orbits) == 1
        return ll
    end

    # TODO: would be nice to find a way to make this non-allocating
    orbits_indexed = sort(collect(pairs(orbits)), by=p->semimajoraxis(p[2]))

    for i in 1:(length(orbits_indexed)-1)
        idx_a, orbit_a = orbits_indexed[i]
        idx_b, orbit_b = orbits_indexed[i+1]
        
        θ_planet_a = θ_system.planets[idx_a]
        θ_planet_b = θ_system.planets[idx_b]

        θ_planet_a = θ_system.planets[idx_a]
        θ_planet_b = θ_system.planets[idx_a]

        # orbit_b has sma > orbit_a by construction
        # sep_farthest_inner = apoapsis(orbit_a)
        # sep_closest_outer = periapsis(orbit_b)
        # closest_approach = sep_closest_outer - sep_farthest_inner
        delta_a = semimajoraxis(orbit_b) - semimajoraxis(orbit_a)

        M_star = max(zero(typeof(θ_planet_b.M)), θ_planet_b.M - θ_planet_a.mass*Octofitter.mjup2msol - θ_planet_b.mass*Octofitter.mjup2msol)
        # calculate Hill radius
        R_H = semimajoraxis(orbit_b) * ((θ_planet_a.mass + θ_planet_b.mass)*Octofitter.mjup2msol/3M_star)^(1/3)

        # hard cutoff
        if delta_a <= 2*sqrt(3)*R_H
            ll = convert(T, -Inf)
            break
        end
    end

    return ll::T
end