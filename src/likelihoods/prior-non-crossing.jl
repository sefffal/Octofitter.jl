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

    n = length(orbits)
    n <= 1 && return ll

    M_star_total = θ_system.M

    # Check all planet pairs for Hill stability (no allocation, Enzyme-compatible)
    for i in 1:n
        for j in (i+1):n
            sma_i = semimajoraxis(orbits[i])
            sma_j = semimajoraxis(orbits[j])

            mass_i = θ_system.planets[i].mass
            mass_j = θ_system.planets[j].mass

            # Determine outer planet's semi-major axis for Hill radius
            sma_outer = max(sma_i, sma_j)
            delta_a = abs(sma_j - sma_i)

            M_star = max(zero(T), M_star_total - mass_i * Octofitter.mjup2msol - mass_j * Octofitter.mjup2msol)
            R_H = sma_outer * ((mass_i + mass_j) * Octofitter.mjup2msol / (3M_star))^(1//3)

            if delta_a <= 2*sqrt(3)*R_H
                return convert(T, -Inf)
            end
        end
    end

    return ll::T
end