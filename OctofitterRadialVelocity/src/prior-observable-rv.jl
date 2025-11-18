
function Octofitter.ln_like(
    like::Union{Octofitter.ObsPriorAstromONeil2019{StarAbsoluteRVObs}, Octofitter.ObsPriorAstromONeil2019{MarginalizedStarAbsoluteRVObs}},
    context::Octofitter.PlanetObservationContext
)
    (; θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start) = context
    T = Octofitter._system_number_type(θ_system)

    # Call the wrapped likelihood's ln_like method
    ln_like_wrapped = Octofitter.ln_like(like.wrapped_like, context)

    orbit = orbits[i_planet]
    # Add period prior
    ln_prior = zero(T)
    P = period(orbit)/365.25
    e = eccentricity(orbit)
    jac = zero(T)
    # We can't use i_orbsol_start for this, since we're actually using the epochs
    # of another likelihood object. Instead loop through each until we find the right epoch.
    # This should normally be faster than solving it again.
    # TODO: I wonder if a Dict would be a better choice here. Something to benchmark.
    # Or maybe a dict of likelihood object to starting index?
    # TODO: problem: this doesn't work with AbsoluteVisual orbits, since the true time is different
    # due to changing light travel time.
    for j in 1:length(like.wrapped_like.table.epoch)
        current_epoch = like.wrapped_like.table.epoch[j]
        local sol
        found = false
        for k in eachindex(orbit_solutions[i_planet])
            sol′ = orbit_solutions[i_planet][k]
            if hasproperty(sol′, :t) && sol′.t == current_epoch
                found = true
                sol = sol′
                break
            elseif  hasproperty(sol′, :sol) && hasproperty(sol′.sol, :t) &&
                    sol′.sol.t == current_epoch
                found = true
                sol = sol′
                break
            end
        end
        # if !found
        #     error("epoch not found in solutions $current_epoch")
        # end
        if !found
            sol = orbitsolve(orbit, current_epoch)
        end
        M = meananom(sol)
        E = eccanom(sol)
        jac += abs(3M*(
            e+cos(E)
        ) + 2*(-2+e^2+e*cos(E)) *sin(E))
    end

    sqrt_eccen = sqrt(1-eccentricity(orbit)^2)
    jac *= cbrt(P) / sqrt_eccen

    ln_prior += 2log(jac)

    return ln_like_wrapped + ln_prior
end