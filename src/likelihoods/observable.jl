#=
This file implements the observable-based priors of
O'Neil 2019
"Improving Orbit Estimates for Incomplete Orbits with a New Approach 
to Priors: with Applications from Black Holes to Planets"

Currently, only the prior for relative astrometry and not radial velocity 
is implemented.
=#


"""
    ObsPriorAstromONeil2019(astrometry_likelihood, period_prior)

Given a an astrometry likelihood (`PlanetRelAstromLikelihood`), apply the
"observable based priors" of K. O'Neil 2019 "Improving Orbit Estimates
for Incomplete Orbits with a New Approach to Priors: with Applications
from Black Holes to Planets".

This prior correction is only correct if you supply Uniform priors on 
all Campbell orbital parameters and a Uniform prior on Period (not
semi-major axis).
This period prior has a significant impact in the
fit and recommendations for its range were not published in the original paper.

## Examples
```julia
astrom_like = PlanetRelAstromLikelihood(astrom_table)

# Apply observable based priors ontop of our uniform Campbell priors:
obs_prior = ObsPriorAstromONeil2019(astrom_like)

# The astrometry lieklihood object is passed as a first parameter
# since the obserable-based priors depend on the observation 
# epochs.

@planet b Visual{KepOrbit} begin
    # Instead of a prior on sma
    # a ~ Uniform(0.001, 10000)

    # Put a prior on period:
	P ~ Uniform(0.001, 2000) # yrs
    a = cbrt(system.M * b.P^2)

    # Keep sine prior on inclination
    i ~ Sine()

    # Rest are uniform
    e ~ Uniform(0.0, 1.0)
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    τ ~ UniformCircular(1.0)
end astrom_like obs_prior
```
"""
struct ObsPriorAstromONeil2019{Likelihood<:AbstractLikelihood} <: AbstractLikelihood
	wrapped_like::Likelihood
	function ObsPriorAstromONeil2019(
        obs::AbstractLikelihood;
    )
		return new{typeof(obs)}(obs)
	end
end
export ObsPriorAstromONeil2019

function Octofitter.ln_like(like::ObsPriorAstromONeil2019{<:PlanetRelAstromLikelihood}, θ_planet, orbit, orbit_solutions, i_orbsol_start) 
# function Octofitter.ln_like(like::ObsPriorAstromONeil2019, θ_planet, orbit,) 

    # Add period prior
    ln_prior = 0.0
    P = period(orbit)/365.25
    e = eccentricity(orbit)
    jac = 0.0
    # We can't use i_orbsol_start for this, since we're actually using the epochs
    # of another likelihood object. Instead loop through each until we find the right epoch.
    # This should normally be faster than solving it again.
    # TODO: I wonder if a Dict would be a better choice here. Something to benchmark.
    # Or maybe a dict of likelihood object to starting index?
    for j in 1:length(like.wrapped_like.table.epoch)
        current_epoch = like.wrapped_like.table.epoch[j]
        local sol
        for k in eachindex(orbit_solutions)
            sol′ = orbit_solutions[k]
            if sol′.sol.t == current_epoch
                sol = sol′
                break
            end
        end
        M = meananom(sol)
        E = eccanom(sol)
        jac += abs(3M*(
            e+cos(E)
        ) + 2*(-2+e^2+e*cos(E)) *sin(E))
    end

    jac *= cbrt(P) / sqrt(1-eccentricity(orbit)^2);

    ln_prior += - 2log(jac)

    return ln_prior
end

# TODO: Add a RadialVelocity correction version in OctofitterRadialVelocity