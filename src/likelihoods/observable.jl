#=
This file implements the observable-based priors of
O'Neil 2019
"Improving Orbit Estimates for Incomplete Orbits with a New Approach 
to Priors: with Applications from Black Holes to Planets"

Currently, only the prior for relative astrometry and not radial velocity 
is implemented.
=#


"""
    ObsPriorAstromONeil2019(astrometry_table)

Given a table of astrometry (in same format as `AstrometryLikelihood`,
but only `epoch` is required), apply the "observable based priors" 
of K. O'Neil 2019 "Improving Orbit Estimates for Incomplete Orbits with  
a New Approach to Priors: with Applications from Black Holes to Planets".

This prior correction is only correct if you supply Uniform priors on 
all orbital parameters.

A Period parameter (in years) must exist in the planet model as `P`,
in addition to the usual semi-major axis `a` (au). Typically
this can be provided as:
```julia
@planet b Visual{KepOrbit} begin
    P ~ Uniform(0.001, 1000)
    a = cbrt(system.M * b.P^2)
    ...
end AstrometryLikelihood(astrom_table) ObsPriorAstromONeil2019(astrom_table)
```
"""
struct ObsPriorAstromONeil2019{TTable<:Table} <: Octofitter.AbstractLikelihood
	table::TTable
	function ObsPriorAstromONeil2019(observations...)
		table = Table(observations...)
		return new{typeof(table)}(table)
	end
end
export ObsPriorAstromONeil2019

function Octofitter.ln_like(like::ObsPriorAstromONeil2019, θ_planet, orbit,)

    jac = 0.0
    for j in 1:length(like.table.epoch)
        current_epoch = like.table.epoch[j]
        M = meananom(orbit, current_epoch)
        E = eccanom(orbit, current_epoch)
        jac += abs(3M*(
            eccentricity(orbit)+cos(E)
        ) + 2*(-2+eccentricity(orbit)^2+eccentricity(orbit)*cos(E)) *sin(E))
    end

    P = period(orbit)/365.25
    jac *= P^(1/3) / sqrt(1-eccentricity(orbit)^2);

    return -2log(jac);
end