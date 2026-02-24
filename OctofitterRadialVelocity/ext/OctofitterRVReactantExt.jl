"""
    OctofitterRVReactantExt

Minimal extension that activates Reactant compilation for
`MarginalizedStarAbsoluteRVObs` when a Reactant device is set.

All compilation logic is in the generic `OctofitterReactantExt` (in Octofitter).
This extension only overrides `ad_backend` to return `ReactantBackend`.
"""
module OctofitterRVReactantExt

using OctofitterRadialVelocity: MarginalizedStarAbsoluteRVObs
using Octofitter
using Octofitter: ReactantBackend, SystemObservationContext, mjup2msol
using PlanetOrbits: orbitsolve_bulk, radvel
using Reactant

const _ReactantDevice = Union{Reactant.XLA.AbstractClient, Reactant.XLA.AbstractDevice}

# When device is a Reactant client (e.g., Reactant.CPU()) or device, use ReactantBackend
function Octofitter.ad_backend(obs::MarginalizedStarAbsoluteRVObs{<:Any, <:Any, <:Reactant.XLA.AbstractClient})
    ReactantBackend(obs.device)
end
function Octofitter.ad_backend(obs::MarginalizedStarAbsoluteRVObs{<:Any, <:Any, <:Reactant.XLA.AbstractDevice})
    ReactantBackend(obs.device)
end

"""
Absolute radial velocity likelihood — broadcasting fallback for Reactant/XLA.

Uses `orbitsolve_bulk` for vectorized orbit solving. All array operations use
broadcasting for Reactant/XLA traceability. This method is selected when the
observation's device is a Reactant XLA client or device.
"""
function Octofitter.ln_like(
    rvlike::MarginalizedStarAbsoluteRVObs{<:Any, <:Any, <:_ReactantDevice},
    ctx::SystemObservationContext
)
    (; θ_system, θ_obs, orbits) = ctx
    epochs = rvlike.table.epoch
    jitter = θ_obs.jitter

    # Bulk solve + RV for each planet (vectorized)
    # Initialize rv_model from the first planet to establish the correct element type
    # (TracedRNumber during Reactant tracing, Float64 for Enzyme)
    rv_model = let
        sol = orbitsolve_bulk(orbits[1], epochs)
        planet_mass = θ_system.planets[1].mass
        radvel(sol, planet_mass * mjup2msol)
    end
    for planet_i in 2:length(orbits)
        sol = orbitsolve_bulk(orbits[planet_i], epochs)
        planet_mass = θ_system.planets[planet_i].mass
        rv_model = rv_model .+ radvel(sol, planet_mass * mjup2msol)
    end

    # Residuals (vectorized)
    resid = rvlike.table.rv .- rv_model

    # Variance per observation (vectorized)
    var = rvlike.table.σ_rv .^ 2 .+ jitter^2

    # Marginalize out the instrument zero point (Orvara paper)
    A = sum(1.0 ./ var)
    B = -2.0 * sum(resid ./ var)
    C = sum(resid .^ 2 ./ var)
    ll = -sum(log.(2π .* var))
    ll -= -B^2 / (4A) + C + log(A)

    return ll
end

end # module OctofitterRVReactantExt
