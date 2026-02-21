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
using Octofitter: ReactantBackend
using Reactant

# When device is a Reactant client (e.g., Reactant.CPU()) or device, use ReactantBackend
function Octofitter.ad_backend(obs::MarginalizedStarAbsoluteRVObs{<:Any, <:Any, <:Reactant.XLA.AbstractClient})
    ReactantBackend(obs.device)
end
function Octofitter.ad_backend(obs::MarginalizedStarAbsoluteRVObs{<:Any, <:Any, <:Reactant.XLA.AbstractDevice})
    ReactantBackend(obs.device)
end

end # module OctofitterRVReactantExt
