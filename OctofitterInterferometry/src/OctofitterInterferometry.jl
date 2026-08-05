module OctofitterInterferometry

using Octofitter
using PlanetOrbits
using Tables, TypedTables

using FITSIO

using LinearAlgebra
using Interpolations
using BlockArrays
using Distributions
using PDMats
using Bumper

"""
    AbstractInterferometryObs <: Octofitter.AbstractObs

Supertype for the interferometric observation types in this package. There is
currently one, [`InterferometryObs`](@ref); the v1 `GRAVITYWideKPObs` is now a
preset of it (see [`GRAVITYWideKPObs`](@ref)).
"""
abstract type AbstractInterferometryObs <: Octofitter.AbstractObs end

# v1 spelled every observation type `…Likelihood`. The aliases are kept so
# that old scripts at least fail on the *arguments* that changed rather than
# on the name.
const AbstractInterferometryLikelihood = AbstractInterferometryObs
export AbstractInterferometryObs, AbstractInterferometryLikelihood

include("oifits.jl")
include("visibility.jl")
include("GRAVITY-correlation.jl")
include("kernel-phases.jl")
include("interferometry.jl")

end
