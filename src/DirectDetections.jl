module DirectDetections
using ComponentArrays
using Distributions

import KissMCMC
using AdvancedHMC
using NamedTupleTools
using ForwardDiff
using Logging

using Statistics
using StatsBase
using NamedTupleTools
using DirectImages
using DirectOrbits
using Base.Threads: @threads
using StaticArrays

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

using RecipesBase

const mjup2msol = 0.0009543

# We use a sonora model grid to tie fluxes to physical properties
include("sonora.jl")
include("types.jl")
include("models.jl")
include("sampling.jl")
include("analysis.jl")

function __init__()
    init_models()
    init_plots()
    return
end

end
