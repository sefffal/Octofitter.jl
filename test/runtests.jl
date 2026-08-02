using Test
using Octofitter
using Distributions
using LinearAlgebra
using Random
using StaticArrays
using DelimitedFiles
using FiniteDiff
using PlanetOrbits

# `test/legacy/` holds the v1 test suite, which exercises the model surface
# this branch replaced. It is kept beside the v1 sources in `src/legacy/` so
# each likelihood's port can bring its tests with it.
const TEST_MODE = get(ENV, "OCTOFITTER_TEST_MODE", "all")

@info "Running Octofitter tests" mode = TEST_MODE

if TEST_MODE in ("all", "unit")
    @testset "Model surface" begin
        include("v2/model.jl")
    end
    @testset "Observations" begin
        include("v2/likelihoods.jl")
    end
    @testset "Hierarchies and propagators" begin
        include("v2/hierarchy.jl")
    end
    @testset "Plotting API" begin
        include("v2/plotting.jl")
    end
    @testset "Gradients and allocations" begin
        include("v2/gradients.jl")
    end
    @testset "v1 numerical regression" begin
        include("v2/v1-regression.jl")
    end
end

if TEST_MODE in ("all", "integration")
    @testset "Sampling" begin
        include("v2/sampling.jl")
    end
end

@info "Tests completed!"
