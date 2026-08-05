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
    # After likelihoods.jl: it reuses that file's reference system.
    @testset "Shared observation front-ends" begin
        include("v2/sky-offset.jl")
    end
    # Also after likelihoods.jl: `reference_system`/`model_system` come from there.
    @testset "Radial velocity" begin
        include("v2/radial-velocity.jl")
    end
    @testset "Photometry and configuration priors" begin
        include("v2/photometry-and-priors.jl")
    end
    @testset "Sky-path helpers" begin
        include("v2/skypath.jl")
    end
    @testset "G23H joint Gaia/Hipparcos astrometry" begin
        include("v2/g23h.jl")
    end
    @testset "Hipparcos IAD" begin
        include("v2/hipparcos.jl")
    end
    @testset "Hierarchies and propagators" begin
        include("v2/hierarchy.jl")
    end
    @testset "IO, interop, flux tables" begin
        include("v2/io.jl")
    end
    # Self-gates its own sampling testsets on OCTOFITTER_TEST_MODE, so it runs
    # in full under "all" and unit-only under "unit".
    @testset "Analysis machinery" begin
        include("v2/analysis-machinery.jl")
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
    # Error-message contracts and silent-wrongness fixes from the
    # consolidation pass. Self-contained; order-independent.
    @testset "Consolidation regressions" begin
        include("v2/consolidation.jl")
    end
    # Loads Pigeons, so `octofit_pigeons`'s package extension is exercised.
    @testset "Parallel tempering (Pigeons extension)" begin
        include("v2/pigeons.jl")
    end
    # One worked model per observation type and public entry point. The
    # subpackage sections skip here — `Pkg.test()` on core resolves only
    # core's deps — and run in the shared `testenv/` sandbox.
    @testset "API smoke" begin
        include("v2/api-smoke.jl")
    end
end

if TEST_MODE in ("all", "integration")
    @testset "Sampling" begin
        include("v2/sampling.jl")
    end
end

@info "Tests completed!"
