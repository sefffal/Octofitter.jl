using Test
using Octofitter
using Distributions
using TypedTables
using CSV
using HDF5

# Check for test mode environment variable
# Run with: OCTOFITTER_TEST_MODE=unit julia --project=. -e 'using Pkg; Pkg.test()'
# Run with: OCTOFITTER_TEST_MODE=integration julia --project=. -e 'using Pkg; Pkg.test()'
# Default: run all tests
const TEST_MODE = get(ENV, "OCTOFITTER_TEST_MODE", "all")

@info "Running Octofitter tests" mode=TEST_MODE

# =============================================================================
# Unit Tests (Fast - no MCMC sampling required)
# =============================================================================
if TEST_MODE in ("all", "unit")
    @testset "Unit Tests" begin
        @info "Running unit tests..."

        @testset "Constructors" begin
            include("unit/constructors.jl")
        end

        @testset "Likelihoods" begin
            include("unit/likelihoods.jl")
        end

        @testset "Distributions" begin
            include("unit/distributions.jl")
        end

        @testset "Priors" begin
            include("unit/priors.jl")
        end

        @testset "I/O" begin
            include("unit/io.jl")
        end
    end
end

# =============================================================================
# Integration Tests (Slower - requires MCMC sampling)
# =============================================================================
if TEST_MODE in ("all", "integration")
    @testset "Integration Tests" begin
        @info "Running integration tests..."

        @testset "Sampling" begin
            include("integration/sampling.jl")
        end

        @testset "Multi-Planet" begin
            include("integration/multi_planet.jl")
        end

        @testset "Joint Fitting" begin
            include("integration/joint_fitting.jl")
        end

        @testset "Cross Validation" begin
            include("integration/cross_validation.jl")
        end

        @testset "Plotting" begin
            include("integration/plotting.jl")
        end
    end
end

@info "Tests completed!"
