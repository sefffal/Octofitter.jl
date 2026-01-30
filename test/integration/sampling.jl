using FiniteDiff
using Random

# Reduced iterations for faster tests
const TEST_ITERATIONS = 50
const TEST_ADAPTATION = 50

@testset "Basic MCMC Sampling" begin
    rng = Random.Xoshiro(1)

    astrom_like = PlanetRelAstromLikelihood(
        Table(
            epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
            ra = [-505.76, -502.57, -498.21, -492.68, -485.98, -478.11, -469.08, -458.90],
            dec = [-66.93, -37.47, -7.93, 21.64, 51.15, 80.54, 109.73, 138.65],
            σ_ra = fill(10.0, 8),
            σ_dec = fill(10.0, 8),
            cor = fill(0.0, 8)
        ),
        name = "sampling_test"
    )

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom_like],
        variables=@variables begin
            a ~ Uniform(0, 100)
            e ~ Uniform(0.0, 0.99)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    TestSys = System(
        name="TestSys",
        companions=[b],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        end
    )

    # Create a single model for all tests to avoid redundant compilation
    model = Octofitter.LogDensityModel(TestSys)

    @testset "Model Construction" begin
        @test model.D == 11
    end

    @testset "Initialization" begin
        Octofitter.default_initializer!(rng, model)
        @test length(model.starting_points) > 0
        @test model.ℓπcallback(model.starting_points[1]) > -1000
    end

    @testset "Sampling" begin
        chain = octofit(rng, model, iterations=TEST_ITERATIONS, adaptation=TEST_ADAPTATION)

        @test all(chain[:logpost] .> -1000)
        # Allow some numerical errors but not too many
        @test mean(chain[:numerical_error]) < 0.15
    end

    @testset "Chain I/O" begin
        chain = octofit(rng, model, iterations=TEST_ITERATIONS, adaptation=TEST_ADAPTATION)

        fname = tempname() * ".fits"
        Octofitter.savechain(fname, chain)
        chain_loaded = Octofitter.loadchain(fname)

        @test chain[:b_a] ≈ chain_loaded[:b_a]
    end

    @testset "HDF5 Chain I/O (Orbitize format)" begin
        chain = octofit(rng, model, iterations=TEST_ITERATIONS, adaptation=TEST_ADAPTATION)

        fname = tempname() * ".hdf5"
        Octofitter.savehdf5(fname, model, chain)
        chain_loaded = Octofitter.loadhdf5(fname)

        # Test that key orbital parameters are preserved through the round-trip
        # Note: HDF5 format uses Orbitize! conventions so some parameters have different names
        @test median(chain[:b_a]) ≈ median(chain_loaded[:b_a]) rtol=1e-4
        @test median(chain[:b_e]) ≈ median(chain_loaded[:b_e]) rtol=1e-4
        @test median(chain[:b_i]) ≈ median(chain_loaded[:b_i]) rtol=1e-4
        @test median(chain[:b_ω]) ≈ median(chain_loaded[:b_ω]) rtol=1e-4
        @test median(chain[:b_Ω]) ≈ median(chain_loaded[:b_Ω]) rtol=1e-4
        @test median(chain[:M]) ≈ median(chain_loaded[:M]) rtol=1e-4
        @test median(chain[:plx]) ≈ median(chain_loaded[:plx]) rtol=1e-4
    end
end

@testset "Autodiff Gradient Comparison" begin
    # Test that FiniteDiff and ForwardDiff produce consistent gradients
    # This verifies the correctness of the autodiff implementation
    rng = Random.Xoshiro(42)

    astrom_like = PlanetRelAstromLikelihood(
        Table(
            epoch = [50000, 50120, 50240, 50360],
            ra = [-505.76, -502.57, -498.21, -492.68],
            dec = [-66.93, -37.47, -7.93, 21.64],
            σ_ra = fill(10.0, 4),
            σ_dec = fill(10.0, 4),
            cor = fill(0.0, 4)
        ),
        name = "gradient_test"
    )

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom_like],
        variables=@variables begin
            a ~ Uniform(0, 100)
            e ~ Uniform(0.0, 0.99)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    GradTestSys = System(
        name="GradTestSys",
        companions=[b],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        end
    )

    # Create models with different autodiff backends
    model_fd = Octofitter.LogDensityModel(GradTestSys, autodiff=AutoFiniteDiff(), verbosity=0)
    model_ad = Octofitter.LogDensityModel(GradTestSys, autodiff=AutoForwardDiff(), verbosity=0)

    # Get a valid starting point by sampling from the prior and transforming
    theta_natural = collect(model_ad.sample_priors(rng))
    theta = model_ad.link(theta_natural)

    # Compare gradients from both backends
    _, grad_fd = model_fd.∇ℓπcallback(theta)
    _, grad_ad = model_ad.∇ℓπcallback(theta)

    # Use relaxed tolerances since finite differences have inherent numerical error
    @test grad_fd ≈ grad_ad atol=1e-3 rtol=1e-4
end
