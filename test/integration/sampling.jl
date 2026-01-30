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
end
