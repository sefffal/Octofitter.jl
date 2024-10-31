using FiniteDiff
using Random
@testset "FIT: Relative Astrometry" begin
    rng = Random.Xoshiro(1)
    astrom_like = PlanetRelAstromLikelihood(
        # Your data here:
        # units are MJD, mas, mas, mas, mas, and correlation.
        (epoch = 50000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50120, ra = -502.570356287689, dec = -37.47217527025044, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50240, ra = -498.2089148883798, dec = -7.927548139010479, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50360, ra = -492.67768482682357, dec = 21.63557115669823, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50480, ra = -485.9770335870402, dec = 51.147204404903704, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50600, ra = -478.1095526888573, dec = 80.53589069730698, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50720, ra = -469.0801731788123, dec = 109.72870493064629, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 50840, ra = -458.89628893460525, dec = 138.65128697876773, σ_ra = 10, σ_dec = 10, cor=0),
    )
    @planet b Visual{KepOrbit} begin
        a ~ Uniform(0, 100) # AU
        e ~ Uniform(0.0, 0.99)
        i ~ Sine() # radians
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(system,b,50000) # use MJD epoch of your data here!!
    end astrom_like
    @system TestSys begin # replace TestSys with the name of your planetary system
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end b
    model_a = Octofitter.LogDensityModel(TestSys, autodiff=:FiniteDiff)
    # Correct dimensionality of problem
    @test model.D == 11
    # Test auto-diff
    model_b = Octofitter.LogDensityModel(TestSys, autodiff=:ForwardDiff)
    theta = collect(model.link(Octofitter.guess_starting_position(rng, model, 100)[1]))
    @test (model_a.∇ℓπcallback(theta)[2]) ≈ (model_b.∇ℓπcallback(theta)[2]) atol=1e-5 rtol=1e-5

    Octofitter.default_initializer!(rng, model_b)
    @test length(model_b.starting_points) > 0
    @test model_a.ℓπcallback(model_b.starting_points[1]) > -1000

    # Test fitting a chain
    chain = octofit(rng, model_b)
    @test all(chain[:logpost] .> -1000)

    # Test no numeical divergences
    @test mean(chain[:numerical_error]) < 0.05

    # Test saving and loading a chain
    Octofitter.savechain("mychain.fits", chain)
    chain_2 = Octofitter.loadchain("mychain.fits")

    @test chain[:b_a] ≈ chain_2[:b_a]
end