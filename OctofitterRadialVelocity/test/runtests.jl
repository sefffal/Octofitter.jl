using Test
using Octofitter
using Distributions
using OctofitterRadialVelocity
using CairoMakie
using Random
@testset "RV Likelihoods" begin
    @testset "StarAbsoluteRVLikelihood" begin
        # Create test data
        vars = @variables begin
            offset ~ Uniform(-1000, 1000)
            jitter ~ LogUniform(0.001, 100)
        end
        rvlike = StarAbsoluteRVLikelihood(
            Table(
                epoch=[50000.0, 50100.0],
                rv=[100.0, 110.0],
                σ_rv=[10.0, 10.0]
            ),
            name="RVDat",
            variables=vars,
        )

        @test rvlike isa StarAbsoluteRVLikelihood
        @test length(rvlike.table) == 2
        @test all([:epoch, :rv, :σ_rv] .∈ Ref(propertynames(rvlike.table)))

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(rvlike, 1:1)
        @test length(subset.table) == 1
    end

    @testset "MarginalizedStarAbsoluteRVLikelihood" begin
        # Create test data
        vars = @variables begin
            jitter ~ LogUniform(0.001, 100)
        end
        rvlike = MarginalizedStarAbsoluteRVLikelihood(
            Table(
                epoch=[50000.0, 50100.0],
                rv=[100.0, 110.0],
                σ_rv=[10.0, 10.0]
            ),
            name="MargRV",
            variables=vars,
        )

        @test rvlike isa MarginalizedStarAbsoluteRVLikelihood
        @test length(rvlike.table) == 2
        @test hasproperty(rvlike.table, :epoch)

        # Test that subsetting is not supported (analytical marginalization couples all obs)
        @test_throws ErrorException Octofitter.likeobj_from_epoch_subset(rvlike, 1:1)
    end

    @testset "PlanetRelativeRVLikelihood" begin
        # Create test data
        vars = @variables begin
            jitter ~ LogUniform(0.001, 100)
        end
        rvlike = PlanetRelativeRVLikelihood(
            (epoch=50000.0, rv=100.0, σ_rv=10.0),
            (epoch=50100.0, rv=110.0, σ_rv=10.0),
            name="test",
            variables=vars,
        )

        @test rvlike isa PlanetRelativeRVLikelihood
        @test length(rvlike.table) == 2
        @test all([:epoch, :rv, :σ_rv] .∈ Ref(propertynames(rvlike.table)))

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(rvlike, 1:1)
        @test length(subset.table) == 1
    end

    @testset "PlanetRelativeRVLikelihood with offset and trend" begin
        # Test construction with offset and trend_function
        trend_fn = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 50000)
        vars = @variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.001, 100)
            trend_slope ~ Normal(0, 1)
        end
        rvlike = PlanetRelativeRVObs(
            (epoch=50000.0, rv=100.0, σ_rv=10.0),
            (epoch=50050.0, rv=115.0, σ_rv=10.0),
            (epoch=50100.0, rv=130.0, σ_rv=10.0),
            name="test_trend",
            trend_function=trend_fn,
            variables=vars,
        )

        @test rvlike isa PlanetRelativeRVObs
        @test length(rvlike.table) == 3
        @test rvlike.trend_function === trend_fn
        @test all([:epoch, :rv, :σ_rv] .∈ Ref(propertynames(rvlike.table)))

        # Test subsetting preserves trend_function
        subset = Octofitter.likeobj_from_epoch_subset(rvlike, 1:2)
        @test length(subset.table) == 2
        @test subset.trend_function === trend_fn

        # Test construction without trend_function uses zero default
        vars2 = @variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.001, 100)
        end
        rvlike_no_trend = PlanetRelativeRVObs(
            (epoch=50000.0, rv=100.0, σ_rv=10.0),
            (epoch=50100.0, rv=110.0, σ_rv=10.0),
            name="test_no_trend",
            variables=vars2,
        )
        @test rvlike_no_trend.trend_function !== nothing
    end
end

@testset "RadialVelocityOrbit Models" begin
    # Create test RV data
    rv_vars = @variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
    end
    rvlike = StarAbsoluteRVLikelihood(
        Table(epoch=[50000.0], rv=[100.0], σ_rv=[10.0]),
        name="test",
        variables=rv_vars,
    )

    b = Planet(
        name=:b,
        basis=RadialVelocityOrbit,
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            e ~ Uniform(0.0, 0.5)
            ω ~ UniformCircular()
            τ ~ Uniform(0,1)
            P ~ Uniform(0.001, 10)
            a = ∛(P^2*M)
            tp = τ*P*365.25 + 50000
            mass ~ Uniform(0,100)
        end
    )

    RVSystem = System(
        name=:RVSystem,
        companions=[b],
        observations=[rvlike],
        variables=@variables begin
        end
    )

    # Test orbit type extraction
    @test Octofitter.orbittype(b) == RadialVelocityOrbit

    model = Octofitter.LogDensityModel(RVSystem)
    chain = octofit(model, adaptation=100,iterations=100)

    fig = octoplot(model, chain)
    @test fig isa Makie.Figure

    fig = Octofitter.rvpostplot(model,chain)
    @test fig isa Makie.Figure
end



@testset "PlanetRelativeRV with offset and trend" begin
    using PlanetOrbits

    # Generate synthetic relative RV data: known planet orbit + offset + linear trend
    true_offset = 50.0    # m/s
    true_slope = 0.1      # m/s/day
    ref_epoch = 50000.0
    true_P = 80.0         # days
    true_M = 1.0
    true_a = ∛(true_P^2 * true_M)

    epochs = collect(range(50000.0, 50200.0, length=20))

    # Generate planet RV signal from known orbit
    true_orbit = RadialVelocityOrbit(a=true_a, e=0.0, ω=0.0, tp=ref_epoch, M=true_M)
    planet_rv = [radvel(orbitsolve(true_orbit, ep)) for ep in epochs]

    # Observed = planet + offset + trend + noise
    Random.seed!(42)
    observed_rv = planet_rv .+ true_offset .+ true_slope .* (epochs .- ref_epoch) .+ randn(length(epochs)) .* 1.0

    vars = @variables begin
        offset ~ Normal(0, 200)
        jitter ~ LogUniform(0.01, 50)
        trend_slope ~ Normal(0, 1)
    end
    rvlike = PlanetRelativeRVObs(
        [(epoch=e, rv=r, σ_rv=1.0) for (e, r) in zip(epochs, observed_rv)]...,
        name="RelRV",
        trend_function=(θ_obs, epoch) -> θ_obs.trend_slope * (epoch - ref_epoch),
        variables=vars,
    )

    # Fix planet orbit so only offset + trend are free
    b = Planet(
        name=:b,
        basis=RadialVelocityOrbit,
        observations=[rvlike],
        variables=@variables begin
            M = $true_M
            e = 0.0
            ω = 0.0
            τ = 0.0
            P = $true_P
            a = ∛(P^2 * M)
            tp = $ref_epoch
            mass = 0.0
        end
    )

    RelRVSys = System(
        name=:RelRVSys,
        companions=[b],
        variables=@variables begin
        end
    )

    model = Octofitter.LogDensityModel(RelRVSys)
    chain = octofit(model, adaptation=200, iterations=200)

    # Observation-level variables are prefixed with planet_observation in the chain
    @test 30 < median(chain[:b_RelRV_offset]) < 70
    @test 0.0 < median(chain[:b_RelRV_trend_slope]) < 0.3
end

using AbstractGPs
@testset "Gaussian Process RV Modeling" begin
    # Create test RV data with stellar activity
    function quasi_periodic_kernel(θ_system)
        η₁ = θ_system.gp_η₁ # amplitude
        η₂ = θ_system.gp_η₂ # decay timescale
        η₃ = θ_system.gp_η₃ # periodicity
        η₄ = θ_system.gp_η₄ # coherence scale
        kernel = η₁^2 * 
            (SqExponentialKernel() ∘ ScaleTransform(1/η₂)) *
            (PeriodicKernel(r=[η₄]) ∘ ScaleTransform(1/η₃))
        return GP(kernel)
    end

    # Create RV data with signal, with frequency that the 
    # RV orbit is constrained against matching
    dat = [
        (epoch=50000.0+t, rv=100.0+10sin(t/25*2pi), σ_rv=0.5)
        for t in 0:3.9:100
    ]
    gp_vars = @variables begin
        offset = 0.0
        jitter = 0.0
        # GP hyperparameters
        gp_η₁ ~ LogUniform(0.1, 100)  # RV amplitude
        gp_η₂ ~ LogUniform(5, 100)    # decay timescale
        gp_η₃ ~ Normal(25, 1)         # rotation period
        gp_η₄ ~ LogUniform(0.2, 10)   # evolution timescale
    end
    rvlike = StarAbsoluteRVLikelihood(
        Table(dat),
        name="HARPS",
        variables=gp_vars,
        gaussian_process=quasi_periodic_kernel
    )

    b = Planet(
        name=:b,
        basis=RadialVelocityOrbit,
        variables=@variables begin
            M = 1.0
            e = 0.0
            ω = 0.0
            mass = 0.0
            τ = 0.0
            a = 1.0
            tp = 50000.0
        end
    )

    GPSystem = System(
        name=:GPSystem,
        companions=[b],
        observations=[rvlike],
        variables=@variables begin
        end
    )

    model = Octofitter.LogDensityModel(GPSystem)
    chain = octofit(model, adaptation=50, iterations=50)

    # Test GP hyperparameters are recovered (obs variables are prefixed with obs name)
    @test 20 < median(chain[:HARPS_gp_η₃]) < 30

    fig = Octofitter.rvpostplot(model,chain)
    @test fig isa Makie.Figure
end
