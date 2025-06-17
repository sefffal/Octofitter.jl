using Test
using Octofitter
using Distributions
using OctofitterRadialVelocity
using CairoMakie
using Random
@testset "RV Likelihoods" begin
    @testset "StarAbsoluteRVLikelihood" begin
        # Create test data
        rvlike = StarAbsoluteRVLikelihood(
            (epoch=50000.0, rv=100.0, σ_rv=10.0),
            (epoch=50100.0, rv=110.0, σ_rv=10.0),
            name="RVDat",
            offset=:rv0,
            jitter=:jitter
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
        rvlike = MarginalizedStarAbsoluteRVLikelihood(
            Table(
                epoch=[50000.0, 50100.0],
                rv=[100.0, 110.0],
                σ_rv=[10.0, 10.0]
            ),
            jitter=:jitter
        )
        
        @test rvlike isa MarginalizedStarAbsoluteRVLikelihood
        @test length(rvlike.table) == 2
        @test hasproperty(rvlike.table, :epoch)
        
        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(rvlike, 1:1)
        @test length(subset.table) == 1
    end

    @testset "PlanetRelativeRVLikelihood" begin
        # Create test data
        rvlike = PlanetRelativeRVLikelihood(
            (epoch=50000.0, rv=100.0, σ_rv=10.0),
            (epoch=50100.0, rv=110.0, σ_rv=10.0),
            name="test",
            jitter=:gamma
        )
        
        @test rvlike isa PlanetRelativeRVLikelihood
        @test length(rvlike.table) == 2
        @test all([:epoch, :rv, :σ_rv] .∈ Ref(propertynames(rvlike.table)))
        
        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(rvlike, 1:1)
        @test length(subset.table) == 1
    end
end

@testset "RadialVelocityOrbit Models" begin
    # Create test RV data
    rvlike = StarAbsoluteRVLikelihood(
        (epoch=50000.0, rv=100.0, σ_rv=10.0),
        name="test",
        offset=:rv0,
        jitter=:jitter
    )

    @planet b RadialVelocityOrbit begin
        e ~ Uniform(0.0, 0.5)
        ω ~ UniformCircular()
        τ ~ Uniform(0,1)
        P ~ Uniform(0.001, 10)
        a = ∛(b.P^2*system.M)
        tp =  b.τ*b.P*365.25 + 50000
        mass ~ Uniform(0,100)
    end

    @system RVSystem begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        jitter ~ LogUniform(0.1, 100)
        rv0 ~ Normal(0, 100)
    end rvlike b

    # Test orbit type extraction
    @test Octofitter.orbittype(b) == RadialVelocityOrbit

    # Test element construction
    θ_system = (
        M = 1.2,
        jitter = 0.2,
        rv0 = 0,
        planets = (
            b = (
                P = 5.0,
                a = 3.1072325059538586,
                e = 0.1,
                ω = 1.0,
                τ = 0.0,
                tp = 50000.0,
                mass=1.0
            ),
        )
    )
    
    orbit = Octofitter.construct_elements(RadialVelocityOrbit, θ_system, θ_system.planets.b)
    @test orbit isa RadialVelocityOrbit
    @test semimajoraxis(orbit) ≈ 3.1072325059538586
    @test period(orbit) ≈ 5*365.25 rtol=1e-3

    model = Octofitter.LogDensityModel(RVSystem)
    chain = octofit(model, adaptation=100,iterations=100)

    fig = octoplot(model, chain)
    @test fig isa Makie.Figure

    fig = Octofitter.rvpostplot(model,chain)
    @test fig isa Makie.Figure
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
    rvlike = StarAbsoluteRVLikelihood(
        dat...,
        name="HARPS",
        offset=:rv0,
        jitter=:jitter,
        gaussian_process=quasi_periodic_kernel
    )

    @planet b RadialVelocityOrbit begin
        e = 0.0
        ω = 0.0
        mass = 0.0
        τ = 0.0
        a = 1.0
        tp =  50000.0
    end

    @system GPSystem begin
        M = 1.0
        jitter = 0.0
        rv0 = 0.0
        
        # GP hyperparameters
        gp_η₁ ~ LogUniform(0.1, 100)  # RV amplitude
        gp_η₂ ~ LogUniform(5, 100)    # decay timescale
        gp_η₃ ~ Normal(25, 1)         # rotation period
        gp_η₄ ~ LogUniform(0.2, 10)   # evolution timescale
    end rvlike b

    model = Octofitter.LogDensityModel(GPSystem)
    chain = octofit(model, adaptation=50, iterations=50)

    # Test GP hyperparameters are recovered
    @test 20 < median(chain[:gp_η₃]) < 30 

    fig = Octofitter.rvpostplot(model,chain)
    @test fig isa Makie.Figure
end
