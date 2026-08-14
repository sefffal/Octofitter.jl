@testset "Prior Specifications" begin
    @testset "Standard Distribution Priors" begin
        # Test normal prior construction
        @test_nowarn Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ truncated(Normal(10, 2), lower=0.1)
                e ~ Beta(1, 2)
                i ~ Normal(1.0, 0.1)
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ UniformCircular()
            end
        )

        # Test LogUniform prior
        @test_nowarn Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ LogUniform(0.1, 100)
                e ~ Uniform(0, 0.99)
                i ~ Normal(1.0, 0.1)
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ UniformCircular()
            end
        )
    end
end

@testset "Orbital Constraint Priors" begin
    # Create test planets for multi-planet systems
    astrom_b = PlanetRelAstromLikelihood(
        Table(epoch=[50000.0], ra=[100.0], dec=[50.0], σ_ra=[1.0], σ_dec=[1.0]),
        name="astrom_b"
    )
    astrom_c = PlanetRelAstromLikelihood(
        Table(epoch=[50000.0], ra=[-200.0], dec=[-100.0], σ_ra=[1.0], σ_dec=[1.0]),
        name="astrom_c"
    )

    @testset "LimitClosestApproachAUPrior" begin
        # Test constructors
        prior_full = LimitClosestApproachAUPrior(0.5, 2.0)
        @test prior_full isa LimitClosestApproachAUPrior
        @test prior_full.hard_closest_approach_au == 0.5
        @test prior_full.soft_closest_approach_au == 2.0

        # Test convenience constructor (hard=0)
        prior_soft_only = LimitClosestApproachAUPrior(1.0)
        @test prior_soft_only.hard_closest_approach_au == 0
        @test prior_soft_only.soft_closest_approach_au == 1.0

        # Test NonCrossingPrior alias
        prior_non_crossing = NonCrossingPrior()
        @test prior_non_crossing isa LimitClosestApproachAUPrior
        @test prior_non_crossing.hard_closest_approach_au == 0.0
        @test prior_non_crossing.soft_closest_approach_au == 0.0

        # Test with single planet system (should return 0)
        b_single = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 100)
                e ~ Uniform(0, 0.5)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
        sys_single = System(
            name="SinglePlanet",
            companions=[b_single],
            observations=[prior_full],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
            end
        )
        model_single = Octofitter.LogDensityModel(sys_single, verbosity=0)
        @test model_single isa Octofitter.LogDensityModel

        # Test with multi-planet system
        b_inner = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 10)  # Inner planet
                e ~ Uniform(0, 0.3)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
        c_outer = Planet(
            name="c",
            basis=Visual{KepOrbit},
            observations=[astrom_c],
            variables=@variables begin
                a ~ LogUniform(15, 50)  # Outer planet (well-separated)
                e ~ Uniform(0, 0.3)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
        sys_multi = System(
            name="MultiPlanet",
            companions=[b_inner, c_outer],
            observations=[prior_full],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
            end
        )
        model_multi = Octofitter.LogDensityModel(sys_multi, verbosity=0)
        @test model_multi isa Octofitter.LogDensityModel

        # Test that we can draw from priors and evaluate
        params = Octofitter.drawfrompriors(model_multi.system)
        @test haskey(params.planets, :b)
        @test haskey(params.planets, :c)
    end

    @testset "HillStabilityPrior" begin
        # Test constructor
        prior = HillStabilityPrior()
        @test prior isa HillStabilityPrior

        # Test with single planet (should return 0)
        b_single = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 100)
                e ~ Uniform(0, 0.5)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.1, 10)  # Mass needed for Hill stability
            end
        )
        sys_single = System(
            name="SinglePlanetHill",
            companions=[b_single],
            observations=[prior],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
            end
        )
        model_single = Octofitter.LogDensityModel(sys_single, verbosity=0)
        @test model_single isa Octofitter.LogDensityModel

        # Test with multi-planet system (well-separated for stability)
        b_inner = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 5)
                e ~ Uniform(0, 0.1)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.1, 5)
                M = system.M  # Needed for Hill stability calculation
            end
        )
        c_outer = Planet(
            name="c",
            basis=Visual{KepOrbit},
            observations=[astrom_c],
            variables=@variables begin
                a ~ LogUniform(20, 50)  # Well-separated
                e ~ Uniform(0, 0.1)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.1, 5)
                M = system.M
            end
        )
        sys_multi = System(
            name="MultiPlanetHill",
            companions=[b_inner, c_outer],
            observations=[prior],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
            end
        )
        model_multi = Octofitter.LogDensityModel(sys_multi, verbosity=0)
        @test model_multi isa Octofitter.LogDensityModel
    end

    @testset "PlanetOrderPrior" begin
        # Create planets for ordering test
        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 10)
                e ~ Uniform(0, 0.3)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
        c = Planet(
            name="c",
            basis=Visual{KepOrbit},
            observations=[astrom_c],
            variables=@variables begin
                a ~ LogUniform(15, 50)
                e ~ Uniform(0, 0.3)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )

        # Test constructor
        prior = PlanetOrderPrior(b, c)
        @test prior isa PlanetOrderPrior
        @test length(prior.planets) == 2
        @test prior.planets[1].name == :b
        @test prior.planets[2].name == :c

        # Test likelihoodname
        name = Octofitter.likelihoodname(prior)
        @test occursin("b", name)
        @test occursin("c", name)

        # Test with system
        sys = System(
            name="OrderedSystem",
            companions=[b, c],
            observations=[prior],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
            end
        )
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        @test model isa Octofitter.LogDensityModel

        # Test that we can draw from priors
        params = Octofitter.drawfrompriors(model.system)
        @test haskey(params.planets, :b)
        @test haskey(params.planets, :c)
    end

    @testset "UnitLengthPrior" begin
        # UnitLengthPrior is automatically created by UniformCircular()
        # Test that UniformCircular creates a UnitLengthPrior

        # Test planet-level UnitLengthPrior (via UniformCircular)
        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 100)
                e ~ Uniform(0, 0.5)
                i ~ Sine()
                ω ~ UniformCircular()  # Creates UnitLengthPrior
                Ω ~ UniformCircular()  # Creates UnitLengthPrior
                θ ~ UniformCircular()  # Creates UnitLengthPrior
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
        sys = System(
            name="UnitLengthTest",
            companions=[b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
            end
        )
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        @test model isa Octofitter.LogDensityModel

        # Verify UnitLengthPriors were created (check planet has x/y components)
        params = Octofitter.drawfrompriors(model.system)
        @test haskey(params.planets.b, :ωx)
        @test haskey(params.planets.b, :ωy)
        @test haskey(params.planets.b, :Ωx)
        @test haskey(params.planets.b, :Ωy)
        @test haskey(params.planets.b, :θx)
        @test haskey(params.planets.b, :θy)

        # Test system-level UnitLengthPrior (via UniformCircular)
        b_shared = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[astrom_b],
            variables=@variables begin
                a ~ LogUniform(1, 100)
                e ~ Uniform(0, 0.5)
                i ~ Sine()
                ω ~ Uniform(0, 2π)
                Ω = system.Ω  # Shared from system
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
        sys_shared = System(
            name="SharedOmega",
            companions=[b_shared],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.0, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
                Ω ~ UniformCircular()  # System-level UniformCircular
            end
        )
        model_shared = Octofitter.LogDensityModel(sys_shared, verbosity=0)
        @test model_shared isa Octofitter.LogDensityModel

        # Verify system has Ωx and Ωy
        params_shared = Octofitter.drawfrompriors(model_shared.system)
        @test haskey(params_shared, :Ωx)
        @test haskey(params_shared, :Ωy)
    end
end

@testset "Prior Sampling" begin
    astrom_data = PlanetRelAstromLikelihood(
        Table(epoch = [50000], ra = [100.0], dec = [50.0], σ_ra = [1.0], σ_dec = [1.0]),
        name="prior_sampling_test"
    )

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom_data],
        variables=@variables begin
            a ~ LogUniform(0.1, 100)
            e ~ Uniform(0, 0.99)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    TestSystem = System(
        name="TestSystem",
        companions=[b],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        end
    )

    # Test drawing from priors
    samples = Octofitter.sample_priors(TestSystem)
    @test length(samples) == 11  # Should match number of free parameters

    # Test multiple samples
    n_samples = 100
    samples_n = Octofitter.sample_priors(TestSystem, n_samples)
    @test length(samples_n) == n_samples
    @test all(length.(samples_n) .== 11)

    # Verify samples are within bounds
    model = Octofitter.LogDensityModel(TestSystem)
    arr2nt = Octofitter.make_arr2nt(TestSystem)

    for sample in samples_n
        params = arr2nt(sample)
        @test 0.1 ≤ params.planets.b.a ≤ 100
        @test 0 ≤ params.planets.b.e ≤ 0.99
        @test 0 ≤ params.planets.b.i ≤ π
        @test params.M ≥ 0.1
        @test params.plx ≥ 0.1
    end
end
