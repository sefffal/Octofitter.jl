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
