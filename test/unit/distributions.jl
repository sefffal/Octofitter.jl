@testset "Special Prior Distributions" begin
    @testset "Sine Distribution" begin
        sine = Sine()
        @test Distributions.minimum(sine) ≈ 0.0 + eps()
        @test Distributions.maximum(sine) ≈ π - eps()
        @test Distributions.insupport(sine, π/2)
        @test !Distributions.insupport(sine, -0.1)
        @test !Distributions.insupport(sine, π + 0.1)

        # Test PDF shape - peak should be at π/2
        @test Distributions.pdf(sine, π/2) > Distributions.pdf(sine, 0.1)
        @test Distributions.pdf(sine, π/2) > Distributions.pdf(sine, π-0.1)
    end

    @testset "UniformCircular" begin
        uc = UniformCircular()
        @test uc.domain ≈ 2π

        uc_custom = UniformCircular(1.0)
        @test uc_custom.domain ≈ 1.0

        # Test in model context
        @test_nowarn Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ LogUniform(0.1, 100)
                e ~ Uniform(0, 0.99)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular(π)  # Custom domain
                θ ~ UniformCircular()
            end
        )
    end

    @testset "KDE Priors" begin
        samples = randn(1000) .+ 10  # Normal(10,1) samples
        kde = Octofitter.KDEDist(samples)

        @test kde isa Octofitter.KDEDist
        @test Distributions.minimum(kde) ≈ minimum(samples)
        @test Distributions.maximum(kde) ≈ maximum(samples)
        @test Distributions.insupport(kde, 10.0)
        @test !Distributions.insupport(kde, minimum(samples) - 1)

        # Test in model context
        @test_nowarn Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ kde
                e ~ Uniform(0, 0.99)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )
    end
end

@testset "Observable-based Priors" begin
    astrom_data = PlanetRelAstromLikelihood(
        Table(epoch = [50000, 50100], ra = [100.0, 110.0], dec = [50.0, 55.0], σ_ra = [1.0, 1.0], σ_dec = [1.0, 1.0]),
        name="obs_prior_test"
    )

    obs_prior = ObsPriorAstromONeil2019(astrom_data)
    @test obs_prior isa ObsPriorAstromONeil2019
    @test obs_prior.wrapped_like === astrom_data

    # Test in model context with period prior
    @test_nowarn Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom_data, obs_prior],
        variables=@variables begin
            e ~ Uniform(0.0, 0.5)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            P ~ LogUniform(0.1, 150)
            a = ∛(system.M * P^2)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    # Test subsetting
    subset = Octofitter.likeobj_from_epoch_subset(obs_prior, 1:1)
    @test subset isa ObsPriorAstromONeil2019
end
