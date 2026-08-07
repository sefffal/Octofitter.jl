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

    # Regression: observation-level variables of the *wrapped* likelihood must
    # still reach it once it is wrapped. `ObsPriorAstromONeil2019` has no `name`
    # field (it supplies its name via `likelihoodname`), and `make_ln_like` used
    # to gate the θ_obs lookup on `hasproperty(like, :name)` — so jitter,
    # platescale and northangle were silently replaced by their defaults.
    @testset "wrapped likelihood receives θ_obs" begin
        astrom_tbl = Table(
            epoch = [58000.0, 58200.0, 58400.0],
            ra    = [  100.0,   110.0,   120.0],
            dec   = [  100.0,    95.0,    90.0],
            σ_ra  = [    5.0,     5.0,     5.0],
            σ_dec = [    5.0,     5.0,     5.0],
        )

        function jitter_sensitivity(; wrap::Bool)
            obs_vars = @variables begin
                jitter = planet.shared_astrom_jitter
            end
            radec = PlanetRelAstromLikelihood(astrom_tbl, name="d_radec", variables=obs_vars)
            obs = wrap ? ObsPriorAstromONeil2019(radec) : radec

            # Everything fixed except the jitter, so ln_like is a function of it alone.
            planet_vars = @variables begin
                shared_astrom_jitter ~ Uniform(0.0, 500.0)
                a = 10.0
                e = 0.2
                i = 0.5
                ω = 0.3
                Ω = 0.4
                tp = 58000.0
            end
            b = Planet(name=:b, basis=Visual{KepOrbit},
                       variables=planet_vars, observations=(obs,))
            sys_vars = @variables begin
                M = 1.0
                plx = 50.0
            end
            sys = System(name=:jitter_prop_test, variables=sys_vars, companions=(b,))

            model = Octofitter.LogDensityModel(sys)
            θ_lo = model.arr2nt([0.001])
            θ_hi = model.arr2nt([300.0])
            lnlike = Octofitter.make_ln_like(model.system, θ_lo)
            return lnlike(model.system, θ_hi) - lnlike(model.system, θ_lo)
        end

        Δ_plain = jitter_sensitivity(wrap=false)
        Δ_wrapped = jitter_sensitivity(wrap=true)

        # Sanity check: jitter matters at all on the unwrapped likelihood.
        @test Δ_plain > 1

        # The O'Neil jacobian term the wrapper adds does not depend on jitter,
        # so the sensitivity to jitter must be identical with and without it.
        @test Δ_wrapped ≈ Δ_plain
    end
end
