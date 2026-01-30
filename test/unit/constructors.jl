@testset "Variables and Priors" begin
    @testset "@variables macro" begin
        # Test that @variables creates a valid tuple of (Priors, Derived)
        vars = @variables begin
            a ~ Uniform(0,1)
            b ~ Normal(0,1)
            c = a + b
        end
        @test vars isa Tuple
        @test length(vars) == 2

        # Test that variables are correctly captured
        priors = first(vars)
        @test haskey(priors.priors, :a)
        @test haskey(priors.priors, :b)
        @test priors.priors[:a] isa Uniform
        @test priors.priors[:b] isa Normal

        derived = vars[2]
        @test haskey(derived.variables, :c)
    end

    @testset "UniformCircular" begin
        uc = UniformCircular()
        @test uc isa UniformCircular
        @test uc.domain ≈ 2π

        uc = UniformCircular(1.0)
        @test uc.domain ≈ 1.0
    end

    @testset "Sine distribution" begin
        sine = Sine()
        @test sine isa Sine
        @test Distributions.minimum(sine) ≈ 0.0 + eps()
        @test Distributions.maximum(sine) ≈ π - eps()
    end
end

@testset "Basic Planet Construction" begin
    astrom = PlanetRelAstromLikelihood(
        Table(epoch=[5000.0], ra=[100.0], dec=[50.0], σ_ra=[1.0], σ_dec=[1.0]),
        name="test_astrom"
    )

    planet = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom],
        variables=@variables begin
            a ~ Uniform(0,1)
            e ~ Beta(1,2)
        end
    )
    @test planet isa Planet
    @test planet.name == :b
    @test length(planet.observations) == 1
    @test first(planet.observations) === astrom

    planet_derived = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom],
        variables=@variables begin
            a ~ Uniform(0,1)
            e ~ Beta(1,2)
            tp = a * e
        end
    )
    @test planet_derived isa Planet
end

@testset "Fixed Position Orbit" begin
    fp = Octofitter.FixedPosition(1.0, 2.0, 3.0)
    @test fp.x == 1.0
    @test fp.y == 2.0
    @test fp.z == 3.0

    fp = Octofitter.FixedPosition(x=1.0, y=2.0, z=3.0)
    @test fp.x == 1.0
    @test fp.y == 2.0
    @test fp.z == 3.0

    @test period(fp) == Inf
    @test meanmotion(fp) == 0.0
    @test eccentricity(fp) == 0.0
    @test totalmass(fp) == 0.0
    @test semimajoraxis(fp) == 0.0
    @test periastron(fp) == 0.0
end

@testset "Orbit Model Types" begin
    @testset "KepOrbit Model" begin
        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ truncated(Normal(10, 4), lower=0.1)
                e ~ Uniform(0.0, 0.5)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            end
        )

        SimpleSystem = System(
            name="SimpleSystem",
            companions=[b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        @test Octofitter.orbittype(b) == Visual{KepOrbit}
        model = Octofitter.LogDensityModel(SimpleSystem)
        @test model isa Octofitter.LogDensityModel
    end

    @testset "ThieleInnes Model" begin
        b = Planet(
            name="b",
            basis=ThieleInnesOrbit,
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                A ~ Normal(0, 10000)
                B ~ Normal(0, 10000)
                F ~ Normal(0, 10000)
                G ~ Normal(0, 10000)
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, A, B, F, G, plx=system.plx)
            end
        )

        TISystem = System(
            name="TISystem",
            companions=[b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        @test Octofitter.orbittype(b) == ThieleInnesOrbit
        model = Octofitter.LogDensityModel(TISystem)
        @test model isa Octofitter.LogDensityModel
    end

    @testset "FixedPosition Model" begin
        b = Planet(
            name="b",
            basis=Visual{Octofitter.FixedPosition},
            observations=[],
            variables=@variables begin
                x ~ Uniform(-2000, 2000)
                y ~ Uniform(-2000, 2000)
            end
        )

        FixedSystem = System(
            name="FixedSystem",
            companions=[b],
            observations=[],
            variables=@variables begin
                plx = 24.4620
            end
        )

        @test Octofitter.orbittype(b) == Visual{Octofitter.FixedPosition}
        model = Octofitter.LogDensityModel(FixedSystem)
        @test model isa Octofitter.LogDensityModel
    end
end
