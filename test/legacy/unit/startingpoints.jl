using Random

# Model with only plain (non-composite) priors, so a starting point can name every
# free variable directly.
function _startingpoints_testmodel()
    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[],
        variables=@variables begin
            a ~ Uniform(1, 20)
            e ~ Uniform(0.0, 0.5)
            i ~ Sine()
            ω ~ Uniform(0, 2pi)
            Ω ~ Uniform(0, 2pi)
            θ ~ Uniform(0, 2pi)
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )
    sys = System(
        name="StartingPointsTestSys",
        companions=[b],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.0, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 1.0), lower=0.1)
        end
    )
    return Octofitter.LogDensityModel(sys; verbosity=0)
end

const _SP_POINT = (;
    M=1.05,
    plx=50.0,
    planets=(; b=(; a=5.0, e=0.3, i=0.6, ω=1.2, Ω=2.4, θ=0.8))
)

@testset "startingpoints!" begin

    @testset "Single point" begin
        model = _startingpoints_testmodel()
        chn = startingpoints!(model, _SP_POINT; verbosity=0)

        @test length(model.starting_points) == 1000
        @test all(==(model.starting_points[1]), model.starting_points)

        # The stored points are in the unconstrained space; they should map back to
        # exactly the values that were requested.
        rt = model.arr2nt(model.invlink(model.starting_points[1]))
        @test rt.M ≈ 1.05
        @test rt.plx ≈ 50.0
        @test rt.planets.b.a ≈ 5.0
        @test rt.planets.b.e ≈ 0.3
        @test rt.planets.b.i ≈ 0.6
        @test rt.planets.b.ω ≈ 1.2
        @test rt.planets.b.Ω ≈ 2.4
        @test rt.planets.b.θ ≈ 0.8

        @test isfinite(model.ℓπcallback(model.starting_points[1]))

        # Returns an inspectable chain of the starting points, as initialize! does.
        @test chn isa Chains
        @test size(chn, 1) == 1000
    end

    @testset "ndraws" begin
        model = _startingpoints_testmodel()
        startingpoints!(model, _SP_POINT; ndraws=25, verbosity=0)
        @test length(model.starting_points) == 25

        # The initial mass matrix is estimated from the spread of the points, so a
        # single stored point is rejected rather than producing a degenerate estimate.
        @test_throws ErrorException startingpoints!(model, _SP_POINT; ndraws=1, verbosity=0)
    end

    @testset "Multiple points" begin
        model = _startingpoints_testmodel()
        point2 = (;
            M=1.0,
            plx=49.0,
            planets=(; b=(; a=6.0, e=0.2, i=0.5, ω=1.0, Ω=2.0, θ=0.5))
        )
        startingpoints!(model, _SP_POINT, point2; verbosity=0)

        @test length(model.starting_points) == 2
        @test model.starting_points[1] != model.starting_points[2]
        @test model.arr2nt(model.invlink(model.starting_points[1])).planets.b.a ≈ 5.0
        @test model.arr2nt(model.invlink(model.starting_points[2])).planets.b.a ≈ 6.0
    end

    @testset "Incomplete point is rejected" begin
        model = _startingpoints_testmodel()
        err = try
            startingpoints!(model, (; M=1.05, planets=(; b=(; a=5.0))); verbosity=0)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        # The error should name what is missing and point at initialize!.
        msg = sprint(showerror, err)
        @test occursin("plx", msg)
        @test occursin("planets.b.e", msg)
        @test occursin("initialize!", msg)

        # A rejected point must not have clobbered the model.
        @test isnothing(model.starting_points)
    end

    @testset "No points is rejected" begin
        model = _startingpoints_testmodel()
        @test_throws ErrorException startingpoints!(model; verbosity=0)
    end

    @testset "Out-of-support point warns" begin
        model = _startingpoints_testmodel()
        # e = 0.9 is outside the Uniform(0, 0.5) prior support.
        bad = (;
            M=1.05,
            plx=50.0,
            planets=(; b=(; a=5.0, e=0.9, i=0.6, ω=1.2, Ω=2.4, θ=0.8))
        )
        @test_logs (:warn,) match_mode=:any startingpoints!(model, bad)
        @test !isfinite(model.ℓπcallback(model.starting_points[1]))
    end

    @testset "No randomness in the result" begin
        # The point of this API (#123): unlike initialize!, the result does not depend
        # on RNG state at all.
        model = _startingpoints_testmodel()
        Random.seed!(1)
        startingpoints!(model, _SP_POINT; verbosity=0)
        first_run = copy(model.starting_points[1])
        Random.seed!(12345)
        startingpoints!(model, _SP_POINT; verbosity=0)
        @test model.starting_points[1] == first_run
    end
end
