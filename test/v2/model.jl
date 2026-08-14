# Model surface: node construction, topology validation, evaluation order.

@testset "node construction" begin
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ Uniform(0.5, 2.0)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 1mjup
        a ~ Uniform(1, 10)
    end)
    @test A.about === nothing
    @test b.about == (:A,)
    # `about=` accepts a node, a Symbol, or a tuple meaning their barycentre
    @test Octofitter.Body(name="c", about=(A, b), variables=@variables begin
        a ~ Uniform(1, 10)
    end).about == (:A, :b)
    @test Octofitter.Body(name="c", about=:A, variables=@variables begin
        a ~ Uniform(1, 10)
    end).about == (:A,)
end

@testset "topology validation" begin
    vars() = @variables begin
        a ~ Uniform(1, 10)
    end
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=vars())
    c = Octofitter.Body(name="c", about=A, variables=vars())

    # No root
    A2 = Octofitter.Body(name="A", about=:b, variables=vars())
    @test_throws r"no root" Octofitter.System(name="s", bodies=[A2, b], variables=@variables begin end)
    # Two roots
    lone = Octofitter.Body(name="z", variables=@variables begin mass = 1.0 end)
    @test_throws r"exactly one root" Octofitter.System(
        name="s", bodies=[A, lone, b], variables=@variables begin end)
    # Wrong number of orbits for the body count: an extra explicit `Orbit`
    # node adds a row without adding a body.
    extra = Octofitter.Orbit(name="extra", exterior=(c,), about=(A, b), variables=vars())
    @test_throws r"exactly 2 orbits" Octofitter.System(
        name="s", bodies=[A, b, c, extra], variables=@variables begin end)
    # Reference to a body that is not in the system
    ghost = Octofitter.Body(name="g", about=:nobody, variables=vars())
    @test_throws r"not a `Body`" Octofitter.System(
        name="s", bodies=[A, ghost], variables=@variables begin end)
    # Duplicate names
    @test_throws r"duplicate node name" Octofitter.System(
        name="s", bodies=[A, b, Octofitter.Body(name="b", about=A, variables=vars())],
        variables=@variables begin end)
end

@testset "frame is chosen by which variables exist" begin
    mk(vars) = Octofitter.System(name="s",
        bodies=[Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end),
                Octofitter.Body(name="b", about=:A, variables=@variables begin a = 1.0 end)],
        variables=vars)
    @test isempty(mk(@variables begin end).framevars)
    @test mk(@variables begin plx = 20.0 end).framevars == [:plx]
    full = mk(@variables begin
        plx = 20.0; ra = 42.0; dec = -31.0; pmra = 5.0; pmdec = -24.0
        rv = 1000.0; ref_epoch = 57388.0
    end)
    @test Set(full.framevars) == Set(Octofitter.FRAME_VARIABLES)
    # A partial absolute frame is an error, not a silent downgrade to parallax.
    @test_throws r"absolute frame needs all of" mk(@variables begin
        plx = 20.0; ra = 42.0; dec = -31.0
    end)
end

@testset "evaluation order: hoisting, and deferred system lines" begin
    # Tier 0: a shared parameter is hoisted to the system block, which is how
    # exact coupling is expressed — one inclination, two users.
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 1mjup
        a = 3.0
        i = system.i_shared
    end)
    c = Octofitter.Body(name="c", about=A, variables=@variables begin
        mass = 1mjup
        a = 9.0
        i = system.i_shared
    end)
    # Tier 1: a system line mentioning a body is deferred until after the
    # body blocks. No new syntax; detected by a static expression walk.
    sys = Octofitter.System(name="s", bodies=[A, b, c], variables=@variables begin
        plx = 20.0
        i_shared ~ Sine()
        mut_inc = b.i - c.i
    end)
    @test sys.deferred == [:mut_inc]
    nt = Octofitter.make_arr2nt(sys)([0.7])
    @test nt.bodies.b.i == nt.i_shared == 0.7
    @test nt.mut_inc == 0.0
    @test haskey(nt, :mut_inc)     # deferred lines still land in the system namespace
end

@testset "a body reading a deferred system variable is a named cycle" begin
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        a = system.derived_from_b
    end)
    err = try
        Octofitter.System(name="s", bodies=[A, b], variables=@variables begin
            derived_from_b = b.a * 2
        end)
        nothing
    catch e
        e
    end
    @test err !== nothing
    msg = sprint(showerror, err)
    @test occursin("circular", msg)
    @test occursin("derived_from_b", msg)
    @test occursin("b", msg)
end

@testset "parameter order is the same everywhere" begin
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ Uniform(0.5, 2.0)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        a ~ Uniform(1, 10)
        e ~ Uniform(0, 0.9)
    end)
    obs = RelAstromObs((epoch=[57000.0], ra=[100.0], dec=[50.0], σ_ra=[1.0], σ_dec=[1.0]);
        target=b, ref=A, name="astrom",
        variables=@variables begin jitter ~ LogUniform(0.1, 10) end)
    sys = Octofitter.System(name="s", bodies=[A, b], observations=[obs],
        variables=@variables begin plx ~ Uniform(1, 100) end)

    names = Octofitter.list_parameter_names(sys)
    @test names == [:plx, :A_mass, :b_a, :b_e, :astrom_jitter]
    @test length(Octofitter._list_priors(sys)) == length(names)

    # arr2nt, the prior density and the prior sampler must agree on the order:
    # feed distinct values and check each lands where its name says.
    nt = Octofitter.make_arr2nt(sys)([11.0, 1.5, 4.0, 0.3, 2.0])
    @test nt.plx == 11.0
    @test nt.bodies.A.mass == 1.5
    @test nt.bodies.b.a == 4.0
    @test nt.bodies.b.e == 0.3
    @test nt.observations.astrom.jitter == 2.0
    @test length(Octofitter.make_prior_sampler(sys)(Random.Xoshiro(0))) == 5
end

@testset "interpolated Symbols are literals, not identifiers" begin
    # Regression: `$x` where x is a Symbol used to splice as an *identifier*,
    # so the generated code did a global lookup and threw an UndefVarError
    # naming a variable the user never wrote.
    key = :plx
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = getproperty(system, $key) / 100
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin a = 1.0 end)
    sys = Octofitter.System(name="s", bodies=[A, b],
        variables=@variables begin plx ~ Uniform(1, 100) end)
    @test Octofitter.make_arr2nt(sys)([50.0]).bodies.A.mass == 0.5
end

@testset "show" begin
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=(A,), variables=@variables begin a = 3.0 end)
    sys = Octofitter.System(name="s", bodies=[A, b],
        variables=@variables begin plx = 20.0 end)
    out = sprint(show, MIME"text/plain"(), sys)
    @test occursin("2 bodies, 1 orbits", out)
    @test occursin("parallax", out)
    @test occursin("Body b  about A", out)
end
