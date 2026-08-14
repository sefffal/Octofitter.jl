# Hierarchies the v1 model could not express, and the propagator seam.

# One physical configuration, spelled two ways.
function moon_system(; convention)
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 10mjup
        a = 5.0; e = 0.1; i = 0.5; ω = 0.3; Ω = 1.4; tp = 55000.0
    end)
    # A moon: astrocentric about its *host planet*, not about the star. Under
    # v1 this was unreachable — topology was inferred from semi-major axis, so
    # the moon would have been treated as an inner planet of the star.
    m = Octofitter.Body(name="m", about=b, variables=@variables begin
        mass = 0.05mjup
        a = 0.02; e = 0.0; i = 0.6; ω = 0.0; Ω = 1.4; tp = 55010.0
    end)
    c = if convention === :jacobi
        Octofitter.Body(name="c", about=(A, b, m), variables=@variables begin
            mass = 3mjup
            a = 20.0; e = 0.2; i = 0.7; ω = 2.0; Ω = 0.5; tp = 54000.0
        end)
    else
        Octofitter.Body(name="c", about=A, variables=@variables begin
            mass = 3mjup
            a = 20.0; e = 0.2; i = 0.7; ω = 2.0; Ω = 0.5; tp = 54000.0
        end)
    end
    return Octofitter.System(name="moon", bodies=[A, b, m, c],
        variables=@variables begin plx = 40.0 end)
end

@testset "a moon is expressible, and mixed conventions are legal" begin
    for conv in (:jacobi, :astrocentric)
        sys = moon_system(convention=conv)
        @test Octofitter.nbodies(sys) == 4
        @test Octofitter.nrows(sys) == 3
        posys = Octofitter.construct_system((; system=sys),
            Octofitter.make_arr2nt(sys)(Float64[]))
        traj = orbitsolve(posys, [56000.0, 56500.0])
        # The moon stays bound to its host: separation from b is ~0.02 AU,
        # nowhere near its 5 AU distance from the star.
        sep_mb = hypot(posx(traj[1], :m, :b), posy(traj[1], :m, :b), posz(traj[1], :m, :b))
        sep_mA = hypot(posx(traj[1], :m, :A), posy(traj[1], :m, :A), posz(traj[1], :m, :A))
        @test 0.015 < sep_mb < 0.025
        @test 4.0 < sep_mA < 6.0
    end

    # Jacobi and astrocentric are *different models* under KeplerianApprox —
    # the rows are the approximation, not a relabelling.
    j = Octofitter.construct_system((; system=moon_system(convention=:jacobi)),
        Octofitter.make_arr2nt(moon_system(convention=:jacobi))(Float64[]))
    a = Octofitter.construct_system((; system=moon_system(convention=:astrocentric)),
        Octofitter.make_arr2nt(moon_system(convention=:astrocentric))(Float64[]))
    tj = orbitsolve(j, [56000.0])
    ta = orbitsolve(a, [56000.0])
    @test raoff(tj[1], :A, PlanetOrbits.barycentre(j)) !=
          raoff(ta[1], :A, PlanetOrbits.barycentre(a))
end

@testset "2+2 quadruple via an explicit Orbit node" begin
    Aa = Octofitter.Body(name="Aa", variables=@variables begin mass = 1.1 end)
    Ab = Octofitter.Body(name="Ab", about=Aa, variables=@variables begin
        mass = 0.9
        a = 0.5; e = 0.1; i = 0.4; ω = 0.2; Ω = 1.0; tp = 55000.0
    end)
    Ba = Octofitter.Body(name="Ba", variables=@variables begin mass = 0.8 end)
    Bb = Octofitter.Body(name="Bb", about=Ba, variables=@variables begin
        mass = 0.7
        a = 0.4; e = 0.05; i = 1.2; ω = 1.1; Ω = 2.0; tp = 55100.0
    end)
    # The wide orbit binds two *pairs*: the case `Body(about=…)` cannot spell,
    # and the reason the explicit `Orbit` node exists. Its name gives it its
    # own chain columns.
    wide = Octofitter.Orbit(name="AB", exterior=(Ba, Bb), about=(Aa, Ab),
        variables=@variables begin
            P ~ Uniform(1e5, 1e6)      # days — the size group's other spelling
            e = 0.3; i = 0.8; ω = 0.1; Ω = 0.6; tp = 50000.0
        end)
    # `Ba` would otherwise be a second root; the wide row is what places it.
    sys = Octofitter.System(name="quad", bodies=[Aa, Ab, Ba, Bb, wide],
        variables=@variables begin plx = 10.0 end)
    @test Octofitter.nbodies(sys) == 4
    @test Octofitter.nrows(sys) == 3
    @test :AB_P in Octofitter.list_parameter_names(sys)

    nt = Octofitter.make_arr2nt(sys)([4e5])
    posys = Octofitter.construct_system((; system=sys), nt)
    # P → a round-trips through the row's gravitating mass (all four bodies).
    @test PlanetOrbits.period(posys, 3) ≈ 4e5 rtol = 1e-10
    traj = orbitsolve(posys, [56000.0])
    @test isfinite(raoff(traj[1], :Ba, :Aa))
end

@testset "AHL21 and KeplerianApprox agree for two bodies" begin
    mkmodel(method) = Octofitter.System(name="two",
        bodies=[Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end),
                Octofitter.Body(name="b", about=:A, variables=@variables begin
                    mass = 30mjup
                    a = 4.0; e = 0.25; i = 0.9; ω = 1.1; Ω = 2.2; tp = 56000.0
                end)],
        variables=(@variables begin plx = 25.0 end),
        method=method)
    epochs = collect(range(56200.0, 57800.0, length=11))
    kep = orbitsolve(Octofitter.construct_system(
            (; system=mkmodel(PlanetOrbits.KeplerianApprox())),
            Octofitter.make_arr2nt(mkmodel(PlanetOrbits.KeplerianApprox()))(Float64[])),
        epochs)
    nb = orbitsolve(Octofitter.construct_system(
            (; system=mkmodel(PlanetOrbits.AHL21(h=20.0, t0=56000.0))),
            Octofitter.make_arr2nt(mkmodel(PlanetOrbits.AHL21(h=20.0, t0=56000.0)))(Float64[])),
        epochs; method=PlanetOrbits.AHL21(h=20.0, t0=56000.0))
    for k in eachindex(epochs)
        @test raoff(kep[k], :b, :A) ≈ raoff(nb[k], :b, :A) rtol = 1e-8
        @test decoff(kep[k], :b, :A) ≈ decoff(nb[k], :b, :A) rtol = 1e-8
    end
end

@testset "the propagator is a model-level choice, not a likelihood's business" begin
    # The same observation object, evaluated under both propagators: nothing
    # in the likelihood knows which one ran.
    epochs = collect(range(56000.0, 58000.0, length=8))
    mk(method) = begin
        A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass = 20mjup
            a = 3.0; e = 0.1; i = 0.8; ω = 0.5; Ω = 1.0; tp = 55500.0
        end)
        obs = RelAstromObs((epoch=epochs, ra=zeros(length(epochs)), dec=zeros(length(epochs)),
                σ_ra=ones(length(epochs)), σ_dec=ones(length(epochs)));
            target=b, ref=A, name="astrom")
        Octofitter.System(name="p", bodies=[A, b], observations=[obs],
            variables=(@variables begin plx = 30.0 end), method=method)
    end
    s1 = mk(PlanetOrbits.KeplerianApprox())
    s2 = mk(PlanetOrbits.AHL21(h=25.0, t0=55500.0))
    l1 = Octofitter.make_ln_like(s1)(s1, Octofitter.make_arr2nt(s1)(Float64[]))
    l2 = Octofitter.make_ln_like(s2)(s2, Octofitter.make_arr2nt(s2)(Float64[]))
    @test isfinite(l1) && isfinite(l2)
    @test l1 ≈ l2 rtol = 1e-6      # exact two-body: the propagators must agree
end
