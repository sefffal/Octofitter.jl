using Test
using PlanetOrbits

@testset "Bulk Orbit Solver" begin

    # Test orbit: moderate eccentricity, non-trivial angles
    orbit = KepOrbit(
        a = 1.0,   # AU
        e = 0.3,
        i = deg2rad(45.0),
        ω = deg2rad(30.0),
        Ω = deg2rad(60.0),
        tp = 0.0,
        M = 1.0,   # M⊙
    )

    # Test epochs spanning ~2 orbits
    times = collect(range(0.0, stop=800.0, length=50))

    @testset "kepler_solver_reactant matches Markley" begin
        # Grid of (M, e) values
        M_vals = range(-π, π, length=100)
        e_vals = [0.0, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999]
        for e in e_vals
            for M in M_vals
                EA_ref = PlanetOrbits.kepler_solver(M, e, PlanetOrbits.Markley())
                EA_new = PlanetOrbits.kepler_solver_reactant(M, e)
                @test EA_new ≈ EA_ref atol=1e-12
            end
        end
    end

    @testset "kepler_solver_reactant edge cases" begin
        # M = 0, e = 0
        @test PlanetOrbits.kepler_solver_reactant(0.0, 0.0) ≈ 0.0 atol=1e-15
        # M = 0, e > 0
        @test PlanetOrbits.kepler_solver_reactant(0.0, 0.5) ≈ 0.0 atol=1e-15
        # M > 0, e = 0
        @test PlanetOrbits.kepler_solver_reactant(1.0, 0.0) ≈ 1.0 atol=1e-12
        # Large M values (should wrap correctly)
        @test PlanetOrbits.kepler_solver_reactant(10π + 0.5, 0.3) ≈
              PlanetOrbits.kepler_solver(10π + 0.5, 0.3, PlanetOrbits.Markley()) atol=1e-12
    end

    @testset "orbitsolve_bulk matches scalar orbitsolve" begin
        bulk = orbitsolve_bulk(orbit, times)

        @test bulk isa PlanetOrbits.OrbitSolutionBulk
        @test length(bulk.ν) == length(times)

        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test bulk.ν[j] ≈ scalar.ν atol=1e-12
            @test bulk.sinν_ω[j] ≈ scalar.sinν_ω atol=1e-12
            @test bulk.cosν_ω[j] ≈ scalar.cosν_ω atol=1e-12
            @test bulk.ecosν[j] ≈ scalar.ecosν atol=1e-12
            @test bulk.r[j] ≈ scalar.r atol=1e-12
        end
    end

    @testset "Bulk observable: radvel" begin
        bulk = orbitsolve_bulk(orbit, times)
        rv_bulk = radvel(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test rv_bulk[j] ≈ radvel(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: posx" begin
        bulk = orbitsolve_bulk(orbit, times)
        x_bulk = posx(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test x_bulk[j] ≈ posx(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: posy" begin
        bulk = orbitsolve_bulk(orbit, times)
        y_bulk = posy(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test y_bulk[j] ≈ posy(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: posz" begin
        bulk = orbitsolve_bulk(orbit, times)
        z_bulk = posz(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test z_bulk[j] ≈ posz(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: velx" begin
        bulk = orbitsolve_bulk(orbit, times)
        vx_bulk = velx(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test vx_bulk[j] ≈ velx(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: vely" begin
        bulk = orbitsolve_bulk(orbit, times)
        vy_bulk = vely(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test vy_bulk[j] ≈ vely(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: velz" begin
        bulk = orbitsolve_bulk(orbit, times)
        vz_bulk = velz(bulk)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test vz_bulk[j] ≈ velz(scalar) atol=1e-10
        end
    end

    @testset "Bulk observable: radvel with mass ratio" begin
        bulk = orbitsolve_bulk(orbit, times)
        M_planet = 0.001  # solar masses
        rv_bulk = radvel(bulk, M_planet)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test rv_bulk[j] ≈ radvel(scalar, M_planet) atol=1e-10
        end
    end

    @testset "Bulk observable: posx with mass ratio" begin
        bulk = orbitsolve_bulk(orbit, times)
        M_planet = 0.001
        x_bulk = posx(bulk, M_planet)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orbit, t)
            @test x_bulk[j] ≈ posx(scalar, M_planet) atol=1e-10
        end
    end

    @testset "Different orbit configurations" begin
        # Near-circular orbit
        orb_circ = KepOrbit(a=2.0, e=0.01, i=deg2rad(10.0), ω=0.0, Ω=0.0, tp=100.0, M=1.5)
        bulk_circ = orbitsolve_bulk(orb_circ, times)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orb_circ, t)
            @test radvel(bulk_circ)[j] ≈ radvel(scalar) atol=1e-10
        end

        # High eccentricity
        orb_ecc = KepOrbit(a=0.5, e=0.95, i=deg2rad(80.0), ω=deg2rad(270.0), Ω=deg2rad(180.0), tp=50.0, M=0.8)
        bulk_ecc = orbitsolve_bulk(orb_ecc, times)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orb_ecc, t)
            # Use rtol for high-eccentricity orbits where RV semiamplitude is large
            @test radvel(bulk_ecc)[j] ≈ radvel(scalar) rtol=1e-12
            @test posx(bulk_ecc)[j] ≈ posx(scalar) rtol=1e-12
        end

        # Edge-on orbit
        orb_edge = KepOrbit(a=1.0, e=0.5, i=deg2rad(90.0), ω=deg2rad(45.0), Ω=deg2rad(90.0), tp=0.0, M=1.0)
        bulk_edge = orbitsolve_bulk(orb_edge, times)
        for (j, t) in enumerate(times)
            scalar = orbitsolve(orb_edge, t)
            @test radvel(bulk_edge)[j] ≈ radvel(scalar) atol=1e-10
        end
    end

    @testset "ForwardDiff through orbitsolve_bulk" begin
        using ForwardDiff

        test_times = collect(range(100.0, stop=500.0, length=10))

        # Derivative w.r.t. semi-major axis
        f_a(a) = sum(radvel(orbitsolve_bulk(
            KepOrbit(a=a, e=0.3, i=deg2rad(45.0), ω=deg2rad(30.0), Ω=deg2rad(60.0), tp=0.0, M=1.0),
            test_times
        )))
        da = ForwardDiff.derivative(f_a, 1.0)
        @test isfinite(da)

        # Derivative w.r.t. eccentricity
        f_e(e) = sum(radvel(orbitsolve_bulk(
            KepOrbit(a=1.0, e=e, i=deg2rad(45.0), ω=deg2rad(30.0), Ω=deg2rad(60.0), tp=0.0, M=1.0),
            test_times
        )))
        de = ForwardDiff.derivative(f_e, 0.3)
        @test isfinite(de)

        # Derivative w.r.t. argument of periapsis
        f_ω(ω) = sum(posx(orbitsolve_bulk(
            KepOrbit(a=1.0, e=0.3, i=deg2rad(45.0), ω=ω, Ω=deg2rad(60.0), tp=0.0, M=1.0),
            test_times
        )))
        dω = ForwardDiff.derivative(f_ω, deg2rad(30.0))
        @test isfinite(dω)

        # Check ForwardDiff derivative against finite difference
        h = 1e-7
        fd_da = (f_a(1.0 + h) - f_a(1.0 - h)) / (2h)
        @test da ≈ fd_da rtol=1e-5

        fd_de = (f_e(0.3 + h) - f_e(0.3 - h)) / (2h)
        @test de ≈ fd_de rtol=1e-5

        fd_dω = (f_ω(deg2rad(30.0) + h) - f_ω(deg2rad(30.0) - h)) / (2h)
        @test dω ≈ fd_dω rtol=1e-5
    end
end
