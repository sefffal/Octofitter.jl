@testset "NSS Catalog Integration" begin

    @testset "TI to Campbell conversion" begin
        # Face-on prograde orbit: A=10, B=0, F=0, G=10
        i, Ω, ω, α = Octofitter._ti_to_campbell(10.0, 0.0, 0.0, 10.0)
        @test α ≈ 10.0
        @test i ≈ 0.0 atol=1e-10  # Face-on

        # Face-on retrograde orbit: A=10, B=0, F=0, G=-10
        i2, _, _, α2 = Octofitter._ti_to_campbell(10.0, 0.0, 0.0, -10.0)
        @test α2 ≈ 10.0
        @test i2 ≈ π atol=1e-10  # Retrograde face-on

        # Round-trip: Campbell → TI → Campbell
        # Create TI constants from known Campbell elements
        a_mas = 50.0
        i_true = 1.2
        Ω_true = 0.8
        ω_true = 2.5
        A = a_mas * (cos(Ω_true)*cos(ω_true) - sin(Ω_true)*sin(ω_true)*cos(i_true))
        B = a_mas * (sin(Ω_true)*cos(ω_true) + cos(Ω_true)*sin(ω_true)*cos(i_true))
        F = a_mas * (-cos(Ω_true)*sin(ω_true) - sin(Ω_true)*cos(ω_true)*cos(i_true))
        G = a_mas * (-sin(Ω_true)*sin(ω_true) + cos(Ω_true)*cos(ω_true)*cos(i_true))

        i_rec, Ω_rec, ω_rec, α_rec = Octofitter._ti_to_campbell(A, B, F, G)
        @test α_rec ≈ a_mas atol=1e-8
        @test i_rec ≈ i_true atol=1e-8
        # (Ω, ω) and (Ω+π, ω+π) are degenerate: same orbit, same TI constants.
        # Check that ω+Ω is preserved (mod 2π), which pins the sky-plane orientation.
        @test cos(Ω_rec + ω_rec) ≈ cos(Ω_true + ω_true) atol=1e-8
        @test sin(Ω_rec + ω_rec) ≈ sin(Ω_true + ω_true) atol=1e-8
        @test cos(Ω_rec - ω_rec) ≈ cos(Ω_true - ω_true) atol=1e-8
        @test sin(Ω_rec - ω_rec) ≈ sin(Ω_true - ω_true) atol=1e-8
    end

    @testset "nss_to_starting_point with ThieleInnes model" begin
        planet_b = Planet(
            name="b",
            basis=ThieleInnesOrbit,
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                A ~ Normal(0, 1000)
                B ~ Normal(0, 1000)
                F ~ Normal(0, 1000)
                G ~ Normal(0, 1000)
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, A, B, F, G, plx=system.plx)
            end
        )
        sys = System(
            name="NSSTest_TI",
            companions=[planet_b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )
        model = Octofitter.LogDensityModel(sys)

        fake_nss = (;
            eccentricity = 0.15,
            period = 1826.25,
            t_periastron = 100.0,
            a_thiele_innes = 5.2,
            b_thiele_innes = 3.1,
            f_thiele_innes = -2.0,
            g_thiele_innes = 4.5,
        )

        result = Octofitter.nss_to_starting_point(fake_nss, model; planet_key=:b)
        @test hasproperty(result, :planets)
        @test hasproperty(result.planets, :b)
        planet_vals = result.planets.b
        @test planet_vals.A ≈ 5.2
        @test planet_vals.B ≈ 3.1
        @test planet_vals.F ≈ -2.0
        @test planet_vals.G ≈ 4.5
        @test planet_vals.e ≈ 0.15
    end

    @testset "nss_to_starting_point with Campbell model (UniformCircular)" begin
        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                a ~ LogUniform(1, 100)
                i ~ Sine()
                Ω ~ UniformCircular()
                ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, a, i, Ω, ω)
            end
        )
        sys = System(
            name="NSSTest_Campbell",
            companions=[planet_b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )
        model = Octofitter.LogDensityModel(sys)

        fake_nss = (;
            eccentricity = 0.2,
            period = 3652.5,
            t_periastron = 200.0,
            a_thiele_innes = -400.0,
            b_thiele_innes = -100.0,
            f_thiele_innes = 50.0,
            g_thiele_innes = -450.0,
        )

        result = Octofitter.nss_to_starting_point(fake_nss, model; planet_key=:b)
        planet_vals = result.planets.b

        # Should have Ωx, Ωy, ωx, ωy (not Ω, ω) due to UniformCircular
        @test hasproperty(planet_vals, :Ωx)
        @test hasproperty(planet_vals, :Ωy)
        @test hasproperty(planet_vals, :ωx)
        @test hasproperty(planet_vals, :ωy)
        @test !hasproperty(planet_vals, :Ω)  # Should NOT have raw Ω
        @test !hasproperty(planet_vals, :ω)

        # Eccentricity should be set
        @test planet_vals.e ≈ 0.2

        # Inclination should be set directly (Sine is not UniformCircular)
        @test hasproperty(planet_vals, :i)
        @test 0 < planet_vals.i < π

        # Semi-major axis should be set from period
        @test hasproperty(planet_vals, :a)
        @test planet_vals.a > 0

        # x,y components should be unit-ish (cos/sin of angle)
        @test planet_vals.Ωx^2 + planet_vals.Ωy^2 ≈ 1.0 atol=1e-10
        @test planet_vals.ωx^2 + planet_vals.ωy^2 ≈ 1.0 atol=1e-10
    end

    @testset "nss_to_starting_point with Period model" begin
        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                P ~ LogUniform(100, 10000)
                i ~ Sine()
                Ω ~ UniformCircular()
                ω ~ UniformCircular()
                M = system.M
                a = ∛(M * (P/PlanetOrbits.year2day_julian)^2)
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000.0; M, e, a, i, Ω, ω)
            end
        )
        sys = System(
            name="NSSTest_Period",
            companions=[planet_b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )
        model = Octofitter.LogDensityModel(sys)

        fake_nss = (;
            eccentricity = 0.1,
            period = 1826.25,
            t_periastron = 50.0,
            a_thiele_innes = 10.0,
            b_thiele_innes = 5.0,
            f_thiele_innes = -3.0,
            g_thiele_innes = 8.0,
        )

        result = Octofitter.nss_to_starting_point(fake_nss, model; planet_key=:b)
        planet_vals = result.planets.b

        # Should map P (free) not a (derived)
        @test hasproperty(planet_vals, :P)
        @test !hasproperty(planet_vals, :a)
        @test planet_vals.P ≈ 1826.25
    end

    @testset "nss_to_starting_point with missing NSS fields" begin
        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                a ~ LogUniform(1, 100)
                i ~ Sine()
                Ω ~ UniformCircular()
                ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, a, i, Ω, ω)
            end
        )
        sys = System(
            name="NSSTest_Missing",
            companions=[planet_b],
            observations=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )
        model = Octofitter.LogDensityModel(sys)

        # NSS solution with only eccentricity and period (no TI constants)
        sparse_nss = (;
            eccentricity = 0.3,
            period = 2000.0,
        )

        result = Octofitter.nss_to_starting_point(sparse_nss, model; planet_key=:b)
        planet_vals = result.planets.b
        @test planet_vals.e ≈ 0.3
        @test hasproperty(planet_vals, :a)  # Should still get a from period
        # Should NOT have angular elements since no TI constants
        @test !hasproperty(planet_vals, :i)
    end

    @testset "nss_to_model_chain" begin
        fake_nss = (;
            eccentricity = 0.15,
            eccentricity_error = 0.02,
            period = 1826.25,
            period_error = 10.0,
            t_periastron = 100.0,
            t_periastron_error = 5.0,
            a_thiele_innes = 5.2,
            a_thiele_innes_error = 0.3,
            b_thiele_innes = 3.1,
            b_thiele_innes_error = 0.2,
            f_thiele_innes = -2.0,
            f_thiele_innes_error = 0.15,
            g_thiele_innes = 4.5,
            g_thiele_innes_error = 0.25,
        )

        nss_model, nss_chain = Octofitter.nss_to_model_chain(fake_nss; plx=50.0, N=500)

        # Check model structure
        @test nss_model isa Octofitter.LogDensityModel
        @test nss_model.system.name == :NSS
        @test length(nss_model.system.planets) == 1

        # Check chain structure
        @test size(nss_chain, 1) == 500
        @test "b_e" in string.(names(nss_chain))
        @test "b_A" in string.(names(nss_chain))
        @test "b_tp" in string.(names(nss_chain))

        # Check values are centered on NSS solution
        using Statistics
        @test median(nss_chain["b_e"][:]) ≈ 0.15 atol=0.01
        @test median(nss_chain["b_A"][:]) ≈ 5.2 atol=0.2

        # Check orbits can be constructed (needed for plotting)
        orbits = Octofitter.construct_elements(nss_model, nss_chain, :b, 1:5)
        @test length(orbits) == 5
        @test first(orbits) isa ThieleInnesOrbit
    end

    @testset "_nss_get_float" begin
        nt = (; a=1.5, b=missing, c=nothing, d=NaN, e=Inf)
        @test Octofitter._nss_get_float(nt, :a) ≈ 1.5
        @test isnothing(Octofitter._nss_get_float(nt, :b))
        @test isnothing(Octofitter._nss_get_float(nt, :c))
        @test isnothing(Octofitter._nss_get_float(nt, :d))
        @test isnothing(Octofitter._nss_get_float(nt, :e))
        @test isnothing(Octofitter._nss_get_float(nt, :nonexistent))
    end
end
