@testset "Gaia DR4" begin

    # Load real test data (from OHP Gaia splinter session)
    test_data_path = joinpath(@__DIR__, "..", "data", "gaia_dr4_test.csv")
    df = CSV.read(test_data_path, DataFrame)

    @testset "GaiaDR4AstromObs construction" begin
        # Test construction with real data format
        ref_epoch_mjd = 57936.375

        # Use a mock gaia_id that won't trigger network lookup in inner constructor
        # We'll use the inner constructor to avoid network dependencies
        priors, derived = @variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 10000)
            dec_offset_mas ~ Normal(0, 10000)
            pmra ~ Uniform(-1000, 1000)
            pmdec ~ Uniform(-1000, 1000)
            plx = system.plx
            ref_epoch = $ref_epoch_mjd
        end

        # Prepare table in expected format
        table = TypedTables.Table(df)
        if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
            table = TypedTables.Table(table; epoch=Octofitter.jd2mjd.(table.obs_time_tcb))
        end
        xyz = TypedTables.Table(Octofitter.geocentre_position_query.(table.epoch))
        table = TypedTables.Table(table; xyz)

        # Mock Gaia solution
        mock_gaia_sol = (
            ra = 180.0,
            dec = 45.0,
            pmra = 5.0,
            pmdec = -10.0,
            parallax = 20.0,
        )

        gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(table), typeof(mock_gaia_sol)}(
            table,
            123456789,
            mock_gaia_sol,
            priors,
            derived,
            "GaiaDR4"
        )

        @test gaia_obs isa GaiaDR4AstromObs
        @test length(gaia_obs.table) == nrow(df)
        @test hasproperty(gaia_obs.table, :epoch)
        @test hasproperty(gaia_obs.table, :xyz)
        @test hasproperty(gaia_obs.table.xyz, :x)
        @test hasproperty(gaia_obs.table.xyz, :y)
        @test hasproperty(gaia_obs.table.xyz, :z)
    end

    # Note: epoch subsetting test skipped because likeobj_from_epoch_subset
    # currently uses the public constructor which requires network access.
    # The function works correctly but can't be tested without network.

    @testset "GaiaDR4AstromObs model building and likelihood" begin
        ref_epoch_mjd = 57936.375

        priors, derived = @variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 10000)
            dec_offset_mas ~ Normal(0, 10000)
            pmra ~ Uniform(-1000, 1000)
            pmdec ~ Uniform(-1000, 1000)
            plx = system.plx
            ref_epoch = $ref_epoch_mjd
        end

        table = TypedTables.Table(df)
        if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
            table = TypedTables.Table(table; epoch=Octofitter.jd2mjd.(table.obs_time_tcb))
        end
        xyz = TypedTables.Table(Octofitter.geocentre_position_query.(table.epoch))
        table = TypedTables.Table(table; xyz)

        mock_gaia_sol = (
            ra = 180.0,
            dec = 45.0,
            pmra = 5.0,
            pmdec = -10.0,
            parallax = 20.0,
        )

        gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(table), typeof(mock_gaia_sol)}(
            table,
            123456789,
            mock_gaia_sol,
            priors,
            derived,
            "GaiaDR4"
        )

        orbit_ref_epoch = mean(gaia_obs.table.epoch)

        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ LogUniform(0.1, 50)
                e ~ Uniform(0, 0.9)
                ω ~ Uniform(0, 2π)
                i ~ Sine()
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.1, 500)
            end
        )

        sys = System(
            name="gaia_dr4_test",
            companions=[b],
            observations=[gaia_obs],
            variables=@variables begin
                M = 1.0
                plx ~ Uniform(1, 100)
            end
        )

        model = Octofitter.LogDensityModel(sys, verbosity=0)
        @test model isa Octofitter.LogDensityModel

        # Test that we can evaluate the likelihood
        # Draw a random sample from priors and evaluate
        arr = model.sample_priors(Random.default_rng())
        ll = model.ℓπcallback(arr)
        @test isfinite(ll)
    end

    @testset "GaiaDR4AstromObs simulate function" begin
        ref_epoch_mjd = 57936.375

        priors, derived = @variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 10000)
            dec_offset_mas ~ Normal(0, 10000)
            pmra ~ Uniform(-1000, 1000)
            pmdec ~ Uniform(-1000, 1000)
            plx = system.plx
            ref_epoch = $ref_epoch_mjd
        end

        table = TypedTables.Table(df)
        if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
            table = TypedTables.Table(table; epoch=Octofitter.jd2mjd.(table.obs_time_tcb))
        end
        xyz = TypedTables.Table(Octofitter.geocentre_position_query.(table.epoch))
        table = TypedTables.Table(table; xyz)

        mock_gaia_sol = (
            ra = 180.0,
            dec = 45.0,
            pmra = 5.0,
            pmdec = -10.0,
            parallax = 20.0,
        )

        gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(table), typeof(mock_gaia_sol)}(
            table,
            123456789,
            mock_gaia_sol,
            priors,
            derived,
            "GaiaDR4"
        )

        orbit_ref_epoch = mean(gaia_obs.table.epoch)

        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ LogUniform(0.1, 50)
                e ~ Uniform(0, 0.9)
                ω ~ Uniform(0, 2π)
                i ~ Sine()
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.1, 500)
            end
        )

        sys = System(
            name="gaia_dr4_test",
            companions=[b],
            observations=[gaia_obs],
            variables=@variables begin
                M = 1.0
                plx ~ Uniform(1, 100)
            end
        )

        model = Octofitter.LogDensityModel(sys, verbosity=0)

        # Test generate_from_params
        params = Octofitter.drawfrompriors(model.system)

        # Test without noise
        sim_system = Octofitter.generate_from_params(model.system, params; add_noise=false)
        @test sim_system isa System
        sim_model = Octofitter.LogDensityModel(sim_system, verbosity=0)
        @test sim_model isa Octofitter.LogDensityModel

        # Test with noise
        sim_system_noisy = Octofitter.generate_from_params(model.system, params; add_noise=true)
        @test sim_system_noisy isa System

        # Verify simulated data is different when noise is added
        sim_obs = first(filter(o -> o isa GaiaDR4AstromObs, sim_system.observations))
        sim_obs_noisy = first(filter(o -> o isa GaiaDR4AstromObs, sim_system_noisy.observations))
        @test sim_obs.table.centroid_pos_al != sim_obs_noisy.table.centroid_pos_al
    end

    @testset "GaiaDR4AstromObs simulate returns NamedTuple" begin
        # This test verifies the fix for the broadcasting error
        # The simulate function should return a NamedTuple with specific fields

        ref_epoch_mjd = 57936.375

        priors, derived = @variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 10000)
            dec_offset_mas ~ Normal(0, 10000)
            pmra ~ Uniform(-1000, 1000)
            pmdec ~ Uniform(-1000, 1000)
            plx = system.plx
            ref_epoch = $ref_epoch_mjd
        end

        table = TypedTables.Table(df)
        if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
            table = TypedTables.Table(table; epoch=Octofitter.jd2mjd.(table.obs_time_tcb))
        end
        xyz = TypedTables.Table(Octofitter.geocentre_position_query.(table.epoch))
        table = TypedTables.Table(table; xyz)

        mock_gaia_sol = (
            ra = 180.0,
            dec = 45.0,
            pmra = 5.0,
            pmdec = -10.0,
            parallax = 20.0,
        )

        gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(table), typeof(mock_gaia_sol)}(
            table,
            123456789,
            mock_gaia_sol,
            priors,
            derived,
            "GaiaDR4"
        )

        orbit_ref_epoch = mean(gaia_obs.table.epoch)

        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ LogUniform(0.1, 50)
                e ~ Uniform(0, 0.9)
                ω ~ Uniform(0, 2π)
                i ~ Sine()
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.1, 500)
            end
        )

        sys = System(
            name="gaia_dr4_test",
            companions=[b],
            observations=[gaia_obs],
            variables=@variables begin
                M = 1.0
                plx ~ Uniform(1, 100)
            end
        )

        model = Octofitter.LogDensityModel(sys, verbosity=0)
        params = Octofitter.drawfrompriors(model.system)

        # Get the observation variables
        θ_obs = params.observations.GaiaDR4

        # Construct orbits
        orbits = map(keys(model.system.planets)) do planet_key
            planet_params = params.planets[planet_key]
            PlanetOrbits.Visual{PlanetOrbits.KepOrbit}(;
                a = planet_params.a,
                e = planet_params.e,
                ω = planet_params.ω,
                i = planet_params.i,
                Ω = planet_params.Ω,
                tp = planet_params.tp,
                M = params.M,
                plx = params.plx,
                ra = mock_gaia_sol.ra,
                dec = mock_gaia_sol.dec,
                pmra = θ_obs.pmra,
                pmdec = θ_obs.pmdec,
                rv = 0.0,
            )
        end

        orbit_solutions = map(orbits) do orbit
            PlanetOrbits.orbitsolve.(orbit, gaia_obs.table.epoch)
        end

        # Call simulate directly
        sim_result = Octofitter.simulate(
            gaia_obs,
            params,
            θ_obs,
            orbits,
            orbit_solutions,
            0
        )

        # Verify it returns a NamedTuple with the expected fields
        @test sim_result isa NamedTuple
        @test haskey(sim_result, :along_scan_residuals_buffer)
        @test haskey(sim_result, :ra_offset_buffer)
        @test haskey(sim_result, :dec_offset_buffer)

        # Verify the arrays have the correct length
        @test length(sim_result.along_scan_residuals_buffer) == length(gaia_obs.table.epoch)
        @test length(sim_result.ra_offset_buffer) == length(gaia_obs.table.epoch)
        @test length(sim_result.dec_offset_buffer) == length(gaia_obs.table.epoch)
    end
end
