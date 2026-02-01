@testset "Likelihood Objects" begin
    @testset "PlanetRelAstromLikelihood" begin
        # Test RA/Dec format
        data_radec = PlanetRelAstromLikelihood(
            Table(epoch=[5000.0, 5100.0], ra=[100.0, 110.0], dec=[50.0, 55.0], σ_ra=[1.0, 1.0], σ_dec=[1.0, 1.0]),
            name="test_radec"
        )
        @test data_radec isa PlanetRelAstromLikelihood
        @test length(data_radec.table) == 2
        @test hasproperty(data_radec.table, :ra)

        # Test sep/PA format
        data_seppa = PlanetRelAstromLikelihood(
            Table(epoch=[5000.0, 5100.0], sep=[100.0, 110.0], pa=[1.0, 1.1], σ_sep=[1.0, 1.0], σ_pa=[0.1, 0.1]),
            name="test_seppa"
        )
        @test data_seppa isa PlanetRelAstromLikelihood
        @test length(data_seppa.table) == 2
        @test hasproperty(data_seppa.table, :sep)

        # Test that invalid column combinations throw errors
        @test_throws Exception PlanetRelAstromLikelihood(
            Table(epoch=[5000.0], ra=[100.0], pa=[1.0], σ_ra=[1.0], σ_pa=[0.1]),
            name="test_invalid"
        )

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(data_seppa, 1:1)
        @test length(subset.table) == 1
    end

    @testset "PhotometryLikelihood" begin
        phot = PhotometryLikelihood(
            Table(band=[:Z, :J], phot=[15.0, 14.0], σ_phot=[0.1, 0.2]),
            name="test_phot"
        )
        @test phot isa PhotometryLikelihood
        @test length(phot.table) == 2
        @test all([:band, :phot, :σ_phot] .∈ Ref(propertynames(phot.table)))

        subset = Octofitter.likeobj_from_epoch_subset(phot, 1:1)
        @test length(subset.table) == 1
    end

    @testset "HGCALikelihood" begin
        gaia_id = 756291174721509376
        hgca = HGCALikelihood(;gaia_id=gaia_id)
        @test hgca isa HGCALikelihood
    end

    @testset "GaiaDR4AstromObs" begin
        # Create mock observation data (avoiding network dependencies)
        N_epochs = 10
        mock_epochs = collect(range(58000.0, 59000.0, length=N_epochs))
        mock_scan_angles = range(0.0, 2π, length=N_epochs)

        # Mock Earth position data (simplified)
        mock_xyz = Table(
            x = fill(0.5, N_epochs),
            y = fill(0.5, N_epochs),
            z = fill(0.2, N_epochs),
            vx = fill(0.0, N_epochs),
            vy = fill(0.0, N_epochs),
            vz = fill(0.0, N_epochs),
        )

        mock_table = Table(
            epoch = mock_epochs,
            scan_pos_angle = collect(mock_scan_angles),
            centroid_pos_al = zeros(N_epochs),
            centroid_pos_error_al = fill(0.1, N_epochs),
            parallax_factor_al = fill(0.5, N_epochs),
            outlier_flag = fill(false, N_epochs),
            xyz = collect(eachrow(mock_xyz)),
        )

        # Mock Gaia solution (minimal fields needed by simulate)
        mock_gaia_sol = (
            ra = 180.0,
            dec = 45.0,
            pmra = 5.0,
            pmdec = -10.0,
            parallax = 20.0,
        )

        # Create observation object using inner constructor (bypasses network)
        ref_epoch_mjd = 57936.375  # DR3 reference epoch
        priors, derived = @variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 100)
            dec_offset_mas ~ Normal(0, 100)
            pmra ~ Normal(5, 10)
            pmdec ~ Normal(-10, 10)
            plx = system.plx
            ref_epoch = $ref_epoch_mjd
        end

        gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(mock_table), typeof(mock_gaia_sol)}(
            mock_table,
            123456789,  # mock gaia_id
            mock_gaia_sol,
            priors,
            derived,
            "GaiaDR4"
        )

        @test gaia_obs isa GaiaDR4AstromObs
        @test length(gaia_obs.table) == N_epochs

        # Test that we can build a full model with this observation
        # This tests the ln_like function compiles correctly
        orbit_ref_epoch = mean(mock_epochs)
        b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=[],
            variables=@variables begin
                a ~ LogUniform(0.5, 20)
                e ~ Uniform(0, 0.5)
                ω ~ Uniform(0, 2π)
                i ~ Sine()
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(1, 100)
            end
        )

        sys = System(
            name="gaia_test",
            companions=[b],
            observations=[gaia_obs],
            variables=@variables begin
                M = 1.0
                plx ~ Uniform(10, 30)
            end
        )

        # This would have failed before the fix due to NamedTuple/Vector mismatch
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        @test model isa Octofitter.LogDensityModel

        # Test that we can draw from priors
        params = Octofitter.drawfrompriors(model.system)
        @test haskey(params, :plx)
        @test haskey(params.planets, :b)

        # Test generate_from_params - this also would have failed before the fix
        sim_system = Octofitter.generate_from_params(model.system, params; add_noise=false)
        @test sim_system isa System

        # Verify the simulated system can also create a valid model
        sim_model = Octofitter.LogDensityModel(sim_system, verbosity=0)
        @test sim_model isa Octofitter.LogDensityModel
    end
end
