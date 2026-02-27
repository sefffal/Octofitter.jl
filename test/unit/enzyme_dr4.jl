# Test Enzyme integration with GaiaDR4AstromObs
# Verifies that the per-term AD machinery correctly uses Enzyme
# for the DR4 astrometry observation type, and that gradients
# match finite differences.

using ADTypes: AutoEnzyme, AutoFiniteDiff
using Random
using Statistics

@testset "GaiaDR4AstromObs Enzyme" begin
    N_epochs = 10
    mock_epochs = collect(range(58000.0, 59000.0, length=N_epochs))
    mock_scan_angles = collect(range(0.0, 2π, length=N_epochs))

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

    mock_gaia_sol = (
        ra = 180.0,
        dec = 45.0,
        pmra = 5.0,
        pmdec = -10.0,
        parallax = 20.0,
    )

    ref_epoch_mjd = 57936.375
    priors, derived = @variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas ~ Normal(0, 100)
        dec_offset_mas ~ Normal(0, 100)
        pmra ~ Normal(5, 10)
        pmdec ~ Normal(-10, 10)
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end

    mean_epoch = sum(mock_epochs) / length(mock_epochs)
    detrend_Δt = collect((mock_epochs .- mean_epoch) ./ 365.25)
    detrend_inv_N = 1.0 / length(mock_epochs)
    detrend_inv_sum_Δt² = 1.0 / sum(detrend_Δt .^ 2)

    gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(mock_table), typeof(mock_gaia_sol)}(
        mock_table,
        123456789,
        mock_gaia_sol,
        priors,
        derived,
        "GaiaDR4",
        false,
        detrend_Δt,
        detrend_inv_N,
        detrend_inv_sum_Δt²,
    )

    orbit_ref_epoch = mean_epoch
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
        name="DR4EnzymeTest",
        companions=[b],
        observations=[gaia_obs],
        variables=@variables begin
            M = 1.0
            plx ~ Uniform(10, 30)
        end
    )

    @testset "ad_backend dispatch" begin
        @test Octofitter.ad_backend(gaia_obs) isa AutoEnzyme
    end

    @testset "Enzyme vs FiniteDiff gradients" begin
        # Build models: default (Enzyme for DR4 term) and FiniteDiff override
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        model_fd = Octofitter.LogDensityModel(sys, autodiff=AutoFiniteDiff(), verbosity=0)

        # Find a finite starting point
        rng = Random.Xoshiro(42)
        theta = nothing
        ll = -Inf
        for attempt in 0:20
            theta_natural = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
            theta = model.link(theta_natural)
            ll = model.ℓπcallback(theta)
            isfinite(ll) && break
        end
        @test isfinite(ll)

        # Compare gradients
        ll_enz, grad_enz = model.∇ℓπcallback(theta)
        ll_fd, grad_fd = model_fd.∇ℓπcallback(theta)

        @test ll_enz ≈ ll_fd
        @test all(isfinite, grad_enz)
        @test all(isfinite, grad_fd)
        @test isapprox(grad_enz, grad_fd, atol=1e-2, rtol=1e-3)
    end

    @testset "primary_star_perturbation mode" begin
        gaia_obs_psp = Octofitter.GaiaDR4AstromObs{typeof(mock_table), typeof(mock_gaia_sol)}(
            mock_table,
            123456789,
            mock_gaia_sol,
            priors,
            derived,
            "GaiaDR4",
            true,  # primary_star_perturbation=true
            detrend_Δt,
            detrend_inv_N,
            detrend_inv_sum_Δt²,
        )

        sys_psp = System(
            name="DR4EnzymeTestPSP",
            companions=[b],
            observations=[gaia_obs_psp],
            variables=@variables begin
                M = 1.0
                plx ~ Uniform(10, 30)
            end
        )

        model_psp = Octofitter.LogDensityModel(sys_psp, verbosity=0)
        model_psp_fd = Octofitter.LogDensityModel(sys_psp, autodiff=AutoFiniteDiff(), verbosity=0)

        theta_psp = nothing
        ll_psp = -Inf
        for attempt in 0:20
            theta_natural = collect(model_psp.sample_priors(Random.Xoshiro(42 + attempt)))
            theta_psp = model_psp.link(theta_natural)
            ll_psp = model_psp.ℓπcallback(theta_psp)
            isfinite(ll_psp) && break
        end
        @test isfinite(ll_psp)

        ll_enz_psp, grad_enz_psp = model_psp.∇ℓπcallback(theta_psp)
        ll_fd_psp, grad_fd_psp = model_psp_fd.∇ℓπcallback(theta_psp)

        @test ll_enz_psp ≈ ll_fd_psp
        @test all(isfinite, grad_enz_psp)
        @test isapprox(grad_enz_psp, grad_fd_psp, atol=1e-2, rtol=1e-3)
    end
end
