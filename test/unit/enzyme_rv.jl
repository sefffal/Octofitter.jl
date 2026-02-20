# Test Enzyme integration with MarginalizedStarAbsoluteRVObs
# Verifies that the per-term AD machinery correctly uses Enzyme
# for the marginalized RV observation type, and that gradients
# match finite differences.
#
# NOTE: Enzyme has a known issue with --check-bounds=yes (which Pkg.test()
# enables by default): the @boundscheck in arr2nt gets garbage length()
# values during Enzyme's forward pass. These tests are skipped when
# --check-bounds=yes is active.

using OctofitterRadialVelocity
using ADTypes: AutoEnzyme
using Random

# Detect --check-bounds=yes (set by Pkg.test / CI)
const _check_bounds_forced = Base.JLOptions().check_bounds == 1

@testset "ad_backend dispatch" begin
    rv_table = Table(
        epoch=[50000.0, 50100.0, 50200.0, 50300.0, 50400.0],
        rv=[100.0, 110.0, 105.0, 95.0, 108.0],
        σ_rv=[10.0, 10.0, 10.0, 10.0, 10.0],
    )
    test_obs = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV")
    if _check_bounds_forced
        @test_skip Octofitter.ad_backend(test_obs) isa AutoEnzyme
    else
        @test Octofitter.ad_backend(test_obs) isa AutoEnzyme
    end
end

@testset "Enzyme vs FiniteDiff gradients" begin
    if _check_bounds_forced
        @test_skip false  # Enzyme + --check-bounds=yes causes @boundscheck issues in arr2nt
    else
        rv_table = Table(
            epoch=[50000.0, 50100.0, 50200.0, 50300.0, 50400.0],
            rv=[100.0, 110.0, 105.0, 95.0, 108.0],
            σ_rv=[10.0, 10.0, 10.0, 10.0, 10.0],
        )
        rvlike = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV")

        b = Planet(
            name="b",
            basis=RadialVelocityOrbit,
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                ω ~ UniformCircular()
                τ ~ Uniform(0, 1)
                P ~ Uniform(0.001, 10)
                a = ∛(P^2 * system.M)
                tp = τ * P * 365.25 + 50000
                mass ~ Uniform(0, 100)
            end
        )

        sys = System(
            name="EnzymeTestSys",
            companions=[b],
            observations=[rvlike],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            end
        )

        # Build models: default (Enzyme for RV term) and FiniteDiff override
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
end
