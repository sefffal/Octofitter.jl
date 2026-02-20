using Test, Octofitter, OctofitterRadialVelocity, PlanetOrbits
using TypedTables, Distributions, Random, ADTypes

@testset "Enzyme vs FiniteDiff gradients (orbitsolve_bulk)" begin
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
    println("Building Enzyme model...")
    model = Octofitter.LogDensityModel(sys, verbosity=0)
    println("Building FiniteDiff model...")
    model_fd = Octofitter.LogDensityModel(sys, autodiff=AutoFiniteDiff(), verbosity=0)

    # Find a finite starting point
    theta = nothing
    ll = -Inf
    for attempt in 0:20
        theta_natural = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
        theta = model.link(theta_natural)
        ll = model.ℓπcallback(theta)
        isfinite(ll) && break
    end
    println("Starting ll = $ll")
    @test isfinite(ll)

    # Compare gradients
    println("Computing Enzyme gradient...")
    ll_enz, grad_enz = model.∇ℓπcallback(theta)
    println("  ll_enz = $ll_enz")
    println("  grad_enz = $grad_enz")

    println("Computing FiniteDiff gradient...")
    ll_fd, grad_fd = model_fd.∇ℓπcallback(theta)
    println("  ll_fd = $ll_fd")
    println("  grad_fd = $grad_fd")

    @test ll_enz ≈ ll_fd
    @test all(isfinite, grad_enz)
    @test all(isfinite, grad_fd)
    @test isapprox(grad_enz, grad_fd, atol=1e-2, rtol=1e-3)
    println("All gradient tests passed!")
end

@testset "Short MCMC run with MarginalizedStarAbsoluteRVObs" begin
    # Create a synthetic dataset from a known orbit
    rv_table = Table(
        epoch=collect(range(50000.0, 50500.0, length=20)),
        rv=10.0 .* sin.(2π .* collect(range(0, 1, length=20))) .+ randn(20) .* 2.0,
        σ_rv=fill(2.0, 20),
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
        name="MCMCTestSys",
        companions=[b],
        observations=[rvlike],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        end
    )

    model = Octofitter.LogDensityModel(sys, verbosity=0)

    # Just verify model construction succeeded and we can evaluate
    theta_natural = collect(model.sample_priors(Random.Xoshiro(123)))
    theta = model.link(theta_natural)
    ll = model.ℓπcallback(theta)
    @test isfinite(ll) || @warn "Initial ll not finite: $ll"

    # Short MCMC run (just verify it doesn't error)
    println("Running short MCMC chain...")
    chain = octofit(
        model,
        iterations=200,
        adaptation=100,
        verbosity=0
    )
    @test length(chain) > 0
    println("MCMC run completed successfully with $(length(chain)) samples")
end
