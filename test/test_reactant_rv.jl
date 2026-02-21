using Test, Octofitter, OctofitterRadialVelocity, PlanetOrbits
using TypedTables, Distributions, Random, ADTypes
using Reactant

@testset "Reactant vs Enzyme gradients (MarginalizedStarAbsoluteRVObs)" begin
    rv_table = Table(
        epoch=[50000.0, 50100.0, 50200.0, 50300.0, 50400.0],
        rv=[100.0, 110.0, 105.0, 95.0, 108.0],
        σ_rv=[10.0, 10.0, 10.0, 10.0, 10.0],
    )

    # Reactant-accelerated version
    rvlike_reactant = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV", device=Reactant.CPU())

    # Enzyme baseline (device=nothing)
    rvlike_enzyme = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV")

    b_vars = @variables begin
        e ~ Uniform(0.0, 0.5)
        ω ~ UniformCircular()
        τ ~ Uniform(0, 1)
        P ~ Uniform(0.001, 10)
        a = ∛(P^2 * system.M)
        tp = τ * P * 365.25 + 50000
        mass ~ Uniform(0, 100)
    end

    sys_vars = @variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    end

    # Build Reactant model
    b_reactant = Planet(
        name="b",
        basis=RadialVelocityOrbit,
        observations=[],
        variables=b_vars,
    )
    sys_reactant = System(
        name="ReactantTestSys",
        companions=[b_reactant],
        observations=[rvlike_reactant],
        variables=sys_vars,
    )
    println("Building Reactant model...")
    model_reactant = Octofitter.LogDensityModel(sys_reactant, verbosity=0)

    # Build Enzyme reference model
    b_enzyme = Planet(
        name="b",
        basis=RadialVelocityOrbit,
        observations=[],
        variables=b_vars,
    )
    sys_enzyme = System(
        name="EnzymeTestSys",
        companions=[b_enzyme],
        observations=[rvlike_enzyme],
        variables=sys_vars,
    )
    println("Building Enzyme model...")
    model_enzyme = Octofitter.LogDensityModel(sys_enzyme, verbosity=0)

    # Find a finite starting point
    theta = nothing
    ll = -Inf
    for attempt in 0:20
        theta_natural = collect(model_enzyme.sample_priors(Random.Xoshiro(42 + attempt)))
        theta = model_enzyme.link(theta_natural)
        ll = model_enzyme.ℓπcallback(theta)
        isfinite(ll) && break
    end
    println("Starting ll = $ll")
    @test isfinite(ll)

    # Compare primal values
    println("Computing primal values...")
    ll_reactant = model_reactant.ℓπcallback(theta)
    ll_enzyme = model_enzyme.ℓπcallback(theta)
    println("  ll_reactant = $ll_reactant")
    println("  ll_enzyme   = $ll_enzyme")
    @test isapprox(ll_reactant, ll_enzyme, rtol=1e-10)

    # Compare gradients
    println("Computing Reactant gradient...")
    ll_r, grad_r = model_reactant.∇ℓπcallback(theta)
    println("  ll_r = $ll_r, grad_r = $grad_r")

    println("Computing Enzyme gradient...")
    ll_e, grad_e = model_enzyme.∇ℓπcallback(theta)
    println("  ll_e = $ll_e, grad_e = $grad_e")

    @test isapprox(ll_r, ll_e, rtol=1e-8)
    @test all(isfinite, grad_r)
    @test all(isfinite, grad_e)
    @test isapprox(grad_r, grad_e, rtol=1e-6)
    println("All Reactant vs Enzyme gradient tests passed!")
end

@testset "Reactant compiled thunk reuse with different θ" begin
    rv_table = Table(
        epoch=[50000.0, 50100.0, 50200.0, 50300.0, 50400.0],
        rv=[100.0, 110.0, 105.0, 95.0, 108.0],
        σ_rv=[10.0, 10.0, 10.0, 10.0, 10.0],
    )
    rvlike = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV", device=Reactant.CPU())

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
        name="ReuseTestSys",
        companions=[b],
        observations=[rvlike],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        end
    )
    model = Octofitter.LogDensityModel(sys, verbosity=0)

    # Evaluate gradient with multiple different θ values
    rng = Random.Xoshiro(42)
    results = []
    for i in 1:5
        theta_natural = collect(model.sample_priors(rng))
        theta = model.link(theta_natural)
        ll = model.ℓπcallback(theta)
        if isfinite(ll)
            ll_g, grad = model.∇ℓπcallback(theta)
            push!(results, (ll=ll_g, grad=copy(grad), theta=copy(theta)))
        end
    end

    # Verify we got at least a few finite results
    @test length(results) >= 2

    # Verify different θ values produce different results (thunk not stuck)
    if length(results) >= 2
        @test results[1].ll != results[2].ll || results[1].theta == results[2].theta
        @test results[1].grad != results[2].grad || results[1].theta == results[2].theta
    end
    println("Compiled thunk reuse test passed with $(length(results)) evaluations")
end

@testset "Reactant short MCMC run" begin
    rv_table = Table(
        epoch=collect(range(50000.0, 50500.0, length=20)),
        rv=10.0 .* sin.(2π .* collect(range(0, 1, length=20))) .+ randn(Random.Xoshiro(99), 20) .* 2.0,
        σ_rv=fill(2.0, 20),
    )
    rvlike = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV", device=Reactant.CPU())

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
        name="MCMCReactantSys",
        companions=[b],
        observations=[rvlike],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        end
    )
    model = Octofitter.LogDensityModel(sys, verbosity=0)

    println("Running short MCMC chain with Reactant...")
    chain = octofit(
        model,
        iterations=100,
        adaptation=50,
        verbosity=0
    )
    @test length(chain) > 0
    println("Reactant MCMC run completed with $(length(chain)) samples")
end
