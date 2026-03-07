#!/usr/bin/env julia
"""
Test whether AbstractGPs works with Enzyme AD through Octofitter's pipeline.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Octofitter
using OctofitterRadialVelocity
using .OctofitterRadialVelocity.AbstractGPs: GP, SqExponentialKernel, ScaleTransform, PeriodicKernel
using Distributions
using PlanetOrbits

# ──────────────────────────────────────────────
# 1. Build a minimal RV + AbstractGPs model
# ──────────────────────────────────────────────

epochs = collect(range(50000.0, 50200.0, length=30))
rv_values = 10.0 .* sin.(2π .* epochs ./ 25.0) .+ randn(30) .* 2.0
σ_rv = fill(2.0, 30)

rv_data = Table(epoch=epochs, rv=rv_values, σ_rv=σ_rv)

rvlike = StarAbsoluteRVObs(
    rv_data,
    name="test_inst",
    gaussian_process = θ_obs -> GP(
        θ_obs.gp_η₁^2 * SqExponentialKernel() ∘ ScaleTransform(1/θ_obs.gp_η₂)
    ),
    variables = @variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
        gp_η₁ ~ LogUniform(1, 100)
        gp_η₂ ~ LogUniform(1, 100)
    end
)

planet_b = Planet(
    name = "b",
    basis = RadialVelocityOrbit,
    observations = [],
    variables = @variables begin
        e = 0.0
        ω = 0.0
        P ~ LogUniform(1, 1000)
        M = system.M
        a = cbrt(M * P^2)
        τ ~ Uniform(0, 1)
        tp = τ * P * 365.25 + 50000
        mass ~ LogUniform(0.001, 10)
    end
)

sys = System(
    name = "test_abstractgps",
    companions = [planet_b],
    observations = [rvlike],
    variables = @variables begin
        M ~ truncated(Normal(1.0, 0.1), lower=0.1)
    end
)

# ──────────────────────────────────────────────
# 2. Build model and test
# ──────────────────────────────────────────────

println("═══ Building LogDensityModel ═══")
model = Octofitter.LogDensityModel(sys)

θ = Octofitter.sample_priors(sys)
println("Parameter vector length: ", length(θ))

println("\n═══ Test 1: Primal log-likelihood ═══")
ll = model.ℓπcallback(θ)
println("  log-likelihood = $ll")
println("  isfinite = $(isfinite(ll))")

println("\n═══ Test 2: Gradient via ∇ℓπcallback (Enzyme) ═══")
# Try multiple draws to find one with finite ll and nonzero gradient
for attempt in 1:20
    global θ = Octofitter.sample_priors(sys)
    ll_test = model.ℓπcallback(θ)
    if isfinite(ll_test) && ll_test > -1e100
        println("  Found good θ on attempt $attempt (ll=$ll_test)")
        break
    end
end
try
    ll_grad, grad = model.∇ℓπcallback(θ)
    println("  log-likelihood = $ll_grad")
    println("  gradient norm  = $(sqrt(sum(grad.^2)))")
    println("  all finite?    = $(all(isfinite, grad))")
    println("  any nonzero?   = $(any(!iszero, grad))")
    println("\n  ✅ AbstractGPs works with Enzyme!")
catch err
    println("  ❌ Enzyme gradient failed:")
    println("  Error type: $(typeof(err))")
    showerror(stdout, err, catch_backtrace())
    println()
end

println("\nDone!")
