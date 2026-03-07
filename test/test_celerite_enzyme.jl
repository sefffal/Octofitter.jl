#!/usr/bin/env julia
"""
Test whether Celerite GP works with Enzyme AD through Octofitter's per-term pipeline.

Key insight: Enzyme operates on LLVM IR with Float64 directly, so Celerite's
hardcoded Float64 types are actually FINE (unlike ForwardDiff which needs Dual propagation).
The main concerns are: try-catch blocks, type instability, and dynamic dispatch.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Octofitter
using OctofitterRadialVelocity
using OctofitterRadialVelocity.Celerite
using Distributions
using PlanetOrbits

# ──────────────────────────────────────────────
# 1. Build a minimal RV + Celerite GP model
# ──────────────────────────────────────────────

# Generate some synthetic RV data
epochs = collect(range(50000.0, 50200.0, length=30))
rv_values = 10.0 .* sin.(2π .* epochs ./ 25.0) .+ randn(30) .* 2.0
σ_rv = fill(2.0, 30)

rv_data = Table(epoch=epochs, rv=rv_values, σ_rv=σ_rv)

rvlike = StarAbsoluteRVObs(
    rv_data,
    name="test_inst",
    gaussian_process = θ_obs -> Celerite.CeleriteGP(
        Celerite.RealTerm(
            log(θ_obs.gp_B * (1 + θ_obs.gp_C) / (2 + θ_obs.gp_C)),
            log(1 / θ_obs.gp_L)
        ) + Celerite.ComplexTerm(
            log(θ_obs.gp_B / (2 + θ_obs.gp_C)),
            -Inf,
            log(1 / θ_obs.gp_L),
            log(2π / θ_obs.gp_Prot)
        )
    ),
    variables = @variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
        gp_B ~ LogUniform(1, 1000)
        gp_C ~ Uniform(0.01, 10)
        gp_L ~ LogUniform(2, 200)
        gp_Prot ~ Uniform(5, 50)
    end
)

planet_b = Planet(
    name = "b",
    basis = RadialVelocityOrbit,
    observations = [],
    variables = @variables begin
        e = 0.0
        ω = 0.0
        P ~ LogUniform(1, 1000)     # days → years
        M = system.M
        a = cbrt(M * P^2)
        τ ~ Uniform(0, 1)
        tp = τ * P * 365.25 + 50000
        mass ~ LogUniform(0.001, 10)
    end
)

sys = System(
    name = "test_celerite",
    companions = [planet_b],
    observations = [rvlike],
    variables = @variables begin
        M ~ truncated(Normal(1.0, 0.1), lower=0.1)
    end
)

# ──────────────────────────────────────────────
# 2. Test primal evaluation
# ──────────────────────────────────────────────

println("═══ Building LogDensityModel ═══")
model = Octofitter.LogDensityModel(sys)

# Draw valid parameters
θ = Octofitter.sample_priors(sys)
println("Parameter vector length: ", length(θ))

println("\n═══ Test 1: Primal log-likelihood ═══")
ll = model.ℓπcallback(θ)
println("  log-likelihood = $ll")
println("  isfinite = $(isfinite(ll))")

# ──────────────────────────────────────────────
# 3. Test gradient (the key question!)
# ──────────────────────────────────────────────

println("\n═══ Test 2: Gradient via ∇ℓπcallback (Enzyme) ═══")
try
    ll_grad, grad = model.∇ℓπcallback(θ)
    println("  log-likelihood = $ll_grad")
    println("  gradient norm  = $(sqrt(sum(grad.^2)))")
    println("  gradient       = $grad")
    println("  all finite?    = $(all(isfinite, grad))")
    println("\n  ✅ Celerite GP works with Enzyme!")
catch err
    println("  ❌ Enzyme gradient failed:")
    println("  Error type: $(typeof(err))")
    showerror(stdout, err, catch_backtrace())
    println()
end

# ──────────────────────────────────────────────
# 4. Allocation check
# ──────────────────────────────────────────────

println("\n═══ Test 3: Allocation check ═══")
# Warm up
model.ℓπcallback(θ)
alloc_primal = @allocated model.ℓπcallback(θ)
println("  Primal allocations: $alloc_primal bytes")

if @isdefined(ll_grad)
    model.∇ℓπcallback(θ)
    alloc_grad = @allocated model.∇ℓπcallback(θ)
    println("  Gradient allocations: $alloc_grad bytes")
end

println("\nDone!")
