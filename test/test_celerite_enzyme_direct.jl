#!/usr/bin/env julia
"""
Direct test of Celerite GP with Enzyme, bypassing the full Octofitter pipeline.
Tests whether the core Celerite operations (compute!, log_likelihood) are
Enzyme-differentiable.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using OctofitterRadialVelocity
using OctofitterRadialVelocity.Celerite
import Enzyme
using Enzyme: autodiff, Reverse, Active, Const, Duplicated

# ──────────────────────────────────────────────
# 1. Set up test data
# ──────────────────────────────────────────────

N = 30
epochs = collect(range(50000.0, 50200.0, length=N))
rv_residuals = randn(N) .* 5.0
σ_rv = fill(2.0, N)

# GP hyperparameters
gp_B = 50.0
gp_C = 2.0
gp_L = 30.0
gp_Prot = 15.0

# ──────────────────────────────────────────────
# 2. Test: basic Celerite forward pass
# ──────────────────────────────────────────────

println("═══ Test 1: Celerite forward pass ═══")
gp = Celerite.CeleriteGP(
    Celerite.RealTerm(
        log(gp_B * (1 + gp_C) / (2 + gp_C)),
        log(1 / gp_L)
    ) + Celerite.ComplexTerm(
        log(gp_B / (2 + gp_C)),
        -Inf,
        log(1 / gp_L),
        log(2π / gp_Prot)
    )
)
Celerite.compute!(gp, epochs, σ_rv)
ll = Celerite.log_likelihood(gp, rv_residuals)
println("  log-likelihood = $ll")
println("  isfinite = $(isfinite(ll))")

# ──────────────────────────────────────────────
# 3. Test: Enzyme gradient of a simple function
#    that wraps Celerite compute! + log_likelihood
# ──────────────────────────────────────────────

println("\n═══ Test 2: Enzyme gradient of Celerite ═══")

# Simple function: hyperparams -> log-likelihood
function celerite_ll(params::Vector{Float64}, epochs::Vector{Float64},
                      σ_rv::Vector{Float64}, rv_residuals::Vector{Float64})
    B, C, L, Prot = params[1], params[2], params[3], params[4]

    kernel = Celerite.RealTerm(
        log(B * (1 + C) / (2 + C)),
        log(1 / L)
    ) + Celerite.ComplexTerm(
        log(B / (2 + C)),
        -Inf,
        log(1 / L),
        log(2π / Prot)
    )

    gp = Celerite.CeleriteGP(kernel)
    Celerite.compute!(gp, epochs, σ_rv)
    return Celerite.log_likelihood(gp, rv_residuals)
end

# Test primal
params = [gp_B, gp_C, gp_L, gp_Prot]
ll_primal = celerite_ll(params, epochs, σ_rv, rv_residuals)
println("  Primal log-likelihood = $ll_primal")

# Test Enzyme gradient
println("  Attempting Enzyme.autodiff...")
try
    dparams = zeros(4)
    autodiff(
        Enzyme.set_runtime_activity(Reverse),
        celerite_ll,
        Active,
        Duplicated(copy(params), dparams),
        Const(epochs),
        Const(σ_rv),
        Const(rv_residuals)
    )
    println("  ✅ Enzyme gradient succeeded!")
    println("  d(ll)/d(B)    = $(dparams[1])")
    println("  d(ll)/d(C)    = $(dparams[2])")
    println("  d(ll)/d(L)    = $(dparams[3])")
    println("  d(ll)/d(Prot) = $(dparams[4])")
    println("  all finite?   = $(all(isfinite, dparams))")
catch err
    println("  ❌ Enzyme gradient failed:")
    println("  Error type: $(typeof(err))")
    showerror(stdout, err, catch_backtrace())
    println()
end

# ──────────────────────────────────────────────
# 4. Test: finite difference comparison
# ──────────────────────────────────────────────

println("\n═══ Test 3: Finite difference comparison ═══")
ε = 1e-6
fd_grad = similar(params)
for i in eachindex(params)
    p_plus = copy(params); p_plus[i] += ε
    p_minus = copy(params); p_minus[i] -= ε
    fd_grad[i] = (celerite_ll(p_plus, epochs, σ_rv, rv_residuals) -
                  celerite_ll(p_minus, epochs, σ_rv, rv_residuals)) / (2ε)
end
println("  FD gradient: $fd_grad")

println("\nDone!")
