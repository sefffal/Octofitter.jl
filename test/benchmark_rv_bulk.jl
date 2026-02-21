#!/usr/bin/env julia
#=
Benchmark: MarginalizedStarAbsoluteRVObs primal and gradient evaluation
as a function of the number of RV observations.

Compares Enzyme (baseline) vs Reactant/XLA compiled gradients.
Run with:
    julia --project=. test/benchmark_rv_bulk.jl
=#

using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using Distributions
using TypedTables
using BenchmarkTools
using Random
using Printf

const USE_REACTANT = "--reactant" in ARGS || "-r" in ARGS
if USE_REACTANT
    using Reactant
end

# ─── Build a model with N_rv radial velocity observations ────────────
function build_rv_model(N_rv; seed=42, device=nothing)
    rng = Random.Xoshiro(seed)
    epochs = collect(range(50_000.0, 50_000.0 + 2000.0, length=N_rv))
    rv_table = Table(
        epoch = epochs,
        rv    = 10.0 .* sin.(2π .* range(0, 2, length=N_rv)) .+ randn(rng, N_rv) .* 3.0,
        σ_rv  = fill(3.0, N_rv),
    )

    rvlike = MarginalizedStarAbsoluteRVObs(rv_table; name="RV",
        device=device,
        variables=@variables begin
            jitter ~ LogUniform(0.001, 100)
        end
    )

    b = Planet(
        name = "b",
        basis = RadialVelocityOrbit,
        observations = [],
        variables = @variables begin
            e  ~ Uniform(0.0, 0.5)
            ω  ~ UniformCircular()
            τ  ~ Uniform(0, 1)
            P  ~ Uniform(0.001, 10)
            a  = ∛(P^2 * system.M)
            tp = τ * P * 365.25 + 50000
            mass ~ Uniform(0, 100)
        end
    )

    sys = System(
        name = "BenchSys",
        companions = [b],
        observations = [rvlike],
        variables = @variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        end
    )

    model = Octofitter.LogDensityModel(sys; verbosity=0)
    return model
end

# ─── Find a finite starting point ────────────────────────────────────
function find_finite_theta(model; max_attempts=50)
    for attempt in 0:max_attempts
        θ_natural = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
        θ_t = model.link(θ_natural)
        ll = model.ℓπcallback(θ_t)
        if isfinite(ll)
            return θ_t, ll
        end
    end
    error("Could not find a finite starting point after $max_attempts attempts")
end

# ─── Benchmark a single configuration ────────────────────────────────
function bench_config(N_rv; samples=200, evals=3, device=nothing)
    model = build_rv_model(N_rv; device=device)
    θ_t, ll = find_finite_theta(model)

    # Warm up both paths
    model.ℓπcallback(θ_t)
    model.∇ℓπcallback(θ_t)

    b_primal = @benchmark $(model.ℓπcallback)($θ_t) samples=samples evals=evals
    b_grad   = @benchmark $(model.∇ℓπcallback)($θ_t) samples=samples evals=evals

    t_primal = median(b_primal).time / 1e3  # μs
    t_grad   = median(b_grad).time / 1e3    # μs
    ratio    = t_grad / t_primal

    allocs_primal = median(b_primal).allocs
    allocs_grad   = median(b_grad).allocs
    mem_primal    = median(b_primal).memory
    mem_grad      = median(b_grad).memory

    return (; N_rv, t_primal, t_grad, ratio,
              allocs_primal, allocs_grad, mem_primal, mem_grad,
              b_primal, b_grad)
end

# ═══════════════════════════════════════════════════════════════════════
# Run benchmarks
# ═══════════════════════════════════════════════════════════════════════
N_values = [5, 10, 20, 50, 100, 200, 500, 1000]

println("="^78)
println("  Benchmark: MarginalizedStarAbsoluteRVObs (orbitsolve_bulk + Enzyme)")
println("  1 planet, RadialVelocityOrbit basis")
println("="^78)

results = []
for N in N_values
    print("  N_rv = $N ... ")
    r = bench_config(N)
    @printf("primal: %8.1f μs | grad: %8.1f μs | ratio: %.1fx\n",
            r.t_primal, r.t_grad, r.ratio)
    push!(results, r)
end

# ─── Summary table ────────────────────────────────────────────────────
println("\n", "="^78)
println("  SUMMARY")
println("="^78)
@printf("  %6s  %10s  %10s  %8s  %12s  %12s\n",
        "N_rv", "primal(μs)", "grad(μs)", "grad/pri", "alloc_pri", "alloc_grad")
println("  ", "-"^72)
for r in results
    @printf("  %6d  %10.1f  %10.1f  %8.1fx  %10d B  %10d B\n",
            r.N_rv, r.t_primal, r.t_grad, r.ratio,
            r.mem_primal, r.mem_grad)
end
println("="^78)

# ─── Scaling analysis ─────────────────────────────────────────────────
println("\n  Scaling (relative to N_rv=$(results[1].N_rv)):")
t0_p = results[1].t_primal
t0_g = results[1].t_grad
N0   = results[1].N_rv
for r in results
    scale_p = r.t_primal / t0_p
    scale_g = r.t_grad / t0_g
    n_ratio = r.N_rv / N0
    @printf("    N=%4d: primal %.1fx (N ratio %.0fx) | grad %.1fx (N ratio %.0fx)\n",
            r.N_rv, scale_p, n_ratio, scale_g, n_ratio)
end
println()

# ═══════════════════════════════════════════════════════════════════════
# Reactant benchmarks (when --reactant flag is passed)
# ═══════════════════════════════════════════════════════════════════════
if USE_REACTANT
    println("\n", "="^78)
    println("  Benchmark: Reactant/XLA compiled gradients")
    println("  1 planet, RadialVelocityOrbit basis, device=Reactant.CPU()")
    println("="^78)

    results_reactant = []
    for N in N_values
        print("  N_rv = $N ... ")
        r_reactant = bench_config(N; device=Reactant.CPU())
        @printf("primal: %8.1f μs | grad: %8.1f μs | ratio: %.1fx\n",
                r_reactant.t_primal, r_reactant.t_grad, r_reactant.ratio)
        push!(results_reactant, r_reactant)
    end

    # ─── Comparison table ─────────────────────────────────────────────
    println("\n", "="^78)
    println("  COMPARISON: Enzyme vs Reactant")
    println("="^78)
    @printf("  %6s  %12s  %12s  %10s\n",
            "N_rv", "Enzyme(μs)", "Reactant(μs)", "Speedup")
    println("  ", "-"^46)
    for (re, rr) in zip(results, results_reactant)
        speedup = re.t_grad / rr.t_grad
        @printf("  %6d  %12.1f  %12.1f  %9.1fx\n",
                re.N_rv, re.t_grad, rr.t_grad, speedup)
    end
    println("="^78)
end
