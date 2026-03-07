#!/usr/bin/env julia
#=
Benchmark suite for PR #91: Per-observation autodiff refactor.

Three test cases:
  1. Simple model (1 PlanetRelAstromLikelihood) — regression check
  2. Multi-instrument RV (5 StarAbsoluteRVObservation, each with jitter+offset) — sparsity win
  3. MarginalizedStarAbsoluteRVObs — Enzyme (new) vs ForwardDiff (old)

Usage:
  julia --project=@octo-per-obs-diff test/benchmark_pr91.jl

Results are written to test/benchmark_results_<branchname>.txt
=#

using Octofitter
using OctofitterRadialVelocity
using Distributions
using TypedTables
using BenchmarkTools
using Printf
using Random
using Dates

# Identify which branch we're on
branch = strip(read(`git rev-parse --abbrev-ref HEAD`, String))
commit = strip(read(`git rev-parse --short HEAD`, String))

results_file = joinpath(@__DIR__, "benchmark_results_$(branch).txt")
io = open(results_file, "w")

function log_result(msg)
    println(msg)
    println(io, msg)
end

log_result("=" ^ 72)
log_result("  Octofitter PR #91 Benchmarks")
log_result("  Branch: $branch  Commit: $commit")
log_result("  Julia: $(VERSION)  Date: $(Dates.now())")
log_result("=" ^ 72)

# ═══════════════════════════════════════════════════════════════════════
# TEST 1: Simple relative astrometry model (regression check)
# ═══════════════════════════════════════════════════════════════════════
log_result("\n" * "─" ^ 72)
log_result("  TEST 1: Simple relative astrometry model")
log_result("─" ^ 72)

astrom_table = Table(
    epoch = collect(range(50_000.0, 50_500.0, length=10)),
    ra    = randn(MersenneTwister(1), 10) .* 100,
    dec   = randn(MersenneTwister(2), 10) .* 100,
    σ_ra  = fill(5.0, 10),
    σ_dec = fill(5.0, 10),
)

planet_b_astrom = Planet(
    name = "b",
    basis = Visual{KepOrbit},
    observations = [
        PlanetRelAstromLikelihood(astrom_table; name="astrom1")
    ],
    variables = @variables begin
        a  ~ LogUniform(0.5, 50)
        e  ~ Uniform(0, 0.5)
        ω  ~ Uniform(0, 2π)
        i  ~ Sine()
        Ω  ~ Uniform(0, 2π)
        θ  ~ Uniform(0, 2π)
        tp = θ_at_epoch_to_tperi(θ, 50_000.0; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(1, 1000)
    end
)

sys_astrom = System(
    name = "AstromBench",
    companions = [planet_b_astrom],
    observations = [],
    variables = @variables begin
        M   ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 5.0), lower=1.0)
    end
)

model1 = Octofitter.LogDensityModel(sys_astrom; verbosity=0)
Octofitter.initialize!(model1; verbosity=0)
t1 = model1.starting_points[1]

# Warm up
model1.ℓπcallback(t1)
model1.∇ℓπcallback(t1)

log_result("  D = $(model1.D)")

b1_ll = @benchmark $(model1.ℓπcallback)($t1) samples=200 evals=10
b1_grad = @benchmark $(model1.∇ℓπcallback)($t1) samples=200 evals=10

log_result(@sprintf("  logdensity          median: %8.1f μs  allocs: %d",
    median(b1_ll).time / 1e3, median(b1_ll).allocs))
log_result(@sprintf("  logdensity+gradient  median: %8.1f μs  allocs: %d",
    median(b1_grad).time / 1e3, median(b1_grad).allocs))

# ═══════════════════════════════════════════════════════════════════════
# TEST 2: Multi-instrument RV (sparsity test)
# ═══════════════════════════════════════════════════════════════════════
log_result("\n" * "─" ^ 72)
log_result("  TEST 2: 5-instrument StarAbsoluteRVObs (sparsity test)")
log_result("─" ^ 72)

N_INST = 5
N_EPOCHS_PER_INST = 20

# Generate synthetic RV data for each instrument
rng = MersenneTwister(42)
rv_obs_list = [
    StarAbsoluteRVObs(
        Table(
            epoch = collect(range(50_000.0 + (k-1)*100, length=N_EPOCHS_PER_INST, step=30.0)),
            rv    = randn(rng, N_EPOCHS_PER_INST) .* 50,
            σ_rv  = fill(5.0, N_EPOCHS_PER_INST),
        );
        name = "inst$k",
        variables = @variables begin
            jitter ~ LogUniform(0.1, 100)
            offset ~ Normal(0, 1000)
        end
    )
    for k in 1:N_INST
]

planet_b_rv = Planet(
    name = "b",
    basis = RadialVelocityOrbit,
    observations = [],
    variables = @variables begin
        e ~ Uniform(0.0, 0.5)
        ω ~ UniformCircular()
        τ ~ Uniform(0, 1)
        P ~ Uniform(0.001, 10)
        a = ∛(P^2 * system.M)
        tp = τ * P * 365.25 + 50000
        mass ~ Uniform(0, 100)
    end
)

sys_rv = System(
    name = "RVSparsityBench",
    companions = [planet_b_rv],
    observations = rv_obs_list,
    variables = @variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    end
)

model2 = Octofitter.LogDensityModel(sys_rv; verbosity=0)
Octofitter.initialize!(model2; verbosity=0)
t2 = model2.starting_points[1]

# Warm up
model2.ℓπcallback(t2)
model2.∇ℓπcallback(t2)

log_result("  D = $(model2.D)")

# Show per-term active parameter info (only available on per-obs-diff branch)
if isdefined(Octofitter, :_compute_active_indices)
    active_all, _ = Octofitter._compute_active_indices(model2.system)
    for (i, active) in enumerate(active_all)
        log_result(@sprintf("    Term %d: D_active = %d / %d  (%.0f%%)",
            i, length(active), model2.D, 100*length(active)/model2.D))
    end
end

b2_ll = @benchmark $(model2.ℓπcallback)($t2) samples=200 evals=10
b2_grad = @benchmark $(model2.∇ℓπcallback)($t2) samples=200 evals=10

log_result(@sprintf("  logdensity          median: %8.1f μs  allocs: %d",
    median(b2_ll).time / 1e3, median(b2_ll).allocs))
log_result(@sprintf("  logdensity+gradient  median: %8.1f μs  allocs: %d",
    median(b2_grad).time / 1e3, median(b2_grad).allocs))

# ═══════════════════════════════════════════════════════════════════════
# TEST 3: MarginalizedStarAbsoluteRVObs (Enzyme vs ForwardDiff)
# ═══════════════════════════════════════════════════════════════════════
log_result("\n" * "─" ^ 72)
log_result("  TEST 3: MarginalizedStarAbsoluteRVObs")
log_result("─" ^ 72)

rv_table_margin = Table(
    epoch = collect(range(50_000.0, 50_500.0, length=30)),
    rv    = randn(MersenneTwister(99), 30) .* 50,
    σ_rv  = fill(5.0, 30),
)

margin_obs = MarginalizedStarAbsoluteRVObs(rv_table_margin, name="MargRV")

planet_b_margin = Planet(
    name = "b",
    basis = RadialVelocityOrbit,
    observations = [],
    variables = @variables begin
        e ~ Uniform(0.0, 0.5)
        ω ~ UniformCircular()
        τ ~ Uniform(0, 1)
        P ~ Uniform(0.001, 10)
        a = ∛(P^2 * system.M)
        tp = τ * P * 365.25 + 50000
        mass ~ Uniform(0, 100)
    end
)

sys_margin = System(
    name = "MargRVBench",
    companions = [planet_b_margin],
    observations = [margin_obs],
    variables = @variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    end
)

model3 = Octofitter.LogDensityModel(sys_margin; verbosity=0)
Octofitter.initialize!(model3; verbosity=0)
t3 = model3.starting_points[1]

# Warm up
model3.ℓπcallback(t3)
model3.∇ℓπcallback(t3)

log_result("  D = $(model3.D)")
if isdefined(Octofitter, :ad_backend)
    log_result("  AD backend: $(Octofitter.ad_backend(margin_obs))")
end

b3_ll = @benchmark $(model3.ℓπcallback)($t3) samples=200 evals=10
b3_grad = @benchmark $(model3.∇ℓπcallback)($t3) samples=200 evals=10

log_result(@sprintf("  logdensity          median: %8.1f μs  allocs: %d",
    median(b3_ll).time / 1e3, median(b3_ll).allocs))
log_result(@sprintf("  logdensity+gradient  median: %8.1f μs  allocs: %d",
    median(b3_grad).time / 1e3, median(b3_grad).allocs))

# ═══════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════
log_result("\n" * "=" ^ 72)
log_result("  SUMMARY — $branch ($commit)")
log_result("=" ^ 72)
log_result(@sprintf("  %-40s %10s %10s %8s", "Test", "ℓπ (μs)", "∇ℓπ (μs)", "allocs"))
log_result("  " * "-" ^ 68)
log_result(@sprintf("  %-40s %10.1f %10.1f %8d",
    "1. Simple astrometry (D=$(model1.D))",
    median(b1_ll).time/1e3, median(b1_grad).time/1e3, median(b1_grad).allocs))
log_result(@sprintf("  %-40s %10.1f %10.1f %8d",
    "2. 5-inst RV sparsity (D=$(model2.D))",
    median(b2_ll).time/1e3, median(b2_grad).time/1e3, median(b2_grad).allocs))
log_result(@sprintf("  %-40s %10.1f %10.1f %8d",
    "3. MarginalizedRV (D=$(model3.D))",
    median(b3_ll).time/1e3, median(b3_grad).time/1e3, median(b3_grad).allocs))
log_result("=" ^ 72)

close(io)
println("\nResults written to: $results_file")
