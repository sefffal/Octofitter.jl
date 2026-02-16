#!/usr/bin/env julia
#=
Benchmark: Per-term active parameter subsetting (gradient sparsity).

Creates two models with identical likelihood structure but different
parameter layouts to demonstrate the speedup from exploiting gradient
sparsity:

  SPARSE: 10 nuisance parameters in each of 2 observation objects
          → each term differentiates through D_active ≈ D - 10 params

  DENSE:  all 20 nuisance parameters at the system level
          → every term differentiates through the full D params

Both models have the same total D, but the sparse layout lets each
per-term ForwardDiff pass use a smaller chunk size.
=#

using Octofitter
using Distributions
using TypedTables
using BenchmarkTools
using OrderedCollections
using Printf

# ─── Fake astrometry data (shared by both models) ─────────────────────
astrom_table = Table(
    epoch = collect(range(50_000.0, 50_500.0, length=10)),
    ra    = randn(10) .* 100,
    dec   = randn(10) .* 100,
    σ_ra  = fill(5.0, 10),
    σ_dec = fill(5.0, 10),
)

N_NUISANCE = 10   # nuisance parameters per observation

# ─── Helper: build @variables block with N nuisance priors ────────────
function _make_nuisance_variables(prefix::Symbol, n::Int)
    priors_dict = OrderedCollections.OrderedDict{Symbol,Any}()
    for i in 1:n
        priors_dict[Symbol(prefix, :_, i)] = Normal(0, 10)
    end
    priors = Octofitter.Priors(priors_dict)
    derived = Octofitter.Derived(
        OrderedCollections.OrderedDict{Symbol,Any}(),
        (), ()
    )
    return (priors, derived)
end

# ═══════════════════════════════════════════════════════════════════════
# MODEL A — "SPARSE": nuisance params live in the observation objects
# ═══════════════════════════════════════════════════════════════════════
function build_sparse_model(; n_obs=2, n_nuisance=N_NUISANCE)
    obs_list = [
        PlanetRelAstromLikelihood(
            astrom_table;
            name="astrom$k",
            variables=_make_nuisance_variables(Symbol(:obs, k), n_nuisance)
        )
        for k in 1:n_obs
    ]

    planet_b = Planet(
        name = "b",
        basis = Visual{KepOrbit},
        observations = obs_list,
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

    sys = System(
        name = "SparseSystem",
        companions = [planet_b],
        observations = [],
        variables = @variables begin
            M   ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 5.0), lower=1.0)
        end
    )

    return Octofitter.LogDensityModel(sys; verbosity=1)
end

# ═══════════════════════════════════════════════════════════════════════
# MODEL B — "DENSE": all nuisance params moved to system level
# ═══════════════════════════════════════════════════════════════════════
function build_dense_model(; n_obs=2, n_nuisance=N_NUISANCE)
    # No per-observation variables
    obs_list = [
        PlanetRelAstromLikelihood(astrom_table; name="astrom$k")
        for k in 1:n_obs
    ]

    planet_b = Planet(
        name = "b",
        basis = Visual{KepOrbit},
        observations = obs_list,
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

    # Build system-level nuisance priors programmatically
    sys_priors_dict = OrderedCollections.OrderedDict{Symbol,Any}(
        :M   => truncated(Normal(1.2, 0.1), lower=0.1),
        :plx => truncated(Normal(50.0, 5.0), lower=1.0),
    )
    for k in 1:n_obs
        for i in 1:n_nuisance
            sys_priors_dict[Symbol(:obs, k, :_, i)] = Normal(0, 10)
        end
    end
    sys_priors = Octofitter.Priors(sys_priors_dict)
    sys_derived = Octofitter.Derived(
        OrderedCollections.OrderedDict{Symbol,Any}(),
        (), ()
    )

    sys = System(
        name = "DenseSystem",
        companions = [planet_b],
        observations = [],
        variables = (sys_priors, sys_derived)
    )

    return Octofitter.LogDensityModel(sys; verbosity=1)
end

# ═══════════════════════════════════════════════════════════════════════
# Benchmark helper
# ═══════════════════════════════════════════════════════════════════════
function run_benchmark(; n_obs, n_nuisance)
    println("\n", "="^70)
    println("  Config: $n_obs observations × $n_nuisance nuisance params each")
    println("="^70)

    println("\nBuilding SPARSE model (nuisance in observations)...")
    model_sparse = build_sparse_model(; n_obs, n_nuisance)
    D_sparse = model_sparse.D

    println("Building DENSE model (nuisance at system level)...")
    model_dense = build_dense_model(; n_obs, n_nuisance)
    D_dense = model_dense.D

    # Show active indices info
    active_all, _ = Octofitter._compute_active_indices(model_sparse.system)
    println("\n  Total D: $D_sparse (both models)")
    for (i, active) in enumerate(active_all)
        println("  Sparse term $i: D_active = $(length(active)) / $D_sparse  ($(round(100*length(active)/D_sparse, digits=1))%)")
    end

    # Sample starting points
    θ_sparse = Octofitter.sample_priors(model_sparse.system)
    θ_sparse_t = model_sparse.link(θ_sparse)
    θ_dense = Octofitter.sample_priors(model_dense.system)
    θ_dense_t = model_dense.link(θ_dense)

    # Warm up
    model_sparse.∇ℓπcallback(θ_sparse_t)
    model_dense.∇ℓπcallback(θ_dense_t)

    println("\n  Benchmarking gradient evaluations...")

    b_sparse = @benchmark $(model_sparse.∇ℓπcallback)($θ_sparse_t) samples=100 evals=5
    b_dense  = @benchmark $(model_dense.∇ℓπcallback)($θ_dense_t)  samples=100 evals=5

    t_sparse = median(b_sparse).time / 1e3  # μs
    t_dense  = median(b_dense).time / 1e3
    speedup  = t_dense / t_sparse

    println("  SPARSE median: $(round(t_sparse, digits=1)) μs")
    println("  DENSE  median: $(round(t_dense, digits=1)) μs")
    println("  Speedup: $(round(speedup, digits=2))x")

    return (; n_obs, n_nuisance, D=D_sparse, t_sparse, t_dense, speedup,
              b_sparse, b_dense)
end

# ═══════════════════════════════════════════════════════════════════════
# Run benchmarks across several configurations
# ═══════════════════════════════════════════════════════════════════════
configs = [
    (n_obs=2,  n_nuisance=10),   # 2×10 = 20 obs params, D=29
    (n_obs=4,  n_nuisance=10),   # 4×10 = 40 obs params, D=49
    (n_obs=2,  n_nuisance=20),   # 2×20 = 40 obs params, D=49
    (n_obs=4,  n_nuisance=20),   # 4×20 = 80 obs params, D=89
]

results = []
for cfg in configs
    r = run_benchmark(; cfg...)
    push!(results, r)
end

# Summary table
println("\n\n", "="^70)
println("  SUMMARY")
println("="^70)
println("  n_obs  n_nuis   D   sparse(μs)  dense(μs)  speedup")
println("  " * "-"^60)
for r in results
    @printf("  %3d    %3d    %3d   %8.1f    %8.1f     %.2fx\n",
            r.n_obs, r.n_nuisance, r.D, r.t_sparse, r.t_dense, r.speedup)
end
println("="^70)
