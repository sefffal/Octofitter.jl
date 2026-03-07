#!/usr/bin/env julia
#=
Benchmark: GaiaDR4AstromObs — Enzyme Reverse vs ForwardDiff

Compares gradient performance (runtime and allocations) for the Gaia DR4
astrometry observation type using different AD backends.

Usage:
  julia --project=@octo-per-obs-diff test/benchmark_dr4_enzyme.jl
=#

using Octofitter
using Distributions
using TypedTables
using BenchmarkTools
using Printf
using Random
using Statistics

using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoEnzyme
using Enzyme: Reverse, set_runtime_activity, Const

# ═══════════════════════════════════════════════════════════════════════
#  Build a GaiaDR4AstromObs model (bypass network with inner constructor)
# ═══════════════════════════════════════════════════════════════════════

function build_dr4_system(; N_epochs=50, primary_star_perturbation=false)
    mock_epochs = collect(range(58000.0, 59500.0, length=N_epochs))
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
        scan_pos_angle = mock_scan_angles,
        centroid_pos_al = randn(MersenneTwister(42), N_epochs) .* 0.5,
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

    ref_epoch_mjd = 57936.375  # DR3 reference epoch
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
        primary_star_perturbation,
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
        name="DR4Bench",
        companions=[b],
        observations=[gaia_obs],
        variables=@variables begin
            M = 1.0
            plx ~ Uniform(10, 30)
        end
    )
    return sys
end

# ═══════════════════════════════════════════════════════════════════════
#  Helper: find a finite starting point
# ═══════════════════════════════════════════════════════════════════════
function find_finite_start(model; max_attempts=50)
    for attempt in 0:max_attempts
        theta_natural = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
        theta = model.link(theta_natural)
        ll = model.ℓπcallback(theta)
        if isfinite(ll)
            return theta
        end
    end
    error("Could not find finite starting point in $max_attempts attempts")
end

# ═══════════════════════════════════════════════════════════════════════
#  Main benchmark
# ═══════════════════════════════════════════════════════════════════════

println("=" ^ 72)
println("  GaiaDR4AstromObs — Enzyme vs ForwardDiff Benchmark")
println("=" ^ 72)

sys = build_dr4_system(N_epochs=50)

# Verify ad_backend dispatch
gaia_obs = sys.observations[1]
backend = Octofitter.ad_backend(gaia_obs)
println("\n  ad_backend(GaiaDR4AstromObs) = $(typeof(backend))")

# ─── TEST A: Default (Enzyme Reverse for DR4 term) ───
println("\n─── TEST A: Default backend (Enzyme Reverse for DR4) ───")
model_default = Octofitter.LogDensityModel(sys; verbosity=0)
Octofitter.initialize!(model_default; verbosity=0)
t_default = model_default.starting_points[1]

# Warm-up
model_default.∇ℓπcallback(t_default)

bA = @benchmark $(model_default.∇ℓπcallback)($t_default) samples=200 evals=5
@printf("  ∇ℓπ  median: %8.1f μs  allocs: %d  memory: %d bytes\n",
    median(bA).time/1e3, median(bA).allocs, median(bA).memory)

# ─── TEST B: ForwardDiff override ───
println("\n─── TEST B: ForwardDiff override ───")
model_fwd = Octofitter.LogDensityModel(sys; autodiff=AutoForwardDiff(), verbosity=0)
Octofitter.initialize!(model_fwd; verbosity=0)
t_fwd = model_fwd.starting_points[1]

model_fwd.∇ℓπcallback(t_fwd)

bB = @benchmark $(model_fwd.∇ℓπcallback)($t_fwd) samples=200 evals=5
@printf("  ∇ℓπ  median: %8.1f μs  allocs: %d  memory: %d bytes\n",
    median(bB).time/1e3, median(bB).allocs, median(bB).memory)

# ─── TEST C: FiniteDiff reference ───
println("\n─── TEST C: FiniteDiff reference (correctness check) ───")
model_fd = Octofitter.LogDensityModel(sys; autodiff=AutoFiniteDiff(), verbosity=0)
Octofitter.initialize!(model_fd; verbosity=0)
t_fd = model_fd.starting_points[1]

model_fd.∇ℓπcallback(t_fd)

bC = @benchmark $(model_fd.∇ℓπcallback)($t_fd) samples=50 evals=2
@printf("  ∇ℓπ  median: %8.1f μs  allocs: %d  memory: %d bytes\n",
    median(bC).time/1e3, median(bC).allocs, median(bC).memory)

# ─── Gradient correctness (same theta for both) ───
println("\n─── Gradient correctness: Enzyme vs FiniteDiff ───")
ll_enz, grad_enz = model_default.∇ℓπcallback(t_default)
ll_fd, grad_fd = model_fd.∇ℓπcallback(t_default)
@printf("  ll (Enzyme):    %.8f\n", ll_enz)
@printf("  ll (FiniteDiff): %.8f\n", ll_fd)
@printf("  ll match:        %s\n", isapprox(ll_enz, ll_fd, atol=1e-6) ? "YES" : "NO ($(abs(ll_enz - ll_fd)))")
grad_match = isapprox(grad_enz, grad_fd, atol=1e-2, rtol=1e-3)
@printf("  grad match:      %s\n", grad_match ? "YES" : "NO")
if !grad_match
    max_diff = maximum(abs.(grad_enz .- grad_fd))
    @printf("  max |diff|:      %.6e\n", max_diff)
    for i in eachindex(grad_enz)
        d = abs(grad_enz[i] - grad_fd[i])
        if d > 1e-2
            @printf("    param %2d: enzyme=%.6f  fd=%.6f  diff=%.6e\n",
                i, grad_enz[i], grad_fd[i], d)
        end
    end
end

# ─── TEST D: Allocation breakdown ───
println("\n─── TEST D: Allocation breakdown ───")
allocs_grad_enz = @allocated model_default.∇ℓπcallback(t_default)
allocs_ll_enz = @allocated model_default.ℓπcallback(t_default)
allocs_grad_fwd = @allocated model_fwd.∇ℓπcallback(t_fwd)
allocs_ll_fwd = @allocated model_fwd.ℓπcallback(t_fwd)

@printf("  Enzyme:     ∇ℓπ = %d bytes,  ℓπ = %d bytes\n", allocs_grad_enz, allocs_ll_enz)
@printf("  ForwardDiff: ∇ℓπ = %d bytes,  ℓπ = %d bytes\n", allocs_grad_fwd, allocs_ll_fwd)

# ─── TEST E: Scaling with number of epochs ───
println("\n─── TEST E: Scaling with number of epochs ───")
for n_ep in [10, 25, 50, 100, 200]
    sys_scaled = build_dr4_system(N_epochs=n_ep)

    # Enzyme (default)
    m_enz = Octofitter.LogDensityModel(sys_scaled; verbosity=0)
    Octofitter.initialize!(m_enz; verbosity=0)
    t_e = m_enz.starting_points[1]
    m_enz.∇ℓπcallback(t_e)  # warmup
    b_enz = @benchmark $(m_enz.∇ℓπcallback)($t_e) samples=100 evals=3

    # ForwardDiff
    m_fwd_s = Octofitter.LogDensityModel(sys_scaled; autodiff=AutoForwardDiff(), verbosity=0)
    Octofitter.initialize!(m_fwd_s; verbosity=0)
    t_f = m_fwd_s.starting_points[1]
    m_fwd_s.∇ℓπcallback(t_f)  # warmup
    b_fwd_s = @benchmark $(m_fwd_s.∇ℓπcallback)($t_f) samples=100 evals=3

    speedup = median(b_fwd_s).time / median(b_enz).time
    @printf("  N=%3d  Enzyme: %7.1f μs (%3d allocs)  FwdDiff: %7.1f μs (%3d allocs)  speedup: %.2fx\n",
        n_ep,
        median(b_enz).time/1e3, median(b_enz).allocs,
        median(b_fwd_s).time/1e3, median(b_fwd_s).allocs,
        speedup)
end

# ─── SUMMARY ───
println("\n" * "=" ^ 72)
println("  SUMMARY — GaiaDR4AstromObs (N=50 epochs)")
println("=" ^ 72)
@printf("  %-30s %10s %10s %12s\n", "Backend", "∇ℓπ (μs)", "allocs", "memory (B)")
println("  " * "-" ^ 62)
@printf("  %-30s %10.1f %10d %12d\n", "Enzyme Reverse (default)",
    median(bA).time/1e3, median(bA).allocs, median(bA).memory)
@printf("  %-30s %10.1f %10d %12d\n", "ForwardDiff (override)",
    median(bB).time/1e3, median(bB).allocs, median(bB).memory)
@printf("  %-30s %10.1f %10d %12d\n", "FiniteDiff (reference)",
    median(bC).time/1e3, median(bC).allocs, median(bC).memory)
speedup_main = median(bB).time / median(bA).time
@printf("\n  Enzyme/ForwardDiff speedup: %.2fx\n", speedup_main)
@printf("  Enzyme/ForwardDiff alloc ratio: %.2fx\n", median(bA).allocs / max(median(bB).allocs, 1))
println("=" ^ 72)
