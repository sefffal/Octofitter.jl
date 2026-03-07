#!/usr/bin/env julia
#=
Diagnostic benchmark: isolate the Enzyme reverse-mode bottleneck for
MarginalizedStarAbsoluteRVObs.

Hypotheses to test:
  H1. Allocations inside ln_like (Vector{T}(undef, L) ×2) are expensive for Enzyme
  H2. D=8 is too small for reverse mode to beat forward mode
  H3. Orbit-solving allocations (_solve_all_orbits) add overhead
  H4. DifferentiationInterface overhead vs raw Enzyme

Usage:
  julia --project=@octo-per-obs-diff test/benchmark_enzyme_diag.jl
=#

using Octofitter
using OctofitterRadialVelocity
using Distributions
using TypedTables
using BenchmarkTools
using Printf
using Random

using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoEnzyme
using Enzyme: Reverse, Forward, set_runtime_activity, Const, Duplicated, autodiff
using DifferentiationInterface: value_and_gradient!, prepare_gradient

# ═══════════════════════════════════════════════════════════════════════
#  Build the model (same as benchmark_pr91.jl Test 3)
# ═══════════════════════════════════════════════════════════════════════

rv_table = Table(
    epoch = collect(range(50_000.0, 50_500.0, length=30)),
    rv    = randn(MersenneTwister(99), 30) .* 50,
    σ_rv  = fill(5.0, 30),
)

margin_obs = MarginalizedStarAbsoluteRVObs(rv_table, name="MargRV")

planet_b = Planet(
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

sys = System(
    name = "EnzymeDiag",
    companions = [planet_b],
    observations = [margin_obs],
    variables = @variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    end
)

println("=" ^ 72)
println("  Enzyme Diagnostic Benchmark")
println("=" ^ 72)

# ═══════════════════════════════════════════════════════════════════════
#  TEST A: Full pipeline — Enzyme Reverse (current default)
# ═══════════════════════════════════════════════════════════════════════

println("\n─── TEST A: Full pipeline — Enzyme Reverse (current default) ───")
model_enz_rev = Octofitter.LogDensityModel(sys; verbosity=0)
Octofitter.initialize!(model_enz_rev; verbosity=0)
t_rev = model_enz_rev.starting_points[1]

# Warm-up
model_enz_rev.∇ℓπcallback(t_rev)

bA = @benchmark $(model_enz_rev.∇ℓπcallback)($t_rev) samples=200 evals=10
@printf("  ∇ℓπ  median: %8.1f μs  allocs: %d\n",
    median(bA).time/1e3, median(bA).allocs)

# ═══════════════════════════════════════════════════════════════════════
#  TEST B: Full pipeline — ForwardDiff override
# ═══════════════════════════════════════════════════════════════════════

println("\n─── TEST B: Full pipeline — ForwardDiff override (autodiff=AutoForwardDiff) ───")
model_fwd = Octofitter.LogDensityModel(sys; autodiff=AutoForwardDiff(), verbosity=0)
Octofitter.initialize!(model_fwd; verbosity=0)
t_fwd = model_fwd.starting_points[1]

model_fwd.∇ℓπcallback(t_fwd)

bB = @benchmark $(model_fwd.∇ℓπcallback)($t_fwd) samples=200 evals=10
@printf("  ∇ℓπ  median: %8.1f μs  allocs: %d\n",
    median(bB).time/1e3, median(bB).allocs)

# ═══════════════════════════════════════════════════════════════════════
#  TEST C: Full pipeline — Enzyme Forward mode
# ═══════════════════════════════════════════════════════════════════════

println("\n─── TEST C: Full pipeline — Enzyme Forward override ───")
bC = nothing
try
    model_enz_fwd = Octofitter.LogDensityModel(sys;
        autodiff=AutoEnzyme(mode=set_runtime_activity(Forward), function_annotation=Const),
        verbosity=0)
    Octofitter.initialize!(model_enz_fwd; verbosity=0)
    t_enz_fwd = model_enz_fwd.starting_points[1]

    model_enz_fwd.∇ℓπcallback(t_enz_fwd)

    bC = @benchmark $(model_enz_fwd.∇ℓπcallback)($t_enz_fwd) samples=200 evals=10
    @printf("  ∇ℓπ  median: %8.1f μs  allocs: %d\n",
        median(bC).time/1e3, median(bC).allocs)
catch e
    println("  SKIPPED: Enzyme Forward through full pipeline not supported")
    println("  Reason: $(sprint(showerror, e))")
    println("  (task_local_storage / IdDict not yet supported in Enzyme Forward mode)")
    println("  → Will test Enzyme Forward on isolated closures in Test D instead")
end

# ═══════════════════════════════════════════════════════════════════════
#  TEST D: Allocation profile — where are the 557 allocs?
# ═══════════════════════════════════════════════════════════════════════

println("\n─── TEST D: Allocation profile (Enzyme Reverse) ───")
println("  Profiling with @allocated and --track-allocation equivalent...")

# Profile: full gradient call
allocs_full = @allocated model_enz_rev.∇ℓπcallback(t_rev)
@printf("  Full ∇ℓπ:  %d bytes allocated\n", allocs_full)

# Profile: just logdensity (no AD)
allocs_ll = @allocated model_enz_rev.ℓπcallback(t_rev)
@printf("  Just ℓπ:   %d bytes allocated\n", allocs_ll)

# Now let's isolate the Enzyme term closure vs the prior gradient.
# We can extract the inner closures from the model construction.
println("\n  Per-component allocation breakdown:")

# Use the DI interface directly to isolate the Enzyme term
# Extract the system observation's closure via reconstructing the same setup
println("  (Constructing isolated closures...)")

# Build the same closure manually to test raw Enzyme vs DI
let
    invlink = model_enz_rev.invlink
    arr2nt = model_enz_rev.arr2nt
    n_planets = length(sys.planets)
    OT = Octofitter._planet_orbit_type(sys.planets[1])
    orbit_constructors = (θ_system -> OT(;merge(θ_system, θ_system.planets[1])...),)
    obs = margin_obs
    obs_name = Octofitter.normalizename(Octofitter.likelihoodname(obs))
    epochs = collect(Float64, obs.table.epoch)

    # This is the inner closure that gets differentiated
    function test_closure(θ_transformed)
        θ_natural = invlink(θ_transformed)
        θ_system = arr2nt(θ_natural)
        orbits = ntuple(i -> orbit_constructors[i](θ_system), n_planets)
        θ_obs = getproperty(θ_system.observations, Symbol(obs_name))
        solutions = Octofitter._solve_all_orbits(orbits, epochs)
        ctx = Octofitter.SystemObservationContext(θ_system, θ_obs, orbits, solutions)
        Octofitter.ln_like(obs, ctx)
    end

    # Test the closure value
    val = test_closure(t_rev)
    @printf("  Closure value check: %.6f\n", val)
    allocs_closure = @allocated test_closure(t_rev)
    @printf("  Closure (no AD):  %d bytes allocated\n", allocs_closure)

    # Test DI value_and_gradient! with Enzyme Reverse
    enz_rev_backend = AutoEnzyme(mode=set_runtime_activity(Reverse), function_annotation=Const)
    grad_buf = similar(t_rev)
    prep_rev = prepare_gradient(test_closure, enz_rev_backend, zero(t_rev))
    value_and_gradient!(test_closure, grad_buf, prep_rev, enz_rev_backend, t_rev) # warmup
    allocs_enz_rev = @allocated value_and_gradient!(test_closure, grad_buf, prep_rev, enz_rev_backend, t_rev)
    @printf("  DI + Enzyme Reverse:  %d bytes allocated\n", allocs_enz_rev)

    b_enz_rev = @benchmark value_and_gradient!($test_closure, $grad_buf, $prep_rev, $enz_rev_backend, $t_rev) samples=200 evals=10
    @printf("  DI + Enzyme Reverse:  %8.1f μs  allocs: %d\n",
        median(b_enz_rev).time/1e3, median(b_enz_rev).allocs)

    # Test DI value_and_gradient! with ForwardDiff
    fwd_backend = AutoForwardDiff()
    prep_fwd = prepare_gradient(test_closure, fwd_backend, zero(t_rev))
    value_and_gradient!(test_closure, grad_buf, prep_fwd, fwd_backend, t_rev) # warmup
    allocs_fwd = @allocated value_and_gradient!(test_closure, grad_buf, prep_fwd, fwd_backend, t_rev)
    @printf("  DI + ForwardDiff:     %d bytes allocated\n", allocs_fwd)

    b_fwd_di = @benchmark value_and_gradient!($test_closure, $grad_buf, $prep_fwd, $fwd_backend, $t_rev) samples=200 evals=10
    @printf("  DI + ForwardDiff:     %8.1f μs  allocs: %d\n",
        median(b_fwd_di).time/1e3, median(b_fwd_di).allocs)

    # Test DI value_and_gradient! with Enzyme Forward
    enz_fwd_backend = AutoEnzyme(mode=set_runtime_activity(Forward), function_annotation=Const)
    prep_enz_fwd = prepare_gradient(test_closure, enz_fwd_backend, zero(t_rev))
    value_and_gradient!(test_closure, grad_buf, prep_enz_fwd, enz_fwd_backend, t_rev) # warmup
    allocs_enz_fwd = @allocated value_and_gradient!(test_closure, grad_buf, prep_enz_fwd, enz_fwd_backend, t_rev)
    @printf("  DI + Enzyme Forward:  %d bytes allocated\n", allocs_enz_fwd)

    b_enz_fwd = @benchmark value_and_gradient!($test_closure, $grad_buf, $prep_enz_fwd, $enz_fwd_backend, $t_rev) samples=200 evals=10
    @printf("  DI + Enzyme Forward:  %8.1f μs  allocs: %d\n",
        median(b_enz_fwd).time/1e3, median(b_enz_fwd).allocs)

    # ═══════════════════════════════════════════════════════════════════════
    #  TEST E: Isolate allocation sources within Enzyme Reverse
    # ═══════════════════════════════════════════════════════════════════════

    println("\n─── TEST E: Isolate allocation sources ───")

    # E1: Just invlink + arr2nt (no orbit solving, no ln_like)
    function closure_just_transform(θ_transformed)
        θ_natural = invlink(θ_transformed)
        θ_system = arr2nt(θ_natural)
        zero(eltype(θ_transformed))
    end
    allocs_transform = @allocated closure_just_transform(t_rev)
    @printf("  Just invlink+arr2nt (no AD): %d bytes\n", allocs_transform)

    prep_transform = prepare_gradient(closure_just_transform, enz_rev_backend, zero(t_rev))
    grad_t = similar(t_rev)
    value_and_gradient!(closure_just_transform, grad_t, prep_transform, enz_rev_backend, t_rev)
    allocs_transform_enz = @allocated value_and_gradient!(closure_just_transform, grad_t, prep_transform, enz_rev_backend, t_rev)
    @printf("  Just invlink+arr2nt (Enzyme Rev): %d bytes\n", allocs_transform_enz)

    # E2: Transform + orbit construction (no solving, no ln_like)
    function closure_with_orbits(θ_transformed)
        θ_natural = invlink(θ_transformed)
        θ_system = arr2nt(θ_natural)
        orbits = ntuple(i -> orbit_constructors[i](θ_system), n_planets)
        zero(eltype(θ_transformed))
    end
    allocs_orbits = @allocated closure_with_orbits(t_rev)
    @printf("  + orbit construction (no AD): %d bytes\n", allocs_orbits)

    prep_orbits = prepare_gradient(closure_with_orbits, enz_rev_backend, zero(t_rev))
    value_and_gradient!(closure_with_orbits, grad_t, prep_orbits, enz_rev_backend, t_rev)
    allocs_orbits_enz = @allocated value_and_gradient!(closure_with_orbits, grad_t, prep_orbits, enz_rev_backend, t_rev)
    @printf("  + orbit construction (Enzyme Rev): %d bytes\n", allocs_orbits_enz)

    # E3: Transform + orbit construction + orbit solving (no ln_like)
    function closure_with_solve(θ_transformed)
        θ_natural = invlink(θ_transformed)
        θ_system = arr2nt(θ_natural)
        orbits = ntuple(i -> orbit_constructors[i](θ_system), n_planets)
        solutions = Octofitter._solve_all_orbits(orbits, epochs)
        zero(eltype(θ_transformed))
    end
    allocs_solve = @allocated closure_with_solve(t_rev)
    @printf("  + orbit solving (no AD): %d bytes\n", allocs_solve)

    prep_solve = prepare_gradient(closure_with_solve, enz_rev_backend, zero(t_rev))
    value_and_gradient!(closure_with_solve, grad_t, prep_solve, enz_rev_backend, t_rev)
    allocs_solve_enz = @allocated value_and_gradient!(closure_with_solve, grad_t, prep_solve, enz_rev_backend, t_rev)
    @printf("  + orbit solving (Enzyme Rev): %d bytes\n", allocs_solve_enz)

    # E4: Full closure (transform + orbits + solve + ln_like)
    allocs_full_enz = @allocated value_and_gradient!(test_closure, grad_buf, prep_rev, enz_rev_backend, t_rev)
    @printf("  Full closure (Enzyme Rev): %d bytes\n", allocs_full_enz)

    println("\n  Incremental allocation cost (Enzyme Reverse):")
    @printf("    invlink+arr2nt:       %d bytes\n", allocs_transform_enz)
    @printf("    + orbit construction: %+d bytes\n", allocs_orbits_enz - allocs_transform_enz)
    @printf("    + orbit solving:      %+d bytes\n", allocs_solve_enz - allocs_orbits_enz)
    @printf("    + ln_like:            %+d bytes\n", allocs_full_enz - allocs_solve_enz)

    # ═══════════════════════════════════════════════════════════════════════
    #  TEST F: Scaling with number of epochs (L)
    # ═══════════════════════════════════════════════════════════════════════

    println("\n─── TEST F: Scaling with number of epochs ───")
    for n_ep in [5, 10, 30, 100]
        rv_tab_scaled = Table(
            epoch = collect(range(50_000.0, 50_500.0, length=n_ep)),
            rv    = randn(MersenneTwister(99), n_ep) .* 50,
            σ_rv  = fill(5.0, n_ep),
        )
        obs_scaled = MarginalizedStarAbsoluteRVObs(rv_tab_scaled, name="ScaledRV")
        sys_scaled = System(
            name = "ScaledDiag",
            companions = [planet_b],
            observations = [obs_scaled],
            variables = @variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            end
        )

        model_scaled = Octofitter.LogDensityModel(sys_scaled; verbosity=0)
        Octofitter.initialize!(model_scaled; verbosity=0)
        t_s = model_scaled.starting_points[1]
        model_scaled.∇ℓπcallback(t_s) # warmup

        b_scaled = @benchmark $(model_scaled.∇ℓπcallback)($t_s) samples=100 evals=5
        @printf("  L=%3d (Enzyme Rev): %8.1f μs  allocs: %d\n",
            n_ep, median(b_scaled).time/1e3, median(b_scaled).allocs)

        # Also test ForwardDiff for same epoch count
        model_scaled_fwd = Octofitter.LogDensityModel(sys_scaled; autodiff=AutoForwardDiff(), verbosity=0)
        Octofitter.initialize!(model_scaled_fwd; verbosity=0)
        t_sf = model_scaled_fwd.starting_points[1]
        model_scaled_fwd.∇ℓπcallback(t_sf) # warmup

        b_scaled_fwd = @benchmark $(model_scaled_fwd.∇ℓπcallback)($t_sf) samples=100 evals=5
        @printf("  L=%3d (ForwardDiff): %8.1f μs  allocs: %d\n",
            n_ep, median(b_scaled_fwd).time/1e3, median(b_scaled_fwd).allocs)
    end
end

# ═══════════════════════════════════════════════════════════════════════
#  SUMMARY
# ═══════════════════════════════════════════════════════════════════════
println("\n" * "=" ^ 72)
println("  SUMMARY — Full pipeline comparison (D=8, L=30)")
println("=" ^ 72)
@printf("  %-35s %10s %10s\n", "Backend", "∇ℓπ (μs)", "allocs")
println("  " * "-" ^ 55)
@printf("  %-35s %10.1f %10d\n", "Enzyme Reverse (current)", median(bA).time/1e3, median(bA).allocs)
@printf("  %-35s %10.1f %10d\n", "ForwardDiff (override)", median(bB).time/1e3, median(bB).allocs)
if bC !== nothing
    @printf("  %-35s %10.1f %10d\n", "Enzyme Forward (override)", median(bC).time/1e3, median(bC).allocs)
else
    @printf("  %-35s %10s %10s\n", "Enzyme Forward (override)", "N/A", "N/A")
end
println("=" ^ 72)

println("\nDiagnostic notes:")
println("  • If Enzyme Reverse >> ForwardDiff: reverse mode overhead dominates at D=8")
println("  • If allocs scale with L: allocations inside ln_like are the bottleneck")
println("  • If Enzyme Forward ≈ ForwardDiff: Enzyme itself is fine, just wrong mode")
println("  • Incremental alloc breakdown shows exactly where the bytes come from")
