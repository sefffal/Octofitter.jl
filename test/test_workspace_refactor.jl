# Test the workspace refactor: model construction, evaluation, and allocation checks
# for already-migrated observation types (PlanetRelAstromObs, GaiaDR4AstromObs)

using Octofitter
using Distributions
using Random
using Test

println("="^60)
println("Test 1: PlanetRelAstromObs model construction + evaluation")
println("="^60)

# Simple relative astrometry data
astrom = PlanetRelAstromLikelihood(
    Table(
        epoch=[50000.0, 50120.0, 50240.0, 50360.0],
        ra=[-505.764, -502.57, -498.209, -492.678],
        dec=[-66.93, -37.47, -7.93, 21.64],
        σ_ra=[10.0, 10.0, 10.0, 10.0],
        σ_dec=[10.0, 10.0, 10.0, 10.0],
    ),
    name="relastrom"
)

b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[astrom],
    variables=@variables begin
        a ~ LogUniform(1, 50)
        e ~ Uniform(0, 0.5)
        ω ~ Uniform(0, 2π)
        i ~ Sine()
        Ω ~ Uniform(0, 2π)
        θ ~ Uniform(0, 2π)
        tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(0.1, 100)
    end
)

sys = System(
    name="WorkspaceTest",
    companions=[b],
    observations=[],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(20, 1), lower=1)
    end
)

println("\nConstructing LogDensityModel...")
model = Octofitter.LogDensityModel(sys, verbosity=2)
println("Model dimension: D = $(model.D)")

# Find a finite starting point
theta, ll = let theta = nothing, ll = -Inf
    for attempt in 0:50
        theta_natural = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
        theta = model.link(theta_natural)
        ll = model.ℓπcallback(theta)
        isfinite(ll) && break
    end
    theta, ll
end
@assert isfinite(ll) "Could not find finite starting point"
println("Found finite starting point: ll = $ll")

# Test gradient evaluation
println("\nTesting gradient evaluation...")
ll_grad, grad = model.∇ℓπcallback(theta)
@assert isfinite(ll_grad) "Gradient evaluation returned non-finite ll"
@assert all(isfinite, grad) "Gradient contains non-finite values"
println("Gradient ll = $ll_grad")
println("Gradient finite: $(all(isfinite, grad))")
println("Max |grad|: $(maximum(abs, grad))")

# Verify primal and gradient agree on ll value
@assert ll ≈ ll_grad "Primal ($ll) and gradient ($ll_grad) disagree on ll value"
println("Primal/gradient ll agreement: $(ll ≈ ll_grad) ✓")

println("\n" * "="^60)
println("Test 2: Allocation check — primal evaluation")
println("="^60)

# Measure allocations inside a function to avoid top-level scope effects
function _measure_allocs(f, args, n_warmup=5)
    for _ in 1:n_warmup; f(args...); end
    @allocated f(args...)
end

allocs_primal = _measure_allocs(model.ℓπcallback, (theta,))
println("Primal allocations: $allocs_primal bytes")

println("\n" * "="^60)
println("Test 3: Allocation check — gradient evaluation")
println("="^60)

allocs_grad = _measure_allocs(model.∇ℓπcallback, (theta,))
println("Gradient allocations: $allocs_grad bytes")

println("\n" * "="^60)
println("Test 4: GaiaDR4AstromObs (Enzyme) model construction + evaluation")
println("="^60)

N_epochs = 10
mock_epochs = collect(range(58000.0, 59000.0, length=N_epochs))
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
    scan_pos_angle = collect(mock_scan_angles),
    centroid_pos_al = zeros(N_epochs),
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

ref_epoch_mjd = 57936.375
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
    false,
    detrend_Δt,
    detrend_inv_N,
    detrend_inv_sum_Δt²,
)

orbit_ref_epoch = mean_epoch
b_dr4 = Planet(
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

sys_dr4 = System(
    name="DR4WorkspaceTest",
    companions=[b_dr4],
    observations=[gaia_obs],
    variables=@variables begin
        M = 1.0
        plx ~ Uniform(10, 30)
    end
)

println("\nConstructing LogDensityModel with GaiaDR4 (Enzyme)...")
model_dr4 = Octofitter.LogDensityModel(sys_dr4, verbosity=2)
println("Model dimension: D = $(model_dr4.D)")

theta_dr4, ll_dr4 = let theta_dr4 = nothing, ll_dr4 = -Inf
    for attempt in 0:50
        theta_natural = collect(model_dr4.sample_priors(Random.Xoshiro(42 + attempt)))
        theta_dr4 = model_dr4.link(theta_natural)
        ll_dr4 = model_dr4.ℓπcallback(theta_dr4)
        isfinite(ll_dr4) && break
    end
    theta_dr4, ll_dr4
end
@assert isfinite(ll_dr4) "Could not find finite starting point for DR4 model"
println("Found finite starting point: ll = $ll_dr4")

println("\nTesting gradient evaluation (Enzyme)...")
ll_grad_dr4, grad_dr4 = model_dr4.∇ℓπcallback(theta_dr4)
@assert isfinite(ll_grad_dr4) "Gradient evaluation returned non-finite ll"
@assert all(isfinite, grad_dr4) "Gradient contains non-finite values"
println("Gradient ll = $ll_grad_dr4")
println("Gradient finite: $(all(isfinite, grad_dr4))")

@assert ll_dr4 ≈ ll_grad_dr4 "Primal/gradient ll disagree"
println("Primal/gradient ll agreement: $(ll_dr4 ≈ ll_grad_dr4) ✓")

println("\n" * "="^60)
println("Test 5: Allocation check — DR4 primal evaluation")
println("="^60)

allocs_dr4_primal = _measure_allocs(model_dr4.ℓπcallback, (theta_dr4,))
println("DR4 primal allocations: $allocs_dr4_primal bytes")

println("\n" * "="^60)
println("Test 6: Allocation check — DR4 gradient evaluation (Enzyme)")
println("="^60)

allocs_dr4_grad = _measure_allocs(model_dr4.∇ℓπcallback, (theta_dr4,))
println("DR4 gradient allocations: $allocs_dr4_grad bytes")

println("\n" * "="^60)
println("Test 7: Isolate allocation sources in primal path")
println("="^60)

# Break down the primal callback to find allocation sources
let model = model, theta = theta
    invlink = model.invlink
    arr2nt = model.arr2nt

    # Test invlink allocation
    invlink(theta)  # warm up
    allocs_invlink = @allocated invlink(theta)
    println("  invlink(θ):        $allocs_invlink bytes")

    θ_natural = invlink(theta)
    arr2nt(θ_natural)  # warm up
    allocs_arr2nt = @allocated arr2nt(θ_natural)
    println("  arr2nt(θ_natural): $allocs_arr2nt bytes")

    println("  (These are the parameter transformation functions outside the AD boundary)")
end

println("\n" * "="^60)
println("Test 8: Isolate eval_term allocation (the inner AD function)")
println("="^60)

# Test that eval_term itself is allocation-free for the RelAstrom model
# We need to extract the workspace and configs
let model = model, theta = theta
    invlink = model.invlink
    arr2nt = model.arr2nt

    θ_natural = invlink(theta)
    θ_sys = arr2nt(θ_natural)

    # Access the primal workspace term configs from the closure
    # Instead, just test _sum_term_likelihoods directly by building what we need
    n_planets = length(model.system.planets)
    orbit_constructors = ntuple(n_planets) do i
        OT = Octofitter._planet_orbit_type(model.system.planets[i])
        let OT=OT, i=i
            (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
        end
    end
    orbits = ntuple(i -> orbit_constructors[i](θ_sys), n_planets)

    # _solve_all_orbits! into pre-allocated buffers (simulating workspace)
    epochs = collect(Float64, model.system.planets[1].observations[1].table.epoch)
    sol0 = PlanetOrbits.orbitsolve(orbits[1], first(epochs))
    sol_bufs = (Vector{typeof(sol0)}(undef, length(epochs)),)

    # Warm up
    Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)
    Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)
    allocs_solve = @allocated Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)
    println("  _solve_all_orbits!:  $allocs_solve bytes")

    # Test ln_like with workspace
    obs = model.system.planets[1].observations[1]
    obs_ws = Octofitter.alloc_obs_workspace(obs, Float64)
    ctx = Octofitter.PlanetObservationContext(θ_sys, θ_sys.planets[1], (;), orbits, sol_bufs, 1, obs_ws)

    Octofitter.ln_like(obs, ctx)
    Octofitter.ln_like(obs, ctx)
    allocs_lnlike = @allocated Octofitter.ln_like(obs, ctx)
    println("  ln_like (with ws):   $allocs_lnlike bytes")
end

println("\n" * "="^60)
println("Test 9: Isolate eval_term allocation for DR4 (Enzyme)")
println("="^60)

let model = model_dr4, theta = theta_dr4
    invlink = model.invlink
    arr2nt = model.arr2nt

    θ_natural = invlink(theta)
    θ_sys = arr2nt(θ_natural)

    n_planets = length(model.system.planets)
    orbit_constructors = ntuple(n_planets) do i
        OT = Octofitter._planet_orbit_type(model.system.planets[i])
        let OT=OT, i=i
            (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
        end
    end
    orbits = ntuple(i -> orbit_constructors[i](θ_sys), n_planets)

    epochs = collect(Float64, gaia_obs.table.epoch)
    sol0 = PlanetOrbits.orbitsolve(orbits[1], first(epochs))
    sol_bufs = (Vector{typeof(sol0)}(undef, length(epochs)),)

    Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)
    Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)
    allocs_solve = @allocated Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)
    println("  _solve_all_orbits!:  $allocs_solve bytes")

    obs = gaia_obs
    obs_ws = Octofitter.alloc_obs_workspace(obs, Float64)
    θ_obs_nt = (;
        astrometric_jitter=0.1,
        ra_offset_mas=0.0, dec_offset_mas=0.0,
        pmra=5.0, pmdec=-10.0,
        plx=θ_sys.plx, ref_epoch=57936.375,
    )
    ctx = Octofitter.SystemObservationContext(θ_sys, θ_obs_nt, orbits, sol_bufs, obs_ws)

    Octofitter.ln_like(obs, ctx)
    Octofitter.ln_like(obs, ctx)
    allocs_lnlike = @allocated Octofitter.ln_like(obs, ctx)
    println("  ln_like (with ws):   $allocs_lnlike bytes")
end

println("\n" * "="^60)
println("SUMMARY")
println("="^60)
println("RelAstrom primal allocations:  $allocs_primal bytes")
println("RelAstrom gradient allocations: $allocs_grad bytes")
println("DR4 primal allocations:         $allocs_dr4_primal bytes")
println("DR4 gradient allocations:       $allocs_dr4_grad bytes")
println()
if allocs_primal == 0
    println("✓ RelAstrom primal: ZERO allocations")
else
    println("✗ RelAstrom primal: $allocs_primal bytes allocated")
end
if allocs_grad == 0
    println("✓ RelAstrom gradient: ZERO allocations")
else
    println("  RelAstrom gradient: $allocs_grad bytes (ForwardDiff DI internals expected)")
end
if allocs_dr4_primal == 0
    println("✓ DR4 primal: ZERO allocations")
else
    println("✗ DR4 primal: $allocs_dr4_primal bytes allocated")
end
if allocs_dr4_grad == 0
    println("✓ DR4 gradient: ZERO allocations")
else
    println("✗ DR4 gradient: $allocs_dr4_grad bytes allocated")
end
