# Diagnostic part 4: Test if avoiding merge/splat in orbit construction reduces Enzyme allocations
#
# Hypothesis: Enzyme allocates because merge(nt1, nt2)... creates complex NamedTuple
# operations that Enzyme can't optimize. Direct field access should be zero-alloc.

using Octofitter
using Distributions
using Random
using Enzyme
using PlanetOrbits

println("="^70)
println("Building test model")
println("="^70)

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
    name="AllocTest",
    companions=[b],
    observations=[],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(20, 1), lower=1)
    end
)

model = Octofitter.LogDensityModel(sys, verbosity=2)
D = model.D
invlink = model.invlink
arr2nt = model.arr2nt

theta, ll = let theta = nothing, ll = -Inf
    for attempt in 0:50
        theta_natural = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
        theta = model.link(theta_natural)
        ll = model.ℓπcallback(theta)
        isfinite(ll) && break
    end
    theta, ll
end
@assert isfinite(ll)

θ_natural = invlink(theta)
θ_sys = arr2nt(θ_natural)

println("\n" * "="^70)
println("TEST 1: Orbit construction approaches")
println("="^70)

# Current approach: merge + splat
f_merge(x) = let
    θn = invlink(x)
    θs = arr2nt(θn)
    orbit = Visual{KepOrbit}(;merge(θs, θs.planets.b)...)
    PlanetOrbits.semimajoraxis(orbit)
end

# Direct field access: no merge/splat
f_direct(x) = let
    θn = invlink(x)
    θs = arr2nt(θn)
    θp = θs.planets.b
    orbit = Visual{KepOrbit}(
        M=θs.M, plx=θs.plx,
        a=θp.a, e=θp.e, ω=θp.ω,
        i=θp.i, Ω=θp.Ω, tp=θp.tp,
        mass=θp.mass
    )
    PlanetOrbits.semimajoraxis(orbit)
end

# Test forward pass
@assert f_merge(theta) ≈ f_direct(theta)

# Enzyme: merge+splat
grad_merge = zeros(D)
for _ in 1:3
    fill!(grad_merge, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_merge), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_merge))
end
fill!(grad_merge, 0.0)
a_merge = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_merge), Enzyme.Active, Enzyme.Duplicated(theta, grad_merge))
println("  merge+splat orbit:    $a_merge bytes")

# Enzyme: direct field access
grad_direct = zeros(D)
for _ in 1:3
    fill!(grad_direct, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_direct), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_direct))
end
fill!(grad_direct, 0.0)
a_direct = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_direct), Enzyme.Active, Enzyme.Duplicated(theta, grad_direct))
println("  direct field access:  $a_direct bytes")

@assert grad_merge ≈ grad_direct "Gradients disagree between merge and direct approaches"
println("  Gradients agree: ✓")

println("\n" * "="^70)
println("TEST 2: Just ln_like, no invlink/arr2nt in Enzyme boundary")
println("="^70)

# Build context objects outside Enzyme boundary, then differentiate ln_like only
epochs = collect(Float64, astrom.table.epoch)
n_planets = 1
orbit_constructors = ntuple(n_planets) do i
    OT = Octofitter._planet_orbit_type(sys.planets[i])
    let OT=OT, i=i
        (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
    end
end
orbits = ntuple(i -> orbit_constructors[i](θ_sys), n_planets)
sol0 = PlanetOrbits.orbitsolve(orbits[1], first(epochs))
sol_bufs = (Vector{typeof(sol0)}(undef, length(epochs)),)
obs_ws = Octofitter.alloc_obs_workspace(astrom, Float64)

# Pre-solve orbits
Octofitter._solve_all_orbits!(sol_bufs, orbits, epochs)

# Build context
θ_planet = θ_sys.planets.b
θ_obs = (;)
ctx = Octofitter.PlanetObservationContext(θ_sys, θ_planet, θ_obs, orbits, sol_bufs, 1, obs_ws)

# ln_like in isolation
for _ in 1:3; Octofitter.ln_like(astrom, ctx); end
a_lnlike = @allocated Octofitter.ln_like(astrom, ctx)
println("  ln_like (primal, no Enzyme): $a_lnlike bytes")

# Now test: can we differentiate ln_like wrt a flat parameter vector
# by having a thin wrapper that constructs orbits and context?

# Wrapper: takes flat params, constructs everything, calls ln_like
struct LnLikeOnly{TObs, TEpochs, TSolBufs, TObsWS}
    obs::TObs
    epochs::TEpochs
    sol_bufs::TSolBufs
    obs_workspace::TObsWS
end
function (f::LnLikeOnly)(x)
    # Manual parameter destructuring (no invlink/arr2nt)
    M = x[1]; plx = x[2]
    a = x[3]; e = x[4]; ω = x[5]; i = x[6]; Ω = x[7]; tp = x[8]; mass = x[9]
    orbit = Visual{KepOrbit}(M=M, plx=plx, a=a, e=e, ω=ω, i=i, Ω=Ω, tp=tp, mass=mass)
    for j in eachindex(f.epochs)
        f.sol_bufs[1][j] = PlanetOrbits.orbitsolve(orbit, f.epochs[j])
    end
    θ_system = (M=M, plx=plx, planets=(b=(a=a, e=e, ω=ω, i=i, Ω=Ω, tp=tp, mass=mass, observations=(;)),),
                observations=(;))
    θ_planet = (a=a, e=e, ω=ω, i=i, Ω=Ω, tp=tp, mass=mass, observations=(;))
    θ_obs = (;)
    ctx = Octofitter.PlanetObservationContext(θ_system, θ_planet, θ_obs, (orbit,), f.sol_bufs, 1, f.obs_workspace)
    return Octofitter.ln_like(f.obs, ctx)
end

# Natural params for the manual evaluator
θ_natural_vec = collect(Float64, θ_natural)

lnlike_only = LnLikeOnly(astrom, epochs, sol_bufs, obs_ws)
sol_bufs_shadow = (Vector{typeof(sol0)}(undef, length(epochs)),)
obs_ws_shadow = Octofitter.alloc_obs_workspace(astrom, Float64)
lnlike_only_shadow = LnLikeOnly(astrom, epochs, sol_bufs_shadow, obs_ws_shadow)

# Verify it gives the same answer
@assert lnlike_only(θ_natural_vec) ≈ Octofitter.ln_like(astrom, ctx) "Forward pass disagrees"
println("  Forward pass agreement: ✓")

grad_lnlike = zeros(D)
for _ in 1:3
    fill!(grad_lnlike, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Duplicated(lnlike_only, lnlike_only_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(copy(θ_natural_vec), grad_lnlike))
end
fill!(grad_lnlike, 0.0)
a_lnlike_enzyme = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
    Enzyme.Duplicated(lnlike_only, lnlike_only_shadow),
    Enzyme.Active,
    Enzyme.Duplicated(θ_natural_vec, grad_lnlike))
println("  ln_like only (Enzyme):       $a_lnlike_enzyme bytes")

# Test without runtime_activity
try
    grad_lnlike_nora = zeros(D)
    for _ in 1:3
        fill!(grad_lnlike_nora, 0.0)
        Enzyme.autodiff(Enzyme.Reverse,
            Enzyme.Duplicated(lnlike_only, lnlike_only_shadow),
            Enzyme.Active,
            Enzyme.Duplicated(copy(θ_natural_vec), grad_lnlike_nora))
    end
    fill!(grad_lnlike_nora, 0.0)
    a_lnlike_nora = @allocated Enzyme.autodiff(Enzyme.Reverse,
        Enzyme.Duplicated(lnlike_only, lnlike_only_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(θ_natural_vec, grad_lnlike_nora))
    println("  ln_like only (no RTA):       $a_lnlike_nora bytes")
catch e
    println("  ln_like only (no RTA):       FAILED: $(typeof(e))")
end

println("\n" * "="^70)
println("TEST 3: Enzyme through invlink (RGF) — does it inherently allocate?")
println("="^70)

# Simple test: sum(invlink(x)), Const
f_sum_invlink(x) = sum(invlink(x))
grad_si = zeros(D)
for _ in 1:3
    fill!(grad_si, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_sum_invlink), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_si))
end
fill!(grad_si, 0.0)
a_si = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_sum_invlink), Enzyme.Active, Enzyme.Duplicated(theta, grad_si))
println("  sum(invlink(x)):             $a_si bytes")

# Without runtime_activity
try
    grad_si2 = zeros(D)
    for _ in 1:3
        fill!(grad_si2, 0.0)
        Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(f_sum_invlink), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_si2))
    end
    fill!(grad_si2, 0.0)
    a_si2 = @allocated Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(f_sum_invlink), Enzyme.Active, Enzyme.Duplicated(theta, grad_si2))
    println("  sum(invlink(x)) no RTA:      $a_si2 bytes")
catch e
    println("  sum(invlink(x)) no RTA:      FAILED: $(typeof(e))")
end

println("\n" * "="^70)
println("SUMMARY: Can we get zero-alloc by narrowing the Enzyme boundary?")
println("="^70)
println("  merge+splat orbit:         $a_merge bytes")
println("  direct field orbit:        $a_direct bytes")
println("  ln_like-only (Enzyme):     $a_lnlike_enzyme bytes")
println("  Full eval_term:            ~9312 bytes (from previous test)")
