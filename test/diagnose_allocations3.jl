# Diagnostic part 3: Identify what in Enzyme is allocating
#
# Test hypothesis: runtime_activity mode causes allocations
# Test hypothesis: RuntimeGeneratedFunctions cause allocations
# Test hypothesis: ntuple inside Enzyme causes allocations

using Octofitter
using Distributions
using Random
using Enzyme

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

invlink = model.invlink
arr2nt = model.arr2nt

println("\n" * "="^70)
println("TEST 1: Enzyme on a trivial function (baseline)")
println("="^70)

# Trivial: just sum the vector
f_trivial(x) = sum(x)
grad_trivial = zeros(D)
for _ in 1:3
    fill!(grad_trivial, 0.0)
    Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(f_trivial), Enzyme.Active, Enzyme.Duplicated(theta, grad_trivial))
end
fill!(grad_trivial, 0.0)
a1 = @allocated Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(f_trivial), Enzyme.Active, Enzyme.Duplicated(theta, grad_trivial))
println("  Trivial (no runtime_activity): $a1 bytes")

for _ in 1:3
    fill!(grad_trivial, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_trivial), Enzyme.Active, Enzyme.Duplicated(theta, grad_trivial))
end
fill!(grad_trivial, 0.0)
a1b = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_trivial), Enzyme.Active, Enzyme.Duplicated(theta, grad_trivial))
println("  Trivial (runtime_activity):    $a1b bytes")

println("\n" * "="^70)
println("TEST 2: Enzyme on invlink (RuntimeGeneratedFunction)")
println("="^70)

# Test invlink through Enzyme
f_invlink(x) = sum(invlink(x))
grad_invlink = zeros(D)
for _ in 1:3
    fill!(grad_invlink, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_invlink), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_invlink))
end
fill!(grad_invlink, 0.0)
a2 = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_invlink), Enzyme.Active, Enzyme.Duplicated(theta, grad_invlink))
println("  sum(invlink(x)):               $a2 bytes")

println("\n" * "="^70)
println("TEST 3: Enzyme on invlink+arr2nt")
println("="^70)

f_invlink_arr2nt(x) = let nt = arr2nt(invlink(x)); nt.M + nt.plx; end
grad_ia = zeros(D)
for _ in 1:3
    fill!(grad_ia, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_invlink_arr2nt), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_ia))
end
fill!(grad_ia, 0.0)
a3 = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_invlink_arr2nt), Enzyme.Active, Enzyme.Duplicated(theta, grad_ia))
println("  invlink+arr2nt (scalars):      $a3 bytes")

println("\n" * "="^70)
println("TEST 4: Enzyme on orbit construction")
println("="^70)

n_planets = length(sys.planets)
orbit_constructors = ntuple(n_planets) do i
    OT = Octofitter._planet_orbit_type(sys.planets[i])
    let OT=OT, i=i
        (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
    end
end

f_orbit(x) = let
    θn = invlink(x)
    θs = arr2nt(θn)
    orbit = orbit_constructors[1](θs)
    PlanetOrbits.semimajoraxis(orbit)
end

grad_orbit = zeros(D)
for _ in 1:3
    fill!(grad_orbit, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_orbit), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_orbit))
end
fill!(grad_orbit, 0.0)
a4 = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(f_orbit), Enzyme.Active, Enzyme.Duplicated(theta, grad_orbit))
println("  orbit construction:            $a4 bytes")

println("\n" * "="^70)
println("TEST 5: Enzyme on orbit construction + solve")
println("="^70)

θ_natural = invlink(theta)
θ_sys = arr2nt(θ_natural)
epochs = collect(Float64, astrom.table.epoch)
orbits_ex = ntuple(i -> orbit_constructors[i](θ_sys), n_planets)
sol0 = PlanetOrbits.orbitsolve(orbits_ex[1], first(epochs))
sol_bufs = (Vector{typeof(sol0)}(undef, length(epochs)),)

# Function that does invlink + arr2nt + orbit + solve, writes to pre-allocated buffers
struct SolveEval{TSolBufs, TOC, TEpochs, TInv, TArr}
    sol_bufs::TSolBufs
    orbit_constructors::TOC
    epochs::TEpochs
    invlink::TInv
    arr2nt::TArr
end
function (se::SolveEval)(x)
    θn = se.invlink(x)
    θs = se.arr2nt(θn)
    orbit = se.orbit_constructors[1](θs)
    for j in eachindex(se.epochs)
        se.sol_bufs[1][j] = PlanetOrbits.orbitsolve(orbit, se.epochs[j])
    end
    return PlanetOrbits.raoff(se.sol_bufs[1][1])
end

solve_eval = SolveEval(sol_bufs, orbit_constructors, epochs, invlink, arr2nt)
sol_bufs_shadow = (Vector{typeof(sol0)}(undef, length(epochs)),)
solve_eval_shadow = SolveEval(sol_bufs_shadow, orbit_constructors, epochs, invlink, arr2nt)

grad_solve = zeros(D)
for _ in 1:3
    fill!(grad_solve, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Duplicated(solve_eval, solve_eval_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(copy(theta), grad_solve))
end
fill!(grad_solve, 0.0)
a5 = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
    Enzyme.Duplicated(solve_eval, solve_eval_shadow),
    Enzyme.Active,
    Enzyme.Duplicated(theta, grad_solve))
println("  orbit + solve:                 $a5 bytes")

println("\n" * "="^70)
println("TEST 6: Full eval_term via Enzyme")
println("="^70)

active_indices_all, _ = Octofitter._compute_active_indices(sys)
tcfgs = let tcfgs = ()
    term_idx = 0
    for i in 1:n_planets
        planet = sys.planets[i]
        for (_, obs) in enumerate(planet.observations)
            term_idx += 1
            active = active_indices_all[term_idx]
            obs_name = hasproperty(obs, :name) ? Octofitter.normalizename(Octofitter.likelihoodname(obs)) : nothing
            obs_epochs = Octofitter._get_obs_epochs(obs)
            tcfg = Octofitter.TermEvalConfig(obs, obs_name, obs_epochs, active, i)
            tcfgs = (tcfgs..., tcfg)
        end
    end
    tcfgs
end

mcfg = Octofitter.ModelEvalConfig(invlink, arr2nt, orbit_constructors, n_planets)
tcfg1 = first(tcfgs)
tw1 = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
evaluator = Octofitter.TermEvaluator(mcfg, tcfg1, tw1)

tw1_shadow = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
evaluator_shadow = Octofitter.TermEvaluator(mcfg, tcfg1, tw1_shadow)

grad_full = zeros(D)
for _ in 1:3
    fill!(grad_full, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Duplicated(evaluator, evaluator_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(copy(theta), grad_full))
end
fill!(grad_full, 0.0)
a6 = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
    Enzyme.Duplicated(evaluator, evaluator_shadow),
    Enzyme.Active,
    Enzyme.Duplicated(theta, grad_full))
println("  Full eval_term:                $a6 bytes")

println("\n" * "="^70)
println("TEST 7: Try without runtime_activity where possible")
println("="^70)

# Prior without runtime_activity
prior_evaluator = Octofitter.PriorEvaluator(invlink, Octofitter.make_ln_prior_transformed(sys))
grad_prior_nora = zeros(D)
try
    for _ in 1:3
        fill!(grad_prior_nora, 0.0)
        Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(prior_evaluator), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_prior_nora))
    end
    fill!(grad_prior_nora, 0.0)
    a7a = @allocated Enzyme.autodiff(Enzyme.Reverse, Enzyme.Const(prior_evaluator), Enzyme.Active, Enzyme.Duplicated(theta, grad_prior_nora))
    println("  Prior (no runtime_activity):   $a7a bytes")
catch e
    println("  Prior (no runtime_activity):   FAILED: $e")
end

# Term without runtime_activity
tw_nora = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
eval_nora = Octofitter.TermEvaluator(mcfg, tcfg1, tw_nora)
tw_nora_shadow = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
eval_nora_shadow = Octofitter.TermEvaluator(mcfg, tcfg1, tw_nora_shadow)
grad_nora = zeros(D)
try
    for _ in 1:3
        fill!(grad_nora, 0.0)
        Enzyme.autodiff(Enzyme.Reverse,
            Enzyme.Duplicated(eval_nora, eval_nora_shadow),
            Enzyme.Active,
            Enzyme.Duplicated(copy(theta), grad_nora))
    end
    fill!(grad_nora, 0.0)
    a7b = @allocated Enzyme.autodiff(Enzyme.Reverse,
        Enzyme.Duplicated(eval_nora, eval_nora_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(theta, grad_nora))
    println("  Term (no runtime_activity):    $a7b bytes")
catch e
    println("  Term (no runtime_activity):    FAILED: $(typeof(e))")
end

println("\n" * "="^70)
println("ALLOCATION BUILDUP SUMMARY")
println("="^70)
println("  Trivial sum:           ~$a1 bytes (baseline)")
println("  + invlink (RGF):       ~$a2 bytes")
println("  + arr2nt:              ~$a3 bytes")
println("  + orbit construction:  ~$a4 bytes")
println("  + orbit solve:         ~$a5 bytes")
println("  Full eval_term:        ~$a6 bytes")
