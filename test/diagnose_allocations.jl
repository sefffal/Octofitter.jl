# Diagnostic: pinpoint allocation sources in primal and gradient callbacks
#
# Tests each step independently to find where type instabilities cause allocations.

using Octofitter
using Distributions
using Random
using Test

println("="^70)
println("Building test model (PlanetRelAstromObs)")
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

# Warm up both callbacks
model.ℓπcallback(theta)
model.ℓπcallback(theta)
model.∇ℓπcallback(theta)
model.∇ℓπcallback(theta)

println("\n" * "="^70)
println("PART 1: Full callback allocations (baseline)")
println("="^70)

allocs_primal = @allocated model.ℓπcallback(theta)
allocs_grad = @allocated model.∇ℓπcallback(theta)
println("  ℓπcallback:   $allocs_primal bytes")
println("  ∇ℓπcallback:  $allocs_grad bytes")

println("\n" * "="^70)
println("PART 2: Primal callback — step-by-step breakdown")
println("="^70)

# Extract the closure captures via the model's stored functions
invlink = model.invlink
arr2nt = model.arr2nt

# Step 1: any(!isfinite, θ)
for _ in 1:3; any(!isfinite, theta); end
a1 = @allocated any(!isfinite, theta)
println("  1. any(!isfinite, θ):       $a1 bytes")

# Step 2: invlink
for _ in 1:3; invlink(theta); end
a2 = @allocated invlink(theta)
println("  2. invlink(θ):              $a2 bytes")

# Step 3: arr2nt
θ_natural = invlink(theta)
for _ in 1:3; arr2nt(θ_natural); end
a3 = @allocated arr2nt(θ_natural)
println("  3. arr2nt(θ_natural):       $a3 bytes")

# Step 4: ln_prior_transformed
# We need to access this via the closure internals. Instead, rebuild it.
ln_prior_transformed = Octofitter.make_ln_prior_transformed(sys)
for _ in 1:3; ln_prior_transformed(θ_natural, true); end
a4 = @allocated ln_prior_transformed(θ_natural, true)
println("  4. ln_prior(θ_natural):     $a4 bytes")

# Step 5: orbit construction
θ_sys = arr2nt(θ_natural)
n_planets = length(sys.planets)
orbit_constructors = ntuple(n_planets) do i
    OT = Octofitter._planet_orbit_type(sys.planets[i])
    let OT=OT, i=i
        (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
    end
end
for _ in 1:3; ntuple(i -> orbit_constructors[i](θ_sys), n_planets); end
a5 = @allocated ntuple(i -> orbit_constructors[i](θ_sys), n_planets)
println("  5. orbit construction:      $a5 bytes")

# Step 6: _sum_term_likelihoods (need to extract tcfgs and tws)
# Build configs matching the model's
active_indices_all, _ = Octofitter._compute_active_indices(sys)
tcfgs = let tcfgs = ()
    term_idx = 0
    for i in 1:n_planets
        planet = sys.planets[i]
        for (i_like, obs) in enumerate(planet.observations)
            term_idx += 1
            active = active_indices_all[term_idx]
            obs_name = hasproperty(obs, :name) ? Octofitter.normalizename(Octofitter.likelihoodname(obs)) : nothing
            epochs = Octofitter._get_obs_epochs(obs)
            tcfg = Octofitter.TermEvalConfig(obs, obs_name, epochs, active, i)
            tcfgs = (tcfgs..., tcfg)
        end
    end
    for (i_obs, obs) in enumerate(sys.observations)
        term_idx += 1
        active = active_indices_all[term_idx]
        obs_name = Octofitter.normalizename(Octofitter.likelihoodname(obs))
        epochs = Octofitter._get_obs_epochs(obs)
        tcfg = Octofitter.TermEvalConfig(obs, obs_name, epochs, active, 0)
        tcfgs = (tcfgs..., tcfg)
    end
    tcfgs
end

mcfg = Octofitter.ModelEvalConfig(invlink, arr2nt, orbit_constructors, n_planets)
primal_tws = Octofitter._build_term_workspaces(mcfg, tcfgs, θ_sys, Float64, D)

orbits = ntuple(i -> orbit_constructors[i](θ_sys), n_planets)
for _ in 1:3; Octofitter._sum_term_likelihoods(θ_sys, orbits, tcfgs, primal_tws); end
a6 = @allocated Octofitter._sum_term_likelihoods(θ_sys, orbits, tcfgs, primal_tws)
println("  6. _sum_term_likelihoods:   $a6 bytes")

# Step 6b: Break down _sum_term_likelihoods components
tcfg1 = first(tcfgs)
tw1 = first(primal_tws)
θ_planet1, θ_obs1 = Octofitter._get_obs_params(tcfg1, θ_sys)
for _ in 1:3; Octofitter._get_obs_params(tcfg1, θ_sys); end
a6a = @allocated Octofitter._get_obs_params(tcfg1, θ_sys)
println("    6a. _get_obs_params:      $a6a bytes")

for _ in 1:3; Octofitter._solve_all_orbits!(tw1.sol_bufs, orbits, tcfg1.epochs); end
a6b = @allocated Octofitter._solve_all_orbits!(tw1.sol_bufs, orbits, tcfg1.epochs)
println("    6b. _solve_all_orbits!:   $a6b bytes")

for _ in 1:3; Octofitter._make_obs_context(tcfg1, θ_sys, θ_planet1, θ_obs1, orbits, tw1.sol_bufs, tw1.obs_workspace); end
a6c = @allocated Octofitter._make_obs_context(tcfg1, θ_sys, θ_planet1, θ_obs1, orbits, tw1.sol_bufs, tw1.obs_workspace)
println("    6c. _make_obs_context:    $a6c bytes")

ctx = Octofitter._make_obs_context(tcfg1, θ_sys, θ_planet1, θ_obs1, orbits, tw1.sol_bufs, tw1.obs_workspace)
for _ in 1:3; Octofitter.ln_like(tcfg1.obs, ctx); end
a6d = @allocated Octofitter.ln_like(tcfg1.obs, ctx)
println("    6d. ln_like:              $a6d bytes")

println("\n  Sum of steps: $(a1+a2+a3+a4+a5+a6) bytes  vs  full: $allocs_primal bytes")

println("\n" * "="^70)
println("PART 3: Type stability checks")
println("="^70)

# Check return types
out_invlink = Core.Compiler.return_type(invlink, typeof((theta,)))
println("  invlink return type:         $(isconcretetype(out_invlink) ? "✓" : "✗") $out_invlink")

out_arr2nt = Core.Compiler.return_type(arr2nt, typeof((θ_natural,)))
println("  arr2nt return type:          $(isconcretetype(out_arr2nt) ? "✓" : "✗") $out_arr2nt")

out_prior = Core.Compiler.return_type(ln_prior_transformed, typeof((θ_natural, true)))
println("  ln_prior return type:        $(isconcretetype(out_prior) ? "✓" : "✗") $out_prior")

out_primal = Core.Compiler.return_type(model.ℓπcallback, typeof((theta,)))
println("  ℓπcallback return type:      $(isconcretetype(out_primal) ? "✓" : "✗") $out_primal")

out_grad = Core.Compiler.return_type(model.∇ℓπcallback, typeof((theta,)))
println("  ∇ℓπcallback return type:     $(isconcretetype(out_grad) ? "✓" : "✗") $out_grad")

# Check _sum_term_likelihoods
out_sum = Core.Compiler.return_type(Octofitter._sum_term_likelihoods, typeof((θ_sys, orbits, tcfgs, primal_tws)))
println("  _sum_term_likelihoods type:  $(isconcretetype(out_sum) ? "✓" : "✗") $out_sum")

println("\n" * "="^70)
println("PART 4: Gradient callback — step-by-step breakdown")
println("="^70)

# We need to access the gradient internals. The easiest way is to reconstruct them.
# Rebuild prior evaluator
prior_evaluator = Octofitter.PriorEvaluator(invlink, ln_prior_transformed)
prior_backend = Octofitter._default_enzyme_backend

using DifferentiationInterface
θ_zero = zero(theta)
prior_prep = prepare_gradient(prior_evaluator, prior_backend, θ_zero)
prior_grad = similar(theta)
total_grad = similar(theta)

# Test prior evaluator
for _ in 1:3; prior_evaluator(theta); end
a_pe = @allocated prior_evaluator(theta)
println("  PriorEvaluator(θ):          $a_pe bytes")

# Test prior gradient
for _ in 1:3; value_and_gradient!(prior_evaluator, prior_grad, prior_prep, prior_backend, theta); end
a_pg = @allocated value_and_gradient!(prior_evaluator, prior_grad, prior_prep, prior_backend, theta)
println("  prior value_and_gradient!:  $a_pg bytes")

# Build term evaluators and specs (matching what logdensitymodel.jl does)
term_backends_tup = ()
term_evaluators_tup = ()
for term_idx in 1:length(tcfgs)
    active = active_indices_all[term_idx]
    D_active = length(active)
    obs_backend = Octofitter.resolve_ad_backend(tcfgs[term_idx].obs, nothing, D_active)
    all_active = D_active == D
    tw = Octofitter._make_term_workspace(mcfg, tcfgs[term_idx], θ_sys, Float64, D)
    obs_backend = Octofitter._ensure_duplicated(obs_backend)
    if all_active
        evaluator = Octofitter.TermEvaluator(mcfg, tcfgs[term_idx], tw)
    else
        evaluator = Octofitter.SparseTermEvaluator(mcfg, tcfgs[term_idx], tw)
    end
    global term_evaluators_tup = (term_evaluators_tup..., evaluator)
    global term_backends_tup = (term_backends_tup..., obs_backend)
end

term_specs = Octofitter._prepare_term_specs(term_evaluators_tup, term_backends_tup, active_indices_all, theta, theta)

# Test individual term evaluator
eval1 = first(term_evaluators_tup)
for _ in 1:3; eval1(theta); end
a_te = @allocated eval1(theta)
println("  TermEvaluator(θ):           $a_te bytes")

# Test individual term gradient
spec1 = first(term_specs)
θ_active1 = theta[spec1.active_indices]
for _ in 1:3
    if spec1.all_active
        value_and_gradient!(spec1.evaluator, spec1.grad_buf, spec1.prep, spec1.backend, theta)
    else
        copyto!(spec1.θ_active_buf, θ_active1)
        value_and_gradient!(spec1.evaluator, spec1.grad_buf, spec1.prep, spec1.backend,
                           spec1.θ_active_buf, DifferentiationInterface.Constant(theta))
    end
end
a_tg = if spec1.all_active
    @allocated value_and_gradient!(spec1.evaluator, spec1.grad_buf, spec1.prep, spec1.backend, theta)
else
    copyto!(spec1.θ_active_buf, θ_active1)
    @allocated value_and_gradient!(spec1.evaluator, spec1.grad_buf, spec1.prep, spec1.backend,
                                   spec1.θ_active_buf, DifferentiationInterface.Constant(theta))
end
println("  term value_and_gradient!:   $a_tg bytes")

# Test _accumulate_term_gradients!
θ_full_base = copy(theta)
fill!(total_grad, 0.0)
for _ in 1:3; Octofitter._accumulate_term_gradients!(total_grad, 0.0, term_specs, theta, θ_full_base); end
a_acc = @allocated Octofitter._accumulate_term_gradients!(total_grad, 0.0, term_specs, theta, θ_full_base)
println("  _accumulate_term_grads!:    $a_acc bytes")

# Check type stability of term evaluator
out_te = Core.Compiler.return_type(eval1, typeof((theta,)))
println("\n  TermEvaluator return type:   $(isconcretetype(out_te) ? "✓" : "✗") $out_te")

out_pe = Core.Compiler.return_type(prior_evaluator, typeof((theta,)))
println("  PriorEvaluator return type:  $(isconcretetype(out_pe) ? "✓" : "✗") $out_pe")

println("\n" * "="^70)
println("SUMMARY")
println("="^70)
println("Full primal:    $allocs_primal bytes")
println("Full gradient:  $allocs_grad bytes")
