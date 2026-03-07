# Diagnostic part 2: Test direct Enzyme calls vs DI to isolate allocation source
#
# If direct Enzyme calls are allocation-free, the DI wrapper is the culprit.

using Octofitter
using Distributions
using Random
using DifferentiationInterface
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
@assert isfinite(ll)

println("\n" * "="^70)
println("TEST 1: Direct Enzyme call for PriorEvaluator")
println("="^70)

invlink = model.invlink
arr2nt = model.arr2nt
ln_prior_transformed = Octofitter.make_ln_prior_transformed(sys)
prior_evaluator = Octofitter.PriorEvaluator(invlink, ln_prior_transformed)

# Direct Enzyme: PriorEvaluator has no mutable state → Const
grad_prior = zeros(D)

# Warm up the direct Enzyme call
for _ in 1:3
    fill!(grad_prior, 0.0)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(prior_evaluator), Enzyme.Active, Enzyme.Duplicated(copy(theta), grad_prior))
end

# Measure
fill!(grad_prior, 0.0)
a_enzyme_prior = @allocated Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), Enzyme.Const(prior_evaluator), Enzyme.Active, Enzyme.Duplicated(theta, grad_prior))
println("  Direct Enzyme prior:    $a_enzyme_prior bytes")

# Compare with DI
prior_backend = Octofitter._default_enzyme_backend
prior_prep = prepare_gradient(prior_evaluator, prior_backend, zero(theta))
prior_grad_di = similar(theta)
for _ in 1:3; value_and_gradient!(prior_evaluator, prior_grad_di, prior_prep, prior_backend, theta); end
a_di_prior = @allocated value_and_gradient!(prior_evaluator, prior_grad_di, prior_prep, prior_backend, theta)
println("  DI prior:               $a_di_prior bytes")

println("\n" * "="^70)
println("TEST 2: Direct Enzyme call for TermEvaluator")
println("="^70)

# Build the term evaluator
n_planets = length(sys.planets)
orbit_constructors = ntuple(n_planets) do i
    OT = Octofitter._planet_orbit_type(sys.planets[i])
    let OT=OT, i=i
        (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
    end
end
mcfg = Octofitter.ModelEvalConfig(invlink, arr2nt, orbit_constructors, n_planets)

active_indices_all, _ = Octofitter._compute_active_indices(sys)
θ_natural = invlink(theta)
θ_sys = arr2nt(θ_natural)

tcfgs = let tcfgs = ()
    term_idx = 0
    for i in 1:n_planets
        planet = sys.planets[i]
        for (_, obs) in enumerate(planet.observations)
            term_idx += 1
            active = active_indices_all[term_idx]
            obs_name = hasproperty(obs, :name) ? Octofitter.normalizename(Octofitter.likelihoodname(obs)) : nothing
            epochs = Octofitter._get_obs_epochs(obs)
            tcfg = Octofitter.TermEvalConfig(obs, obs_name, epochs, active, i)
            tcfgs = (tcfgs..., tcfg)
        end
    end
    tcfgs
end

tcfg1 = first(tcfgs)
tw1 = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
evaluator = Octofitter.TermEvaluator(mcfg, tcfg1, tw1)

# Check: is this term all_active?
D_active = length(active_indices_all[1])
all_active = D_active == D
println("  D = $D, D_active = $D_active, all_active = $all_active")

# Direct Enzyme: TermEvaluator has mutable state → Duplicated
# Create the shadow evaluator (all mutable buffers zeroed)
tw1_shadow = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
evaluator_shadow = Octofitter.TermEvaluator(mcfg, tcfg1, tw1_shadow)

grad_term = zeros(D)

# Warm up
for _ in 1:3
    fill!(grad_term, 0.0)
    # Zero the shadow workspace manually
    for buf in tw1_shadow.sol_bufs
        fill!(buf, first(buf))  # Reset to avoid stale shadow state
    end
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Duplicated(evaluator, evaluator_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(copy(theta), grad_term)
    )
end

# Measure direct Enzyme
fill!(grad_term, 0.0)
a_enzyme_term = @allocated Enzyme.autodiff(
    Enzyme.set_runtime_activity(Enzyme.Reverse),
    Enzyme.Duplicated(evaluator, evaluator_shadow),
    Enzyme.Active,
    Enzyme.Duplicated(theta, grad_term)
)
println("  Direct Enzyme term:     $a_enzyme_term bytes")

# Compare with DI
obs_backend_dup = Octofitter._ensure_duplicated(Octofitter._default_enzyme_backend)
tw1_di = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
evaluator_di = Octofitter.TermEvaluator(mcfg, tcfg1, tw1_di)
prep_di = prepare_gradient(evaluator_di, obs_backend_dup, zero(theta))
grad_di = similar(theta)

for _ in 1:3; value_and_gradient!(evaluator_di, grad_di, prep_di, obs_backend_dup, theta); end
a_di_term = @allocated value_and_gradient!(evaluator_di, grad_di, prep_di, obs_backend_dup, theta)
println("  DI term:                $a_di_term bytes")

println("\n" * "="^70)
println("TEST 3: Check if sparse term has different behavior")
println("="^70)

# Build sparse evaluator (even if all_active, let's test the sparse path)
tw_sparse = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
sparse_eval = Octofitter.SparseTermEvaluator(mcfg, tcfg1, tw_sparse)
active = active_indices_all[1]
θ_active = theta[active]

# Direct Enzyme for sparse
tw_sparse_shadow = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
sparse_eval_shadow = Octofitter.SparseTermEvaluator(mcfg, tcfg1, tw_sparse_shadow)
grad_sparse = zeros(length(active))

for _ in 1:3
    fill!(grad_sparse, 0.0)
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Duplicated(sparse_eval, sparse_eval_shadow),
        Enzyme.Active,
        Enzyme.Duplicated(copy(θ_active), grad_sparse),
        Enzyme.Const(theta)
    )
end

fill!(grad_sparse, 0.0)
a_enzyme_sparse = @allocated Enzyme.autodiff(
    Enzyme.set_runtime_activity(Enzyme.Reverse),
    Enzyme.Duplicated(sparse_eval, sparse_eval_shadow),
    Enzyme.Active,
    Enzyme.Duplicated(θ_active, grad_sparse),
    Enzyme.Const(theta)
)
println("  Direct Enzyme sparse:   $a_enzyme_sparse bytes")

# DI sparse
tw_sparse_di = Octofitter._make_term_workspace(mcfg, tcfg1, θ_sys, Float64, D)
sparse_eval_di = Octofitter.SparseTermEvaluator(mcfg, tcfg1, tw_sparse_di)
prep_sparse_di = prepare_gradient(sparse_eval_di, obs_backend_dup, zero(θ_active),
                                   DifferentiationInterface.Constant(zero(theta)))
grad_sparse_di = similar(θ_active)
θ_active_buf = similar(θ_active)

for _ in 1:3
    copyto!(θ_active_buf, θ_active)
    value_and_gradient!(sparse_eval_di, grad_sparse_di, prep_sparse_di, obs_backend_dup,
                        θ_active_buf, DifferentiationInterface.Constant(theta))
end
copyto!(θ_active_buf, θ_active)
a_di_sparse = @allocated value_and_gradient!(sparse_eval_di, grad_sparse_di, prep_sparse_di, obs_backend_dup,
                                              θ_active_buf, DifferentiationInterface.Constant(theta))
println("  DI sparse:              $a_di_sparse bytes")

println("\n" * "="^70)
println("TEST 4: Primal in a function wrapper")
println("="^70)

# Test if primal allocations are due to top-level scope
function test_primal_allocs(model, theta, n)
    for _ in 1:n
        model.ℓπcallback(theta)
    end
    return @allocated model.ℓπcallback(theta)
end
test_primal_allocs(model, theta, 3)
a_fn = test_primal_allocs(model, theta, 3)
println("  ℓπcallback in function: $a_fn bytes")

function test_grad_allocs(model, theta, n)
    for _ in 1:n
        model.∇ℓπcallback(theta)
    end
    return @allocated model.∇ℓπcallback(theta)
end
test_grad_allocs(model, theta, 3)
a_fn_g = test_grad_allocs(model, theta, 3)
println("  ∇ℓπcallback in function: $a_fn_g bytes")

println("\n" * "="^70)
println("SUMMARY: Direct Enzyme vs DI")
println("="^70)
println("  Prior:  Enzyme=$a_enzyme_prior  DI=$a_di_prior")
println("  Term:   Enzyme=$a_enzyme_term   DI=$a_di_term")
println("  Sparse: Enzyme=$a_enzyme_sparse DI=$a_di_sparse")
