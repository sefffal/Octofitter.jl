# Diagnostic: Find type instabilities that prevent the Julia compiler from
# optimizing merge/splat and other patterns before they hit LLVM/Enzyme.
#
# Uses @code_warntype to check each step of eval_term for type instabilities.

using Octofitter
using Distributions
using Random
using InteractiveUtils

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

n_planets = length(sys.planets)
construct_orbits = Octofitter._make_construct_orbits(sys)
mcfg = Octofitter.ModelEvalConfig(invlink, arr2nt, construct_orbits, Val(n_planets))

active_indices_all, _ = Octofitter._compute_active_indices(sys)
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

println("\n" * "="^70)
println("CHECK 1: @code_warntype eval_term")
println("="^70)
println("\n--- eval_term(θ, mcfg, tcfg, tw) ---\n")
@code_warntype Octofitter.eval_term(theta, mcfg, tcfg1, tw1)

println("\n" * "="^70)
println("CHECK 2: @code_warntype on individual steps")
println("="^70)

println("\n--- invlink(θ) ---\n")
@code_warntype invlink(theta)

println("\n--- arr2nt(θ_natural) ---\n")
@code_warntype arr2nt(θ_natural)

println("\n--- construct_orbits(θ_sys) ---\n")
@code_warntype construct_orbits(θ_sys)

println("\n--- merge(θ_sys, θ_sys.planets.b) ---\n")
@code_warntype merge(θ_sys, θ_sys.planets.b)

println("\n--- _get_obs_params ---\n")
@code_warntype Octofitter._get_obs_params(tcfg1, θ_sys)

orbits = construct_orbits(θ_sys)
println("\n--- _solve_all_orbits! ---\n")
@code_warntype Octofitter._solve_all_orbits!(tw1.sol_bufs, orbits, tcfg1.epochs)

θ_planet1, θ_obs1 = Octofitter._get_obs_params(tcfg1, θ_sys)
Octofitter._solve_all_orbits!(tw1.sol_bufs, orbits, tcfg1.epochs)

println("\n--- _make_obs_context ---\n")
@code_warntype Octofitter._make_obs_context(tcfg1, θ_sys, θ_planet1, θ_obs1, orbits, tw1.sol_bufs, tw1.obs_workspace)

ctx = Octofitter._make_obs_context(tcfg1, θ_sys, θ_planet1, θ_obs1, orbits, tw1.sol_bufs, tw1.obs_workspace)
println("\n--- ln_like ---\n")
@code_warntype Octofitter.ln_like(tcfg1.obs, ctx)

println("\n" * "="^70)
println("CHECK 3: construct_orbits RuntimeGeneratedFunction")
println("="^70)
println("\n--- construct_orbits(θ_sys) ---\n")
@code_warntype construct_orbits(θ_sys)

println("\n" * "="^70)
println("CHECK 4: TermEvaluator call")
println("="^70)
evaluator = Octofitter.TermEvaluator(mcfg)
println("\n--- evaluator(θ) ---\n")
@code_warntype evaluator(theta, tcfg1)

println("\n" * "="^70)
println("CHECK 5: PriorEvaluator call")
println("="^70)
ln_prior_transformed = Octofitter.make_ln_prior_transformed(sys)
prior_eval = Octofitter.PriorEvaluator(invlink, ln_prior_transformed)
println("\n--- prior_evaluator(θ) ---\n")
@code_warntype prior_eval(theta)

println("\n" * "="^70)
println("CHECK 6: Full primal callback")
println("="^70)
println("\n--- ℓπcallback(θ) ---\n")
@code_warntype model.ℓπcallback(theta)

println("\n" * "="^70)
println("CHECK 7: _sum_term_likelihoods")
println("="^70)
primal_tws = Octofitter._build_term_workspaces(mcfg, tcfgs, θ_sys, Float64, D)
println("\n--- _sum_term_likelihoods ---\n")
@code_warntype Octofitter._sum_term_likelihoods(θ_sys, orbits, tcfgs, primal_tws)
