# [Samplers](@id samplers)

Octofitter provides three built-in samplers:
* No U-turn Hamiltonian Monte Carlo (via `octofit`)
* Non-reversible parallel tempered Monte Carlo (via `octofit_pigeons`)
* Rejection sampling (via `octofit_rejection`)

Many additional samplers can be used through the LogDensityProblems.jl interface, but they are not tested.

## Workflow
If the posterior is unimodal (even if it has a complicated shape), go ahead and use AdvancedHMC (`chain = octofit(model)`). This uses a single computer core and is in many cases very efficient.

If the posterior is multimodal, and the modes are quite separated, then use Pigeons (`chain, pt = octofit_pigeons(model, n_rounds=12)`).

For very low-dimensional problems (1--3 parameters), or when you need independent samples, use rejection sampling (`chain = octofit_rejection(model, draws=1_000_000)`).

Read more about these samplers below.


## Hamiltonian Monte Carlo (NUTS)

The recommended choice for almost all problems is Hamiltonian Monte Carlo. It can be run using the `octofit` function:


```julia
chain = octofit(model)
```

!!! note
    Start julia with `julia --threads=auto` to make sure you have multiple threads available. `octofit` is single-threaded, but may calculate the likelihood of your model in parallel if you have many data points (100s or more).

This sampling  method makes use of derivative information, and is much more efficient. This package by default uses the No U-Turn sampler, as implemented in AdvancedHMC.jl.

Derviatives for a complex model are usualy tedious to code, but Octofitter uses ForwardDiff.jl to generate them automatically.

When using HMC, only a few chains are necessary. This is in contrast to Affine Invariant MCMC based packages where hundreds or thousands of walkers are required.
One chain should be enough to cover the whole posterior, but you can run a few different chains to make sure each has converged to the same distribution.


Similarily, fewer samples are required. This is because unlike Affine Invariant MCMC, HMC produces samples that are much less correlated after each step (i.e. the autocorrelation time is much shorter).

Before sampling, run [`initialize!`](@ref) to find good starting points — `octofit` will
do it for you if you haven't, using Pathfinder to warm up the sampler and reduce
convergence times significantly.

The method signature of `octofit` is as follows:
```julia
octofit(
    [rng::Random.AbstractRNG],
    model::Octofitter.LogDensityModel,
    target_accept::Number=0.8;
    adaptation=1000,
    iterations=1000,
    drop_warmup=true,
    max_depth=12,
    verbosity=2,
)
```
The only required argument is `model`.
The two positional arguments are `model`, the model you wish to sample; and `target_accept`, the acceptance rate that should be targeted during windowed adaptation. During this time, the step size and mass matrix will be adapted (see AdvancedHMC.jl for more information). The number of steps taken during adaptation is controlled by `adaptation`, and by default these are dropped; pass `drop_warmup=false` to keep them. The total number of posterior samples produced is given by `iterations`.
`max_depth` controls the maximum tree depth of the sampler. 

!!! note
    Fewer than 1000 adaptation steps will produce a warning. It is there for a reason —
    short adaptation is the most common cause of a badly-scaled mass matrix and a chain
    that looks converged but isn't.

## Pigeons Non-Reversible Parallel Tempering

Pigeons implements non-reversible parallel tempering. You can read more about it here:
[http://pigeons.run](https://pigeons.run/stable/). Pigeons is slower if you only run it on a single (or a few) computer cores, but can scale up very well over many cores or compute nodes. It can reliably sample from multimodal posteriors.

Rather than one chain on the posterior, Pigeons runs a ladder of chains interpolating
between a reference distribution (the prior, or a fitted variational reference) and the
posterior. The reference is easy to move around in, and replicas swap up and down the
ladder, so a mode that HMC could never reach — it can only jump 2–3σ gaps — stays
reachable. That is why the multimodal tutorials in this manual
([likelihood maps](@ref fit-likemap), [proper motion anomaly](@ref fit-pma),
[detection limits](@ref detection-limits), [G23H](@ref fit-g23h)) all sample this way.

Two further reasons to reach for it:

* Its local explorer is a **slice sampler**, which needs only log-density evaluations. So
  it samples models that cannot be differentiated (see below), and models with discrete
  variables, which HMC cannot move at all.
* It estimates the **log evidence ratio** as a by-product, via `Pigeons.stepping_stone(pt)`
  — see [Multi-Planet RV Fits](@ref fit-rv-multi), which uses it to compare a one-planet
  model against a two-planet one.

!!! note "Pigeons must be installed separately"
    `octofit_pigeons`'s methods live in a **package extension**. Run `pkg> add Pigeons`,
    and put `using Pigeons` in scope before calling it — otherwise the function exists
    with no methods and you get a `MethodError`.

Pigeons can be run locally with one or more Julia threads.
!!! note
    Start julia with `julia --threads=auto` to make sure you have multiple threads available for sampling.

You can get started with Pigeons by running:
```julia
using Pigeons
model = Octofitter.LogDensityModel(sys)
chain, pt = octofit_pigeons(model, n_rounds=10)
```

Note the two return values: unlike `octofit`, `octofit_pigeons` returns a named tuple
`(;chain, pt)`. The `chain` is an ordinary `MCMCChains.Chains` you can pass to
[`octoplot`](@ref) and `octocorner`; `pt` is the Pigeons state, which carries the
diagnostics, the evidence estimate, and everything needed to resume.

The method signature of `octofit_pigeons` is as follows:
```julia
chain, pt = octofit_pigeons(
    target::Octofitter.LogDensityModel;
    n_rounds::Int,
    n_chains::Int=16,
    n_chains_variational::Int=16,
    checkpoint::Bool=false,
    multithreaded=true,
    variational=GaussianReference(first_tuning_round=5),
    pigeons_kw... # forwarded to Pigeons.Inputs
)
```

`n_rounds` is the only required keyword. Each round doubles the number of samples, so
`n_rounds=10` yields 1024 draws — increase it until `log(Z₁/Z₀)` in the report stops
drifting.

By default, this will use:
* 16 chains between the posterior and the prior
* 16 chains between the posterior and a variational reference
* the `SliceSampler` local explorer

The number of chains should ideally be set to twice the value of `Λ` in the resulting table.
If you notice `Λ` is not approximately 8, you should adjust `n_chains` and `n_chains_variational` to be approximately twice the value of `Λ` and `Λ_var` respectively.

Pass `explorer=SliceSampler()` (or another Pigeons explorer) to override the local move,
and `n_chains_variational=0, variational=nothing` to drop the variational leg — worth
doing when the posterior is too structured for a Gaussian reference to help. Models with
discrete variables disable the variational reference automatically.

A nice feature of Pigeons is that you can resume sampling for additional rounds without having to start over:
```julia
pt = increment_n_rounds!(pt, 1)
chain, pt = octofit_pigeons(pt)
```

## Rejection Sampling

Rejection sampling is the simplest sampling method. It draws samples from the prior and accepts or rejects each one based on the likelihood. Accepted samples are independent (no autocorrelation), making diagnostics straightforward. However, it can be very inefficient for high-dimensional problems or when the posterior is much narrower than the prior.

Rejection sampling is a good choice when:
* Your model has very few free parameters (1--3 dimensions)
* You want independent, uncorrelated posterior samples
* You want a quick sanity check before running a longer HMC chain
* Gradient-based sampling is not possible (e.g. discrete parameters)

```julia
chain = octofit_rejection(model; draws=1_000_000)
```

The method signature of `octofit_rejection` is as follows:
```julia
octofit_rejection(
    [rng::Random.AbstractRNG],
    model::Octofitter.LogDensityModel;
    draws=100_000,
    verbosity=2,
)
```

The `draws` parameter controls how many prior samples are drawn. The number of accepted posterior samples depends on the acceptance rate, which is reported after sampling. If the acceptance rate is very low, consider using `octofit` (HMC) instead.

## Models that cannot be differentiated

A few likelihood terms are not compatible with forward-mode automatic differentiation,
and a gradient-based sampler cannot move a model containing one:

* The Gaia RV variability channel of [`G23HObs`](@ref) (`Distributions.NoncentralChisq`
  does not accept `ForwardDiff.Dual`). Pass `include_rv=false`, which is what
  [`HGCAObs`](@ref) does for you.
* Gaussian processes evaluated with the Celerite backend, which is `Float64`-only.

The symptom is characteristic: the log density evaluates to a finite number while the
gradient comes back `-Inf` with an all-zero gradient vector, so the sampler cannot take a
step. Use a derivative-free sampler — `octofit_pigeons`, whose slice-sampler explorer only
ever evaluates the log density, or `octofit_rejection` — or disable the offending term.
The [RV + GP tutorial](@ref fit-rv-gp) takes the first route for its Celerite model.

## Sampling across a cluster

See [Distributed Sampling](@ref) for guidance on running `octofit_pigeons` across multiple
nodes of a cluster with MPI.

## Advanced Usage: Additional Samplers
This section is for people interested in developing support for new samplers with Octofitter.

Octofitter converts your model specification into an `Octofitter.LogDensityModel` which implements the [LogDensityProblems.jl interface](https://www.tamaspapp.eu/LogDensityProblems.jl/dev/).

That way, you can sample from your model using a wide variety of Julia based samplers.
These samplers may return results in less convenient formats, and for example, may need you to map their results back to the natural domain of your variables using `model.link` or `model.invlink`.

For convenience, Octofitter bundles special support for the No U-Turn Sampler (NUTS) as implemented by AdvancedHMC.jl (see above).

In order to use the results of most other samplers, you will need a function to map
results from their transformed variables back to their natural domain and reconstruct the chains:

```julia
# Results are in normalized parameter space and need to be mapped back to their constrained support

# Function to map the samples back to their natural domain
function remapchain(mc_samples)
    logpost = map(s->s.lp, mc_samples)
    # Transform samples back to constrained support
    samples = map(mc_samples) do s
        θ_t = s.params
        θ = model.invlink(θ_t)
        return θ
    end
    chain_res = model.arr2nt.(samples)
    chain = Octofitter.result2mcmcchain(chain_res)
    return MCMCChains.setinfo(
        chain,
        (;
            # start_time,
            # stop_time,
            model=model.system,
            logpost=logpost,
        )
    )
end
```

### AdvancedMH
Here is an example of using a separate package to sample from a model---in this case, AdvancedMH. For other packages, see their documentation for full details.

Note: this sampler does not work well and is just provided as a reference for how to use an arbitrary sampling package.

```julia
using AdvancedMH
using MCMCChains: Chains

# Construct model from a system (see elsewhere in docs)
model = Octofitter.LogDensityModel(sys)

# Set up a random walk sampler with a joint multivariate Normal proposal.
using LinearAlgebra
spl = RWMH(MvNormal(zeros(model.D), I))

# Find initial guess by drawing from priors
initial_θ = Octofitter.guess_starting_position(model,50_000)[1]
initial_θ_t = model.link(initial_θ) # Map to unconstrainted parameterization

# Sample from the posterior.
chn_norm = sample(
    model,
    spl,
    1_000_000;
    chain_type=Any,
    init_params=initial_θ_t
)

chn_mh = remapchain(chn_norm)
```

### Emcee (affine invariant sampler)

!!! warning
    This example is under construction

We can use the AdvancedMH package to implement a sampler that is similar to emcee.py. This might be helpful for reproducing the results of packages like orbitize! that are based on this sampler, but is not recommended otherwise.

```julia
using AdvancedMH
using MCMCChains: MCMCChains, Chains, chainscat

# Construct model from a system (see elsewhere in docs)
model = Octofitter.LogDensityModel(sys)

initial_θ = Octofitter.guess_starting_position(model,50_000)[1]
initial_θ_t = model.link(initial_θ) # Map to unconstrainted parameterization


using LinearAlgebra

# Set up our sampler with a joint multivariate Normal proposal.
spl = Ensemble(1_000, StretchProposal(MvNormal(zeros(model.D), I)))

# Sample from the posterior.
start_time = time()
chn_raw = sample(
    model,
    spl,
    1_000;
    chain_type=Any,
    init_params=initial_θ_t
)
stop_time = time()
# Results are in normalized parameter space and need to be mapped back to their constrained support

# Function to map the samples from all walkers back to their natural domain
function remapchain(mc_samples_by_chain)
    chains = map(mc_samples_by_chain) do mc_samples
        logpost = map(s->s.lp, mc_samples)
        # Transform samples back to constrained support
        samples = map(mc_samples) do s
            θ_t = s.params
            θ = model.invlink(θ_t)
            return θ
        end
        chain_res = model.arr2nt.(samples)
        chain = Octofitter.result2mcmcchain(chain_res)
        return MCMCChains.setinfo(
            chain,
            (;
                start_time,
                stop_time,
                model=model.system,
                logpost=logpost,
            )
        )
    end
    chainscat(chains...)
end
# Remap back to natural domain 
chn_all = remapchain(chn_raw)
# Discard some burn in
chn = chn_all[500:end,:,:];
```

### Tempering
The package MCMCTempering can be used to temper most Julia MCMC samplers.
Here is an example with AdvancedMH. 

!!! note
    MCMCTempering is under active development. The API might evolve, and you may
    need to ensure you're using the latest `#main` branch rather than published
    release.


```julia
using MCMCTempering, AdvancedMH, MCMCChains
MCMCTempering.getparams(transition::AdvancedMH.Transition) = transition.params

tempered_sampler = tempered(sampler, 25);

# Sample from the posterior.
chn_norm = sample(
    model, tempered_sampler, 1_000_000;
    discard_initial=100_000, chain_type=Any,
    init_params=initial_θ_t
)


chn = remapchain(chn_norm)
```


### Customized NUTS Sampling
This example shows how to customize different aspects of the default NUTS
sampler.
```julia
using AdvancedHMC
initial_θ = Octofitter.guess_starting_position(model,150_000)[1]
initial_θ_t = model.link(initial_θ)
metric = DenseEuclideanMetric(model.D)
hamiltonian = Hamiltonian(metric, model)
ϵ = find_good_stepsize(hamiltonian, initial_θ_t)

integrator = JitteredLeapfrog(ϵ, 0.1) # 10% normal distribution on step size to help in areas of high curvature. 
# integrator = Leapfrog(ϵ)
# κ = NUTS{MultinomialTS,GeneralisedNoUTurn}(integrator, max_depth=12)
κ = NUTS{SliceTS,GeneralisedNoUTurn}(integrator, max_depth=12)

mma = MassMatrixAdaptor(metric)
ssa = StepSizeAdaptor(0.75, integrator)
adaptor = StanHMCAdaptor(mma, ssa) 
sampler = AdvancedHMC.HMCSampler(κ, metric, adaptor)

# Sample from the posterior.
chn_norm = sample(
    # model, tempered_sampler, 500;
    model, sampler, 500,
    nadapts = 250,
    discard_initial=250, chain_type=Any,
    init_params=initial_θ_t
)

function remapchain(mc_samples::AbstractArray{<:AdvancedHMC.Transition})
    stat = map(s->s.stat, mc_samples)
    logpost = map(s->s.z.ℓπ.value, mc_samples)
    
    mean_accept = mean(getproperty.(stat, :acceptance_rate))
    ratio_divergent_transitions = mean(getproperty.(stat, :numerical_error))
    mean_tree_depth = mean(getproperty.(stat, :tree_depth))

    println("""
    Sampling report for chain:
    mean_accept         = $mean_accept
    ratio_divergent_transitions        = $ratio_divergent_transitions
    mean_tree_depth     = $mean_tree_depth\
    """)

    # Report some warnings if sampling did not work well
    if ratio_divergent_transitions == 1.0
        @error "Numerical errors encountered in ALL iterations. Check model and priors."
    elseif ratio_divergent_transitions > 0.1
        @warn "Numerical errors encountered in more than 10% of iterations" ratio_divergent_transitions
    end
    # Transform samples back to constrained support
    samples = map(mc_samples) do s
        θ_t = s.z.θ
        θ = model.invlink(θ_t)
        return θ
    end
    chain_res = model.arr2nt.(samples)
    chain = Octofitter.result2mcmcchain(chain_res)
    return MCMCChains.setinfo(
        chain,
        (;
            # start_time,
            # stop_time,
            model=model.system,
            logpost=logpost,
        )
    )
end

chn = remapchain(chn_norm)

```
