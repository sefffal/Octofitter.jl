# [Samplers](@id samplers)


## Hamiltonian Monte Carlo (NUTS)

The recommended choice for almost all problems is Hamiltonian Monte Carlo. It can be run using the `Octofitter.advancedhmc` function.
This sampling  method makes use of derivative information, and is much more efficient. This package by default uses the No U-Turn sampler, as implemented in AdvancedHMC.jl.

Derviatives for a complex model are usualy tedious to code, but Octofitter uses ForwardDiff.jl to generate them automatically.

When using HMC, only a few chains are necessary. This is in contrast to Affine Invariant MCMC based packages where hundreds or thousands of walkers are required.
One chain should be enough to cover the whole posterior, but you can run a few different chains to make sure each has converged to the same distribution.

Similarily, fewer samples are required. This is because unlike Affine Invariant MCMC, HMC produces samples that are much less correlated after each step (i.e. the autocorrelation time is much shorter).

### Usage

The method signature of `octofit` is as follows:
```julia
octofit(
    [rng::Random.AbstractRNG],
    model::Octofitter.LogDensityModel,
    target_accept::Number=0.8,
    ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
    adaptation,
    iterations,
    drop_warmup=true,
    max_depth=12,
    initial_samples=50_000,
    initial_parameters=nothing,
    step_size=nothing,
    verbosity=2,
)
```
The only required arguments are `model`, `adaptation`, and `iterations`.
The two positional arguments are `model`, the model you wish to sample; and `target_accept`, the acceptance rate that should be targeted during windowed adaptation. During this time, the step size and mass matrix will be adapted (see AdvancedHMC.jl for more information). The number of steps taken during adaptation is controlled by `adaptation`. You can prevent these samples from being dropped by pasing `include_adaptation=false`. The total number of posterior samples produced are given by `iterations`. These include the adaptation steps that may be discarded.
`tree_depth` controls the maximum tree depth of the sampler. `initial_parameters` is an optional way to pass a starting point for the chain. If you don't pass a default position, one will be selected by drawing `initial_samples` from the priors. The sample with the highest posterior value will be used as the starting point.


## Other Samplers

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

### Nested Sampling
The image likelihood that is implemented in this package is not suitable for use with nested sampling. Other observation types should work fine.


### AdvancedMH
Here is an example of using a separate package to sample from a model---in this case, AdvancedHM. For other packages, see their documentation for full details.

Note: this sampler does not work well and is just provided as a reference for how to use an arbitrary sampling package.

```julia
using AdvancedMH
using MCMCChains: Chains

# Construct model from a system (see elsewhere in docs)
model = Octofitter.LogDensityModel(system)

# Set up a random walk sampler with a joint multivariate Normal proposal.
using LinearAlgebra
spl = RWMH(MvNormal(zeros(model.D), I))

# Find initial guess by drawing from priors
initial_θ = Octofitter.guess_starting_position(model.system,50_000)[1]
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
using MCMCChains: Chains

# Construct model from a system (see elsewhere in docs)
model = Octofitter.LogDensityModel(system)

initial_θ = Octofitter.guess_starting_position(model.system,50_000)[1]
initial_θ_t = model.link(initial_θ) # Map to unconstrainted parameterization


using LinearAlgebra

# Set up our sampler with a joint multivariate Normal proposal.
spl = Ensemble(1_000, StretchProposal(MvNormal(zeros(model.D), I)))

# Sample from the posterior.
chn_norm = sample(
    model,
    spl,
    1_000;
    chain_type=Any,
    init_params=initial_θ_t
)
# Results are in normalized parameter space and need to be mapped back to their constrained support

chn = remapchain(chn_norm)
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
initial_θ = Octofitter.guess_starting_position(model.system,150_000)[1]
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

### Tempered NUTS
Combine the above AdvancedHMC example with the following:

!!! note
    For this sampler to work well, a more careful adapation scheme is likely necessary.

```julia
using AdvancedHMC
using MCMCTempering
MCMCTempering.getparams(transition::AdvancedHMC.Transition) = transition.z.θ

sampler = AdvancedHMC.HMCSampler(κ, metric, adaptor)
tempered_sampler = tempered(sampler, 10);

# Sample from the posterior.
chn_norm = sample(
    model, tempered_sampler, 1000;
    nadapts = 500,
    discard_initial=0, chain_type=Any,
    init_params=initial_θ_t
)
chn_tempered = remapchain(chn_norm)
```

## Hamiltonian Monte Carlo on the GPU
This example shows how one might deploy Octofitter with AdvancedHMC on a GPU.
The No U-Turn Sampler is not well suited since it includes a dynamically chosen path length, but we can use a variant of static HMC.

Further work is likely necessary to ensure all observations are also stored on the GPU.

```julia
using AdvancedHMC, CUDA
initial_θ = Octofitter.guess_starting_position(model.system,150_000)[1]
initial_θ_t = model.link(initial_θ)
metric = DenseEuclideanMetric(model.D)
hamiltonian = Hamiltonian(metric, model)
ϵ = find_good_stepsize(hamiltonian, initial_θ_t)

integrator = JitteredLeapfrog(ϵ, 0.1)
N_steps = 15
κ = HMCKernel(Trajectory{EndPointTS}(integrator, FixedNSteps(N_steps)))

mma = MassMatrixAdaptor(metric)
ssa = StepSizeAdaptor(0.75, integrator)
adaptor = StanHMCAdaptor(mma, ssa) 
sampler = AdvancedHMC.HMCSampler(κ, metric, adaptor)

# Sample from the posterior.
chn_norm = sample(
    # model, tempered_sampler, 500;
    model, sampler, 50_000,
    nadapts = 10_000,
    discard_initial=10_000, chain_type=Any,
    init_params=initial_θ_t
)

function remapchain(mc_samples::AbstractArray{<:AdvancedHMC.Transition})
    stat = map(s->s.stat, mc_samples)
    logpost = map(s->s.z.ℓπ.value, mc_samples)
    
    mean_accept = mean(getproperty.(stat, :acceptance_rate))
    if hasproperty(first(mc_samples), :numerical_error)
        ratio_divergent_transitions = mean(getproperty.(stat, :numerical_error))
    else
        ratio_divergent_transitions = 0
    end
    if hasproperty(first(mc_samples), :tree_depth)
        mean_tree_depth = mean(getproperty.(stat, :tree_depth))
    else
        mean_tree_depth = N_steps
    end

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

```julia
using Plots
plotchains(chn, :B, kind=:astrometry)
```

```julia
plotchains(chn, :B, kind=:radvel)
```


