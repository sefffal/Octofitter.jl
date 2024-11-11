# [Samplers](@id samplers)

We recommend using one of the following MCMC samplers:
* No U-turn Hamiltonian Monte Carlo (via `octofit`)
* Non-reversible parallel tempered Monte Carlo  (via `octofit_pigeons`)

Many additional samplers can be used through the LogDensityProblems.jl interface, but they are not tested.

## Workflow
When you're testing a new model and/or data, we recommend you test it quickly with Pathfinder (`chains = octoquick(model)`). This will return a rough approximation of the posterior and will pick up if it contains multiple modes. 

If the posterior is unimodal (even if it has a complicated shape), go ahead and use AdvancedHMC (`chains = octofit(model)`). This uses a single computer core and is in many cases very efficient.

If the posterior is multimodal, and the modes are quite separated, then use Pigeons (`chains, pt = octofit_pigeons(model, n_rounds=12)`).

Read mode about these samplers below.


## Pathfinder
You can use the function `octoquick` to generate a very rough approximation of the posterior. This uses the multi-pathfinder approximate inference algorithm.

The useage of `octoquick` is similar to `octofit`:
```julia
chain = octoquick(model)
```

These results are not statistically meaningful, but should give you some very rough idea of how the model fits the data in just a few seconds.


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

`octofit` will internally use Pathfinder to warm up the sampler, reducing convergence times signficantly. 


The method signature of `octofit` is as follows:
```julia
octofit(
    [rng::Random.AbstractRNG],
    model::Octofitter.LogDensityModel,
    target_accept::Number=0.8,
    adaptation=1000,
    iterations=1000,
    drop_warmup=true,
    max_depth=12,
    verbosity=2,
)
```
The only required arguments are `model`, `adaptation`, and `iterations`.
The two positional arguments are `model`, the model you wish to sample; and `target_accept`, the acceptance rate that should be targeted during windowed adaptation. During this time, the step size and mass matrix will be adapted (see AdvancedHMC.jl for more information). The number of steps taken during adaptation is controlled by `adaptation`. You can prevent these samples from being dropped by pasing `include_adaptation=false`. The total number of posterior samples produced are given by `iterations`. These include the adaptation steps that may be discarded.
`max_depth` controls the maximum tree depth of the sampler. 

## Pigeons Non-Reversible Parallel Tempering
Pigeons implements non-reversible parallel tempering. You can read more about it here:
[http://pigeons.run](https://pigeons.run/stable/). Pigeons is slower if you only run it on a single (or a few) computer cores, but can scale up very well over many cores or compute nodes. It can reliably sample from multimodal posteriors.

!!! note
   Pigeons must be installed as a separate package install it, run 
    `pkg> add Pigeons`


Pigeons can be run locally with one or more Julia threads.
!!! note
    Start julia with `julia --threads=auto` to make sure you have multiple threads available for sampling.

You can get started with Pigeons by running:
```julia
using Pigeons
model = Octofitter.LogDensityModel(System)
chain, pt = octofit_pigeons(model)
```

The method signature of `octofit_pigeons` is as follows:
```julia
chain, pt = octofit_pigeons(
    target::Octofitter.LogDensityModel;
    n_rounds::Int,
    pigeons_kw... # forwarded to Pigeons.Inputs
)
```

By default, this will use:
* 16 chains between the posterior and the prior
* 16 chains between the posterior and a variational reference
* the SliceSampler local explorer

The number of chains should ideally be set to twice the value of `Λ` in the resulting table.
If you notice `Λ` is not approximately 8, you should adjust `n_chains` and `n_chains_variational` to be approximately twice the value of `Λ` and `Λ_var` respectively.


A nice feature of Pigeons is that you can resume sampler for additional rounds without having to start over:
```julia
pt = increment_n_rounds!(pt, 1)
chain, pt = octofit_pigeons(pt)
```

## Distributed Sampling



This guide shows how you can sample from Octofitter models using a cluster.
If you just want to sample across multiple cores on the same computer, start julia with multiple threads (`julia --threads=auto`) and use `octofit_pigeons`.

If your problem is challenging enough to benefit from parallel sampling across multiple nodes in a cluster, you might consider using Pigeons with MPI by following this guide. 

## MPI Launcher Script

We will use a Julia script to submit the batch job to the cluster. The script will define the model and start the sampling process. The sampler can then run in the background, and you can periodically load the results in from the checkpoint file to examine them after each round of sampling.

Here is an example:
```julia
using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using CairoMakie
using PairPlots
using DataFrames
using Distributions

# Specify your data as usual
astrom_like = PlanetRelAstromLikelihood(
    # Your data here:
    (epoch = 50000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50120, ra = -502.570356287689, dec = -37.47217527025044, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50240, ra = -498.2089148883798, dec = -7.927548139010479, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50360, ra = -492.67768482682357, dec = 21.63557115669823, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50480, ra = -485.9770335870402, dec = 51.147204404903704, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50600, ra = -478.1095526888573, dec = 80.53589069730698, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50720, ra = -469.0801731788123, dec = 109.72870493064629, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50840, ra = -458.89628893460525, dec = 138.65128697876773, σ_ra = 10, σ_dec = 10, cor=0),
)

# build your model as usual
@planet b Visual{KepOrbit} begin
    a ~ Uniform(0, 100) # AU
    e ~ Uniform(0.0, 0.99)
    i ~ Sine() # radians
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50000) # use MJD epoch of your data here!!
end astrom_like
@system Tutoria begin # replace Tutoria with the name of your planetary system
    M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
end b
model = Octofitter.LogDensityModel(Tutoria)
```


## Launcher Script
Use this script to launch your MPI job.


```julia
include("distributed-model.jl")
pt = pigeons(
    target = Pigeons.LazyTarget(MyLazyTarget()),
    record = [traces; round_trip; record_default()],
    on = Pigeons.MPIProcesses(
        n_mpi_processes = n_chains,
        n_threads = 1,
        dependencies = [abspath("distributed-model.jl")]
    ),
    # Pass additional flags to the HPC scheduler here
    # See here for more details: https://pigeons.run/stable/reference/#Pigeons.MPIProcesses
    # add_to_submission = ["#PBS -A my_user_allocation_code"] # pbs
    add_to_submission = [ # slurm
        "#SBATCH --account=my_user_name",
        "#SBATCH --time=24:00:00",
    ],
     # HPC modules to load on each worker
    environment_modules: ["StdEnv/2023", "intel", "openmpi", "julia/1.10", "hdf5"]
)
```


!!! info 
    Don't submit this script to your cluster. Run it on a login node and it will submit the job for you.

## Troubleshooting

If you run into library issues with MPI and/or HDF5, you may need to tell Julia
to use the system provided versions. 

Here is an example that works on AllianceCanada clusters, and may be adaptable to other slurm-based systems:
```julia
using Preferences, HDF5

set_preferences!(
    HDF5,
    "libhdf5" => ENV["EBROOTHDF5"]*"/lib/libhdf5_hl.so",
    "libhdf5_hl" => ENV["EBROOTHDF5"]*"/lib/libhdf5_hl.so",
    force = true
)

modelfname = ARGS[1]
n_proc = parse(Int, ARGS[2])

Pigeons.setup_mpi(
    submission_system = :slurm,
    environment_modules = ["StdEnv/2023", "intel", "openmpi", "julia/1.10", "hdf5"],
    library_name = ENV["EBROOTOPENMPI"]*"/lib/libmpi",
    add_to_submission = [
        "#SBATCH --time=24:00:00",
        "#SBATCH --account=def-account-name",
        "#SBATCH --mem-per-cpu=8g"
    ]
)
println("Setup MPIProcesses")
```

## Examine Results
After one or more sampling rounds have completed, you can run this command to load the results so far for analysis.

```julia

# If still in current session, just pass the `pt` object:
results = Chains(model, pt)

# Else, if the sampling has been running in the background, run:
pt = PT(mpi_run)
model = pt.inputs.target
results = Chains(model, pt)


octocorner(model, results, small=true)
```


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
model = Octofitter.LogDensityModel(system)

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

