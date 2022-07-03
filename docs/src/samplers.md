# [Samplers](@id samplers)

DirectDetections.jl includes support for three Monte Carlo samplers: Affine Invariant MCMC (KissMCMC.jl), the No U-Turn Sampler (NUTS) which is Hamiltonian Monte Carlo (AdvancedHMC.jl), and finally, nested sampling (NestedSamplers.jl).

## Hamiltonian Monte Carlo

The recommended choice for almost all problems is Hamiltonian Monte Carlo. It can be run using the `DirectDetections.hmc` function.
This sampling
 method makes use of derivative information, and is much more efficient. This package by default uses the No U-Turn sampler, as implemented in AdvancedHMC.jl.

Derviatives for a complex model are usualy tedious to code, but DirectDetections uses ForwardDiff.jl to generate them automatically.

When using HMC, only a few chains are necessary. This is in contrast to Affine Invariant MCMC based packages where hundreds or thousands of walkers are required.
One chain should be enough to cover the whole posterior, but you can run a few different chains to make sure each has converged to the same distribution.

Similarily, fewer samples are required. This is because unlike Affine Invariant MCMC, HMC produces samples that are much less correlated after each step (i.e. the autocorrelation time is much shorter).

### Usage

The method signature of `DirectDetections.hmc` is as follows:
```julia
function hmc(
    rng::Random.AbstractRNG,
    system::System,
    target_accept::Number=0.8,
    ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
    num_chains=1,
    adaptation,
    iterations,
    thinning=1,
    discard_initial=adaptation,
    tree_depth=10,
    initial_samples=50_000,
    initial_parameters=nothing,
    step_size=nothing,
    verbosity=2,
    autodiff=ForwardDiff
)
```
The only required arguments are `system`, `adaptation`, and `iterations`.
The two positional arguments  are `system`, the model you wish to sample; and `target_accept`, the acceptance rate that should be targeted during windowed adaptation. During this time, the step size and mass matrix will be adapted (see AdvancedHMC.jl for more information). The number of steps taken during adaptation is controlled by `adaptation`. You can prevent these samples from being dropped by pasing `include_adaptation=false`. The total number of posterior samples produced are given by `iterations`. These include the adaptation steps that may be discarded.
`tree_depth` controls the maximum tree depth of the sampler. `initial_parameters` is an optional way to pass a starting point for the chain. If you don't pass a default position, one will be selected by drawing `initial_samples` from the priors. The sample with the highest posterior value will be used as the starting point.


## Nested Sampling
Support for nested sampling is still work in progress.