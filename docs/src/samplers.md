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
hmc(
    system::System,
    target_accept=0.65;
    numwalkers=1,
    burnin,
    numsamples_perwalker,
    initial_samples=100_000,
    initial_parameters=nothing,
    tree_depth=10,
    include_warmup=false
)
```
The two positional arguments  are `system`, the model you wish to sample; and `target_accept`, the acceptance rate that should be targeted during windowed adaptation. During this time, the step size and mass matrix will be adapted (see AdvancedHMC.jl for more information). The number of steps taken during adaptation is controlled by `burnin`. You can prevent these samples from being dropped by pasing `include_warmup=false`.
`tree_depth` controls the maximum tree depth of the sampler. `initial_parameters` is an optional way to pass a starting point for the chain. If you don't pass a default position, one will be selected by drawing `initial_samples` from the priors. The sample with the highest posterior value will be used as the starting point.

If you want to sample more than one chain in parallel, you can adjust `numwalkers`. This will use multiple threads. Due to garbage collection, it might be more efficient to instead sample one chain in multiple separate Julia processes.

## Affine Invariant MCMC
Affine Invariant MCMC is the sampling method used in the ever-popular `emcee` Python package. It is included here to aid in reproducing results from papers using that method.

This method is very robust, and does not require gradient information to sample the posterior. However, DirectDetections calculates gradients automatically using ForwardDiff.jl.

One limitation of Affine Invariant MCMC is that it performs very poorly on posteriors with high curvature. Unfortuantely, fitting orbits of planets using Keplerian elements often produces just such a posterior. In these situations, the sampling will reduce to a random walk around the mode. This is hard to diagnose, and results in inflated uncertainty estimates.

This can be partially mitigated using tempering, as in the `ptemcee` Python package. 

One application where Affine Invariant MCMC is useful, is if you require completely flat/uniform priors and/or priors that have a sharp discontinuity. These are very challenging for HMC.

The method signature of `DirectDetections.mcmc` is:
```julia
mcmc(
    system::System;
    burnin,
    numwalkers,
    numsamples_perwalker,
    thinning = 1,
    squash = true,
)
```


## Nested Sampling
Support for nested sampling is still work in progress.