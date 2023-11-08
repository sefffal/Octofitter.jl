# Parallel Sampling

!!! warn
    Parallel sampling broke after an update to AdvancedHMC.jl. This page will be updated with new instructions in the future.

You can sample from multiple chains in parallel using either threads (useful for quick tests) processes (high throughput), or simply running multiple copies of the program and merging the results afterwards (recommended for clusters).

## Threads

You can sample from multiple threads by passing `Octofitter.MCMCThreads` as the ensemble type.
Specify `num_chains` greater than 1.
You may also want to reduce the verbosity (perhaps to 0) as the overlapping log messages can be very confusing.
```julia
chain = Octofitter.advancedhmc(
    model, 0.85,
    Octofitter.MCMCThreads();
    num_chains=4,
    adaptation =   500,
    iterations =  1000,
    verbosity = 0,
    tree_depth = 15
)
```

!!! note
    Chains sampled in parallel using threads by default begin with the same initial parameter values. In practice, they  rapidly diffuse duing the adapataion phase; hoever, it may be worth running at least one separate chain to ensure your all your results aren't stuck in one tight local maximum.

## Independent Copies
The best way to scale up Octofitter to run on a larger cluster is to write out your model as a script that samples a single chain. Then, write out the chain using e.g. `CSV.write(Table(chain))` at the end.
Simply start the script N times to sample from N independent chains.
Finally at the end, concatenate your chains together.

This is efficient because it allows chains that finish early to release resources back to the cluster. In many cluster environments, you can for example prepare 500 copies of the script and allow only 50 of them to run simultaneously.

## Processes

!!! danger
    Distributed.jl based multi-processing is not currently functional within Octofitter due to an interaction with its runtime code generation.


You can sample from multiple worker processes using Distributed.jl and by passing `Octofitter.MCMCDistributed` as the ensemble type.
Make sure Octofitter is loaded on all worker processes, and the model is also defined all on processes using e.g. the `@everywhere` macro.
Specify `num_chains` greater than 1.
You may also want to reduce the verbosity (perhaps to 0) as the overlapping log messages can be very confusing.

AdvancedHMC.jl has some internal array allocations. With a high number of workers, this can lead to contention between threads due to Julia's stop-the-world garbage collector. On the other hand, using worker processes via Distributed.jl prevents this overhead at the expense of more communcation overhead at the start and end of sampling and more memory utilization.

```julia
chain = octofit(
    model, 0.85,
    Octofitter.MCMCDistributed();
    num_chains=8,
    adaptation =   500,
    iterations =  1000,
    verbosity = 2,
    tree_depth = 15
)
```