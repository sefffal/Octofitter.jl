# [Chains](@id chains)

This page describes the format of the Monte Carlo chains created by DirectDetections.jl


The output of the samplers in DirectDetections is a vector of chains. For example:
```julia
chains, stats = DirectDetections.hmc(system, numwalkers=5, ...)
```

In this example, 5 chains will be returned. `chains[1]` accesses the first chain, `chains[2]` the second, and so on.


Each chain is a ComponentArray from the [ComponentArrays.jl](https://jonniedie.github.io/ComponentArrays.jl/stable/) package. The stucture follows the definition of your model.

You can access them by name or index as follows.
In all of these cases, you will get an array of samples from the posterior. If you used MCMC via `DirectDetections.mcmc` this will be a matrix of N samples × N walkers, unless `squash_=true` was passed to the sampler. For HMC via `DirectDetections.hmc`, it will be a vector.

## Properties of the System
You can access parameters belonging to the system via:
```julia
chains[1].μ
chains[1].plx
```

If you defined a co-planar model, or a model with hierarchical inclination, you can access inclination and longitude of the ascending node in a similar way:
```julia
chains[1].i
chains[1].Ω
```



## Propeties of Planets
You can access parameters belonging to planets by index via:
```julia
chains[1].planets[1].a
chains[1].planets[1].i
chains[1].planets[1].e
# etc
```

Planets are indexed in the order they were passed to create the `System` model.