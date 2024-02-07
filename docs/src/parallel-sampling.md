# Parallel Sampling


!!! note
    Octofitter's default sampler (Hamiltonian Monte Carlo) is not easily parallelizable; however, it performs excellently on a single core. Give it a try before assuming you need to sample with multiple cores or nodes.


This guide shows how you can sample from Octofitter models using a cluster.
If you just want to sample across multiple cores on the same computer, start julia with multiple threads (`julia --threads=auto`) and use `octofit_pigeons`.


If your problem is challenging enough to benefit from parallel sampling across multiple nodes in a cluster, you might consider using Pigeons with MPI. 

## Model Script
Start by creating a script that only loads your data and defines your model. At the end of the script, add some boilerplate shown below. Here is an example:
```


using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using CairoMakie
using PairPlots
using CSV
using DataFrames
using Distributions

# load your data
astrom_like = PlanetRelAstromLikelihood(
    # Your data here:
    # units are MJD, mas, mas, mas, mas, and correlation.
    (epoch = 50000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50120, ra = -502.570356287689, dec = -37.47217527025044, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50240, ra = -498.2089148883798, dec = -7.927548139010479, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50360, ra = -492.67768482682357, dec = 21.63557115669823, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50480, ra = -485.9770335870402, dec = 51.147204404903704, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50600, ra = -478.1095526888573, dec = 80.53589069730698, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50720, ra = -469.0801731788123, dec = 109.72870493064629, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50840, ra = -458.89628893460525, dec = 138.65128697876773, σ_ra = 10, σ_dec = 10, cor=0),
)

# define your model
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
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b
model = Octofitter.LogDensityModel(Tutoria)

# Copy this boilerplate
struct MyTargetFlag end 
using Pigeons
Pigeons.instantiate_target(flag::MyTargetFlag) = model
```


## Launcher Script

```julia
include("distributed-model.jl")
pt = pigeons(
    target = Pigeons.LazyTarget(MyTargetFlag()),
    record = [traces; round_trip; record_default()],
    on = ChildProcess(
        n_local_mpi_processes = 4,
        dependencies = ["distributed-model.jl"]
    )
)
results = Chains(model, pt)
octocorner(model, results,small=true)
```