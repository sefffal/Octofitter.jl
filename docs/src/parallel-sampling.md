# Distributed Sampling


!!! note
    Octofitter's default sampler (Hamiltonian Monte Carlo) is not easily parallelizable; however, it performs excellently on a single core. Give it a try before assuming you need to sample with multiple cores or nodes.


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
