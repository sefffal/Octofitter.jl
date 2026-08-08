# Distributed Sampling


!!! note
    Octofitter's default sampler (Hamiltonian Monte Carlo) is not easily parallelizable; however, it performs excellently on a single core. Give it a try before assuming you need to sample with multiple cores or nodes.

!!! note "Requires Pigeons"
    Everything on this page depends on [`octofit_pigeons`](@ref), whose methods live in a
    package extension: `pkg> add Pigeons` and `using Pigeons`, or the function exists with
    no methods. The MPI and cluster plumbing shown here has not been re-verified on a
    real cluster recently — the model definition and the API calls are current, but treat
    the launcher scripts as a starting point. See [Samplers](@ref samplers).


This guide shows how you can sample from Octofitter models using many cores, or a cluster.

## Multiple cores on your own computer

The easiest way to use every core of your own machine is the `cores` keyword of
[`octofit_pigeons`](@ref):

```julia
using Octofitter, Pigeons
model = Octofitter.LogDensityModel(sys)
chain, pt = octofit_pigeons(model; n_rounds=12, cores=8)
```

This runs the sampler in `cores` separate worker processes, which for expensive
models (many RV epochs, absolute astrometry, images) is often about twice as
fast as sampling with threads. Each run spends a minute or two launching the
workers and compiling the model in them before sampling begins, so it pays off
for long fits rather than quick tests — `octofit_pigeons` prints a hint when
one path or the other looks clearly better. Progress appears in your session as
usual, results checkpoint each round under `results/` in the current directory,
and the same `(chain, pt)` comes back.

For models with thousands of epochs, `threads_per_process=2` (or 4) splits the
same core budget into fewer workers that each also thread their Kepler /
trajectory solve: `cores=16, threads_per_process=2` runs 8 workers × 2 threads.

Alternatively, start Julia with multiple threads (`julia --threads=auto`) and
`octofit_pigeons` will sample one chain per thread with no startup cost — the
better choice for quick runs and cheap models.

If your problem is challenging enough to benefit from parallel sampling across multiple nodes in a cluster, you might consider using Pigeons with MPI by following the rest of this guide.

## MPI Launcher Script

We will use a Julia script to submit the batch job to the cluster. The script will define the model and start the sampling process. The sampler can then run in the background, and you can periodically load the results in from the checkpoint file to examine them after each round of sampling.

Here is an example (save it as `distributed-model.jl`):
```julia
using Octofitter
using PlanetOrbits
using CairoMakie
using PairPlots
using DataFrames
using Distributions

# Specify your data as usual
astrom_dat = Table(
    epoch = [50000.0, 50120.0, 50240.0, 50360.0, 50480.0, 50600.0, 50720.0, 50840.0],
    ra    = [-505.7637580573554, -502.570356287689, -498.2089148883798, -492.67768482682357,
             -485.9770335870402, -478.1095526888573, -469.0801731788123, -458.89628893460525],
    dec   = [-66.92982418533026, -37.47217527025044, -7.927548139010479, 21.63557115669823,
             51.147204404903704, 80.53589069730698, 109.72870493064629, 138.65128697876773],
    σ_ra  = fill(10.0, 8),
    σ_dec = fill(10.0, 8),
    cor   = zeros(8),
)

# build your model as usual
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.2, 0.1), lower=0.1)   # [M⊙]
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0
        a ~ Uniform(0, 100)     # AU
        e ~ Uniform(0.0, 0.99)
        i ~ Sine()              # radians
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 50000.0         # use an MJD epoch near your data here!!
    end
)

astrom_obs = RelAstromObs(astrom_dat; target=b, ref=A, name="astrom")

sys = System(   # replace Tutoria with the name of your planetary system
    name="Tutoria",
    bodies=[A, b],
    observations=[astrom_obs],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
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
    environment_modules = ["StdEnv/2023", "intel", "openmpi", "julia/1.10", "hdf5"]
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
