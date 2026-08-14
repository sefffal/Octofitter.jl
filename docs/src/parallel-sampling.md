# Distributed Sampling


Octofitter's default sampler (Hamiltonian Monte Carlo) is not easily parallelizable (only a single chain sampled at a time); however, it performs excellently on a single core. Give it a try before assuming you need to sample with multiple cores or nodes. If you have a lot of data (e.g. RV models with 1000 epochs) you can start julia with multiple threads, and `octofit` will multi-thread the kepler solve giving a moderate speedup. 

This guide shows how you can sample from Octofitter models with the Pigeons parallel tempered sampler using many cores, or a cluster.

## Multiple cores on your own computer

The easiest way to use every core of your own machine is the `cores` keyword of
[`octofit_pigeons`](@ref):

```julia
using Octofitter, Pigeons
model = Octofitter.LogDensityModel(sys)
chain, pt = octofit_pigeons(model; n_rounds=12, cores=8)
```

This runs the sampler in `cores` separate worker processes, which for expensive
models (many RV epochs, absolute astrometry, images) is often faster than sampling with threads
Each run spends a minute or two launching the
workers and compiling the model in them before sampling begins, so it pays off
for long fits rather than quick tests — `octofit_pigeons` prints a hint when
one path or the other looks clearly better. Progress appears in your session as
usual, results checkpoint each round under `results/` in the current directory,
and the same `(chain, pt)` comes back.

For models with thousands of epochs, `threads_per_process=2` (or 4) splits the
same core budget into fewer workers that each also thread their Kepler /
trajectory solve: `cores=16, threads_per_process=2` runs 8 workers × 2 threads.
Chains within a worker sample one after another, so only reach for this once
there is already one worker per chain (e.g. `cores=32,
threads_per_process=2` with the default 16 chains).

Alternatively, start Julia with multiple threads (`julia --threads=auto`) and
`octofit_pigeons` will sample one chain per thread with no startup cost — the
better choice for quick runs and cheap models.

!!! note "cores= inside a cluster job"
    `cores=` targets your own machine. It also works inside a single-node
    cluster allocation (`salloc` / `sbatch`), with one caveat: the bundled
    MPI may try to launch the workers through the cluster's own job
    launcher and fail silently. Set `ENV["HYDRA_BOOTSTRAP"] = "fork"`
    (or `export HYDRA_BOOTSTRAP=fork`) first. For sampling *across* nodes,
    use the MPI launcher below instead.

If your problem is challenging enough to benefit from parallel sampling across multiple nodes in a cluster, you might consider using Pigeons with MPI by following the rest of this guide.

## [Threading the trajectory solve](@id kepsolve-threads)

Independently of how many chains you run, the *likelihood itself* can use more
than one thread. PlanetOrbits' `orbitsolve!` takes a `threads=` keyword that
splits the epoch union into contiguous chunks solved on concurrent tasks; the
chunk boundaries are aligned to the SIMD block, so the chunked result is
bit-for-bit identical to the serial one, values and ForwardDiff gradients both.

Octofitter's generated likelihood passes that keyword through, driven by a
single flag:

```julia
Octofitter._kepsolve_use_threads[] = true
```

Most of the time you never touch it, because whoever starts the sampler already
knows whether the threads are spoken for:

* [`octofit`](@ref) (HMC, one chain) turns it **on** after initialization when
  Julia has more than one thread and the model has at least
  `2 * PlanetOrbits.MIN_EPOCHS_PER_TASK` (currently 1024) epochs.
* [`octofit_pigeons`](@ref) sampling with threads turns it **off** — the chains
  are already one per thread.
* `octofit_pigeons(...; cores=N, threads_per_process=T)` with `T > 1` turns it
  **on** inside each worker process.

You have to set it yourself in only one situation: when you launch the workers
yourself, as in the MPI launcher below. It is safe to set unconditionally —
`orbitsolve!` silently runs serial when there are fewer than 512 epochs per
task, and always for the `AHL21` N-body propagator, which has to march through
time in order. So it does nothing for N-body models, and nothing for small
ones.

!!! note
    `_kepsolve_use_threads` is spelled with a leading underscore because it is
    internal: it is a switch on how the likelihood is evaluated, not part of the
    modelling API, and it may move. Setting it never changes the numbers you
    get — only how many cores compute them.

Save this as `threads_setup.jl` next to your launcher; the MPI launcher below
hands it to every rank:

```julia
# threads_setup.jl — included in each MPI rank before the model is deserialized.
using Octofitter

# Thread the Kepler / trajectory solve inside the likelihood. A no-op below
# PlanetOrbits.MIN_EPOCHS_PER_TASK epochs per task and for the AHL21 N-body
# propagator, and bit-identical to the serial solve when it does engage — so
# there is nothing to guard here beyond "did we actually get threads?".
Octofitter._kepsolve_use_threads[] = Threads.nthreads() > 1
```

## The model script

We split this into two files: one that builds the model, and one that submits the
batch job. The sampler then runs in the background, and you can periodically load
the results in from the checkpoint file to examine them after each round of
sampling.

Here is the model script (save it as `distributed-model.jl`):
```julia
using Octofitter
using PlanetOrbits
using CairoMakie
using PairPlots
using DataFrames
using Distributions

# Specify your data as usual
# In practice, this model is way to fast to benefit from MPI -- but we show it as a quick example.
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


## The launcher script

Use this script to launch your MPI job. Save it as `submit.jl` next to
`distributed-model.jl` and `threads_setup.jl`.

Run `Pigeons.setup_mpi(...)` once on the login node first (see
[Troubleshooting](@ref) below): it is what tells Pigeons which scheduler you are
on, which environment modules to load, and which lines to add to the submission
script. Those settings are saved permanently under `~/.pigeons`, and are *not*
arguments of `MPIProcesses`.

```julia
using Octofitter
using OctofitterRadialVelocity   # only if your model uses RV data
using Pigeons

include("distributed-model.jl")  # defines `model`

# --- job knobs ---------------------------------------------------------------
const NROUNDS   = 12          # doubling schedule: the last round is 2^NROUNDS scans
const NCHAINS   = 32          # parallel tempering chains
const NMPI      = NCHAINS     # one MPI rank per chain
const NTHREADS  = 1           # Julia threads per rank (--cpus-per-task)
const WALLTIME  = "24:00:00"
const MEMPERCPU = "8gb"       # --mem-per-cpu

output_dir = @__DIR__

# The explorer Octofitter uses by default. Swap in another Pigeons explorer
# (AutoMALA, AAPS, ...) or your own here.
explorer = Pigeons.SliceSampler()

# Mirrors octofit_pigeons, but with multithreaded=false so that MPI parallelizes
# across PT chains (one chain per rank => its own Julia heap and its own GC) and
# the per-process threads are left to the Kepler solver — see threads_setup.jl.
inputs = Pigeons.Inputs(;
    target = model,
    record = [Pigeons.traces; Pigeons.round_trip; Pigeons.record_default(); Pigeons.index_process],
    multithreaded = false,
    show_report = true,
    n_rounds = NROUNDS,
    n_chains = NCHAINS,
    n_chains_variational = 0,
    variational = nothing,
    checkpoint = true,
    explorer = explorer,
)

@info "Submitting Pigeons MPI job" NROUNDS NCHAINS NMPI NTHREADS WALLTIME MEMPERCPU

# With >1 thread per rank, tell OpenMPI's mpiexec NOT to bind each rank to a
# single core (its default), or all of a rank's threads would be pinned to one
# core. `--bind-to none` lets a rank's threads roam over its --cpus-per-task
# cores. Override with MPIEXEC_BIND if needed.
bind_mode = get(ENV, "MPIEXEC_BIND", NTHREADS > 1 ? "none" : "")
mpiexec_args = isempty(bind_mode) ? `` : `--bind-to $bind_mode`

result = pigeons(
    inputs,
    MPIProcesses(
        n_mpi_processes = NMPI,
        n_threads = NTHREADS,
        walltime = WALLTIME,
        memory = MEMPERCPU,
        mpiexec_args = mpiexec_args,
        # Children need Octofitter loaded to deserialize the model; any file
        # defining types the model refers to (e.g. a custom explorer) must be
        # included so those definitions exist BEFORE the PT arguments are
        # deserialized; threads_setup.jl flips on threaded Kepler solving
        # inside each rank.
        dependencies = [
            Octofitter,
            OctofitterRadialVelocity,
            abspath(joinpath(@__DIR__, "threads_setup.jl")),
        ],
    ),
)

open(joinpath(output_dir, "exec_folder.txt"), "w") do io
    println(io, result.exec_folder)
end
@info "Submitted MPI job" exec_folder=result.exec_folder
println("EXEC_FOLDER=", result.exec_folder)
println("After the job finishes, load with:")
println("  julia --project=. load-results.jl ", result.exec_folder)
```

A few things worth calling out about this pattern:

* **No lazy target.** Earlier versions of this page wrapped the model in a
  `Pigeons.LazyTarget` so each rank rebuilt it. That is no longer needed — pass
  `target = model` directly and Pigeons serializes it to the execution folder,
  as long as every package that defines a piece of the model is listed in
  `dependencies` (`Octofitter` plus whichever `Octofitter*` package supplies
  your observation types). This is exactly what `octofit_pigeons(...; cores=N)`
  does locally.
* **`multithreaded = false` is deliberate.** MPI is doing the parallelism across
  chains; leaving Pigeons' own threading on would have the ranks fight for the
  same cores. Threads within a rank go to the trajectory solve instead — see
  [Threading the trajectory solve](@ref kepsolve-threads).
* **`checkpoint = true` is how results come back.** The job writes a checkpoint
  every round under `results/`; nothing is returned in memory.
* `dependencies` entries that are `Module`s become `using` lines and entries
  that are paths become `include` calls, emitted in order at the top of each
  rank's launch script, before the model is deserialized.
* **The knobs map straight onto the scheduler.** Under SLURM, `MPIProcesses`
  emits `-t WALLTIME`, `--ntasks=NMPI`, `--cpus-per-task=NTHREADS` and
  `--mem-per-cpu=MEMPERCPU`, so the job asks for `NMPI * NTHREADS` cores in
  total. Don't repeat any of those four in `add_to_submission` — they would be
  duplicated directives. Keep `NTHREADS = 1` unless your model has thousands of
  epochs; then 2 or 4 threads per rank buys a faster likelihood, and
  `threads_setup.jl` is what makes the extra threads useful.

While the job runs, `Pigeons.queue_status(result)` and `Pigeons.watch(result)`
report the scheduler's view and the job's stdout/stderr.

!!! info 
    Don't submit this script to your cluster. Run it on a login node and it will submit the job for you.

## Troubleshooting

`Pigeons.setup_mpi` has to be run once, on the login node, before any
`MPIProcesses` submission — otherwise the submission errors with
`call setup_mpi(..) first`. This is also where the scheduler directives
(`add_to_submission`) and `environment_modules` live: they are *not* keywords of
`MPIProcesses` or of `pigeons`.

If you run into library issues with MPI and/or HDF5, you may need to tell Julia
to use the system provided versions. 

Here is an example that works on AllianceCanada clusters, and may be adaptable to other slurm-based systems:
```julia
using Pigeons
using Preferences, HDF5

set_preferences!(
    HDF5,
    "libhdf5" => ENV["EBROOTHDF5"]*"/lib/libhdf5_hl.so",
    "libhdf5_hl" => ENV["EBROOTHDF5"]*"/lib/libhdf5_hl.so",
    force = true
)

Pigeons.setup_mpi(
    submission_system = :slurm,
    environment_modules = ["StdEnv/2023", "intel", "openmpi", "julia/1.10", "hdf5"],
    library_name = ENV["EBROOTOPENMPI"]*"/lib/libmpi",
    # Walltime, ntasks, cpus-per-task and mem-per-cpu are emitted by
    # MPIProcesses from WALLTIME / NMPI / NTHREADS / MEMPERCPU — put only the
    # directives it does not cover here.
    add_to_submission = [
        "#SBATCH --account=def-account-name",
    ]
)
println("Setup MPIProcesses")
```

## Examine Results
After one or more sampling rounds have completed, you can load the results so far
for analysis. Save this as `load-results.jl`:

```julia
using Octofitter
using OctofitterRadialVelocity   # the same packages the model needs
using Pigeons
using CairoMakie, PairPlots

# The launcher printed this and wrote it to exec_folder.txt; it is a path of the
# form results/all/<something>.
exec_folder = isempty(ARGS) ? strip(read("exec_folder.txt", String)) : ARGS[1]

pt = PT(exec_folder)          # loads the latest completed round's checkpoint
model = pt.inputs.target
results = Chains(model, pt)

octocorner(model, results, small=true)
```

If you are still in the session that ran the sampler (e.g. the `cores=` path
above), you already have `pt` and only need the last two lines.
