#=
Gaia DR4 completeness trial — one array-job element.

Runs a single injection-recovery trial from the completeness grid defined in
`common.jl`. The job index comes from SLURM_ARRAY_TASK_ID (1-based).

No detection decision is made here: the trial serializes the whole posterior
chain and the true injected parameters, and `assemble_results.jl` applies the
criterion afterwards — so a threshold can be changed without re-sampling.

Usage:
    julia --project=. --threads=auto completeness_trial.jl

Environment variables:
    SLURM_ARRAY_TASK_ID       1-based job index into the completeness grid
    OCTO_COMPLETENESS_QUICK   1 → the reduced grid and short chains (see common.jl)
=#

using Serialization

include(joinpath(@__DIR__, "common.jl"))

# Reads the cache `setup.jl` wrote; queries the archive and GOST only if it is
# missing, which on a compute node without internet will fail loudly rather
# than silently substituting different data.
target = load_target(@__DIR__)
sys = build_system(target)

jobs = completeness_jobs(; masses=MASSES, separations=SEPARATIONS, n_trials=N_TRIALS)

job_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
if job_id > length(jobs)
    @warn "Job index $job_id exceeds total jobs $(length(jobs)); exiting."
    exit(0)
end
job = jobs[job_id]

@info "Starting completeness trial" job_id mass_msol = job.mass mass_mjup = job.mass / mjup separation_au = job.separation trial = job.i_trial

result = run_completeness_trial(job, sys, sampler; inject, add_noise=true, verbosity=1)

outdir = joinpath(@__DIR__, "results")
mkpath(outdir)
outfile = joinpath(outdir, "trial_$(lpad(job_id, 4, '0')).jls")
serialize(outfile, result)

@info "Trial complete" detected = detection(result.chain, result.θ_true) outfile
