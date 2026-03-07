#=
Completeness / sensitivity mapping via injection-recovery.

This module provides a framework for computing detection completeness as a
function of companion mass and separation (or any two parameters). The workflow
supports both local execution and cluster-scale parallelism:

  Phase 1: `completeness_jobs()`       — generate lightweight job descriptions
  Phase 2: `run_completeness_trial()`  — execute a single injection-recovery trial
  Phase 3: `assemble_completeness()`   — combine results into a CompletenessMap

For convenience, `completeness_map()` runs all three phases locally.

NOTE: For efficiency, each trial initializes the sampler at the true injected
parameters rather than running full blind initialization. This dramatically
reduces sampling cost but means the completeness estimate is optimistic about
convergence. The results reflect the *statistical* detectability of a signal,
not the ability to blindly discover it.
=#

# ──────────────────────────────────────────────────────────────────────
# Data structures
# ──────────────────────────────────────────────────────────────────────

"""
    CompletenessJob

A lightweight, serializable description of a single injection-recovery trial.
Contains the grid indices, physical values, and RNG seed — everything needed
to reproduce the trial deterministically.

# Fields
- `i_mass::Int` — index into the mass grid
- `i_sep::Int` — index into the separation grid
- `i_trial::Int` — trial number within this grid cell
- `mass::Float64` — companion mass [Mjup]
- `separation::Float64` — semi-major axis [AU] (or period, depending on usage)
- `seed::UInt64` — RNG seed for reproducibility
"""
struct CompletenessJob
    i_mass::Int
    i_sep::Int
    i_trial::Int
    mass::Float64
    separation::Float64
    seed::UInt64
end
export CompletenessJob

"""
    CompletenessResult

The outcome of a single injection-recovery trial.

# Fields
- `job::CompletenessJob` — the job that produced this result
- `detected::Bool` — whether the detection criterion was satisfied
"""
struct CompletenessResult
    job::CompletenessJob
    detected::Bool
end
export CompletenessResult

"""
    CompletenessMap

Assembled completeness results on a 2D grid of mass × separation.

# Fields
- `masses::Vector{Float64}` — mass grid values [Mjup]
- `separations::Vector{Float64}` — separation grid values [AU]
- `completeness::Matrix{Float64}` — fraction of trials detected (mass × sep)
- `n_detected::Matrix{Int}` — number of detections per cell
- `n_total::Matrix{Int}` — number of trials per cell
"""
struct CompletenessMap
    masses::Vector{Float64}
    separations::Vector{Float64}
    completeness::Matrix{Float64}
    n_detected::Matrix{Int}
    n_total::Matrix{Int}
end
export CompletenessMap


# ──────────────────────────────────────────────────────────────────────
# Phase 1: Generate jobs
# ──────────────────────────────────────────────────────────────────────

"""
    completeness_jobs(; masses, separations, n_trials=5)

Generate a list of [`CompletenessJob`](@ref) descriptions for a completeness grid.

Each job specifies a (mass, separation) grid point and a trial index. Jobs are
independent and can be dispatched to separate processes or cluster nodes.

# Arguments
- `masses` — iterable of companion masses [Mjup]
- `separations` — iterable of semi-major axes [AU]
- `n_trials::Int=5` — number of independent trials per grid cell

# Returns
`Vector{CompletenessJob}` — one job per (mass, separation, trial) combination.

# Example: cluster dispatch
```julia
jobs = completeness_jobs(masses=10 .^ range(-1, 2, 15), separations=10 .^ range(-0.3, 1.7, 15), n_trials=10)
# Write job index from SLURM_ARRAY_TASK_ID
job = jobs[parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
result = run_completeness_trial(job, system, sampler, detection; inject=my_inject)
# Save result...
```
"""
function completeness_jobs(;
    masses,
    separations,
    n_trials::Int=5,
)
    masses_vec = collect(Float64, masses)
    seps_vec = collect(Float64, separations)
    jobs = CompletenessJob[]
    sizehint!(jobs, length(masses_vec) * length(seps_vec) * n_trials)
    for (im, m) in enumerate(masses_vec)
        for (is, s) in enumerate(seps_vec)
            for it in 1:n_trials
                seed = hash((im, is, it, :completeness))
                push!(jobs, CompletenessJob(im, is, it, m, s, seed))
            end
        end
    end
    return jobs
end
export completeness_jobs


# ──────────────────────────────────────────────────────────────────────
# Phase 2: Run a single trial
# ──────────────────────────────────────────────────────────────────────

"""
    run_completeness_trial(job, system, sampler, detection_criterion; inject, add_noise=true, verbosity=0)

Execute a single injection-recovery trial.

# Arguments
- `job::CompletenessJob` — job description (grid point + seed)
- `system::System` — template system with priors, observations, and planets
- `sampler` — callable `(model) -> chain`; e.g. `m -> octofit(m, iterations=5000)`
- `detection_criterion` — callable `(chain, θ_true) -> Bool`; returns whether
  the injected companion was recovered

# Keyword Arguments
- `inject` — callable `(mass, separation) -> NamedTuple`; maps grid values to
  parameter overrides applied to the drawn prior sample. Must return overrides
  for *free* (prior) parameters only, not derived parameters.
  Example: `(m, s) -> (; planets=(; b=(; mass=m, a=s)))`
- `add_noise::Bool=true` — whether to add measurement noise to simulated data
- `verbosity::Int=0` — logging verbosity (0=silent, 1=info, 2=debug)

# Returns
[`CompletenessResult`](@ref) — contains the job and whether detection succeeded.

# Details
1. Draws parameters from `system`'s priors using a seeded RNG
2. Overrides parameters using `inject(job.mass, job.separation)`
3. Simulates observations via `generate_from_params`
4. Builds a `LogDensityModel` from the simulated system
5. **Initializes the sampler at the true parameters** (see module note)
6. Calls `sampler(model)` to obtain a posterior chain
7. Calls `detection_criterion(chain, θ_true)` to determine recovery

# Example
```julia
inject = (mass, sep) -> (; planets=(; b=(; planet_present=1.0, mass_prime=mass, a=sep)))

detection = function(chain, θ_true)
    p = mean(chain["b_planet_present"])
    BF = p / (1 - p)
    return BF > 3
end

result = run_completeness_trial(job, system,
    model -> octofit(model, iterations=5000, verbosity=0),
    detection;
    inject=inject,
)
```
"""
function run_completeness_trial(
    job::CompletenessJob,
    system::System,
    sampler,
    detection_criterion;
    inject,
    add_noise::Bool=true,
    verbosity::Int=0,
)
    rng = Xoshiro(job.seed)

    # 1. Draw parameters from priors
    θ_flat = sample_priors(rng, system)
    arr2nt = make_arr2nt(system)
    θ_nt = arr2nt(θ_flat)

    # 2. Apply parameter overrides from inject function
    overrides = inject(job.mass, job.separation)
    θ_flat = _apply_overrides!(θ_flat, arr2nt, overrides)
    θ_nt = arr2nt(θ_flat)

    verbosity >= 2 && @info "Trial $(job.i_trial) at mass=$(job.mass), sep=$(job.separation)" θ_nt

    # 3. Simulate observations from true parameters
    sim_system = generate_from_params(system, θ_nt; add_noise)

    # 4. Build model from simulated system
    model = LogDensityModel(sim_system; verbosity=0)

    # 5. Initialize at true parameters (the "cheat" for efficiency)
    _initialize_at_truth!(model, θ_flat)

    # 6. Run sampler
    verbosity >= 1 && @info "Sampling trial $(job.i_trial) at mass=$(round(job.mass, sigdigits=3)) Mjup, sep=$(round(job.separation, sigdigits=3)) AU"
    chain = sampler(model)

    # 7. Check detection
    detected = detection_criterion(chain, θ_nt)::Bool

    verbosity >= 1 && @info "Trial $(job.i_trial): detected=$(detected)"

    return CompletenessResult(job, detected)
end
export run_completeness_trial


# ──────────────────────────────────────────────────────────────────────
# Phase 3: Assemble results
# ──────────────────────────────────────────────────────────────────────

"""
    assemble_completeness(results; masses, separations)

Combine a collection of [`CompletenessResult`](@ref)s into a [`CompletenessMap`](@ref).

# Arguments
- `results` — iterable of `CompletenessResult`
- `masses` — the mass grid used to generate the jobs
- `separations` — the separation grid used to generate the jobs

# Returns
[`CompletenessMap`](@ref) with completeness fractions on the mass × separation grid.
"""
function assemble_completeness(
    results;
    masses,
    separations,
)
    masses_vec = collect(Float64, masses)
    seps_vec = collect(Float64, separations)
    n_mass = length(masses_vec)
    n_sep = length(seps_vec)

    n_detected = zeros(Int, n_mass, n_sep)
    n_total = zeros(Int, n_mass, n_sep)

    for r in results
        j = r.job
        n_total[j.i_mass, j.i_sep] += 1
        n_detected[j.i_mass, j.i_sep] += r.detected
    end

    completeness = n_detected ./ max.(n_total, 1)

    return CompletenessMap(masses_vec, seps_vec, completeness, n_detected, n_total)
end
export assemble_completeness


# ──────────────────────────────────────────────────────────────────────
# Convenience: run everything locally
# ──────────────────────────────────────────────────────────────────────

"""
    completeness_map(system, sampler, detection_criterion; inject, masses, separations, n_trials=5, add_noise=true, verbosity=1)

Compute a completeness map by running injection-recovery trials locally.

This is a convenience wrapper that calls [`completeness_jobs`](@ref),
[`run_completeness_trial`](@ref), and [`assemble_completeness`](@ref) in sequence.

For cluster-scale work, use the three-phase API directly.

# Arguments
- `system::System` — template system
- `sampler` — callable `(model) -> chain`
- `detection_criterion` — callable `(chain, θ_true) -> Bool`

# Keyword Arguments
- `inject` — callable `(mass, separation) -> NamedTuple` of parameter overrides
- `masses` — grid of companion masses [Mjup]
- `separations` — grid of semi-major axes [AU]
- `n_trials::Int=5` — trials per grid cell
- `add_noise::Bool=true` — add measurement noise to simulated data
- `verbosity::Int=1` — logging level

# Returns
[`CompletenessMap`](@ref)

# Example
```julia
using Octofitter, Distributions

# Define system with spike-and-slab detection model
planet_b = Planet(
    name = "b",
    basis = Visual{KepOrbit},
    observations = [my_rv_data, my_gaia_data],
    variables = @variables begin
        planet_present ~ Bernoulli(0.5)
        mass_prime ~ LogUniform(0.01, 100)  # Mjup
        mass = planet_present * mass_prime
        a ~ LogUniform(0.1, 100)            # AU
        e ~ Uniform(0, 0.5)
        i ~ Sine()
        Ω ~ UniformCircular()
        ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, a, i, ω, Ω)
    end,
)

sys = System(name="target", companions=[planet_b], observations=[...],
    variables = @variables begin
        M ~ truncated(Normal(1.0, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.5), lower=0.1)
    end,
)

# Run completeness map
cmap = completeness_map(
    sys,
    model -> octofit(model, iterations=5000, verbosity=0),
    (chain, θ) -> begin
        p = mean(chain["b_planet_present"])
        return p / (1 - p) > 3  # Bayes factor > 3
    end;
    inject = (mass, sep) -> (; planets=(; b=(; planet_present=1.0, mass_prime=mass, a=sep))),
    masses = 10 .^ range(-1, 2, length=12),
    separations = 10 .^ range(-0.3, 1.7, length=12),
    n_trials = 5,
)

# Plot
using CairoMakie
completenessplot(cmap)
```
"""
function completeness_map(
    system::System,
    sampler,
    detection_criterion;
    inject,
    masses,
    separations,
    n_trials::Int=5,
    add_noise::Bool=true,
    verbosity::Int=1,
)
    jobs = completeness_jobs(; masses, separations, n_trials)
    n_jobs = length(jobs)

    verbosity >= 1 && @info "Running $n_jobs completeness trials" n_masses=length(masses) n_separations=length(separations) n_trials

    results = Vector{CompletenessResult}(undef, n_jobs)
    for (i, job) in enumerate(jobs)
        verbosity >= 1 && print("\r  Trial $i / $n_jobs")
        results[i] = run_completeness_trial(
            job, system, sampler, detection_criterion;
            inject, add_noise, verbosity=max(0, verbosity - 1),
        )
    end
    verbosity >= 1 && println()

    cmap = assemble_completeness(results; masses, separations)

    if verbosity >= 1
        detected_total = sum(cmap.n_detected)
        trials_total = sum(cmap.n_total)
        @info "Completeness map complete" detected=detected_total total=trials_total overall_rate=round(detected_total/max(trials_total,1), digits=3)
    end

    return cmap
end
export completeness_map


# ──────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────

"""
Apply parameter overrides from a nested NamedTuple to a flat parameter vector.
Uses sentinel-value matching to find parameter indices (same approach as
`extract_fixed_params` in initialization.jl).
"""
function _apply_overrides!(θ_flat, arr2nt, overrides::NamedTuple)
    θ_nt = arr2nt(θ_flat)
    _apply_overrides_recursive!(θ_flat, θ_nt, overrides, "")
    return θ_flat
end
_apply_overrides!(θ_flat, arr2nt, ::Nothing) = θ_flat
_apply_overrides!(θ_flat, arr2nt, overrides::Tuple{}) = θ_flat

function _apply_overrides_recursive!(θ_flat, full_nt, partial_nt::NamedTuple, path)
    for name in propertynames(partial_nt)
        val = getproperty(partial_nt, name)
        if !hasproperty(full_nt, name)
            error("Cannot override parameter '$(path).$(name)': not found in model. " *
                  "Only free (prior) parameters can be overridden, not derived parameters.")
        end
        full_val = getproperty(full_nt, name)
        if val isa NamedTuple
            _apply_overrides_recursive!(θ_flat, full_val, val, "$(path).$(name)")
        else
            # Find index in flat array via sentinel matching
            sentinel = full_val
            idx = findfirst(==(sentinel), θ_flat)
            if isnothing(idx)
                error("Cannot override parameter '$(path).$(name)': could not locate in " *
                      "flat parameter vector. This parameter may be derived (not a free prior variable).")
            end
            θ_flat[idx] = Float64(val)
        end
    end
end

"""
Initialize a model's starting points at the true (injected) parameter values.
This is the "cheat" that avoids expensive blind initialization for completeness trials.
"""
function _initialize_at_truth!(model::LogDensityModel, θ_flat)
    θ_transformed = model.link(θ_flat)
    # Verify the true parameters have finite posterior density
    lp = model.ℓπcallback(θ_transformed)
    if !isfinite(lp)
        @warn "True parameters have non-finite log-posterior ($lp); initialization may fail"
    end
    model.starting_points = fill(θ_transformed, 8)
    return nothing
end
