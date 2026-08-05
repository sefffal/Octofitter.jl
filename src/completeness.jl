# ---------------------------------------------------
# Detection-completeness / sensitivity maps
#
# Injection-recovery on a grid, in three phases so the expensive middle one
# can be a cluster array job:
#
#   Phase 1: `completeness_jobs()`       — lightweight job descriptions
#   Phase 2: `run_completeness_trial()`  — inject, simulate, sample
#   Phase 3: `assemble_completeness()`   — apply a detection criterion
#
# Detection is *purely* phase 3. Each trial keeps the whole posterior chain and
# the true injected parameters, so a threshold can be changed without
# re-running the sampler. `completeness_map()` runs all three locally.
#
# NOTE: each trial initializes the sampler at the true injected parameters
# rather than running blind initialization. That cuts the cost enormously but
# makes the estimate optimistic about convergence: it measures the
# *statistical* detectability of a signal, not the ability to discover it
# blind.
# ---------------------------------------------------

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
- `mass::Float64` — companion mass [M⊙]
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

The output of a single injection-recovery trial. Stores the full posterior chain
and the true injected parameters so that detection criteria can be applied (and
re-applied) after the fact.

# Fields
- `job::CompletenessJob` — the job description that produced this result
- `chain::Chains` — full posterior chain from the sampler
- `θ_true::NamedTuple` — the true injected parameter values
"""
struct CompletenessResult{TChain,TParams}
    job::CompletenessJob
    chain::TChain
    θ_true::TParams
end
export CompletenessResult

"""
    CompletenessMap

Assembled completeness results on a 2D grid of mass × separation.

# Fields
- `masses::Vector{Float64}` — mass grid values [M⊙] (v9 has one mass unit;
  write `5mjup` for a Jupiter-mass grid point)
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
- `masses` — iterable of companion masses [M⊙] (v9 has one mass unit; write `5mjup` for a Jupiter-mass grid point)
- `separations` — iterable of semi-major axes [AU]
- `n_trials::Int=5` — number of independent trials per grid cell

# Returns
`Vector{CompletenessJob}` — one job per (mass, separation, trial) combination.

# Example: cluster dispatch
```julia
jobs = completeness_jobs(masses=10 .^ range(-1, 2, 15), separations=10 .^ range(-0.3, 1.7, 15), n_trials=10)
# Write job index from SLURM_ARRAY_TASK_ID
job = jobs[parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
result = run_completeness_trial(job, system, sampler; inject=my_inject)
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
# Phase 2: Run a single trial (inject → simulate → sample)
# ──────────────────────────────────────────────────────────────────────

"""
    run_completeness_trial(job, system, sampler; inject, add_noise=true, verbosity=0)

Execute a single injection-recovery trial: inject a companion, simulate
observations, and sample the posterior. No detection decision is made here —
that happens in [`assemble_completeness`](@ref).

The returned [`CompletenessResult`](@ref) stores the full posterior chain and
the true injected parameters, allowing detection criteria to be applied and
iterated on after the fact.

# Arguments
- `job::CompletenessJob` — job description (grid point + seed)
- `system::System` — template system with priors, observations, and bodies
- `sampler` — callable `(model) -> chain`; e.g. `m -> octofit(m, iterations=5000)`

# Keyword Arguments
- `inject` — callable `(mass, separation) -> NamedTuple`; maps grid values to
  parameter overrides applied to the drawn prior sample. Must return overrides
  for *free* (prior) parameters only, not derived parameters.
  Example: `(m, s) -> (; bodies=(; b=(; mass=m, a=s)))`
- `add_noise::Bool=true` — whether to add measurement noise to simulated data
- `verbosity::Int=0` — logging verbosity (0=silent, 1=info, 2=debug)

# Returns
[`CompletenessResult`](@ref) containing the job, posterior chain, and true parameters.

# Details
1. Draws parameters from `system`'s priors using a seeded RNG
2. Overrides parameters using `inject(job.mass, job.separation)`
3. Simulates observations via [`generate_from_params`](@ref)
4. Builds a `LogDensityModel` from the simulated system
5. **Initializes the sampler at the true parameters** (see the note at the top
   of this file)
6. Calls `sampler(model)` to obtain a posterior chain
7. Returns the chain and true parameters for later analysis

# Example
```julia
result = run_completeness_trial(job, system,
    model -> octofit(model, iterations=5000, verbosity=0);
    inject = (mass, sep) -> (; bodies=(; b=(; mass=mass, a=sep))),
)

# Inspect the posterior chain
result.chain
result.θ_true
```
"""
function run_completeness_trial(
    job::CompletenessJob,
    system::System,
    sampler;
    inject,
    add_noise::Bool=true,
    verbosity::Int=0,
)
    rng = Random.Xoshiro(job.seed)

    # 1. Draw parameters from priors
    θ_flat = sample_priors(rng, system)

    # 2. Apply parameter overrides from the inject function
    θ_flat = _apply_overrides!(θ_flat, system, inject(job.mass, job.separation))
    θ_nt = make_arr2nt(system)(θ_flat)

    verbosity >= 2 && @info "Trial $(job.i_trial) at mass=$(job.mass), sep=$(job.separation)" θ_nt

    # 3. Simulate observations from the true parameters. Seeded separately from
    #    the parameter draw (see `_with_seed`) so that re-running one array task
    #    reproduces its data as well as its truth.
    sim_system = _with_seed(hash((job.seed, :noise))) do
        generate_from_params(system, θ_nt; add_noise)
    end

    # 4. Build a model from the simulated system
    model = LogDensityModel(sim_system; verbosity=0)

    # 5. Initialize at the true parameters (the "cheat" for efficiency)
    _initialize_at_truth!(model, θ_flat)

    # 6. Run the sampler
    verbosity >= 1 && @info "Sampling trial $(job.i_trial) at mass=$(round(job.mass, sigdigits=3)) M⊙, sep=$(round(job.separation, sigdigits=3)) AU"
    chain = sampler(model)

    verbosity >= 1 && @info "Trial $(job.i_trial) complete"

    return CompletenessResult(job, chain, θ_nt)
end
export run_completeness_trial


# ──────────────────────────────────────────────────────────────────────
# Phase 3: Apply detection criterion and assemble
# ──────────────────────────────────────────────────────────────────────

"""
    assemble_completeness(results, detection_criterion; masses, separations)

Apply a detection criterion to a collection of [`CompletenessResult`](@ref)s
and assemble a [`CompletenessMap`](@ref).

Detection is applied here — not during the trial — so you can call this
function multiple times with different criteria to iterate on thresholds
without re-running the sampler.

# Arguments
- `results` — iterable of `CompletenessResult`
- `detection_criterion` — callable `(chain, θ_true) -> Bool`; returns whether
  the injected companion was recovered in a given trial
- `masses` — the mass grid used to generate the jobs
- `separations` — the separation grid used to generate the jobs

# Returns
[`CompletenessMap`](@ref) with completeness fractions on the mass × separation grid.

# Example
```julia
# Assemble with a Bayes factor threshold:
cmap_bf3 = assemble_completeness(results,
    (chain, θ) -> let p = mean(chain["b_planet_present"]); p/(1-p) > 3 end;
    masses=masses, separations=seps,
)

# Try a stricter threshold on the same results:
cmap_bf10 = assemble_completeness(results,
    (chain, θ) -> let p = mean(chain["b_planet_present"]); p/(1-p) > 10 end;
    masses=masses, separations=seps,
)

# Or a simple mass recovery criterion:
cmap_mass = assemble_completeness(results,
    (chain, θ) -> quantile(vec(chain["b_mass"]), 0.05) > 0.1;
    masses=masses, separations=seps,
)
```
"""
function assemble_completeness(
    results,
    detection_criterion;
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
        detected = detection_criterion(r.chain, r.θ_true)::Bool
        n_total[j.i_mass, j.i_sep] += 1
        n_detected[j.i_mass, j.i_sep] += detected
    end

    completeness = n_detected ./ max.(n_total, 1)

    return CompletenessMap(masses_vec, seps_vec, completeness, n_detected, n_total)
end
export assemble_completeness


# ──────────────────────────────────────────────────────────────────────
# Convenience: run everything locally
# ──────────────────────────────────────────────────────────────────────

"""
    completeness_map(system, sampler, detection_criterion; inject, masses, separations, n_trials=5, add_noise=true, verbosity=1) -> (CompletenessMap, Vector{CompletenessResult})

Compute a completeness map by running injection-recovery trials locally.

This is a convenience wrapper that calls [`completeness_jobs`](@ref),
[`run_completeness_trial`](@ref), and [`assemble_completeness`](@ref) in sequence.
Returns both the map and the full results vector, so you can re-apply different
detection criteria without re-sampling.

For cluster-scale work, use the three-phase API directly.

# Arguments
- `system::System` — template system
- `sampler` — callable `(model) -> chain`
- `detection_criterion` — callable `(chain, θ_true) -> Bool`

# Keyword Arguments
- `inject` — callable `(mass, separation) -> NamedTuple` of parameter overrides
- `masses` — grid of companion masses [M⊙]
- `separations` — grid of semi-major axes [AU]
- `n_trials::Int=5` — trials per grid cell
- `add_noise::Bool=true` — add measurement noise to simulated data
- `verbosity::Int=1` — logging level

# Returns
`(cmap::CompletenessMap, results::Vector{CompletenessResult})` — the assembled
map and the raw results for re-thresholding.

# Example
```julia
using Octofitter, Distributions

# Run completeness map
cmap, results = completeness_map(
    sys,
    model -> octofit(model, iterations=5000, verbosity=0),
    (chain, θ) -> quantile(vec(chain["b_mass"]), 0.05) > 0.1;
    inject = (mass, sep) -> (; bodies=(; b=(; mass=mass, a=sep))),
    masses = 10 .^ range(-1, 2, length=12),
    separations = 10 .^ range(-0.3, 1.7, length=12),
    n_trials = 5,
)

# Plot
using CairoMakie
completenessplot(cmap)

# Try a different threshold without re-running:
cmap_strict = assemble_completeness(results,
    (chain, θ) -> quantile(vec(chain["b_mass"]), 0.05) > 1.0;
    masses = 10 .^ range(-1, 2, length=12),
    separations = 10 .^ range(-0.3, 1.7, length=12),
)
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
            job, system, sampler;
            inject, add_noise, verbosity=max(0, verbosity - 1),
        )
    end
    verbosity >= 1 && println()

    cmap = assemble_completeness(results, detection_criterion; masses, separations)

    if verbosity >= 1
        detected_total = sum(cmap.n_detected)
        trials_total = sum(cmap.n_total)
        @info "Completeness map complete" detected=detected_total total=trials_total overall_rate=round(detected_total/max(trials_total,1), digits=3)
    end

    return cmap, results
end
export completeness_map


# ──────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────

"""
    _prior_slots(system) -> Vector{(group, owner, key, range)}

Where every free variable lives in the flat parameter vector, in the order
`make_arr2nt` consumes it: system priors, then each node's, then each
observation's. `group` is `:system`, `:bodies` or `:observations`; `range`
spans one index for a scalar prior and several for a vector-valued one.

v8 located an override's index by *sentinel matching* — searching the flat
vector for a value equal to the one the nested NamedTuple showed. That is
wrong whenever two parameters happen to draw the same number (constants,
discrete variables, a Dirac prior), and it cannot address one element of a
vector-valued prior at all. The layout is known exactly, so look it up.
"""
function _prior_slots(sys::System)
    out = Tuple{Symbol,Symbol,Symbol,UnitRange{Int}}[]
    i = 0
    for (k, d) in sys.priors.priors
        n = length(d)
        push!(out, (:system, :system, k, (i+1):(i+n)))
        i += n
    end
    for node in sys.nodes, (k, d) in node.priors.priors
        n = length(d)
        push!(out, (:bodies, node.name, k, (i+1):(i+n)))
        i += n
    end
    for o in sys.observations
        hasproperty(o, :priors) || continue
        nm = normalizename(likelihoodname(o))
        for (k, d) in o.priors.priors
            n = length(d)
            push!(out, (:observations, nm, k, (i+1):(i+n)))
            i += n
        end
    end
    return out
end

"""
    _apply_overrides!(θ_flat, system, overrides::NamedTuple)

Write the values in a nested override NamedTuple into the flat parameter
vector. The nesting matches the model's own parameter structure: system
variables at the top level, node variables under `bodies`, observation
variables under `observations`.

Only *free* (prior) variables can be overridden — a derived variable is a
function of them and setting it would be silently discarded, so it is an
error instead.
"""
function _apply_overrides!(θ_flat, sys::System, overrides::NamedTuple)
    if hasproperty(overrides, :planets)
        error("""
        Parameter overrides are nested under `bodies`, not `planets`, in v9:

            inject = (mass, sep) -> (; bodies=(; b=(; mass=mass, a=sep)))

        (`planets` became `bodies` when observations stopped living under a
        companion; the host star is a body like any other.)
        """)
    end
    slots = _prior_slots(sys)
    for name in propertynames(overrides)
        val = getproperty(overrides, name)
        if name === :bodies || name === :observations
            for owner in propertynames(val)
                sub = getproperty(val, owner)
                for key in propertynames(sub)
                    _set_slot!(θ_flat, slots, sys, name, owner, key, getproperty(sub, key))
                end
            end
        else
            _set_slot!(θ_flat, slots, sys, :system, :system, name, val)
        end
    end
    return θ_flat
end
_apply_overrides!(θ_flat, ::System, ::Nothing) = θ_flat
_apply_overrides!(θ_flat, ::System, ::Tuple{}) = θ_flat

function _set_slot!(θ_flat, slots, sys::System, group::Symbol, owner::Symbol, key::Symbol, val)
    idx = findfirst(s -> s[1] === group && s[2] === owner && s[3] === key, slots)
    if isnothing(idx)
        path = group === :system ? string(key) : "$group.$owner.$key"
        error("""
        Cannot override parameter `$path`: it is not a free (prior) variable of this model.

        Only variables declared with `~` can be overridden. Derived variables (declared
        with `=`) are computed from them, so setting one would have no effect. The model's
        free variables are:
        $(join(("  " * (s[1] === :system ? string(s[3]) : "$(s[1]).$(s[2]).$(s[3])") for s in slots), '\n'))
        """)
    end
    rng = slots[idx][4]
    if length(rng) == 1
        θ_flat[first(rng)] = _coerce(θ_flat[first(rng)], val)
    else
        length(val) == length(rng) || error(
            "Override for `$owner.$key` has $(length(val)) elements, but that variable has " *
            "$(length(rng)).")
        for (j, i) in enumerate(rng)
            θ_flat[i] = _coerce(θ_flat[i], val[j])
        end
    end
    return θ_flat
end

# Keep the element type a prior draw produced: a discrete variable stays
# discrete, so `model.link` still sees the type its Bijector expects.
_coerce(old::Bool, v) = !iszero(v)
_coerce(old, v) = convert(typeof(old), v)

"""
    _initialize_at_truth!(model, θ_flat; ndraws=8)

Point the sampler's starting points at the true (injected) parameter values.
This is the shortcut that lets a completeness grid skip blind initialization;
see the note at the top of this file for what it costs.
"""
function _initialize_at_truth!(model::LogDensityModel, θ_flat; ndraws::Int=8)
    # Preserve the element types of a prior draw, so a discrete variable is
    # not promoted to Float64 — the same invariant `startingpoints!` keeps.
    s = model.sample_priors(Random.default_rng())
    θ_transformed = convert.(typeof.(s), model.link(θ_flat))
    lp = model.ℓπcallback(θ_transformed)
    if !isfinite(lp)
        @warn "True parameters have non-finite log-posterior ($lp); initialization may fail"
    end
    model.starting_points = fill(θ_transformed, ndraws)
    return nothing
end

# ──────────────────────────────────────────────────────────────────────
# Makie-extension stubs (methods added by OctofitterMakieExt)
# ──────────────────────────────────────────────────────────────────────

"""
    completenessplot(cmap::CompletenessMap; kwargs...)

Heatmap of a completeness map over mass and separation. Requires a Makie
backend to be loaded (e.g. `using CairoMakie`).
"""
function completenessplot end
completenessplot(args...; kwargs...) = _require_makie("completenessplot")
export completenessplot

"""
    completenessplot!(gridposition, cmap::CompletenessMap; kwargs...)

Draw a completeness heatmap into an existing figure layout. Requires Makie.
See [`completenessplot`](@ref).
"""
function completenessplot! end
completenessplot!(args...; kwargs...) = _require_makie("completenessplot!")
export completenessplot!
