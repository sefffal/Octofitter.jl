# Parallel tempering via Pigeons.jl.
#
# Ported from the v1 extension (parked at `src/legacy/ext/OctofitterPigeonsExt/`).
# Only two things were v1-model-shaped and had to change:
#
#   * `_has_non_sampleable_priors` walked `system.planets` and each planet's
#     own observation list. There is one flat observation list now, plus the
#     `@variables`-emitted prior terms that hang off the system and its nodes
#     (`sys.priorterms`), so both have to be checked.
#   * the `Inputs` method referenced an undefined `model` when deciding whether
#     to thread the likelihood — a latent `UndefVarError` on the non-threaded
#     path. It reads `inputs.target` now.
module OctofitterPigeonsExt

using Random
using Octofitter
using Pigeons
using MCMCChains
using LinearAlgebra
using Logging

function (model::Octofitter.LogDensityModel)(θ)
    return model.ℓπcallback(θ)
end

function Pigeons.initialization(model::Octofitter.LogDensityModel, rng::AbstractRNG, chain_no::Int)

    Octofitter.get_starting_point!!(rng, model)

    # `model.starting_points` is a `Vector`, so validate the chain index with a
    # bounds check. (v1 used `haskey`, which is not defined for `Vector`/`Int`
    # and only worked incidentally when some other loaded package provided it;
    # under MPI the child process loads a minimal set of packages.)
    if isnothing(model.starting_points) || !checkbounds(Bool, model.starting_points, chain_no)
        error("Insufficient starting points provided. Provide at least one per chain (model.starting_points is too short: got $(isnothing(model.starting_points) ? "nothing" : length(model.starting_points)) for chain $chain_no).")
    end
    initial_θ_t = collect(model.starting_points[chain_no])
    initial_θ = model.invlink(initial_θ_t)
    initial_logpost = model.ℓπcallback(initial_θ_t)

    if any(!isfinite, initial_θ_t) || any(!isfinite, initial_θ) || !isfinite(initial_logpost)
        error("Could not find a starting point with finite arguments initial_logpost=$initial_logpost, initial_θ_t=$initial_θ_t, initial_θ=$(model.arr2nt(initial_θ))")
    end

    return initial_θ_t
end

# Valid for the reference model only.
function Pigeons.sample_iid!(model_reference::Octofitter.LogDensityModel, replica, shared)
    if _has_non_sampleable_priors(model_reference)
        Pigeons.step!(shared.explorer, replica, shared)
    else
        θ = model_reference.sample_priors(replica.rng)
        θ_t = model_reference.link(θ)
        replica.state .= θ_t
    end
end

"""
    _has_non_sampleable_priors(model) -> Bool

Whether the reference model still carries terms that reshape the prior, so
that drawing from the declared distributions is *not* a draw from the
reference. Two sources in v9: prior-shaped observations in the system's flat
observation list, and the terms `@variables` emitted for `lhs ~ dist` /
`LL +=` lines and for the `UnitLengthPrior` behind each `UniformCircular`
(`sys.priorterms`, which is a tuple of `owner => term` pairs).

`prior_only_model(sys; exclude_all=true)` is the model for which this is
`false`; the default `exclude_all=false` keeps those terms deliberately.

`BlankLikelihood` is excluded even though `_isprior` is true for it: it
contributes exactly zero, so it does not reshape the prior. (v8 got this by
accident — it never defined `_isprior` for `BlankLikelihood`, so it fell
through to the `false` fallback. Counting it here would send every
`prior_only_model` reference down the explorer path instead of drawing IID.)
"""
_reshapes_prior(@nospecialize term) =
    Octofitter._isprior(term) && !(term isa Octofitter.BlankLikelihood)

function _has_non_sampleable_priors(model)
    sys = model.system
    any(_reshapes_prior, sys.observations) && return true
    for (_, term) in sys.priorterms
        _reshapes_prior(term) && return true
    end
    return false
end

function Pigeons.default_reference(target::Octofitter.LogDensityModel)
    reference_sys = Octofitter.prior_only_model(target.system)
    reference = Octofitter.LogDensityModel(reference_sys; verbosity=0)
    reference.starting_points = target.starting_points
    return reference
end

function Pigeons.default_explorer(::Octofitter.LogDensityModel)
    return Pigeons.SliceSampler()
end

# Decide whether to thread the *likelihood*. Spawning a task costs ~450 ns and
# an orbit solve ~32 ns, so inner threading only pays when there is ~15x more
# data than there are threads — and never when Pigeons is already running one
# chain per thread.
function _set_kepsolve_threading!(sys, sampler_is_threaded::Bool)
    if sampler_is_threaded
        Octofitter._kepsolve_use_threads[] = false
    else
        Octofitter._kepsolve_use_threads[] =
            Octofitter._count_epochs(sys) > 15 * Threads.nthreads()
    end
    return Octofitter._kepsolve_use_threads[]
end

# A single posterior evaluation, timed at a finite point drawn from the prior.
# Used only to decide whether to print the `cores` advisory, so it returns
# `nothing` rather than throwing when no finite draw turns up (e.g. models
# whose prior mass is almost entirely outside the likelihood's domain).
function _estimate_eval_seconds(model::Octofitter.LogDensityModel)
    try
        rng = Random.Xoshiro(0)
        for _ in 1:20
            θ_t = model.link(collect(model.sample_priors(rng)))
            isfinite(model.ℓπcallback(θ_t)) || continue
            # Twice, keeping the faster: the model constructor already ran the
            # callback so this is warm, but a GC pause can still pollute one run.
            t1 = @elapsed model.ℓπcallback(θ_t)
            t2 = @elapsed model.ℓπcallback(θ_t)
            return min(t1, t2)
        end
    catch
    end
    return nothing
end

# Slice sampling costs roughly this many posterior evaluations per dimension
# per scan (doubling plus shrinkage, averaged over the tempering ladder). Only
# used for order-of-magnitude runtime advisories, never for control flow.
const _SLICE_EVALS_PER_DIM = 10

_fmt_duration(s) = s < 120 ? string(round(Int, s), " seconds" ) :
    s < 7200 ? string(round(s/60; digits=1), " minutes") :
    string(round(s/3600; digits=1), " hours")

# The modules a worker process must load before it can deserialize the model:
# Octofitter itself plus whichever package defines each observation's type
# (OctofitterRadialVelocity, OctofitterImages, ...). Users combining data
# kinds would otherwise have to know this list — and the failure mode for
# getting it wrong is a deserialization error inside an MPI child.
function _worker_dependencies(system, extra)
    deps = Module[Octofitter]
    for obs in system.observations
        m = Base.moduleroot(parentmodule(typeof(obs)))
        m in deps || push!(deps, m)
    end
    return vcat(deps, extra)
end

Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    target::Octofitter.LogDensityModel;
    n_rounds::Int,
    n_chains::Int=16,
    n_chains_variational::Int=16,
    checkpoint::Bool=false,
    multithreaded=true,
    variational = GaussianReference(first_tuning_round = 5),
    cores::Union{Nothing,Int}=nothing,
    threads_per_process::Int=1,
    dependencies::AbstractVector=[],
    pigeons_kw...
)
    @nospecialize

    threads_per_process >= 1 || throw(ArgumentError("threads_per_process must be a positive integer"))
    if !isnothing(cores)
        cores >= 1 || throw(ArgumentError("cores must be a positive integer"))
        if cores == 1
            @info "cores=1 samples in the current process; pass cores=2 or more to use separate worker processes."
            cores = nothing
        end
    end

    # Variational rereference often errors on models with discrete variables.
    contains_discrete_variables = any(isa.(Octofitter.sample_priors(Random.default_rng(), target.system), Integer))
    if contains_discrete_variables && (!isnothing(variational) || n_chains_variational != 0)
        @info "Variational reference is not supported with discrete variables; disabling and setting n_chains_variational=0."
        variational = nothing
        n_chains_variational = 0
    end

    if _has_non_sampleable_priors(target)
        @warn "This model has priors that cannot be sampled IID."
    end

    # Order-of-magnitude runtime advisory. Deliberately talks about "worker
    # processes", never MPI: the audience is users who have no reason to know
    # what MPI is, and the mechanism is an implementation detail.
    t_eval = _estimate_eval_seconds(target)
    if !isnothing(t_eval)
        chains_total = n_chains + n_chains_variational
        core_seconds = 2.0^(n_rounds + 1) * _SLICE_EVALS_PER_DIM * target.D * t_eval * chains_total
        if isnothing(cores)
            est = core_seconds / min(Threads.nthreads(), chains_total)
            if est > 600
                suggested = min(chains_total, Sys.CPU_THREADS, 16)
                @info """This model is expensive to evaluate: sampling may take on the order of $(_fmt_duration(est)).
                On a multi-core machine, `octofit_pigeons(model; n_rounds=$n_rounds, cores=$suggested)` runs the sampler
                in separate worker processes, which is often about twice as fast for models like this one
                (each launch spends a minute or two starting workers before sampling begins)."""
            end
        else
            est = core_seconds / cores
            if est < 60
                @info """This model evaluates quickly (sampling itself may take well under a minute), so the
                one-or-two minute startup cost of `cores=$cores` worker processes will likely dominate.
                Dropping the `cores` keyword and sampling with threads may well be faster overall."""
            end
        end
    end

    if !isnothing(cores)
        # `cores` is the total CPU budget: processes × threads_per_process
        # never exceeds it, so asking for cores=8 on an 8-core laptop cannot
        # oversubscribe no matter how the threads are split.
        n_processes, rem = divrem(cores, threads_per_process)
        if n_processes < 1
            throw(ArgumentError("cores=$cores is smaller than threads_per_process=$threads_per_process; there is nothing to launch"))
        end
        if rem != 0
            @warn "cores=$cores is not a multiple of threads_per_process=$threads_per_process; launching $n_processes processes ($(n_processes*threads_per_process) cores used)."
        end
        # Worker processes cannot report back through memory; the checkpoint
        # written each round is how results return, so it is always on here.
        inputs = Pigeons.Inputs(;
            target,
            record = [traces; round_trip; record_default(); index_process],
            multithreaded = false,
            show_report = true,
            n_rounds,
            n_chains,
            n_chains_variational,
            checkpoint = true,
            variational,
            pigeons_kw...
        )
        return _octofit_pigeons_childprocess(inputs, n_processes, threads_per_process,
            _worker_dependencies(target.system, dependencies))
    end

    _set_kepsolve_threading!(target.system, multithreaded && n_chains + n_chains_variational > 1)
    @info "Sampler running with multiple threads     : $multithreaded"
    @info "Likelihood evaluated with multiple threads: $(Octofitter._kepsolve_use_threads[])"
    inputs = Pigeons.Inputs(;
        target,
        record = [traces; round_trip; record_default(); index_process],
        multithreaded=true,
        show_report=true,
        n_rounds,
        n_chains,
        n_chains_variational,
        checkpoint,
        variational,
        pigeons_kw...
    )
    return octofit_pigeons(inputs)
end

Base.@nospecializeinfer function _octofit_pigeons_childprocess(
    inputs::Pigeons.Inputs,
    n_processes::Int,
    threads_per_process::Int,
    worker_dependencies::Vector,
)
    @nospecialize

    # A worker is a fresh process, so the parent's `_kepsolve_use_threads`
    # setting does not reach it — globals are not serialized with the model.
    # Pigeons `include`s path-valued dependencies in the worker after loading
    # the module-valued ones, which is exactly the hook needed to set the flag
    # once Octofitter is loaded there.
    if threads_per_process > 1
        worker_init = joinpath(mktempdir(), "octofitter-worker-init.jl")
        write(worker_init, "Octofitter._kepsolve_use_threads[] = Threads.nthreads() > 1\n")
        worker_dependencies = vcat(worker_dependencies, worker_init)
    end

    cpus = n_processes * threads_per_process
    desc = threads_per_process == 1 ? "$n_processes worker processes" :
        "$n_processes worker processes × $threads_per_process threads ($cpus cores)"
    @info """Starting $desc for sampling.
    The first minute or two goes to loading packages and compiling the model in each
    worker — progress will appear below once sampling begins. Results are checkpointed
    every round under `$(joinpath(pwd(), "results"))`."""

    # MPICH's default binding pins each rank to one core, which is right for
    # single-threaded ranks and disastrous for threaded ones: a rank's
    # threads timeshare its one core (measured ~4x slower end to end).
    mpiexec_args = threads_per_process > 1 ? `-bind-to none` : ``

    start_time = time()
    result = try
        pigeons(inputs, Pigeons.ChildProcess(;
            n_local_mpi_processes = n_processes,
            n_threads = threads_per_process,
            dependencies = worker_dependencies,
            mpiexec_args,
        ))
    catch err
        @error """Sampling in worker processes failed. If the error below points to loading or
        deserialization in a worker, check that every package used to build the model is
        loaded, passing any beyond Octofitter and its companions via the `dependencies`
        keyword. Dropping the `cores` keyword samples in-process with threads instead."""
        rethrow(err)
    end
    pt = Pigeons.load(result)
    stop_time = time()

    mcmcchains = Chains(pt.inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            mcmcchains.info...,
            start_time,
            stop_time,
            model_name=pt.inputs.target.system.name,
            sampler="pigeons"
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end

Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    pt::Pigeons.PT
)
    @nospecialize

    start_time = time()
    pt = pigeons(pt)
    stop_time = time()

    mcmcchains = Chains(pt.inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            mcmcchains.info...,
            start_time,
            stop_time,
            model_name=pt.inputs.target.system.name,
            sampler="pigeons"
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end

Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    inputs::Pigeons.Inputs
)
    @nospecialize

    _set_kepsolve_threading!(inputs.target.system, inputs.multithreaded)

    start_time = time()
    pt = pigeons(inputs)
    stop_time = time()

    mcmcchains = Chains(inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            mcmcchains.info...,
            start_time,
            stop_time,
            sampler="pigeons"
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end

Base.@nospecializeinfer function MCMCChains.Chains(
    model::Octofitter.LogDensityModel,
    pt::Pigeons.PT,
    chain_num::Union{Nothing,Int}=nothing
)
    ln_prior = Octofitter.make_ln_prior_transformed(model.system)
    # `loglike` is the data terms only; the prior-shaped observations it skips
    # are reported with the parameter priors under `logprior`, matching the
    # other samplers.
    θ_example = model.arr2nt(model.sample_priors(Random.default_rng()))
    ln_like = Octofitter.make_ln_like(model.system, θ_example)
    ln_like_data = Octofitter.make_ln_like(model.system, θ_example; include_priors=false)

    # Resolve the array back into the nested named tuple structure used
    # internally, augmented with the sampler's own diagnostics.
    if !isnothing(chain_num)
        samples = get_sample(pt, chain_num)
    else
        samples = get_sample(pt)
    end
    chain_res = map(samples) do sample
        θ_t = @view(sample[begin:begin+model.D-1])
        logpot = sample[model.D+1]
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        ll = ln_like(model.system, resolved_namedtuple)
        ll_data = ln_like_data(model.system, resolved_namedtuple)
        lp = ln_prior(θ, true)
        # `logpot` is the tempered log potential and does not equal ll + lp.
        # `ll - ll_data` is the prior-shaped part of the likelihood sum.
        return merge((;
            loglike=ll_data,
            logprior=lp + (ll - ll_data),
            logpost=ll+lp,
            pigeons_logpotential = logpot
        ), resolved_namedtuple)
    end
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :loglike,
            :logpost,
            :logprior,
            :pigeons_logpotential,
        ])
    )
    global_barrier_variational = pt isa Pigeons.StabilizedPT ? Pigeons.global_barrier(pt) : nothing
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            mcmcchains.info...,
            model_name=try; pt.inputs.target.system.name; catch; nothing; end,
            logevidence_ratio=try; Pigeons.stepping_stone(pt); catch; nothing; end,
            global_barrier_variational,
            global_barrier=try; Pigeons.global_barrier(pt); catch; nothing; end,
            start_time=0,
            stop_time=try; sum(pt.shared.reports.summary.last_round_max_time); catch; nothing; end,
            sampler="pigeons"
        )
    )
    return mcmcchains_with_info
end

end
