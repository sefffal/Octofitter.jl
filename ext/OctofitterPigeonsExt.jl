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

Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    target::Octofitter.LogDensityModel;
    n_rounds::Int,
    n_chains::Int=16,
    n_chains_variational::Int=16,
    checkpoint::Bool=false,
    multithreaded=true,
    variational = GaussianReference(first_tuning_round = 5),
    pigeons_kw...
)
    @nospecialize

    _set_kepsolve_threading!(target.system, multithreaded && n_chains + n_chains_variational > 1)

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
