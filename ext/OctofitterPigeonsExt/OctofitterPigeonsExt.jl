module OctofitterPigeonsExt
using Random
using Octofitter
using Pigeons
using MCMCChains
using LinearAlgebra
using Logging
using Pathfinder

function (model::Octofitter.LogDensityModel)(θ)
    return model.ℓπcallback(θ)
end
function Pigeons.initialization(model::Octofitter.LogDensityModel, rng::AbstractRNG, chain_no::Int)

    Octofitter.get_starting_point!!(rng, model)

    if isnothing(model.starting_points) || !haskey(model.starting_points, chain_no)
        @show model.starting_points chain_no
        error("Insufficient starting points provided. Provide at least one per chain. (model.starting_points is too short)")
    end
    initial_θ_t = collect(model.starting_points[chain_no])
    # initial_θ_t = collect(Octofitter.get_starting_point!!(rng, model))
    initial_θ = model.invlink(initial_θ_t)
    initial_logpost = model.ℓπcallback(initial_θ_t)

    if any(!isfinite, initial_θ_t) || any(!isfinite, initial_θ) || !isfinite(initial_logpost)
        error("Could not find a starting point with finite arguments initial_logpost=$initial_logpost, initial_θ_t=$initial_θ_t, initial_θ=$(model.arr2nt(initial_θ))")
    end

    # @info "Determined initial position" chain_no initial_θ initial_θ_nt=model.arr2nt(initial_θ) initial_logpost
    # @info "Determined initial position" chain_no initial_logpost
    
    return initial_θ_t
end

# Valid for reference model only
function Pigeons.sample_iid!(model_reference::Octofitter.LogDensityModel, replica, shared)
    if _has_non_sampleable_priors(model_reference)
        Pigeons.step!(shared.explorer, replica, shared)
    else
        θ = model_reference.sample_priors(replica.rng)
        θ_t = model_reference.link(θ)
        replica.state .= θ_t
    end
end

function _has_non_sampleable_priors(model)
    ret = false
    ret |= any(Octofitter._isprior, model.system.observations)
    for planet in model.system.planets
        ret |= any(Octofitter._isprior, planet.observations)
    end
    return ret
end

function Pigeons.default_reference(target::Octofitter.LogDensityModel)
    reference_sys = prior_only_model(target.system)
    # TODO: adapt Pigeons to work using DifferentiationInterface
    reference = Octofitter.LogDensityModel(reference_sys; verbosity=0)
    reference.starting_points = target.starting_points
    return reference
end


function Pigeons.default_explorer(::Octofitter.LogDensityModel)
    return Pigeons.SliceSampler()
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

    # Turn off likelihood parallelism if we are sampling one chain per thread
    if multithreaded && n_chains + n_chains_variational > 1
        Octofitter._kepsolve_use_threads[] = false
    # Turn on likelihood parallelism if we have ~15x more data than threads.
    # This is based on some local benchmarks. Spawning tasks takes about 450ns;
    # an orbit solve takes about 32ns, or 1/14 as long.
    else
        Octofitter._kepsolve_use_threads[] = Octofitter._count_epochs(target.system) > 15*Threads.nthreads()
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
            mcmcchains.info,
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

    # Turn off likelihood parallelism if we are sampling one chain per thread
    if inputs.multithreaded
        Octofitter._kepsolve_use_threads[] = false
    # Turn on likelihood parallelism if we have ~15x more data than threads.
    # This is based on some local benchmarks. Spawning tasks takes about 450ns;
    # an orbit solve takes about 32ns, or 1/14 as long.
    else
        Octofitter._kepsolve_use_threads[] = _count_epochs(model.system) > 15*Threads.nthreads()
    end

    start_time = time()
    pt = pigeons(inputs)
    stop_time = time()

    mcmcchains = Chains(inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            mcmcchains.info,
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
    ln_like = Octofitter.make_ln_like(model.system, model.arr2nt(model.sample_priors(Random.default_rng())))

    # Resolve the array back into the nested named tuple structure used internally.
    # Augment with some internal fields
    if !isnothing(chain_num)
        samples = get_sample(pt, chain_num)
    else
        samples = get_sample(pt)
    end
    chain_res = map(samples) do sample 
        θ_t = @view(sample[begin:begin+model.D-1])
        logpot = sample[model.D+1]
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        ll = ln_like(model.system, resolved_namedtuple)
        lp = ln_prior(θ,true)
        # logpot does not equal ll + lp, so I'm not fully sure what it is.
        return merge((;
            loglike=ll,
            logprior=lp,
            logpost=ll+lp,
            pigeons_logpotential = logpot
        ), resolved_namedtuple)
    end
    # Then finally flatten and convert into an MCMCChain object / table.
    # Mark the posterior, likelihood, numerical error flag, and tree depth as internal
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :loglike,
            :logpost,
            :logprior,
            :pigeons_logpotential,
        ])
    )
    global_barrier_variational = pt isa Pigeons.StabilizedPT ?  Pigeons.global_barrier(pt) : nothing
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            model_name=try;pt.inputs.target.system.name; catch; nothing; end,
            logevidence_ratio=try;Pigeons.stepping_stone(pt); catch; nothing; end,
            global_barrier_variational,
            global_barrier=try;Pigeons.global_barrier(pt); catch; nothing; end,
            start_time=0,
            stop_time=try;sum(pt.shared.reports.summary.last_round_max_time); catch; nothing; end,
            sampler="pigeons"
        )
    )
    return mcmcchains_with_info
end

end