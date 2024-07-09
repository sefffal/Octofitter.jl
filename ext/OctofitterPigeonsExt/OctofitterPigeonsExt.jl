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
    verbosity=1
    initial_parameters = nothing
    initial_samples = 1000

    local result_pf = nothing
    ldm_any = Octofitter.LogDensityModelAny(model)
    # Use Pathfinder to initialize HMC. Works but currently disabled.
    verbosity >= 2 && @info "Determining initial positions using pathfinder"
    # It seems to hit a PosDefException sometimes when factoring a matrix.
    # When that happens, the next try usually succeeds.
    # start_time = time()
    for i in 1:8
        try
            if !isnothing(initial_parameters)
                initial_θ = initial_parameters
                # Transform from constrained support to unconstrained support
                initial_θ_t = model.link(initial_θ)
                errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
                result_pf = with_logger(errlogger) do 
                    Pathfinder.multipathfinder(
                        ldm_any, 1000;
                        nruns=8,
                        init=fill(collect(initial_θ_t),8),
                        progress=verbosity > 1,
                        maxiters=25_000,
                        rng=rng,
                        # maxtime=25.0,
                        # reltol=1e-4,
                    )
                end
            else
                init_sampler = function(rng, x) 
                    initial_θ, mapv = Octofitter.guess_starting_position(rng,model.system,initial_samples)
                    initial_θ_t = model.link(initial_θ)
                    x .= initial_θ_t
                end
                errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
                result_pf = with_logger(errlogger) do 
                    Pathfinder.multipathfinder(
                        ldm_any, 1000;
                        nruns=8,
                        init_sampler=Octofitter.CallableAny(init_sampler),
                        progress=verbosity > 1,
                        maxiters=25_000,
                        # maxtime=25.0,
                        reltol=1e-8,
                        rng=rng,
                    ) 
                end
            end
            # Check pareto shape diagnostic to see if we have a good approximation
            # If we don't, just try again
            if result_pf.psis_result.pareto_shape > 3
                verbosity > 3 && display(result_pf)
                verbosity >= 4 && display(result_pf.psis_result)
                verbosity > 2 && @warn "Restarting pathfinder" i
                continue
            end
            verbosity >= 3 && "Pathfinder complete"
            break
        catch ex
            result_pf = nothing
            if ex isa PosDefException
                verbosity > 2 && @warn "Mass matrix failed to factorize. Restarting pathfinder" i
                continue
            end
            if ex isa InterruptException
                rethrow(ex)
            end
            @error "Unexpected error occured running pathfinder" exception=(ex, catch_backtrace())
            break
        end
    end
    if isnothing(result_pf)
        # Guess initial starting positions by drawing from priors a bunch of times
        # and picking the best one (highest likelihood).
        # Then transform that guess into our unconstrained support
        if isnothing(initial_parameters)
            verbosity >= 1 && @info "Guessing a starting location by sampling from prior" initial_samples
            initial_θ, mapv = Octofitter.guess_starting_position(rng,model.system,initial_samples)
            verbosity > 2 && @info "Found starting location" θ=Octofitter.stringify_nested_named_tuple(model.arr2nt(initial_θ))
            # Transform from constrained support to unconstrained support
            initial_θ_t = model.link(initial_θ)
        else
            initial_θ = initial_parameters
            # Transform from constrained support to unconstrained support
            initial_θ_t = model.link(initial_θ)
        end
    else # !isnothing(result_pf)
        # Start using a draw from the typical set as estimated by Pathfinder
        # TODO: ideally we want to run multi-pathfinder once, store the result in the object,
        # and then pull out single draws on each chain.
        if result_pf.psis_result.pareto_shape < 3
            verbosity >= 4 && @info "PSIS result good; starting with sample from typical set"
            initial_θ_t = collect(last(eachcol(result_pf.draws))) # result_pf.fit_distribution.μ
        else
            verbosity >= 4 && @info "PSIS result bad; starting with distribution mean"
            initial_θ_t = collect(mean(result_pf.fit_distribution))
        end
        initial_θ = collect(model.invlink(initial_θ_t))
    end


    if verbosity >= 4
        initial_logpost = model.ℓπcallback(initial_θ_t)
        @info "Determined initial position" initial_θ initial_θ_nt=model.arr2nt(initial_θ) initial_logpost
    end
    
    return initial_θ_t
end

# Valid for reference model only
function Pigeons.sample_iid!(model_reference::Octofitter.LogDensityModel, replica, shared)
    # This could in theory be done without any array allocations
    θ = sample_priors(replica.rng, model_reference.system)
    # t = @elapsed θ, logpost = Octofitter.Octofitter.guess_starting_position(replica.rng, model_reference.system, 100_000)
    # @info "Sampled from reference" logpost t
    θ_t = model_reference.link(θ)
    replica.state .= θ_t
end

function Pigeons.default_reference(target::Octofitter.LogDensityModel)
    reference_sys = prior_only_model(target.system)
    # Note we could run into issues if their priors aren't well handled by the default
    # autodiff backend
    reference = Octofitter.LogDensityModel(reference_sys)
    return reference
end


function Pigeons.default_explorer(target::Octofitter.LogDensityModel)
    return Pigeons.Compose(
        Pigeons.SliceSampler(),
        Pigeons.AutoMALA()
    )
end


"""
octofit_pigeons(model; nrounds, n_chains=[auto])

Use Pigeons.jl to sample from intractable posterior distributions.

```julia
model = Octofitter.LogDensityModel(System, autodiff=:ForwardDiff, verbosity=4)
chain, pt = octofit_pigeons(model)
```
"""
Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    target::Octofitter.LogDensityModel;
    n_rounds::Int,
    n_chains::Int=16,
    n_chains_variational::Int=16,
    checkpoint::Bool=false,
    pigeons_kw...
)
    @nospecialize
    inputs = Pigeons.Inputs(;
        target,
        record = [traces; round_trip; record_default(); index_process],
        multithreaded=true,
        show_report=true,
        n_rounds,
        n_chains,
        n_chains_variational,
        variational = GaussianReference(first_tuning_round = 5),
        checkpoint,
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
            start_time,
            stop_time,
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end
Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    inputs::Pigeons.Inputs
)
    @nospecialize

    start_time = time()
    pt = pigeons(inputs)
    stop_time = time()

    mcmcchains = Chains(inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end
Base.@nospecializeinfer function MCMCChains.Chains(
    model::Octofitter.LogDensityModel,
    pt::Pigeons.PT,
    chain_num::Int=pt.inputs.n_chains
)
    ln_like = Octofitter.make_ln_like(model.system, model.arr2nt(Octofitter.sample_priors(model.system)))

    # Resolve the array back into the nested named tuple structure used internally.
    # Augment with some internal fields
    samples = get_sample(pt, chain_num)
    chain_res = map(samples) do sample 
        θ_t = @view(sample[begin:begin+model.D-1])
        logpost = sample[model.D+1]
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        loglike = ln_like(model.system, resolved_namedtuple)
        return merge((;
            loglike,
            logpost,
        ), resolved_namedtuple)
    end
    # Then finally flatten and convert into an MCMCChain object / table.
    # Mark the posterior, likelihood, numerical error flag, and tree depth as internal
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :loglike,
            :logpost
        ])
    )
    return mcmcchains
end

end