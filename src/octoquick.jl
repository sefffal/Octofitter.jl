
octoquick(model::LogDensityModel; kwargs...) = 
    octoquick(Random.default_rng(), model; kwargs...)
Base.@nospecializeinfer function octoquick(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    model::LogDensityModel;
    verbosity::Int=2,
    initial_samples::Int=10_000,
    nruns=cld(16,Threads.nthreads())*Threads.nthreads(),
)
    @nospecialize

    # start_time = time()
    local result_pf = nothing
    ldm_any = LogDensityModelAny(model)
    
    init_sampler = function(rng, x) 
        initial_θ, mapv = guess_starting_position(rng,model.system,initial_samples)
        initial_θ_t = model.link(initial_θ)
        x .= initial_θ_t
    end
    local start_time, stop_time
    errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
    result_pf = with_logger(errlogger) do 
        start_time = time()
        result_pf = Pathfinder.multipathfinder(
            ldm_any, 1000*nruns;
            # Do at least 16 runs, rounding up to nearest multiple
            # of nthreads
            nruns,
            init_sampler=CallableAny(init_sampler),
            progress=verbosity > 1,
            maxiters=10_000,
            executor=Transducers.ThreadedEx(),
            # reltol=1e-4,
        ) 
        stop_time = time()
        return result_pf
    end

    augmented_resolved_samples = map(eachcol(result_pf.draws)) do θ_t
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        logpost = model.ℓπcallback(θ)
        return merge((;logpost = logpost,), resolved_namedtuple)
    end
    pathfinder_chain =  Octofitter.result2mcmcchain(
        augmented_resolved_samples,
        Dict(:internals => [
            :logpost
        ])
    )
    pathfinder_chain_with_info = MCMCChains.setinfo(
        pathfinder_chain,
        (;
            start_time,
            stop_time,
            result_pf,
        )
    )

    return pathfinder_chain_with_info
end
export octoquick