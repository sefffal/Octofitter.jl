function gridsearch(model, N=4; q=0.65)

    # Get a list of all priors for the system as a flat vector of
    # Distribution objects
    prior_distributions = _list_priors(model.system)

    # We will perform a grid search over these parmeters. Some of them
    # may be unbounded (eg M ~ Normal())
    q_start = (1-q)/2
    q_stop  = 1 - (1-q)/2
    starts = quantile.(prior_distributions, q_start)::Vector{Float64}
    stops = quantile.(prior_distributions, q_stop)::Vector{Float64}
    ranges = range.(starts, stops, length=N)

    θ = zeros(length(ranges))
    # @show prod([length.(ranges)...])
    @showtime LL = zeros(length.(ranges)...)::Array{Float64,length(prior_distributions)}
    @show length(LL)
    @time Threads.@threads for I in CartesianIndices(LL)
        for variable_i in 1:length(prior_distributions)
            θ[variable_i] = ranges[variable_i][I[variable_i]]
        end
        LL[I] = model.ℓπcallback(model.link(θ))
    end

    I_best = argmax(LL)
    @show LL[I_best]
    θ_best = getindex.(ranges, Tuple(I_best))

    # # Reshape each range so that it points in a different dimension along the grid
    # ranges_reshaped = map(enumerate(ranges)) do (i,range)
    #     reshape(range, fill(1,i-1)..., :)
    # end

    # # # Now broadcast over the multidimensional grid
    # @time LL = broadcast(ranges_reshaped...) do θ...
    #     return sum(θ)
    # end
    # Evaluate log post and log like
    logpost = LL[I_best]
    resolved_namedtuple = model.arr2nt(θ_best)
    ln_like = make_ln_like(model.system, resolved_namedtuple)
    loglike = ln_like(model.system, resolved_namedtuple)
    nt = (; logpost, loglike, model.arr2nt(θ_best)...)
    chain_best = result2mcmcchain(
        [nt], 
        Dict(:internals => [:logpost, :loglike])
    )

    return (;LL, chain_best)

end