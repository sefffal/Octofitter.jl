# This file contains functions for initializing a chain
# These functions run on a LogDensity model and store a set of 
# plausible initial positions. These shouldn't usually be the maximum
# a-posteriori estimate, but draws from an approximation of the 
# typical set.

"""
    θ_transformed = get_starting_point!!(rng=Random.default_rng(), model::LogDensityModel)

Retrieve or generate a starting point for a model. Draws samples from `model.starting_points`.
If `model.starting_points` has not yet been set, it first runs `Octofitter.default_initializer!(model)`.
"""
function get_starting_point!!(rng::Random.AbstractRNG, model::LogDensityModel)
    if isnothing(model.starting_points)
        default_initializer!(rng, model)
    end
    return rand(rng, model.starting_points)
end
function get_starting_point!!(model::LogDensityModel; kwargs...)
    return get_starting_point!!(Random.default_rng(), model; kwargs...)
end

"""
    default_initializer!(model::LogDensityModel; initial_samples=100_000)

Prepare a set of possible starting points for a given LogDensityModel.
Use multi-pathfinder (Pathfinder.jl) to optimize and fit a variational approximation.
If this fails repeatedly, simply draw `initial_samples` from the prior and keeping
`ndraws` samples with the highest posterior density.
"""
function default_initializer!(model::LogDensityModel; kwargs...)
    return default_initializer!(Random.default_rng(), model; kwargs...)
end
function default_initializer!(rng::Random.AbstractRNG, model::LogDensityModel; initial_point = nothing, nruns=8, ntries=2, ndraws=1000, initial_samples=10000, verbosity=1)


    local result_pf = nothing
    local metric = nothing
    ldm_any = LogDensityModelAny(model)
    verbosity >= 1 && @info "Determining initial positions and metric using pathfinder"
    # It can sometimes hit a PosDefException sometimes when factoring a matrix.
    # When that happens, the next try usually succeeds.
    try
        for i in 1:ntries
            verbosity >= 3 && @info "Starting multipathfinder run"
            init_sampler = function(rng, x) 
                if isnothing(initial_point) || length(initial_point) < model.D
                    if verbosity > 3
                        @info "drawing new starting guess by sampling IID from priors"
                    end
                    initial_θ, mapv = guess_starting_position(rng,model,initial_samples)
                    if verbosity > 3
                        @info "Starting point drawn" initial_logpost=mapv
                    end
                end
                if !isnothing(initial_point)
                    if length(initial_point) < model.D
                        initial_θ = (initial_point..., initial_θ[length(initial_point)+1:end]...)
                    else
                        initial_θ = initial_point
                    end
                end
                initial_θ_t = model.link(initial_θ)
                x .= initial_θ_t
            end
            errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
            initial_mt = _kepsolve_use_threads[]
            _kepsolve_use_threads[] = false
            result_pf = with_logger(errlogger) do 
                result_pf = Pathfinder.multipathfinder(
                    ldm_any, ndraws;
                    nruns,
                    init_sampler=CallableAny(init_sampler),
                    progress=verbosity > 1,
                    maxiters=25_000,
                    reltol=1e-6,
                    rng=rng,
                    ntries=1,
                    executor=Pathfinder.Transducers.PreferParallel(),
                    optimizer=Pathfinder.Optim.LBFGS(;
                        m=6,
                        linesearch=Pathfinder.Optim.LineSearches.BackTracking(),
                        alphaguess=Pathfinder.Optim.LineSearches.InitialHagerZhang()
                    )
                ) 
                return result_pf
            end
            _kepsolve_use_threads[] = initial_mt
            # Check pareto shape diagnostic to see if we have a good approximation
            # If we don't, just try again
            if result_pf.psis_result.pareto_shape > 3
                verbosity > 3 && display(result_pf)
                verbosity >= 4 && display(result_pf.psis_result)
                i<ntries && verbosity > 2 && @warn "Restarting pathfinder" i
                continue
            end
            
            verbosity >= 3 && "Pathfinder complete"
            verbosity > 2 &&  display(result_pf)
            break
        end
    catch
    end
    
    if !isnothing(result_pf)
        model.starting_points = collect.(eachcol(result_pf.draws))
        logposts = model.ℓπcallback.(model.starting_points)
        initial_logpost_range = extrema(logposts)
    end
    # Occasionally there is a failure mode of pathfinder where, despite starting it at a reasonable spot, it returns garbage
    # starting draws that are orders of magnitude worse.
    # Check for this by ensuring the highest a-posteriori pathfinder draw is better than a random guess
    _, random_guess_logpost = guess_starting_position(rng,model,100)
    if isnothing(result_pf) || maximum(initial_logpost_range) < random_guess_logpost
        verbosity >= 1 && @warn "Falling back to sampling from the prior and keeping the $ndraws samples with highest posterior density."
        samples_t = map(1:1000) do _
            initial_θ, mapv = guess_starting_position(rng,model,max(1,initial_samples÷100))
            initial_θ_t = model.link(initial_θ)
            return initial_θ_t
        end
        # samples = sample_priors(rng, model, ndraws)
        # samples_t = model.link.(samples)
        logposts = model.ℓπcallback.(samples_t)
        II = sortperm(logposts, rev=true)[begin:ndraws]
        model.starting_points = samples_t[II]
        initial_logpost_range = extrema(@view logposts[II])
        logposts = logposts[II]
    end

    if verbosity >= 1
        @info "Found a sample of initial positions" initial_logpost_range
    end

    return model.arr2nt(model.invlink(model.starting_points[argmax(logposts)]))
end