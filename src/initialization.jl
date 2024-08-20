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
function default_initializer!(rng::Random.AbstractRNG, model::LogDensityModel; nruns=8, ntries=2, ndraws=1000, initial_samples=10000, verbosity=1)


    local result_pf = nothing
    local metric = nothing
    ldm_any = LogDensityModelAny(model)
    verbosity >= 1 && @info "Determining initial positions and metric using pathfinder"
    # It can sometimes hit a PosDefException sometimes when factoring a matrix.
    # When that happens, the next try usually succeeds.
    for i in 1:ntries
        try
            verbosity >= 3 && @info "Starting multipathfinder run"
            init_sampler = function(rng, x) 
                if verbosity > 3
                    @info "drawing new starting guess by sampling IID from priors"
                end
                initial_θ, mapv = guess_starting_position(rng,model,initial_samples)
                if verbosity > 3
                    @info "Starting point drawn" initial_logpost=mapv
                end
                initial_θ_t = model.link(initial_θ)
                x .= initial_θ_t
            end
            errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
            result_pf = with_logger(errlogger) do 
                Pathfinder.multipathfinder(
                    ldm_any, ndraws;
                    nruns,
                    init_sampler=CallableAny(init_sampler),
                    progress=verbosity > 1,
                    maxiters=25_000,
                    reltol=1e-6,
                    rng=rng,
                    optimizer=Pathfinder.Optim.LBFGS(;
                        m=6,
                        linesearch=Pathfinder.Optim.LineSearches.BackTracking(),
                        alphaguess=Pathfinder.Optim.LineSearches.InitialHagerZhang()
                    )
                ) 
            end
            # Check pareto shape diagnostic to see if we have a good approximation
            # If we don't, just try again
            if result_pf.psis_result.pareto_shape > 3
                verbosity > 3 && display(result_pf)
                verbosity >= 4 && display(result_pf.psis_result)
                i<ntries && verbosity > 2 && @warn "Restarting pathfinder" i
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
    
    if !isnothing(result_pf)
        model.starting_points = collect.(eachcol(result_pf.draws))
        logposts = model.ℓπcallback.(model.starting_points)
        initial_logpost_range = extrema(logposts)
    else
        verbosity > 1 && @warn("Falling back to sampling from the prior and keeping the $ndraws samples with highest posterior density.")
        samples = sample_priors(rng, model, ndraws)
        samples_t = model.link.(samples)
        logposts = model.ℓπcallback.(samples_t)
        II = sortperm(logposts)[end-ndraws+1:end]
        model.starting_points = samples_t[II]
        initial_logpost_range = extrema(@view logposts[II])
        logposts = logposts[II]
    end

    if verbosity >= 1
        @info "Found a sample of initial positions" initial_logpost_range
    end

    return model.arr2nt(model.invlink(model.starting_points[argmax(logposts)]))
end