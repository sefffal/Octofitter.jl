# This file contains functions for initializing a chain
# These functions run on a LogDensity model and store a set of 
# plausible initial positions. These shouldn't usually be the maximum
# a-posteriori estimate, but draws from an approximation of the 
# typical set.

function guess_starting_position(model::LogDensityModel, args...; kwargs...)
    return guess_starting_position(Random.default_rng(), model, args...; kwargs...)
end
# Sample IID from the prior many times and return the highest posterior density sample.
# If there is applicable planets and data, at each guess attempt a step of 
# Orbits for the Impatient (Blunt et al) to improve the rate at which samples 
# are generated close to the data.
function guess_starting_position(rng::Random.AbstractRNG, model::LogDensityModel, N=500_000, enable_ofti=true)
    if !_kepsolve_use_threads[]
        # Seed RNGs for each thread
        rngs = [
            Xoshiro(reinterpret(UInt64,rand(rng)))
            for i in 1: Threads.nthreads()
        ]
        bestlogposts = fill(-Inf64, 1:Threads.nthreads())
        bestparams = fill(model.sample_priors(rng), 1:Threads.nthreads())

        Threads.@threads for i in 1:N
            tid = Threads.threadid()
            params = model.sample_priors(rngs[tid])
            params_t = model.link(params)
            logpost = model.ℓπcallback(params_t)
            
            params_ofti = ofti_step(rng, model, params)
            params_ofti_t = model.link(params_ofti)
            logpost_ofti = model.ℓπcallback(params_ofti_t)
            if logpost_ofti > logpost
                logpost = logpost_ofti
                params = params_ofti
            end

            logpost = model.ℓπcallback(params_t)
            if logpost > bestlogposts[tid]
                bestparams[tid] = params
                bestlogposts[tid] = logpost
            end
        end
        I_max = argmax(bestlogposts)
        return bestparams[I_max], bestlogposts[I_max]
    else
        bestparams = model.sample_priors(rng)
        bestlogpost = -Inf64
        for _ in 1:N
            params = model.sample_priors(rng)
            params_t = model.link(params)
            logpost = model.ℓπcallback(params_t)
            # println("logpost_initial = ", logpost)
            # Try an OFTI step and only keep it if it helps
            params_ofti = ofti_step(rng, model, params)
            # @show model.arr2nt(params)
            params_ofti_t = model.link(params_ofti)
            logpost_ofti = model.ℓπcallback(params_ofti_t)
                # println("logpost_ofti = ", logpost_ofti)
            if logpost_ofti > logpost
                logpost = logpost_ofti
                params = params_ofti
            end
                # println("logpost_after = ", logpost)
            if logpost > bestlogpost
                bestlogpost = logpost
                bestparams = params
            end        
        end
        return bestparams, bestlogpost
    end
end


# Implementation of a similar algorithm to Orbits for the Impatient, by Blunt et al.
# For each planet (if present) and if that planet has astrometry data, move the planet
# to the right PA and then scale the orbit to match. 
# This is a tiny bit different than OFTI, where the longitude of ascending node is adjusted
# to rotate the orbit. We typically sample from PA-at-epoch variable, so we can just adjust
# that directly. 
function ofti_step(rng, model, params)

    # Check if there are even any planets in the model, and bail out if not
    if length(model.system.planets) == 0
        return params
    end

    # θ_system = model.arr2nt(params)

    for planet_key in keys(model.system.planets)
        planet = model.system.planets[planet_key]

        astrom_likes = filter(obs->obs isa PlanetRelAstromLikelihood, planet.observations)
        if isempty(astrom_likes)
            continue
        end
        astrom_like = rand(rng, astrom_likes)
        data_index = rand(rng, eachindex(astrom_like.table.epoch))
        
        # θ_planet = θ_system.planets[planet_key]
        
        # orbit = Octofitter.construct_elements(Octofitter.orbittype(planet), θ_system, θ_planet)

        ind_a = _get_planet_param_index(model, params, planet_key, :a)
        # TODO: add support for sampling from period as well.
        if isnothing(ind_a)
            continue
        end
        # TODO: add support for sampling from tau as well.
        if haskey(planet.priors.priors, :θ)
            uniform_circular_param = false
            ind_θ = _get_planet_param_index(model, params, planet_key, :θ)
            old_θ = params[ind_θ]
        elseif haskey(planet.priors.priors, :θx) && haskey(planet.priors.priors, :θx)
            uniform_circular_param = true
            ind_θy = _get_planet_param_index(model, params, planet_key, :θy)
            ind_θx = _get_planet_param_index(model, params, planet_key, :θx)
            old_θ = atan(params[ind_θy],params[ind_θx])
        else
            continue
        end

        
        t = astrom_like.table.epoch[data_index]
        # pa_model = posangle(orbit, t)

        # We randomize the guess within uncertainties a bit, but no need to get the correlation
        if hasproperty(astrom_like.table, :sep)
            sep_data = astrom_like.table.sep[data_index] #+ randn(rng)*astrom_like.table.σ_sep[data_index]
            pa_data = astrom_like.table.pa[data_index] #+ randn(rng)*astrom_like.table.σ_pa[data_index]
        else
            x = astrom_like.table.ra[data_index]
            y = astrom_like.table.dec[data_index]
            sep_data = sqrt(x^2 + y^2) #+ randn(rng)*astrom_like.table.σ_ra[data_index] 
            pa_data = rem2pi(atan(x,y), RoundDown) #+ randn(rng)*astrom_like.table.σ_dec[data_index],RoundDown)
        end

        new_θ = pa_data

        # Update parameters
        if uniform_circular_param
            new_θx, new_θy = sincos(new_θ)
            params = Base.setindex(params, new_θy, ind_θy)
            params = Base.setindex(params, new_θx, ind_θx)
        else
            new_θ = clamp(
                new_θ,
                nextfloat(float(minimum(planet.priors.priors[:θ]))),
                prevfloat(float(maximum(planet.priors.priors[:θ]))),
            )
            params = Base.setindex(params, new_θ, ind_θ)
        end

        ##########
        # Now calculate separation at this updated location

        θ_system = model.arr2nt(params)
        θ_planet = θ_system.planets[planet_key]
        
        orbit = Octofitter.construct_elements(Octofitter.orbittype(planet), θ_system, θ_planet)
        sep_model = projectedseparation(orbit, t)

        scale_sep = sep_data/sep_model
        new_a = params[ind_a] * scale_sep

        # Clamp to prior support
        new_a = clamp(
            new_a,
            nextfloat(float(minimum(planet.priors.priors[:a]))),
            prevfloat(float(maximum(planet.priors.priors[:a]))),
        )
        
        params = Base.setindex(params, new_a, ind_a)

        # This doesn't work, need to make it more general.
        # Find the sep at epoch.
        # Scale to right size.
        # Optimize the theta or tau parameter.
        # Done.
    end
    return params
end
function _get_planet_param_index(model, params, planet_key, param_name)
    # TODO: Hack! Find index in flat array corresponding to desired parameters using sentials.
    # This is really gross, we really want an nt2arr function.
    params_search_nt = model.arr2nt(params)
    if !hasproperty(params_search_nt.planets[planet_key], param_name)
        return nothing
    end
    sentinal = params_search_nt.planets[planet_key][param_name]
    return findfirst(==(sentinal), params)
end

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

using OptimizationOptimJL, OptimizationBBO

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
function default_initializer!(rng::Random.AbstractRNG, model::LogDensityModel; nruns=8, ntries=2, ndraws=1000, verbosity=1)
    ldm_any = LogDensityModelAny(model)

    # Pathfinder (and especially multipathfinder) do not work well with global optimization methods.
    # Instead, we do a two-step process. 
    # Find the global MAP point, then initialize multi-pathfinder in Gaussian ball around that point.
    verbosity > 0 && @info "Performing a global optimization to search for a starting position, bounded to the 0.1% to 99.9% percentiles of the priors."

    priors = Octofitter._list_priors(model.system)
    lb = model.link(quantile.(priors,0.001))
    ub = model.link(quantile.(priors,0.999))
    f = OptimizationFunction(
        (u,p)->-p.ℓπcallback(u),
        grad=(G,u,p)->(G .= p.∇ℓπcallback(u)[2])
    )
    prob = Optimization.OptimizationProblem(f, model.link(quantile.(priors,0.5)), model; lb, ub)
    Random.seed!(rand(rng, UInt64))
    sol = solve(prob, BBO_adaptive_de_rand_1_bin(), rel_tol=1e-3, maxiters = 1_000_000, )
    
    model.starting_points = fill(sol.u, 1000)
    initial_logpost_range = (-sol.objective, -sol.objective)
    if verbosity > 1
        @info "Found the global maximum logpost" MAP=-sol.objective
    end

    # TODO: we don't really need to use pathfinder in this case, we should look into
    # more rigourous variational methods
    local result_pf = nothing
    verbosity >= 1 && @info "Determining initial positions using pathfinder, around that location."
    # It can sometimes hit a PosDefException sometimes when factoring a matrix.
    # When that happens, the next try usually succeeds.
    try
        for i in 1:ntries
            verbosity >= 3 && @info "Starting multipathfinder run"
            init_sampler = function(rng, x)
                for _ in 1:10
                    # take a random step away from the MAP value according to the gradient
                    x .= sol.u -0.1 .* rand.(rng) .* model.∇ℓπcallback(sol.u)[2] 
                    if all(lb .< x .< ub) && isfinite(model.ℓπcallback(x))
                        return
                    end
                end
                error("Could not find starting point within 0.001 - 0.999 quantiles of priors.")
            end
            errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
            initial_mt = _kepsolve_use_threads[]
            _kepsolve_use_threads[] = nruns == 1
            result_pf = with_logger(errlogger) do 
                result_pf = Pathfinder.multipathfinder(
                    ldm_any, ndraws;
                    nruns,
                    init_sampler,
                    progress=verbosity > 1,
                    maxiters=25_000,
                    reltol=1e-6,
                    rng=rng,
                    ntries=1,
                    executor=Pathfinder.Transducers.PreferParallel(),
                    optimizer=Pathfinder.Optim.BFGS(;
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
    catch err
        @warn err
    end

    if !isnothing(result_pf)
        model.starting_points = collect.(eachcol(result_pf.draws))
    end
    logposts = model.ℓπcallback.(model.starting_points)
    initial_logpost_range = extrema(logposts)

    if initial_logpost_range[2] < -sol.objective - 10 
        if verbosity >= 1
            @warn "Pathfinder produced samples with log-likelihood 10 worse than global max. Will just initialize at global max."
        end
        model.starting_points = fill(sol.u, 1000)
        logposts = fill(-sol.objective, 1000)
        initial_logpost_range = (-sol.objective, -sol.objective)
    else
        if verbosity >= 1
            @info "Found a sample of initial positions" initial_logpost_range
        end
    end

    return model.arr2nt(model.invlink(model.starting_points[argmax(logposts)]))
end


# Helper function for testing that the pathfinder initialization gives reasonable results
function _startingpoints2chain(model)
    solnts =  [(;logpost=0,model.arr2nt(model.invlink(s))...,) for s in model.starting_points]
    chn = Octofitter.result2mcmcchain(solnts,  Dict(:internals => [:logpost]))
    return chn
end