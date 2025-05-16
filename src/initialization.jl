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
function guess_starting_position(rng::Random.AbstractRNG, model::LogDensityModel, N=500_000, )
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
                # println("logpost_after = ", logpost)
            if logpost > bestlogpost
                bestlogpost = logpost
                bestparams = params
            end        
        end
        return bestparams, bestlogpost
    end
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
_gradorder(::LogDensityProblems.LogDensityOrder{K}) where K = K
# function default_initializer!(rng::Random.AbstractRNG, model::LogDensityModel; nruns=8, ntries=2, ndraws=1000, verbosity=1)
#     order = _gradorder(LogDensityProblems.capabilities(model))
#     if order == 0
#         # Have to just guess by sampling
#         bestparams, bestlogpost = guess_starting_position(rng, model)
#         model.starting_points = fill(model.link(bestparams),ndraws)
#         logposts = fill(bestlogpost,ndraws)
#     else
#         # Can use global optimization followed by pathfinder
#         logposts = optimization_and_pathfinder_based_initializer(rng, model; nruns, ntries, ndraws, verbosity)
#     end
#     return model.arr2nt(model.invlink(model.starting_points[argmax(logposts)]))
# end

# function optimization_and_pathfinder_based_initializer(rng, model; nruns=8, ntries=2, ndraws=1000, verbosity=1)
#     # Pathfinder (and especially multipathfinder) do not work well with global optimization methods.
#     # Instead, we do a two-step process. 
#     # Find the global MAP point, then initialize multi-pathfinder in Gaussian ball around that point.
#     verbosity > 0 && @info "Performing a global optimization to search for a starting position, bounded to the 0.1% to 99.9% percentiles of the priors."

#     priors = Octofitter._list_priors(model.system)
#     lb = model.link(quantile.(priors,0.001))
#     ub = model.link(quantile.(priors,0.999))
#     f = OptimizationFunction(
#         (u,p)->-p.ℓπcallback(u),
#         grad=(G,u,p)->(G .= p.∇ℓπcallback(u)[2])
#     )
#     prob = Optimization.OptimizationProblem(f, model.link(model.sample(rng)), model; lb, ub)
#     Random.seed!(rand(rng, UInt64))
#     sol = solve(prob, BBO_adaptive_de_rand_1_bin(), rel_tol=1e-3, maxiters = 1_000_000, )
    
#     model.starting_points = fill(sol.u, 1000)
#     initial_logpost_range = (-sol.objective, -sol.objective)
#     if verbosity > 1
#         @info "Found the global maximum logpost" MAP=-sol.objective
#     end

#     # TODO: we don't really need to use pathfinder in this case, we should look into
#     # more rigourous variational methods
#     local result_pf = nothing
#     verbosity >= 1 && @info "Determining initial positions using pathfinder, around that location."
#     # It can sometimes hit a PosDefException sometimes when factoring a matrix.
#     # When that happens, the next try usually succeeds.
#     try
#         for i in 1:ntries
#             verbosity >= 3 && @info "Starting multipathfinder run"
#             init_sampler = function(rng, x)
#                 for _ in 1:10
#                     # take a random step away from the MAP value according to the gradient
#                     x .= sol.u -0.1 .* rand.(rng) .* model.∇ℓπcallback(sol.u)[2] 
#                     if all(lb .< x .< ub) && isfinite(model.ℓπcallback(x))
#                         return
#                     end
#                 end
#                 error("Within pathfinder, could not find a finite starting point within 0.001 - 0.999 quantiles of priors.")
#             end
#             errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
#             initial_mt = _kepsolve_use_threads[]
#             _kepsolve_use_threads[] = nruns == 1
#             result_pf = with_logger(errlogger) do 
#                 result_pf = Pathfinder.multipathfinder(
#                     model, ndraws;
#                     nruns,
#                     init_sampler,
#                     progress=verbosity > 1,
#                     maxiters=25_000,
#                     reltol=1e-6,
#                     rng=rng,
#                     ntries=1,
#                     executor=Pathfinder.Transducers.PreferParallel(),
#                     optimizer=Pathfinder.Optim.BFGS(;
#                         linesearch=Pathfinder.Optim.LineSearches.BackTracking(),
#                         alphaguess=Pathfinder.Optim.LineSearches.InitialHagerZhang()
#                     )
#                 ) 
#                 return result_pf
#             end
#             _kepsolve_use_threads[] = initial_mt
#             # Check pareto shape diagnostic to see if we have a good approximation
#             # If we don't, just try again
#             if result_pf.psis_result.pareto_shape > 3
#                 verbosity > 3 && display(result_pf)
#                 verbosity >= 4 && display(result_pf.psis_result)
#                 i<ntries && verbosity > 2 && @warn "Restarting pathfinder" i
#                 continue
#             end
            
#             verbosity >= 3 && "Pathfinder complete"
#             verbosity > 2 &&  display(result_pf)
#             break
#         end
#     catch err
#         @warn err
#     end

#     if !isnothing(result_pf)
#         model.starting_points = collect.(eachcol(result_pf.draws))
#     end
#     logposts = model.ℓπcallback.(model.starting_points)
#     initial_logpost_range = extrema(logposts)

#     if initial_logpost_range[2] < -sol.objective - 10 
#         if verbosity >= 1
#             @warn "Pathfinder produced samples with log-likelihood 10 worse than global max. Will just initialize at global max."
#         end
#         model.starting_points = fill(sol.u, 1000)
#         logposts = fill(-sol.objective, 1000)
#         initial_logpost_range = (-sol.objective, -sol.objective)
#     else
#         if verbosity >= 1
#             @info "Found a sample of initial positions" initial_logpost_range
#         end
#     end

#     return logposts
# end


# Helper function for testing that the pathfinder initialization gives reasonable results
function _startingpoints2chain(model)
    solnts =  [(;logpost=0,model.arr2nt(model.invlink(s))...,) for s in model.starting_points]
    chn = Octofitter.result2mcmcchain(solnts,  Dict(:internals => [:logpost]))
    return chn
end





function identify_param_types(model)
    # Sample from priors to get a representative parameter vector
    sample = model.sample_priors(Random.default_rng())
    
    # Identify discrete (integer) parameters
    discrete_indices = findall(x -> x isa Integer, sample)
    continuous_indices = findall(x -> x isa AbstractFloat, sample)
    
    return discrete_indices, continuous_indices
end

"""
    default_initializer!(model::LogDensityModel, fixed_params=nothing; kwargs...)

Initialize the model with optional fixed parameters provided as a named tuple.
Fixed parameters will be held constant during optimization and sampling.
"""
function default_initializer!(model::LogDensityModel, fixed_params=nothing; kwargs...)
    return default_initializer!(Random.default_rng(), model, fixed_params; kwargs...)
end

function default_initializer!(rng::Random.AbstractRNG, 
                                       model::LogDensityModel, 
                                       fixed_params=nothing; 
                                       nruns=8, 
                                       ntries=2, 
                                       ndraws=1000, 
                                       verbosity=1)
    # # If no fixed parameters, just call the original initializer
    # if isnothing(fixed_params)
    #     return default_initializer!(rng, model; nruns, ntries, ndraws, verbosity)
    # end
    
    # Extract fixed parameters and their indices
    fixed_values, fixed_indices = extract_fixed_params(model, fixed_params)
    
    if verbosity > 0 && !isempty(fixed_indices)
        @info "Initializing with $(length(fixed_indices)) fixed parameters"
    end
    
    # Create masks for fixed and variable parameters
    all_indices = collect(1:model.D)
    variable_indices = setdiff(all_indices, fixed_indices)
    
    # Check if we have mixed data types
    order = _gradorder(LogDensityProblems.capabilities(model))
    
    if order == 0
        # Identify discrete and continuous parameters
        sample = model.sample_priors(rng)
        discrete_indices = findall(x -> x isa Integer, sample)
        continuous_indices = findall(x -> x isa AbstractFloat, sample)
        
        # Adjust for already fixed parameters
        discrete_indices = setdiff(discrete_indices, fixed_indices)
        continuous_indices = setdiff(continuous_indices, fixed_indices)
        
        if verbosity > 0
            @info "Mixed parameter types detected" discrete=length(discrete_indices) continuous=length(continuous_indices) fixed=length(fixed_indices)
        end
        
        # 1. Sample discrete parameters and combine with fixed parameters
        verbosity > 0 && @info "Sampling for starting positions of discrete values"
        bestparams, bestlogpost = guess_starting_position_with_fixed(
            rng, model, fixed_values, fixed_indices, discrete_indices
        )
        
        # 2. If we have continuous parameters, optimize them
        if !isempty(continuous_indices)
            # Optimize continuous parameters while keeping discrete and fixed parameters constant
            combined_fixed_indices = union(fixed_indices, discrete_indices)
            combined_fixed_values = vcat(fixed_values, bestparams[discrete_indices])
            
            bestparams, bestlogpost = optimize_continuous_params(
                rng, model, bestparams, combined_fixed_values, combined_fixed_indices, continuous_indices;
                verbosity
            )
        end
        
        # Store the results
        model.starting_points = fill(model.link(bestparams), ndraws)
        logposts = fill(bestlogpost, ndraws)
    else
        # Can use global optimization followed by pathfinder
        logposts = optimization_and_pathfinder_with_fixed(
            rng, model, fixed_values, fixed_indices, variable_indices;
            nruns, ntries, ndraws, verbosity
        )
    end
    
    return model.arr2nt(model.invlink(model.starting_points[argmax(logposts)]))
end


"""
Extract fixed parameters from a partial named tuple and return their values and indices.
"""
function extract_fixed_params(model, partial_nt)
    # Create dummy parameter array to construct full named tuple
    dummy_params = model.sample_priors(Random.default_rng())
    full_nt = model.arr2nt(dummy_params)
    
    fixed_values = Float64[]
    fixed_indices = Int[]
    
    function process_tuple!(current_partial, current_full, path="")
        for name in propertynames(current_partial)
            if !hasproperty(current_full, name)
                @warn "Parameter $name at $path not found in model parameters"
                continue
            end
            
            val = getproperty(current_partial, name)
            full_val = getproperty(current_full, name)
            
            if val isa NamedTuple
                # Recursive case: nested named tuple
                process_tuple!(val, full_val, path * "." * String(name))
            else
                # Find index in the flat array using sentinel values
                sentinel_value = full_val
                sentinal_idx = findfirst(x -> x == sentinel_value, dummy_params)
                
                if sentinal_idx !== nothing
                    push!(fixed_values, val)
                    push!(fixed_indices, sentinal_idx)
                else
                    error("Could not find parameter $name in model priors. Is it a free parameter in your model? Derived parameters (`x = ...`) cannot be set directly, only sampled parameters (`x ~ ...`). `UniformCircular()` values also cannot be specified manually, change them to be `Uniform(0,2pi)`")
                end
            end
        end
    end
    
    # Start processing
    process_tuple!(partial_nt, full_nt)
    
    return fixed_values, fixed_indices
end

"""
Sample from priors with fixed parameters.
"""
function guess_starting_position_with_fixed(
    rng, model, fixed_values, fixed_indices, discrete_indices=Int[];
    N=10_000, enable_ofti=true
)
    # Sample parameters with fixed values inserted
    function sample_with_fixed()
        params = model.sample_priors(rng)
        # Insert fixed values
        for (i, idx) in enumerate(fixed_indices)
            params[idx] = fixed_values[i]
        end
        return params
    end
    
    # If we only care about discrete indices beyond the fixed ones
    function sample_with_fixed_and_discrete()
        params = model.sample_priors(rng)
        # Insert fixed values
        for (i, idx) in enumerate(fixed_indices)
            params[idx] = fixed_values[i]
        end
        # Only keep the discrete parameters we care about
        # (the rest will be overwritten during optimization)
        return params
    end
    
    # Choose the appropriate sampling function
    sample_fn = isempty(discrete_indices) ? sample_with_fixed : sample_with_fixed_and_discrete
    
    # Initialize with a valid sample
    bestparams = sample_fn()
    params_t = model.link(bestparams)
    bestlogpost = model.ℓπcallback(params_t)
    
    # Sample and keep the best
    for _ in 1:N
        params = sample_fn()
        params_t = model.link(params)
        logpost = model.ℓπcallback(params_t)
        
        # Optionally try an OFTI step
        if enable_ofti
            params_ofti = ofti_step(rng, model, params)
            # Make sure fixed values remain fixed
            for (i, idx) in enumerate(fixed_indices)
                params_ofti[idx] = fixed_values[i]
            end
            params_ofti_t = model.link(params_ofti)
            logpost_ofti = model.ℓπcallback(params_ofti_t)
            if logpost_ofti > logpost
                logpost = logpost_ofti
                params = params_ofti
            end
        end
        
        if logpost > bestlogpost
            bestlogpost = logpost
            bestparams = params
        end
    end
    
    return bestparams, bestlogpost
end
"""
Extract fixed parameters from a partial named tuple and return their values and indices.
"""
function extract_fixed_params(model, partial_nt)
    # Create dummy parameter array to construct full named tuple
    dummy_params = model.sample_priors(Random.default_rng())
    full_nt = model.arr2nt(dummy_params)
    
    fixed_values = Float64[]
    fixed_indices = Int[]
    
    function process_tuple!(current_partial, current_full, path="")
        for name in propertynames(current_partial)
            if !hasproperty(current_full, name)
                @warn "Parameter $name at $path not found in model parameters"
                continue
            end
            
            val = getproperty(current_partial, name)
            full_val = getproperty(current_full, name)
            
            if val isa NamedTuple
                # Recursive case: nested named tuple
                process_tuple!(val, full_val, path * "." * String(name))
            else
                # Find index in the flat array using sentinel values
                sentinel_value = full_val
                sentinal_idx = findfirst(x -> x == sentinel_value, dummy_params)
                
                if sentinal_idx !== nothing
                    push!(fixed_values, val)
                    push!(fixed_indices, sentinal_idx)
                else
                    @warn "Could not locate parameter $name at $path in flat array"
                end
            end
        end
    end
    
    # Start processing
    process_tuple!(partial_nt, full_nt)
    
    return fixed_values, fixed_indices
end

"""
Sample from priors with fixed parameters.
"""
function guess_starting_position_with_fixed(
    rng, model, fixed_values, fixed_indices, discrete_indices=Int[];
    N=10_000, 
)
    # Sample parameters with fixed values inserted
    function sample_with_fixed()
        params = collect(model.sample_priors(rng))
        # Insert fixed values
        for (i, idx) in enumerate(fixed_indices)
            params[idx] = fixed_values[i]
        end
        return params
    end
    
    # Initialize with a valid sample
    bestparams = sample_with_fixed()
    params_t = model.link(bestparams)
    bestlogpost = model.ℓπcallback(params_t)
    
    # Sample and keep the best
    for _ in 1:N
        params = sample_with_fixed()
        params_t = model.link(params)
        logpost = model.ℓπcallback(params_t)
        
       
        if logpost > bestlogpost
            bestlogpost = logpost
            bestparams = params
        end
    end
    
    return bestparams, bestlogpost
end

"""
Optimize continuous parameters while keeping discrete and fixed parameters constant.
"""
function optimize_continuous_params(
    rng, model, initial_params, fixed_values, fixed_indices, continuous_indices;
    verbosity=1
)
    # Create a mapping function to go from reduced space to full space
    function reduced_to_full(reduced_params)
        full_params = copy(initial_params)
        # Fill in the continuous parameters
        for (i, idx) in enumerate(continuous_indices)
            full_params[idx] = reduced_params[i]
        end
        # Ensure fixed parameters stay fixed
        for (i, idx) in enumerate(fixed_indices)
            full_params[idx] = fixed_values[i]
        end
        return full_params
    end
    
    # Create wrapper callbacks for reduced parameter space
    function reduced_ℓπcallback(reduced_params)
        full_params = reduced_to_full(reduced_params)
        transformed = model.link(full_params)
        return model.ℓπcallback(transformed)
    end
    
    function reduced_∇ℓπcallback(reduced_params)
        full_params = reduced_to_full(reduced_params)
        transformed = model.link(full_params)
        
        # Get full gradient
        logpost, full_grad = model.∇ℓπcallback(transformed)
        
        # Extract only the components for continuous parameters
        reduced_grad = zeros(length(continuous_indices))
        for (i, idx) in enumerate(continuous_indices)
            reduced_grad[i] = full_grad[idx]
        end
        
        return logpost, reduced_grad
    end
    
    # Extract initial values for continuous parameters
    reduced_initial = initial_params[continuous_indices]
    
    # Create reduced optimization problem
    f = OptimizationFunction(
        (u, p) -> -reduced_ℓπcallback(u),
        grad = (G, u, p) -> begin
            _, grad = reduced_∇ℓπcallback(u)
            G .= -grad  # Negate for minimization
        end
    )
    
    # Get prior bounds for continuous parameters
    priors = Octofitter._list_priors(model.system)
    lb_full = model.link([eltype(p) <: Integer ? NaN : quantile(p, 0.001) for p in priors])
    ub_full = model.link([eltype(p) <: Integer ? NaN : quantile(p, 0.999) for p in priors])
    
    # Extract bounds for continuous parameters
    lb = lb_full[continuous_indices]
    ub = ub_full[continuous_indices]
    
    # Create and solve reduced problem
    prob = Optimization.OptimizationProblem(f, reduced_initial, nothing; lb, ub)
    Random.seed!(rand(rng, UInt64))
    
    if verbosity > 0
        @info "Optimizing $(length(continuous_indices)) continuous parameters"
    end
    
    sol = solve(prob, BBO_adaptive_de_rand_1_bin(), rel_tol=1e-3, maxiters=1_000_000)
    
    # Convert optimal solution back to full parameter space
    optimal_params = reduced_to_full(sol.u)
    optimal_logpost = -sol.objective
    
    if verbosity > 0
        @info "Continuous parameter optimization complete" logpost=optimal_logpost
    end
    
    return optimal_params, optimal_logpost
end

"""
Run optimization and pathfinder with fixed parameters.
"""
function optimization_and_pathfinder_with_fixed(
    rng, model, fixed_values, fixed_indices, variable_indices;
    nruns=8, ntries=2, ndraws=1000, verbosity=1
)
    # First, identify which parameters are discrete vs continuous
    sample = model.sample_priors(rng)
    discrete_indices = findall(x -> x isa Integer, sample)
    continuous_indices = findall(x -> x isa AbstractFloat, sample)
    
    # For BBO: keep all specified parameters fixed
    # For pathfinder: only keep discrete parameters fixed
    discrete_fixed_indices = intersect(fixed_indices, discrete_indices)
    discrete_fixed_values = fixed_values[indexin(discrete_fixed_indices, fixed_indices)]
    
    if verbosity > 2
        n_freed = length(setdiff(fixed_indices, discrete_fixed_indices))
        @info "Parameter fixing strategy" total_fixed=length(fixed_indices) discrete_fixed=length(discrete_fixed_indices) continuous_freed=n_freed
    end
    
    # Define conversion functions between reduced and full parameter spaces for BBO
    fixed_values_linked = zeros(model.D)
    fixed_values_linked[fixed_indices] .= fixed_values
    fixed_values_linked = model.link(fixed_values_linked)
    
    function reduced_to_full(reduced_params)
        full_params = zeros(model.D)
        # Fill in variable parameters
        for (i, idx) in enumerate(variable_indices)
            try
                full_params[idx] = reduced_params[i]
            catch err
                rethrow(err)
            end
        end
        # Fill in fixed parameters
        for (i, idx) in enumerate(fixed_indices)
            full_params[idx] = fixed_values_linked[idx]
        end
        return full_params
    end
    
    # Create wrapper callbacks for reduced parameter space during BBO
    function reduced_ℓπcallback(reduced_params)
        full_params = reduced_to_full(reduced_params)
        return model.ℓπcallback(full_params)
    end
    
    function reduced_∇ℓπcallback(reduced_params)
        full_params = reduced_to_full(reduced_params)
        
        # Get full gradient
        logpost, full_grad = model.∇ℓπcallback(full_params)
        
        # Extract gradient components for variable parameters
        reduced_grad = full_grad[variable_indices]
        
        return logpost, reduced_grad
    end
    
    # Get bounds for the reduced parameters for BBO
    priors = Octofitter._list_priors(model.system)
    lb_full = model.link(quantile.(priors, 0.00001))
    ub_full = model.link(quantile.(priors, 0.99999))
    
    lb = lb_full[variable_indices]
    ub = ub_full[variable_indices]
    
    # Define optimization function for reduced parameters
    f = OptimizationFunction(
        (u, p) -> -reduced_ℓπcallback(u),
        grad = (G, u, p) -> begin
            _, grad = reduced_∇ℓπcallback(u)
            G .= -grad  # Negate for minimization
        end
    )
    
    if verbosity > 0
        @info "Performing global optimization. $(length(fixed_indices)) parameters are fixed by the user."
    end
    
    # Initial reduced parameters (sample until all are in bounds)
    local initial
    for i in 1:100
        initial = model.link(model.sample_priors(rng))
        if any(.!(lb_full .< initial .< ub_full))
            continue
        else
            break
        end
    end
    reduced_initial = initial[variable_indices]
    
    # Set up and solve reduced optimization problem
    prob = Optimization.OptimizationProblem(f, reduced_initial, nothing; lb, ub)
    Random.seed!(rand(rng, UInt64))
    sol = solve(prob, BBO_adaptive_de_rand_1_bin(), rel_tol=1e-3, maxiters=1_000_000)

    verbosity > 2 && display(sol.original)
    
    # Convert solution to full parameter space
    opt_full = reduced_to_full(sol.u)
    
    # Initialize starting points around optimum
    model.starting_points = fill(opt_full, ndraws)
    initial_logpost = -sol.objective

    full_params = reduced_to_full(sol.u)

    if !isfinite(initial_logpost)
        error("Could not find finite starting point")
    end
    
    if verbosity > 1
        @info "Found global maximum logpost" MAP=initial_logpost opt_full
    end
    
    # -------------------------------------------------------------
    # For pathfinder: create new conversion functions that only fix discrete parameters
    # -------------------------------------------------------------
    
    # Determine which indices are free for pathfinder (all except discrete fixed params)
    pathfinder_fixed_indices = discrete_fixed_indices
    pathfinder_variable_indices = setdiff(1:model.D, pathfinder_fixed_indices)
    
    # Define a new reduced to full mapping for pathfinder
    pathfinder_fixed_values_linked = zeros(model.D)
    for (i, idx) in enumerate(discrete_fixed_indices)
        # Find the original index in the fixed_indices array
        orig_idx = findfirst(x -> x == idx, fixed_indices)
        if orig_idx !== nothing
            pathfinder_fixed_values_linked[idx] = fixed_values_linked[idx]
        end
    end
    
    function pathfinder_reduced_to_full(reduced_params)
        full_params = zeros(model.D)
        # Fill in variable parameters (which now include continuous fixed params)
        for (i, idx) in enumerate(pathfinder_variable_indices)
            full_params[idx] = reduced_params[i]
        end
        # Fill in only the discrete fixed parameters
        for idx in pathfinder_fixed_indices
            full_params[idx] = pathfinder_fixed_values_linked[idx]
        end
        return full_params
    end
    
    # Create new callback functions for pathfinder
    function pathfinder_ℓπcallback(reduced_params)
        full_params = pathfinder_reduced_to_full(reduced_params)
        return model.ℓπcallback(full_params)
    end
    
    function pathfinder_∇ℓπcallback(reduced_params)
        full_params = pathfinder_reduced_to_full(reduced_params)
        
        # Get full gradient
        logpost, full_grad = model.∇ℓπcallback(full_params)
        
        # Extract gradient components for pathfinder variable parameters
        reduced_grad = full_grad[pathfinder_variable_indices]
        
        return logpost, reduced_grad
    end
    
    # Define optimization function for pathfinder
    pathfinder_f = OptimizationFunction(
        # (u, p) -> -pathfinder_ℓπcallback(u),
        (u, p) -> begin
            l = -pathfinder_ℓπcallback(u)
            if !isfinite(l)
                @warn "no finite " l maxlog=10
            end
            l
        end,
        grad = (G, u, p) -> begin
            _, grad = pathfinder_∇ℓπcallback(u)
            G .= -grad  # Negate for minimization
        end
    )
    
    # Convert the BBO optimum to pathfinder reduced space
    pathfinder_opt = zeros(length(pathfinder_variable_indices))
    for (i, idx) in enumerate(pathfinder_variable_indices)
        pathfinder_opt[i] = opt_full[idx]
    end
    
    # Run multipathfinder with only discrete parameters fixed
    local result_pf = nothing
    local draws_full = nothing
    
    try
        for i in 1:ntries
            verbosity >= 3 && @info "Starting multipathfinder run with only discrete parameters fixed" free_params=length(pathfinder_variable_indices) fixed_params=length(pathfinder_fixed_indices)
            
            errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Debug : Logging.Error)
            initial_mt = _kepsolve_use_threads[]
            _kepsolve_use_threads[] = nruns == 1

            # Generate initial values for pathfinder
            initial_vals = map(1:nruns) do _
                x = pathfinder_opt

                # take a random step away from the max like solution and aim for a drop in
                # log like of about 10
                target_drop=30.0 + 10randn(rng)
                max_iter=20

                d = randn(rng, length(x)) # Generate random unit vector
                d = d / norm(d)
                L₀ = pathfinder_ℓπcallback(x) # initial value
                α_min, α_max = 0.0, 0.1 # Binary search for step size
                # Expand upper bound if needed
                while (L₀ - pathfinder_ℓπcallback(x + α_max * d) < target_drop) && (α_max < 1e6)
                    α_max *= 2.0
                end
                
                new_x = x
                for _ in 1:max_iter # Binary search - always same number of iterations
                    α_mid = (α_min + α_max) / 2.0
                    new_x = x + α_mid * d
                    drop = L₀ - pathfinder_ℓπcallback(new_x)
                    if drop < target_drop
                        α_min = α_mid
                    else
                        α_max = α_mid
                    end
                    new_x = new_x
                end
                return new_x
            end
                
            result_pf, draws_full = with_logger(errlogger) do 
                # Use pathfinder-specific reduced model
                result_reduced = Pathfinder.multipathfinder(
                    pathfinder_f, ndraws;
                    nruns, init=initial_vals,
                    progress = verbosity > 1,
                    maxiters = 25_000,
                    reltol = 1e-6,
                    rng = rng,
                    ntries = 1,
                    executor = Pathfinder.Transducers.SequentialEx(),
                    # executor = Pathfinder.Transducers.PreferParallel(),
                    optimizer = Pathfinder.Optim.BFGS(;
                        linesearch = Pathfinder.Optim.LineSearches.BackTracking(),
                        alphaguess = Pathfinder.Optim.LineSearches.InitialHagerZhang()
                    )
                )
                verbosity > 1 && @info "Done pathfinder"
                
                # Convert reduced pathfinder result to full space
                draws_full = zeros(ndraws, model.D)
                for j in 1:ndraws
                    draws_full[j, :] .= pathfinder_reduced_to_full(result_reduced.draws[:, j])
                end
                
                return result_reduced, draws_full
            end
            
            _kepsolve_use_threads[] = initial_mt
            
            # Check pareto shape diagnostic
            if result_pf.psis_result.pareto_shape > 3
                verbosity > 3 && display(result_pf)
                verbosity >= 4 && display(result_pf.psis_result)
                i < ntries && verbosity > 2 && @warn "Restarting pathfinder" i
                continue
            end
            
            verbosity >= 3 && @info "Pathfinder complete"
            verbosity > 2 && display(result_pf)
            break
        end
    catch err
        @warn "Pathfinder error: $err"
        @warn "Using BBO result as fallback"
    end
    
    # Use pathfinder results if available
    if !isnothing(result_pf) && !isnothing(draws_full)
        model.starting_points = collect(eachrow(draws_full))
    end
    
    # Calculate log posterior for starting points
    logposts = model.ℓπcallback.(model.starting_points)
    
    # Check if pathfinder produced good results
    if maximum(logposts) < initial_logpost - 10
        if verbosity >= 1
            @warn "Pathfinder produced samples with log-likelihood 10 worse than global max. Will just initialize at global max."
        end
        model.starting_points = fill(opt_full, ndraws)
        logposts = fill(initial_logpost, ndraws)
    else
        if verbosity >= 1
            @info "Found sample of initial positions" logpost_range=extrema(logposts) mean_logpost=mean(logposts)
        end
    end
    
    return logposts
end