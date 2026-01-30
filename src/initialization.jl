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
If `model.starting_points` has not yet been set, it first runs `Octofitter.initialize!(model)`.
"""
function get_starting_point!!(rng::Random.AbstractRNG, model::LogDensityModel)
    if isnothing(model.starting_points)
        initialize!(rng, model)
    end
    return rand(rng, model.starting_points)
end
function get_starting_point!!(model::LogDensityModel; kwargs...)
    return get_starting_point!!(Random.default_rng(), model; kwargs...)
end

using OptimizationOptimJL, OptimizationBBO


_gradorder(::LogDensityProblems.LogDensityOrder{K}) where K = K

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
    initialize!(model::LogDensityModel, fixed_params=nothing; kwargs...)

Initialize the model with optional fixed parameters provided as a named tuple.
Fixed parameters will be held constant during optimization and sampling.

The `fixed_params` can include:
- System-level variables: `(; plx=24.4, pmra=10.2, ...)`
- Planet variables: `(; planets=(; b=(; a=1.5, e=0.1, ...), ...))`  
- Observation variables: `(; observations=(; ObsName=(; var1=val1, var2=val2, ...), ...))`

Available keyword arguments include:
 - `verbosity=1`: control extra logging, can be 0 for silent, up to 4 for debugging info
 - `pathfinder_autodiff=AutoForwardDiff()`: what autodiff backend to use for initialization (not necessarily the same one used for the model in general)
 - `nruns=8`: how many runs of multi-pathfinder to use 
 - `ntries=2`: how many times can pathfinder fail and restart
 - `ndraws=1000`: how many draws to return from the pathfinder approximation

Example:
```julia
init_chain = initialize!(model, (;
    plx=24.4,
    pmra=10.2,
    planets=(;
        b=(;
            a=1.5,
            e=0.1,
        )
    ),
    observations=(;
        GaiaRV=(;
            offset_gaiarv=-50.0,
            jitter_gaiarv=0.1,
        ),
        GaiaDR4=(;
            astrometric_jitter=0.05,
        )
    )
))
```
"""
function initialize!(model::LogDensityModel, fixed_params=nothing; kwargs...)
    return initialize!(Random.default_rng(), model, fixed_params; kwargs...)
end

function initialize!(rng::Random.AbstractRNG, 
    model::LogDensityModel, 
    fixed_params=nothing; 
    nruns=8, 
    ntries=2, 
    ndraws=1000, 
    verbosity=1,
    pathfinder_autodiff=AutoForwardDiff(),
    use_pathfinder=true,
)
    # # If no fixed parameters, just call the original initializer
    # if isnothing(fixed_params)
    #     return initialize!(rng, model; nruns, ntries, ndraws, verbosity)
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
    
    start_time = fill(time(), 1) 
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
        if isempty(continuous_indices)
            # Store the results
            model.starting_points = fill(model.link(bestparams), ndraws)
            logposts = fill(bestlogpost, ndraws)
        else
            # Optimize continuous parameters while keeping discrete and fixed parameters constant
            # De-duplicate in cases where we have *fixed* values provided for discrete variables
            combined_fixed_indices = vcat(discrete_indices, fixed_indices)
            combined_fixed_values = vcat(bestparams[discrete_indices], fixed_values)
            d = Dict(zip(combined_fixed_indices,combined_fixed_values))
            combined_fixed_indices = sort(collect(keys(d)))
            combined_fixed_values = map(k->getindex(d,k), combined_fixed_indices)
            
            # bestparams, bestlogpost = optimize_continuous_params(
            #     rng, model, bestparams, combined_fixed_values, combined_fixed_indices, continuous_indices;
            #     verbosity
            # )

            logposts = optimization_and_pathfinder_with_fixed(
                rng, model, combined_fixed_values, combined_fixed_indices, variable_indices;
                nruns, ntries, ndraws, verbosity, pathfinder_autodiff, use_pathfinder
            )
        end
        
        
    else
        # Can use global optimization followed by pathfinder
        logposts = optimization_and_pathfinder_with_fixed(
            rng, model, fixed_values, fixed_indices, variable_indices;
            nruns, ntries, ndraws, verbosity, pathfinder_autodiff, use_pathfinder
        )
    end
    stop_time = fill(time(), 1)
    
    starting_point_chain = _startingpoints2chain(model)

    # Store meta data, including flag that this was sampled just from pathfinder,
    # not an mcmc
    starting_point_chain_with_info = MCMCChains.setinfo(
        starting_point_chain,
        (;
            start_time,
            stop_time,
            model_name=Symbol("$(model.system.name)-init"),
            sampler="pathfinder"
        )
    )

    return starting_point_chain_with_info
end



"""
Extract fixed parameters from a partial named tuple and return their values and indices.
Supports system variables, planet variables, and observation variables.
"""
function extract_fixed_params(model, partial_nt)
    # Create dummy parameter array to construct full named tuple
    dummy_params = model.sample_priors(Random.default_rng())
    full_nt = model.arr2nt(dummy_params)
    
    fixed_values = []
    fixed_indices = Int[]
    
    function process_tuple!(current_partial, current_full, path="")
        for name in propertynames(current_partial)
            val = getproperty(current_partial, name)
            # Special handling for observations
            if name == :observations && val isa NamedTuple
                process_observations!(val, current_full, path * ".observations")
                continue
            end
            
            if !hasproperty(current_full, name)
                error("Could not find parameter $name in model. You can only provide free parameters (e.g `x ~ ...`) and not derived parameters (e.g. `x = `). You also cannot provide values for variables that are e.g. `Ω ~ UniformCircular()`. You can instead supply `Ωx` and `Ωy`, or replace the distribution with `Uniform(0,2pi)`.")
                continue
            end
            
            full_val = getproperty(current_full, name)
            
            if val isa NamedTuple
                # Recursive case: nested named tuple
                process_tuple!(val, full_val, path * "." * String(name))
            else
                # Find index in the flat array using sentinel values
                sentinel_value = full_val[1]
                sentinal_idx = findfirst(x -> x == sentinel_value, dummy_params)
                
                if sentinal_idx !== nothing
                    push!(fixed_values, val)
                    push!(fixed_indices, sentinal_idx)
                else
                    error("Could not find parameter $name in model. You can only provide free parameters (e.g `x ~ ...`) and not derived parameters (e.g. `x = `). You also cannot provide values for variables that are e.g. `Ω ~ UniformCircular()`. You can instead supply `Ωx` and `Ωy`, or replace the distribution with `Uniform(0,2pi)`.")
                end
            end
        end
    end
    
    function process_observations!(obs_partial, current_full, path)
        # Process observation variables: observations=(; ObsName=(; var1=val1, var2=val2, ...), ...)
        for obs_name in propertynames(obs_partial)
            obs_vars = getproperty(obs_partial, obs_name)
            
            if !(obs_vars isa NamedTuple)
                error("Expected observation variables for $obs_name to be a NamedTuple, got $(typeof(obs_vars))")
            end
            
            # Normalize the observation name to match how it's stored in the model
            normalized_obs_name = Octofitter.normalizename(string(obs_name))
            
            # For observation variables, we need to look them up in the full_nt.observations
            for var_name in propertynames(obs_vars)
                var_value = getproperty(obs_vars, var_name)
                
                # Try to find the variable in observations
                found = false
                if hasproperty(current_full, :observations)
                    # Look for the normalized observation name first
                    if hasproperty(current_full.observations, normalized_obs_name)
                        full_obs_vars = getproperty(current_full.observations, normalized_obs_name)
                        
                        # Look for the exact variable name
                        if hasproperty(full_obs_vars, var_name)
                            sentinel_value = getproperty(full_obs_vars, var_name)[1]
                            sentinal_idx = findfirst(x -> x == sentinel_value, dummy_params)
                            
                            if sentinal_idx !== nothing
                                push!(fixed_values, var_value)
                                push!(fixed_indices, sentinal_idx)
                                found = true
                            end
                        end
                        
                        # Also try prefixed versions (e.g., likelihoodname_varname)
                        if !found
                            var_name_str = string(var_name)
                            prefixed_name = Symbol(string(normalized_obs_name) * "_" * var_name_str)
                            
                            if hasproperty(full_obs_vars, prefixed_name)
                                sentinel_value = getproperty(full_obs_vars, prefixed_name)[1]
                                sentinal_idx = findfirst(x -> x == sentinel_value, dummy_params)
                                
                                if sentinal_idx !== nothing
                                    push!(fixed_values, var_value)
                                    push!(fixed_indices, sentinal_idx)
                                    found = true
                                end
                            end
                        end
                    end
                    
                    # If not found with normalized name, try all observation names
                    if !found
                        for full_obs_name in propertynames(current_full.observations)
                            full_obs_vars = getproperty(current_full.observations, full_obs_name)
                            # Look for the exact variable name
                            if hasproperty(full_obs_vars, var_name)
                                sentinel_value = getproperty(full_obs_vars, var_name)[1]
                                sentinal_idx = findfirst(x -> x == sentinel_value, dummy_params)
                                if sentinal_idx !== nothing
                                    push!(fixed_values, var_value)
                                    push!(fixed_indices, sentinal_idx)
                                    found = true
                                    break
                                end
                            end
                        end
                    end
                end
                
                if !found
                    available_obs = if hasproperty(current_full, :observations)
                        [(obs_name, propertynames(getproperty(current_full.observations, obs_name))) for obs_name in propertynames(current_full.observations)]
                    else
                        "none"
                    end
                    error("Could not find observation variable $var_name for observation $obs_name (normalized: $normalized_obs_name) in model. Available observations and their variables: $available_obs")
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
            val = fixed_values[i]
            if val isa Number
                params[idx] = val
            else
                params[idx:(idx+length(val)-1)] .= val
            end
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
Run optimization and pathfinder with fixed parameters.
"""
function optimization_and_pathfinder_with_fixed(
    rng, model, fixed_values, fixed_indices, variable_indices;
    nruns=8, ntries=2, ndraws=1000, verbosity=1,
    pathfinder_autodiff=AutoForwardDiff(),
    use_pathfinder=true,
)
    # First, identify which parameters are discrete vs continuous
    sample = model.sample_priors(rng)
    discrete_indices = findall(x -> x isa Integer, sample)
    # continuous_indices = findall(x -> x isa AbstractFloat, sample)
    
    # For BBO: keep all specified parameters fixed
    # For pathfinder: only keep discrete parameters fixed
    discrete_fixed_indices = intersect(fixed_indices, discrete_indices)
    # discrete_fixed_values = fixed_values[indexin(discrete_fixed_indices, fixed_indices)]
    
    if verbosity > 2
        n_freed = length(setdiff(fixed_indices, discrete_fixed_indices))
        @info "Parameter summary" total_fixed=length(fixed_indices) discrete_fixed=length(discrete_fixed_indices) continuous_free=n_freed
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
    
    # function reduced_∇ℓπcallback(reduced_params)
    #     full_params = reduced_to_full(reduced_params)
        
    #     # Get full gradient
    #     logpost, full_grad = model.∇ℓπcallback(full_params)
        
    #     # Extract gradient components for variable parameters
    #     reduced_grad = full_grad[variable_indices]
        
    #     return logpost, reduced_grad
    # end

    # Get prior bounds for continuous parameters
    samples = stack(model.link.([model.sample_priors(rng) for _ in 1:10000]))
    lb_full = minimum(samples,dims=2)[:]
    ub_full = maximum(samples,dims=2)[:]
    
    # Extract bounds for continuous parameters
    lb = convert(Vector{Float64},lb_full[variable_indices])
    ub = convert(Vector{Float64},ub_full[variable_indices])

    reduced_initial_prime, bestlogpost = guess_starting_position_with_fixed(
        rng, model, fixed_values, fixed_indices, discrete_indices
    )
    reduced_initial = model.link(reduced_initial_prime)[variable_indices]
    if !isfinite(reduced_ℓπcallback(reduced_initial))
        opt_full = reduced_to_full(reduced_initial)
        @show opt_full model.invlink(opt_full) model.arr2nt(model.invlink(opt_full))
        error("Starting point for global optimization is not finite; did you fix parameters to invalid values?")
    end
    
    if isempty(reduced_initial)
        opt_full = reduced_to_full(reduced_initial)
        model.starting_points = fill(opt_full, ndraws)
        initial_logpost = model.ℓπcallback(opt_full)
    else

        # Define optimization function for reduced parameters
        f = OptimizationFunction(
            (u, p) -> -reduced_ℓπcallback(u),
            # grad = (G, u, p) -> begin
            #     _, grad = reduced_∇ℓπcallback(u)
            #     G .= -grad  # Negate for minimization
            # end
        )
        
        if verbosity > 0
            @info "Starting values not provided for all parameters! Guessing starting point using global optimization:" num_params=length(variable_indices) num_fixed=length(fixed_indices)
        end

        # Set up and solve reduced optimization problem
        reduced_initial = convert(Vector{Float64}, reduced_initial)
        lb = convert(Vector{Float64}, lb)
        ub = convert(Vector{Float64}, ub)
        if !all(lb .< reduced_initial .< ub)
            @warn "The initial guess parameters fell outside the 0.00001 or 0.9999 quantile range of the priors, so global optimization search is disabled. Returning initial guess only, which may be far from the global mode."
            opt_full = reduced_to_full(reduced_initial)
            initial_logpost = model.ℓπcallback(opt_full)
        else
            prob = Optimization.OptimizationProblem(f, reduced_initial, nothing; lb, ub)
            Random.seed!(rand(rng, UInt64))
            sol = solve(prob, BBO_adaptive_de_rand_1_bin(), rel_tol=1e-3, maxiters=1_000_000, show_trace=verbosity>2, show_every=1000)

            verbosity > 2 && display(sol.original)
            
            # Convert solution to full parameter space
            opt_full = reduced_to_full(sol.u)
                    
            initial_logpost = -sol.objective
        end
        # Initialize starting points around optimum
        model.starting_points = fill(opt_full, ndraws)

        full_params = opt_full  # Use opt_full which is set in both branches
    end

    if !isfinite(initial_logpost)
        error("Could not find finite starting point")
    end
    
    if verbosity > 1
        @info "Found global maximum logpost" MAP=initial_logpost
        println(opt_full)
    end
    
    # -------------------------------------------------------------
    # For pathfinder: create new conversion functions that only fix discrete parameters
    # -------------------------------------------------------------
    
    local result_pf = nothing
    local draws_full = nothing
    
    if use_pathfinder
    
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
            full_params = zeros(eltype(reduced_params), model.D)
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
        
        # function pathfinder_∇ℓπcallback(reduced_params)
        #     full_params = pathfinder_reduced_to_full(reduced_params)
            
        #     # Get full gradient
        #     logpost, full_grad = model.∇ℓπcallback(full_params)
            
        #     # Extract gradient components for pathfinder variable parameters
        #     reduced_grad = full_grad[pathfinder_variable_indices]
            
        #     return logpost, reduced_grad
        # end

        autodiff_type(model::LogDensityModel{D,Tℓπ,T∇ℓπ,TSys,TLink,TInvLink,TArr2nt,TPriSamp,ADType}) where {D,Tℓπ,T∇ℓπ,TSys,TLink,TInvLink,TArr2nt,TPriSamp,ADType} = ADType
        
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
            pathfinder_autodiff
            # this is because users might pick finitediff for models with discrete variables
            # but we still want to optimize efficiently here where we've masked out all the discrete variables
        )
        
        # Convert the BBO optimum to pathfinder reduced space
        pathfinder_opt = zeros(length(pathfinder_variable_indices))
        for (i, idx) in enumerate(pathfinder_variable_indices)
            pathfinder_opt[i] = opt_full[idx]
        end
        

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
                    # log like of about 200
                    target_drop=200.0 + 10randn(rng)
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
                    i < ntries && verbosity > 2 && @warn "Restarting pathfinder from new starting points" i
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
    end
    
    # Use pathfinder results if available
    if !isnothing(result_pf) && !isnothing(draws_full)
        verbosity > 2 && @info "Saving starting positions"
        model.starting_points = collect.(eachrow(draws_full))
    end
    
    # Calculate log posterior for starting points
    logposts = model.ℓπcallback.(model.starting_points)
    
    # Check if pathfinder produced good results
    if maximum(logposts) < initial_logpost - 20
        if verbosity >= 1
            @warn "Pathfinder produced samples with log-likelihood 20 worse than global max. Will just initialize at global max." logpost_range=extrema(logposts) mean_logpost=mean(logposts) initial_logpost
        end
        model.starting_points = fill(opt_full, ndraws)
        logposts = fill(initial_logpost, ndraws)
    else
        if verbosity >= 1
            @info "Found sample of initial positions" logpost_range=extrema(logposts) mean_logpost=mean(logposts)
        end
    end

    # Ensure discrete variables are still discrete (not promited to Float64)
    s = model.sample_priors(rng)
    model.starting_points = map(model.starting_points) do θ
        convert.(typeof.(s), θ)
    end
    
    return logposts
end


default_initializer! = initialize!
export initialize!


