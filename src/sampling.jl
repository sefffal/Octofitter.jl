using DiffResults
using LinearAlgebra
using Preferences
using Pathfinder
using Transducers
using CovarianceEstimation
export sample_priors

sample_priors(arg::Union{Planet,System,<:LogDensityModel}, args...; kwargs...) = sample_priors(Random.default_rng(), arg, args...; kwargs...)
# Sample priors from system once
function sample_priors(rng::Random.AbstractRNG, system::System)
    # Collect all priors in the same order as _list_priors and make_arr2nt
    all_priors = []
    
    # System priors
    for prior_distribution in values(system.priors.priors)
        push!(all_priors, prior_distribution)
    end
    
    # System observation priors
    for obs in system.observations
        if hasproperty(obs, :priors)
            for prior_distribution in values(obs.priors.priors)
                push!(all_priors, prior_distribution)
            end
        end
    end
    
    # Planet priors
    for planet in system.planets
        for prior_distribution in values(planet.priors.priors)
            push!(all_priors, prior_distribution)
        end
        
        # Planet observation priors
        for obs in planet.observations
            if hasproperty(obs, :priors)
                for prior_distribution in values(obs.priors.priors)
                    push!(all_priors, prior_distribution)
                end
            end
        end
    end
    
    # Sample from all priors
    priors_flat_sampled = map(prior -> rand(rng, prior), all_priors)
    return priors_flat_sampled
end
# Sample priors from system many times
sample_priors(rng::Random.AbstractRNG, system::System, N::Number) = [sample_priors(rng, system) for _ in 1:N]
sample_priors(rng::Random.AbstractRNG, model::LogDensityModel, N::Number) = [model.sample_priors(rng) for _ in 1:N]




# """
#     construct_elements(θ_system, θ_planet)

# Given a named tuple for of parameters from a System (θ_system) and Planet (θ_planet),
# return a `Visual{KepOrbit} PlanetOrbits.jl.
# """
# function construct_elements(::Type{Visual{KepOrbit}}, θ_system, θ_planet)
#     return Visual{KepOrbit}(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{AbsoluteVisual{KepOrbit}}, θ_system, θ_planet)
#     @show "dec"
#     return AbsoluteVisual{KepOrbit}(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{KepOrbit}, θ_system, θ_planet)
#     @show "abc"
#     return KepOrbit(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{ThieleInnesOrbit}, θ_system, θ_planet)
#     return ThieleInnesOrbit(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{RadialVelocityOrbit}, θ_system, θ_planet)
#     return RadialVelocityOrbit(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{CartesianOrbit}, θ_system, θ_planet)
#     return CartesianOrbit(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{Visual{CartesianOrbit}}, θ_system, θ_planet)
#     return Visual{CartesianOrbit}(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{FixedPosition}, θ_system, θ_planet)
#     return FixedPosition(merge(θ_system,θ_planet))
# end
# function construct_elements(::Type{Visual{FixedPosition}}, θ_system, θ_planet)
#     return Visual{FixedPosition}(merge(θ_system,θ_planet))
# end

"""
    construct_elements(chains, :b, 4)

Given a Chains object, a symbol matching the name of a planet, and an index,
construct a PlanetOrbits.jl orbit object.
"""
function construct_elements(model::LogDensityModel, chain::Chains, planet_key::Union{String,Symbol}, i)
    nts = mcmcchain2result(model,chain,i)
    return map(nts) do nt
        θ_planet = getproperty(nt.planets, planet_key)
        return orbit(;merge(nt,θ_planet)...)
    end
end
function construct_elements(model::LogDensityModel, chain::Chains, planet_key::Union{String,Symbol}, i::Number)
    nt = mcmcchain2result(model,chain,i)
    return construct_elements(nt, planet_key)
end
function construct_elements(result::NamedTuple, planet_key::Union{String,Symbol})
    θ_planet = getproperty(result.planets, planet_key)
    return orbit(;merge(result,θ_planet)...)
end
construct_elements(model::LogDensityModel, chain::Chains, planet::Planet, args...; kwargs...) = construct_elements(model, chain, planet.name, args...; kwargs...) 


# Fallback when no random number generator is provided (as is usually the case)
Base.@nospecializeinfer function advancedhmc(model::LogDensityModel, target_accept::Number=0.8; kwargs...)
    return advancedhmc(Random.default_rng(), model, target_accept; kwargs...)
end

"""
    octofit(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel
        target_accept::Number=0.8,
        ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
        adaptation,
        iterations,
        drop_warmup=true,
        max_depth=12,
        initial_samples= pathfinder ? 500 : 250_000,  # deprecated
        initial_parameters=nothing, # deprecated
        step_size=nothing,
        verbosity=2,
    )

Sample from the posterior defined by `model` using Hamiltonian Monte Carlo with 
the No U-Turn Sampler from AdvancedHMC.jl.
"""
Base.@nospecializeinfer function octofit(args...; kwargs...)
    return advancedhmc(args...; kwargs...)
end
export octofit

# Define some wrapper functions that hide type information
# so that we don't have to recompile pathfinder() with each 
# new model.
# It's worth it for the sampler, but not pathfinder.
struct CallableAny
    func::Function
end
(ca::CallableAny)(args...;kwargs...) = ca.func(args...;kwargs...)

struct LogDensityModelAny
    ldm::LogDensityModel
end
LogDensityProblems.logdensity(ldm_any::LogDensityModelAny, θ) = LogDensityProblems.logdensity(ldm_any.ldm, θ)
LogDensityProblems.logdensity_and_gradient(ldm_any::LogDensityModelAny, θ) = LogDensityProblems.logdensity_and_gradient(ldm_any.ldm, θ)
LogDensityProblems.dimension(ldm_any::LogDensityModelAny) = LogDensityProblems.dimension(ldm_any.ldm)
LogDensityProblems.capabilities(ldm_any::Type{<:LogDensityModelAny}) = LogDensityProblems.capabilities(ldm_any.ldm)



"""
The method signature of Octofitter.hmc is as follows:

    chain = advancedhmc(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel
        target_accept::Number=0.8,
        adaptation=1000,
        iterations=1000,
        drop_warmup=true,
        max_depth=12,
    )

Sample from the posterior defined by `model` using Hamiltonian Monte Carlo with 
the No U-Turn Sampler from AdvancedHMC.jl.
"""
Base.@nospecializeinfer function advancedhmc(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    model::LogDensityModel,
    target_accept::Number=0.8;
    adaptation::Int=1000,
    iterations::Int=1000,
    drop_warmup::Bool=true,
    max_depth::Int=12,
    verbosity::Int=2,
)
    @nospecialize
    if adaptation < 1000
        @warn "At least 1000 steps of adapation are recomended for good sampling"
    end

    # inialize if not already done or set by user
    get_starting_point!!(model)

    local metric = nothing
    for diag_eps in [0; 10.0 .^ range(-8, 0)]
        try
            # This can fail, triggering an exception
            S =  (cov(SimpleCovariance(), stack(model.starting_points)'))
            metric = DenseEuclideanMetric(S .+ Diagonal(diag_eps .* ones(model.D)))
            break
        catch err
            continue
        end
    end
    if isnothing(metric)
        verbosity > 1 && @warn("Falling back to initializing the diagonals with the prior interquartile ranges.")
        # We already sampled from the priors earlier to get the starting positon.
        # Use those variance estimates and transform them into the unconstrainted space.
        # variances_t = (model.link(initial_θ .+ sqrt.(variances)/2) .- model.link(initial_θ .- sqrt.(variances)/2)).^2
        # p = _list_priors(model.system)
        samples = eachrow(stack(Octofitter.sample_priors(model, 1000)))
        variances_t = (model.link(quantile.(samples, 0.85)) .- model.link(quantile.(samples, 0.15))).^2
        # metric = DenseEuclideanMetric(model.D)
        metric = DenseEuclideanMetric(collect(Diagonal(variances_t)))
        if verbosity >= 3
            print("Initial mass matrix M⁻¹ from priors\n")
            display(metric.M⁻¹)
        end
        if any(v->!isfinite(v)||v==0, variances_t)
            error("failed to initialize mass matrix")
        end
    end

    initial_θ_t = rand(rng, model.starting_points)
    initial_θ = model.invlink(initial_θ_t)

    if verbosity >= 4
        @info "flat starting point" initial_θ
        @info "transformed flat starting point" initial_θ_t
    end

    verbosity >= 3 && @info "Creating hamiltonian"
    hamiltonian = Hamiltonian(metric, model)
    verbosity >= 3 && @info "Finding good stepsize"
    ϵ = find_good_stepsize(rng, hamiltonian, initial_θ_t)
    verbosity >= 3 && @info "Found initial stepsize" ϵ 

    # Create integrator
    integrator = Leapfrog(ϵ)
    # integrator = JitteredLeapfrog(ϵ, 0.05) # 5% normal distribution on step size to help in areas of high curvature. 

    verbosity >= 3 && @info "Creating kernel"
    kernel = HMCKernel(Trajectory{MultinomialTS}(integrator, GeneralisedNoUTurn(;max_depth)))

    
    verbosity >= 3 && @info "Creating adaptor"
    # if isnothing(result_pf)
    # if metric isa Pathfinder.RankUpdateEuclideanMetric
        # adaptor = StepSizeAdaptor(target_accept, integrator)
    # else
        mma = MassMatrixAdaptor(metric)
        ssa = StepSizeAdaptor(target_accept, integrator)
        adaptor = StanHMCAdaptor(mma, ssa) 
    # end


    # Turn on likelihood parallelism if we have ~15x more data than threads.
    # This is based on some local benchmarks. Spawning tasks takes about 450ns;
    # an orbit solve takes about 32ns, or 1/14 as long.
    threads_avail = Threads.nthreads()
    n_epochs = _count_epochs(model.system) 
    Octofitter._kepsolve_use_threads[] = threads_avail > 1  && n_epochs > 15
    if verbosity >= 1
        @info "Kepler solver will use multiple threads: $(Octofitter._kepsolve_use_threads[] )" threads_avail n_epochs > 15
    end

    initial_parameters = initial_θ_t

    verbosity >= 1 && @info "Sampling, beginning with adaptation phase..."
    start_time = fill(time(), 1)
    mc_samples, stats = AdvancedHMC.sample(
        rng,
        hamiltonian,
        kernel,
        initial_parameters,
        iterations+adaptation,
        adaptor,
        adaptation;
        verbose=verbosity>=4,
        progress=verbosity>=1,
        drop_warmup,
    )
    stop_time = fill(time(), 1)
    
    verbosity >= 1 && @info "Sampling compete. Building chains."

    # Rebuild just the likelihood function (should already be compiled anyways)
    ln_like = make_ln_like(model.system, model.arr2nt(initial_θ))

    # Go through chain and repackage results
    numerical_error = getproperty.(stats, :numerical_error)
    actual_tree_depth = getproperty.(stats, :tree_depth)
     
    mean_accept = mean(getproperty.(stats, :acceptance_rate))
    ratio_divergent_transitions = mean(numerical_error)
    mean_tree_depth = mean(actual_tree_depth)
    max_tree_depth_frac = mean(==(max_depth), actual_tree_depth)

    total_steps = round(Int,sum(getproperty.(stats, :n_steps))*(iterations+adaptation)/iterations)

    elapsed = only(stop_time) - only(start_time)
    verbosity >= 1 && println("""
    Sampling report for chain:
    mean_accept                 = $mean_accept
    ratio_divergent_transitions = $ratio_divergent_transitions
    mean_tree_depth             = $mean_tree_depth
    max_tree_depth_frac         = $max_tree_depth_frac
    total_steps                 = $total_steps
    μs/step (approx.)           = $(round(elapsed/total_steps*1e6,sigdigits=3))\
    """)

    # Report some warnings if sampling did not work well
    if ratio_divergent_transitions == 1.0
        @error "Numerical errors encountered in ALL iterations. Check model and priors."
    elseif ratio_divergent_transitions > 0.05
        @warn "Numerical errors encountered in more than 5% of iterations. Results are likely incorrect. Your model might have high curvature, and could be improved. Otherwise, increasing the target acceptance rate (second argument to `octofit`) might help" ratio_divergent_transitions
    end
    # if max_tree_depth_frac > 0.1
    #     @warn "Maximum tree depth hit in more than 10% of iterations (reduced efficiency)" max_tree_depth_frac
    # end

    # Resolve the array back into the nested named tuple structure used internally.
    # Augment with some internal fields
    chain_res = map(zip(stats, mc_samples)) do (stat, θ_t)
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        loglike = ln_like(model.system, resolved_namedtuple)
        return merge((;
            stat.n_steps,
            stat.is_accept,
            stat.acceptance_rate,
            stat.hamiltonian_energy,
            stat.hamiltonian_energy_error,
            stat.max_hamiltonian_energy_error,
            stat.tree_depth,
            stat.numerical_error,
            stat.step_size,
            stat.nom_step_size,
            stat.is_adapt,
            loglike = loglike,
            logpost = stat.log_density,
        ), resolved_namedtuple)
    end
    # Then finally flatten and convert into an MCMCChain object / table.
    # Mark the posterior, likelihood, numerical error flag, and tree depth as internal
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :n_steps
            :is_accept
            :acceptance_rate
            :hamiltonian_energy
            :hamiltonian_energy_error
            :max_hamiltonian_energy_error
            :tree_depth
            :numerical_error
            :step_size
            :nom_step_size
            :is_adapt
            :loglike
            :logpost
            :tree_depth
            :numerical_error
        ])
    )

    # Concatenate the log posteriors and make them the same shape as the chains (N_iters,N_vars,N_chains)
    # logposts_mat = reduce(hcat, logposts)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
            samples_transformed=mc_samples,
            adaptor,
            initial_metric=metric,
            model_name=model.system.name,
            sampler="nuts"
        )
    )
    return mcmcchains_with_info
end

# Helper function for displaying nested named tuples in a compact format.
function stringify_nested_named_tuple(val::Any) # fallback
    string(val)
end
function stringify_nested_named_tuple(num::Number)
    string(round(num,digits=1))
end
function stringify_nested_named_tuple(nt::NamedTuple)
    "(;"*join(map(keys(nt)) do k
        "$k="*stringify_nested_named_tuple(nt[k])
    end, ", ")*")"
end

"""
Convert a vector of component arrays returned from sampling into an MCMCChains.Chains
object.

!!! warning
    Currently any nested named tuples must appear in the final position ie.
    `(;a=1,b=2,c=(;d=1,e=2))`.
"""
function result2mcmcchain(chain_in, sectionmap=Dict())
    # There is a specific column name convention used by MCMCChains to indicate
    # that multiple parameters form a group. Instead of planets.X.a, we adapt our to X_a 
    # accordingly (see flatten_named_tuple)
    flattened_labels = keys(flatten_named_tuple(first(chain_in)))
    data = zeros(length(chain_in), length(flattened_labels))
    
    for (i, sample) in enumerate(chain_in)
        # Instead of using nested Iterators.flatten, we need to extract values
        # in the same order as flatten_named_tuple creates keys
        vals = Float64[]
        
        # System-level variables
        for key in keys(sample)
            if key in (:planets, :observations)
                continue
            end
            if sample[key] isa Number
                push!(vals, sample[key])
            else
                for val in sample[key]
                    push!(vals, val)
                end
            end
        end
        
        # System observations
        for obs in keys(get(sample, :observations, (;)))
            for key in keys(sample.observations[obs])
                if sample.observations[obs][key] isa Number
                    push!(vals, sample.observations[obs][key])
                else
                    for val in sample.observations[obs][key]
                        push!(vals, val)
                    end
                end
            end
        end
        
        # Planet-level variables and their observations
        for pl in keys(get(sample, :planets, (;)))
            for key in keys(sample.planets[pl])
                if key == :observations
                    continue
                end
                if sample.planets[pl][key] isa Number
                    push!(vals, sample.planets[pl][key])
                else
                    for val in sample.planets[pl][key]
                        push!(vals, val)
                    end
                end
            end
            
            # Planet observations
            for obs in keys(get(sample.planets[pl], :observations, (;)))
                for key in keys(sample.planets[pl].observations[obs])
                    if sample.planets[pl].observations[obs][key] isa Number
                        push!(vals, sample.planets[pl].observations[obs][key])
                    else
                        for val in sample.planets[pl].observations[obs][key]
                            push!(vals, val)
                        end
                    end
                end
            end
        end
        
        # Fill the data matrix
        for (j, val) in enumerate(vals)
            data[i,j] = val
        end
    end
    
    c = Chains(data, [string(l) for l in flattened_labels], sectionmap)
    return c
end


"""
    mcmcchain2result(model, chain_in,)

Does the opposite of result2mcmcchain: given a model and a chain, return a vector of named tuples.
"""
function mcmcchain2result(model, chain, ii=(:))

    # Quickly construct a named tuple template
    θ = model.sample_priors(Random.default_rng())
    nt = model.arr2nt(θ)

    planetkeys = string.(keys(model.system.planets))
    
    # Get observation keys
    sysobs_keys = string.(keys(get(nt, :observations, (;))))
    planet_obs_keys = Dict{String, Vector{String}}()
    for pl in keys(get(nt, :planets, (;)))
        planet_obs_keys[string(pl)] = collect(string.(keys(get(nt.planets[pl], :observations, (;)))))
    end

    # Map output keys in the named tuple to one or more input keys in the chain
    # Complicated because in the chain representation, array-valued variables get
    # flattened out; e.g. (;a=1,b=(1,2,3)) becomes four columns [a, b_1, b_2, b_3]
    key_mapping = Pair{Symbol,Vector{Symbol}}[]
    
    # System variables
    for key in keys(nt)
        if key in (:planets, :observations)
            continue
        end
        if nt[key] isa Number
            push!(key_mapping, key => [key])
        else
            arr = Symbol[]
            push!(key_mapping, key => arr)
            for i in eachindex(nt[key])
                key_i = Symbol(key, '_', i)
                push!(arr, key_i)
            end
        end
    end
    
    # System observations
    for obs in keys(get(nt, :observations, (;)))
        for key in keys(nt.observations[obs])
            if nt.observations[obs][key] isa Number
                k = Symbol(obs, '_', key)
                push!(key_mapping, k => [k])
            else
                arr = Symbol[]
                push!(key_mapping, Symbol(obs, '_', key) => arr)
                for i in eachindex(nt.observations[obs][key])
                    key_i = Symbol(obs, '_', key, '_', i)
                    push!(arr, key_i)
                end
            end
        end
    end
    
    # Planet variables
    for pl in keys(get(nt, :planets, (;)))
        for key in keys(nt.planets[pl])
            if key == :observations
                continue
            end
            if nt.planets[pl][key] isa Number
                k = Symbol(pl, '_', key)
                push!(key_mapping, k => [k])
            else
                arr = Symbol[]
                push!(key_mapping, Symbol(pl, '_', key) => arr)
                for i in eachindex(nt.planets[pl][key])
                    key_i = Symbol(pl, '_', key, '_', i)
                    push!(arr, key_i)
                end
            end
        end
        
        # Planet observations
        for obs in keys(get(nt.planets[pl], :observations, (;)))
            for key in keys(nt.planets[pl].observations[obs])
                if nt.planets[pl].observations[obs][key] isa Number
                    k = Symbol(pl, '_', obs, '_', key)
                    push!(key_mapping, k => [k])
                else
                    arr = Symbol[]
                    push!(key_mapping, Symbol(pl, '_', obs, '_', key) => arr)
                    for i in eachindex(nt.planets[pl].observations[obs][key])
                        key_i = Symbol(pl, '_', obs, '_', key, '_', i)
                        push!(arr, key_i)
                    end
                end
            end
        end
    end
    
    # These are the labels corresponding to the flattened named tuple without the *planet_key* prepended
    IIs = broadcast(1:size(chain,1),(1:size(chain,3))') do i,j
        return (i,j)
    end
    
    function reform((i,j))
        # System variables
        nt_sys = Dict{Symbol,Any}()
        for (kout,kins) in key_mapping
            # Skip if this is a planet or observation key
            if any(map(pk->startswith(string(kout),pk*"_"), planetkeys)) || 
               any(map(ok->startswith(string(kout),ok*"_"), sysobs_keys))
                continue
            end
            if length(kins) == 1
                if haskey(chain, kins[])
                    nt_sys[kout] = chain[i,kins[],j]
                else
                    nt_sys[kout] = missing
                end
            else
                nt_sys[kout] = [
                    chain[i,kin,j]
                    for kin in kins
                ]
            end
        end
        
        # System observations
        nt_observations = map(collect(sysobs_keys)) do ok
            nt_obs = Dict{Symbol,Any}()
            for (kout,kins) in key_mapping
                if !startswith(string(kout),ok*"_")
                    continue
                end
                # Skip if this is actually a planet observation
                is_planet_obs = false
                for pk in planetkeys
                    if startswith(string(kout),pk*"_"*ok*"_")
                        is_planet_obs = true
                        break
                    end
                end
                if is_planet_obs
                    continue
                end
                
                kout_clean = Symbol(replace(string(kout), r"^"*string(ok)*"_" =>""))
                if length(kins) == 1
                    if haskey(chain, kins[])
                        nt_obs[kout_clean] = chain[i,kins[],j]
                    else
                        nt_obs[kout_clean] = missing
                    end
                else
                    nt_obs[kout_clean] = [
                        chain[i,kin,j]
                        for kin in kins
                    ]
                end
            end
            return Symbol(ok) => namedtuple(nt_obs)
        end
        
        # Planets and their observations
        nt_planets = map(collect(planetkeys)) do pk
            nt_pl = Dict{Symbol,Any}()
            nt_pl_obs = Dict{Symbol,Dict}()
            
            for (kout,kins) in key_mapping
                if !startswith(string(kout),pk*"_")
                    continue
                end
                
                # Check if this is a planet observation
                is_obs = false
                obs_name = ""
                for ok in get(planet_obs_keys, pk, String[])
                    if startswith(string(kout),pk*"_"*ok*"_")
                        is_obs = true
                        obs_name = ok
                        break
                    end
                end
                
                if is_obs
                    # This is a planet observation variable
                    kout_clean = Symbol(replace(string(kout), r"^"*string(pk)*"_"*obs_name*"_" =>""))
                    if !haskey(nt_pl_obs, Symbol(obs_name))
                        nt_pl_obs[Symbol(obs_name)] = Dict{Symbol,Any}()
                    end
                    obs_dict = nt_pl_obs[Symbol(obs_name)]
                    if length(kins) == 1
                        if haskey(chain, kins[])
                            obs_dict[kout_clean] = chain[i,kins[],j]
                        else
                            obs_dict[kout_clean] = missing
                        end
                    else
                    obs_dict[kout_clean] = [
                            chain[i,kin,j]
                            for kin in kins
                        ]
                    end
                else
                    # This is a regular planet variable
                    kout_clean = Symbol(replace(string(kout), r"^"*string(pk)*"_" =>""))
                    if length(kins) == 1
                        if haskey(chain, kins[])
                            nt_pl[kout_clean] = chain[i,kins[],j]
                        else
                            nt_pl[kout_clean] = missing
                        end
                    else
                            nt_pl[kout_clean] = [
                            chain[i,kin,j]
                            for kin in kins
                        ]
                    end
                end
            end
            # Convert observation dicts to named tuples
            if isempty(nt_pl_obs)
                nt_pl[:observations] = (;)
            else
                nt_pl[:observations] = namedtuple(Dict(
                    k => namedtuple(v)
                    for (k,v) in nt_pl_obs
                ))
            end
            
            return Symbol(pk) => namedtuple(nt_pl)
        end
        
        # Construct final result
        result = namedtuple(nt_sys)
        result = merge(result, (;observations=isempty(nt_observations) ? (;) : namedtuple(nt_observations)))
        result = merge(result, (;planets=isempty(nt_planets) ? (;) : namedtuple(nt_planets)))
        return result
    end
    
    if ii isa Number
        return reform(IIs[ii])
    else
        return broadcast(reform, IIs[ii])
    end
end



# Used for flattening a nested named tuple posterior sample into a flat named tuple
# suitable to be used as a Tables.jl table row.
# Used for flattening a nested named tuple posterior sample into a flat named tuple
# suitable to be used as a Tables.jl table row.
function flatten_named_tuple(nt)
    pairs = Pair{Symbol, Float64}[]
    
    # System-level variables
    for key in keys(nt)
        if key in (:planets, :observations)
            continue
        end
        if nt[key] isa Number
            push!(pairs, key => nt[key])
        else
            for i in eachindex(nt[key])
                key_i = Symbol(key, '_', i)
                push!(pairs, key_i => nt[key][i])
            end
        end
    end
    
    # System observations
    for obs in keys(get(nt, :observations, (;)))
        for key in keys(nt.observations[obs])
            if nt.observations[obs][key] isa Number
                push!(pairs, Symbol(obs, '_', key) => nt.observations[obs][key])
            else
                for i in eachindex(nt.observations[obs][key])
                    push!(pairs, Symbol(obs, '_', key, '_', i) => nt.observations[obs][key][i])
                end
            end
        end
    end
    
    # Planet-level variables and their observations
    for pl in keys(get(nt, :planets, (;)))
        for key in keys(nt.planets[pl])
            if key == :observations
                continue
            end
            if nt.planets[pl][key] isa Number
                push!(pairs, Symbol(pl, '_', key) => nt.planets[pl][key])
            else
                for i in eachindex(nt.planets[pl][key])
                    push!(pairs, Symbol(pl, '_', key, '_', i) => nt.planets[pl][key][i])
                end
            end
        end
        
        # Planet observations
        for obs in keys(get(nt.planets[pl], :observations, (;)))
            for key in keys(nt.planets[pl].observations[obs])
                if nt.planets[pl].observations[obs][key] isa Number
                    push!(pairs, Symbol(pl, '_', obs, '_', key) => nt.planets[pl].observations[obs][key])
                else
                    for i in eachindex(nt.planets[pl].observations[obs][key])
                        push!(pairs, Symbol(pl, '_', obs, '_', key, '_', i) => nt.planets[pl].observations[obs][key][i])
                    end
                end
            end
        end
    end
    
    return namedtuple(pairs)
end

include("octoquick.jl")