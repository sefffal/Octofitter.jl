

export sample_priors
sample_priors(planet::Planet) = rand.(ComponentArray(planet.priors.priors))
sample_priors(planet::Planet,N) = [sample_priors(planet) for _ in 1:N]


function sample_priors(system::System)
    sampled = ComponentVector(
        merge(NamedTuple(rand.(ComponentArray(system.priors.priors))),
        # (;planets=[sample_priors(planet) for planet in system.planets])
        (;planets=namedtuple(collect(keys(system.planets)), [
            ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(sample_priors(planet)))]))
            for planet in system.planets
        ]))
    ))
    return getdata(sampled)
end

sample_priors(system::System,N) = [sample_priors(system) for _ in 1:N]





# Instead of just calling mean for the distributions, we sample and then take the mean of the data.
# This does add a little jitter, but some distributions do not directly define the mean function!
# Specifically, truncated(InverseGamma()) does not work, and this is very useful.
# mean_priors(planet::Planet) = Statistics.mean.(Statistics.rand.(planet.priors.priors,1000))
function mean_priors(system::System)
    priors_all = ComponentVector(;
        NamedTuple(system.priors.priors)...,
        planets=[planet.priors.priors for planet in system.planets]
    )
    # return Statistics.mean.(Statistics.rand.(priors_all,1000))
    return Statistics.mean.(priors_all)
end


function guess_starting_position(system, N=500_000)

    # TODO: this shouldn't have to allocate anything, we can just loop keeping the best.
    θ = sample_priors(system, N)
    arr2nt = DirectDetections.make_arr2nt(system) 

    posts = zeros(N)
    ln_prior = make_ln_prior(system)
    Threads.@threads for i in eachindex(posts)
        θ_res = arr2nt(θ[i])
        posts[i] = ln_prior(θ[i]) + ln_like(θ_res, system)
    end
    # posts = map(eachrow(A)) do c
    #     DirectDetections.ln_post(ComponentVector(c, ax), system)
    # end
    mapv,mapi = findmax(posts)
    best = θ[mapi]
    
    # @info "Found good location" mapv best=NamedTuple(best)

    return best
end

using Optim
function optimize_starting_position(ℓπ, initial_θ_t)

    @info "Optimizing starting location "
    results = optimize(θ_t-> -ℓπ(θ_t), initial_θ_t, LBFGS(), autodiff = :forward, Optim.Options(iterations=100_000, time_limit=5, allow_f_increases=true))
    # results = optimize(θ_t-> -ℓπ(θ_t), initial_θ_t, NelderMead(), Optim.Options(iterations=100_000, time_limit=5))
    # results = optimize(θ_t-> -ℓπ(θ_t), initial_θ_t, ParticleSwarm(), Optim.Options(iterations=100_000, time_limit=5))

    intial_objective = ℓπ(initial_θ_t)

    minimizer = Optim.minimizer(results)

    @info "Starting location improved" logimprovement=(ℓπ(minimizer) - intial_objective)
    return minimizer
end


"""
    construct_elements(θ_system, θ_planet)

Given a named tuple for of parameters from a System (θ_system) and Planet (θ_planet),
return a `KeplerianElements` from DirectOrbits.jl.
"""
function construct_elements(θ_system, θ_planet)
    return KeplerianElements((;
        θ_system.μ,
        θ_system.plx,
        θ_planet.i,
        θ_planet.Ω,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    ))
end


"""
    construct_elements(chains, :b, 4)

Given a Chains object, a symbol matching the name of a planet, and an index,
construct a `KeplerianElements` from DirectOrbits of that planet from that
index of the chains.
"""
function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, i::Union{Integer,CartesianIndex})
    pk = string(planet_key)
    return KeplerianElements((;
        μ=chain["μ"][i],
        plx=chain["plx"][i],
        i=chain[pk*"[i]"][i],
        Ω=chain[pk*"[Ω]"][i],
        ω=chain[pk*"[ω]"][i],
        e=chain[pk*"[e]"][i],
        τ=chain[pk*"[τ]"][i],
        a=chain[pk*"[a]"][i],
    ))
end

"""
    construct_elements(chains, :b, [4,5,10])

Given a Chains object, a symbol matching the name of a planet, and an array of indices,
construct a `KeplerianElements` from DirectOrbits of that planet from those indices
of the chains.
"""
function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, ii::AbstractArray{<:Union{Integer,CartesianIndex}})
    pk = string(planet_key)
    μs=chain["μ"]
    plxs=chain["plx"]
    is=chain[pk*"[i]"]
    Ωs=chain[pk*"[Ω]"]
    ωs=chain[pk*"[ω]"]
    es=chain[pk*"[e]"]
    τs=chain[pk*"[τ]"]
    as=chain[pk*"[a]"]
    return map(ii) do i
        KeplerianElements((;
            μ=μs[i],
            plx=plxs[i],
            i=is[i],
            Ω=Ωs[i],
            ω=ωs[i],
            e=es[i],
            τ=τs[i],
            a=as[i],
        ))
    end
end




# Fallback when no random number generator is provided (as is usually the case)
function hmc(system::System, target_accept::Number=0.8, ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial(); kwargs...)
    return hmc(Random.default_rng(), system, target_accept, ensemble; kwargs...)
end

function hmc(
    rng::Random.AbstractRNG,
    system::System, target_accept::Number=0.8,
    ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
    num_chains=1,
    adaptation,
    iterations,
    discard_initial=adaptation,
    tree_depth=10,
    initial_samples=50_000,
    initial_parameters=nothing,
    step_size=nothing,
    progress=true
)

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(system)
    D = length(initial_θ_0)

    ln_prior_transformed = make_ln_prior_transformed(system)
    arr2nt = DirectDetections.make_arr2nt(system) 

    priors_vec = _list_priors(system)
    Bijector_invlinkvec = make_Bijector_invlinkvec(priors_vec)

    # Capture these variables in a let binding to improve performance
    ℓπ = let system=system, ln_prior_transformed=ln_prior_transformed, arr2nt=arr2nt, Bijector_invlinkvec=Bijector_invlinkvec
        function (θ_t)
            # Transform back from the unconstrained support to constrained support for the likelihood function
            θ = Bijector_invlinkvec(θ_t)
            θ_res = arr2nt(θ)
            ll = ln_prior_transformed(θ) + ln_like(θ_res, system)
            return ll
        end
    end

    if isnothing(initial_parameters)
        progress && @info "Guessing a good starting location by sampling from priors" initial_samples
        initial_θ = guess_starting_position(system,initial_samples)
        # Transform from constrained support to unconstrained support
        initial_θ_t = Bijectors.link.(priors_vec, initial_θ)

        # initial_θ_guess = guess_starting_position(system,initial_samples)
        # # Transform from constrained support to unconstrained support
        # initial_θ_guess_t = Bijectors.link.(priors_vec, initial_θ_guess)

        # # @show initial_θ_guess
        # initial_θ_t = optimize_starting_position(ℓπ, initial_θ_guess_t)
        
        # # Just for display
        # initial_θ = Bijectors.invlink.(priors_vec, initial_θ_t)

        # # @show initial_θ_guess initial_θ
        # # @show priors_vec

    else
        initial_θ = initial_parameters
        # Transform from constrained support to unconstrained support
        initial_θ_t = Bijectors.link.(priors_vec, initial_θ)
    end

    # Define a Hamiltonian system
    metric = DenseEuclideanMetric(D)
    hamiltonian = Hamiltonian(metric, ℓπ, ForwardDiff)

    if !isnothing(step_size)
        initial_ϵ = step_size
    else
        initial_ϵ = find_good_stepsize(hamiltonian, initial_θ_t)
        progress && @info "Found initial stepsize" initial_ϵ
    end


    integrator = Leapfrog(initial_ϵ)
    # integrator = TemperedLeapfrog(initial_ϵ, 1.05)


    # # We have change some parameters when running with image data
    if !isnothing(system.images) && target_accept > 0.6 && isnothing(step_size)
        @warn "Sampling from images with target accept greater than 0.6. This may lead to insufficient exploration."
    end

    mma = MassMatrixAdaptor(metric)
    if isnothing(step_size)
        progress && @info "Will adapt step size and mass matrix"
        ssa = StepSizeAdaptor(target_accept, integrator)
        adaptor = StanHMCAdaptor(mma, ssa) 
    else
        progress && @info "Will adapt mass matrix only"
        adaptor = MassMatrixAdaptor(metric)
    end

    model = AdvancedHMC.DifferentiableDensityModel(ℓπ, ForwardDiff)

    # κ = NUTS{MultinomialTS,GeneralisedNoUTurn}(integrator, max_depth=tree_depth) 
    # κ = NUTS{SliceTS, StrictGeneralisedNoUTurn}(integrator, max_depth=tree_depth) 
    

    # Had some good results with this one:
    # κ = NUTS{MultinomialTS, StrictGeneralisedNoUTurn}(integrator, max_depth=tree_depth) 

    κ = NUTS(integrator, max_depth=tree_depth) 
    sampler = AdvancedHMC.HMCSampler(κ, metric, adaptor)


    start_time = fill(time(), num_chains)

    # Neat: it's possible to return a live iterator
    # We could use this to build e.g. live plotting while the code is running
    # once the analysis code is ported to Makie.
    # return  AbstractMCMC.steps(
    #     rng,
    #     model,
    #     sampler,
    #     nadapts = adaptation,
    #     init_params = initial_θ_t,
    #     discard_initial = adaptation,
    #     progress=progress,
    #     verbose=false
    # )


    # function callback(rng, model, sampler, transition, state, iteration; kwargs...)
    #     # nadapts init_params
    #     @show propertynames(kwargs)
    #     @show iteration
    # end

    mc_samples_all_chains = sample(
        rng,
        model,
        sampler,
        ensemble,
        iterations+adaptation, num_chains;
        nadapts = adaptation,
        init_params = initial_θ_t,
        discard_initial,
        progress=progress,
        # callback
    )
    stop_time = fill(time(), num_chains)
    
    @info "Sampling compete. Building chains."
    # Go through each chain and repackage results
    chains = MCMCChains.Chains[]
    logposts = Vector{Float64}[]
    for (i,mc_samples) in enumerate(mc_samples_all_chains)
        stat = map(s->s.stat, mc_samples)
        logpost = map(s->s.z.ℓπ.value, mc_samples)
     
        mean_accept = mean(getproperty.(stat, :acceptance_rate))
        num_err_frac = mean(getproperty.(stat, :numerical_error))
        mean_tree_depth = mean(getproperty.(stat, :tree_depth))
        max_tree_depth_frac = mean(getproperty.(stat, :tree_depth) .== tree_depth)
    
        progress && println("""\
        Sampling report for chain $i:
        mean_accept =         $mean_accept
        num_err_frac =        $num_err_frac
        mean_tree_depth =     $mean_tree_depth
        max_tree_depth_frac = $max_tree_depth_frac
        """)

        # Report some warnings if sampling did not work well
        if num_err_frac == 1.0
            @error "Numerical errors encountered in ALL iterations. Check model and priors."
        elseif num_err_frac > 0.1
            @warn "Numerical errors encountered in more than 10% of iterations" num_err_frac
        end
        if max_tree_depth_frac > 0.1
            @warn "Maximum tree depth hit in more than 10% of iterations (reduced efficiency)" max_tree_depth_frac
        end

        logpost = map(s->s.z.ℓπ.value, mc_samples)
    
        # Transform samples back to constrained support
        samples = map(mc_samples) do s
            θ_t = s.z.θ
            θ = Bijectors.invlink.(priors_vec, θ_t)
            return θ
        end
        chain_res = arr2nt.(samples)
        push!(chains, DirectDetections.result2mcmcchain(system, chain_res))
        push!(logposts, logpost)
    end

    # Concatenate the independent chains now that we have remapped / resolved the variables.
    mcmcchains = AbstractMCMC.chainscat(chains...)

    # Concatenate the log posteriors and make them the same shape as the chains (N_iters,N_vars,N_chains)
    logposts_mat = reduce(hcat, logposts)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
            model=system,
            logpost=logposts_mat,
            _restart=(;
                model,
                sampler,
                adaptor,
                state = last.(mc_samples_all_chains)
            )
        )
    )
    return mcmcchains_with_info
end


include("tempered-sampling.jl")


"""
Convert a vector of component arrays returned from sampling into an MCMCChains.Chains
object.
"""
function result2mcmcchain(system, chains_in_0)
    chains_in = ComponentArray.(chains_in_0)
    # `system` not currently used, but a more efficient/robust mapping in future might require it.

    # There is a specific column name convention used by MCMCChains to indicate
    # that multiple parameters form a group. Instead of planets.X.a, we adapt our to X[a] 
    # accordingly
    flattened_labels = replace.(labels(first(chains_in)), r"planets\.([^\.]+).([^\.]+)" => s"\1[\2]")
    c = Chains(
        reduce(vcat, getdata.(chains_in)'),
        flattened_labels
    )
    return c
end
