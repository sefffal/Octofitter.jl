using MCMCTempering, AbstractMCMC



struct Joint{Tℓprior, Tℓll}
    ℓprior      :: Tℓprior
    ℓlikelihood :: Tℓll
end

function (joint::Joint)(θ)
    return joint.ℓprior(θ) .+ joint.ℓlikelihood(θ)
end


struct TemperedJoint{Tℓprior, Tℓll, T<:Real} # <: Function
    ℓprior      :: Tℓprior
    ℓlikelihood :: Tℓll
    β           :: Real
end
function TemperedJoint(a,b,c)
    return TemperedJoint{typeof(a), typeof(b), typeof(c)}(a,b,c)
end

function (tj::TemperedJoint)(θ)
    return tj.ℓprior(θ) .+ (tj.ℓlikelihood(θ) .* tj.β)
end



function MCMCTempering.make_tempered_model(
    model::DifferentiableDensityModel,
    β::Real
)
    # @show typeof(model.ℓπ.ℓprior) typeof(model.ℓπ.ℓlikelihood) typeof(β)
    # error()
    ℓπ_β = TemperedJoint(model.ℓπ.ℓprior, model.ℓπ.ℓlikelihood, β)
    ∂ℓπ∂θ_β = TemperedJoint(model.∂ℓπ∂θ.ℓprior, model.∂ℓπ∂θ.ℓlikelihood, β)
    model_β = DifferentiableDensityModel(ℓπ_β, ∂ℓπ∂θ_β)
    return model_β
end

function MCMCTempering.make_tempered_loglikelihood(
    model::DifferentiableDensityModel,
    β::Real
)
    function logπ(z)
        return model.ℓπ.ℓlikelihood(z) * β
    end
    return logπ
end

function MCMCTempering.get_params(trans::AdvancedHMC.Transition)
    return trans.z.θ
end

# Fallback when no random number generator is provided (as is usually the case)
function temperedhmc(system::System, target_accept=0.9; kwargs...)
    return temperedhmc(Random.default_rng(), system, target_accept; kwargs...)
end

function temperedhmc(
    rng::Random.AbstractRNG,
    system::System, target_accept=0.9;
    adaptation,
    iterations,
    tree_depth=10,
    include_adapatation=false,
    initial_samples=100_000,
    initial_parameters=nothing,

    step_size=nothing,
    progress=true
)

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(system)
    D = length(initial_θ_0)

    ln_prior = make_ln_prior(system)
    arr2nt = DirectDetections.make_arr2nt(system) 
    priors_vec = _list_priors(system)

    # Capture these variables in a let binding to improve performance
    # ln_post = let system=system, ln_prior=ln_prior, arr2nt=arr2nt
    #     function (θ)
    #         θ_res = arr2nt(θ)
    #         ll = ln_prior(θ) + ln_like(θ_res, system)
    #         return ll
    #     end
    # end
    ℓlikelihood = let system=system, arr2nt=arr2nt, priors_vec=priors_vec
        function (θ_t)
            # Transform back from the unconstrained support to constrained support for the likelihood function
            θ = Bijectors.invlink.(priors_vec, θ_t)
            θ_res = arr2nt(θ)
            return ln_like(θ_res, system)
        end
    end
    
    # model = AdvancedHMC.DifferentiableDensityModel(ℓπ, ForwardDiff)

    # Define the target distribution
    # ℓprior(θ) = ln_prior(θ)
    function ℓprior(θ_t)
        θ = Bijectors.invlink.(priors_vec, θ_t)
        ln_prior(θ)
    end
    ∂ℓprior∂θ(θ) = (ℓprior(θ), ForwardDiff.gradient(ℓprior, θ))
    ∂ℓlikelihood∂θ(θ) = (ℓlikelihood(θ), ForwardDiff.gradient(ℓlikelihood, θ))
    
    model = DifferentiableDensityModel(
        Joint(ℓprior, ℓlikelihood),
        Joint(∂ℓprior∂θ, ∂ℓlikelihood∂θ)
    )

    if isnothing(initial_parameters)
        initial_θ = guess_starting_position(system,initial_samples)
    else
        initial_θ = initial_parameters
    end
    

    # Transform from constrained support to unconstrained support
    initial_θ_t = Bijectors.link.(priors_vec, initial_θ)

    # Define a Hamiltonian system
    metric = DenseEuclideanMetric(D)
    
    hamiltonian = Hamiltonian(metric, model.ℓπ, model.∂ℓπ∂θ)

    if !isnothing(step_size)
        initial_ϵ = step_size
    else
        initial_ϵ = find_good_stepsize(hamiltonian, initial_θ_t)
    end

    integrator = Leapfrog(initial_ϵ)

    # # We have change some parameters when running with image data
    if !isnothing(system.images) && target_accept > 0.6 && isnothing(step_size)
        @warn "Sampling from images with target accept greater than 0.6. This may lead to insufficient exploration."
    end

    mma = MassMatrixAdaptor(metric)
    if isnothing(step_size)
        @info "Adapting step size and mass matrix"
        ssa = StepSizeAdaptor(target_accept, integrator)
        adaptor = StanHMCAdaptor(mma, ssa) 
    else
        @info "Adapting mass matrix only"
        adaptor = MassMatrixAdaptor(metric)
    end



    # κ = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    κ = NUTS(integrator, max_depth=tree_depth) 
    spl = AdvancedHMC.HMCSampler(κ, metric, adaptor)

    tempered_sampler = tempered(spl, 4)
    start_time = time()

    AbstractMCMC.setprogress!(true)
    # mc_samples = sample(
    #     model,
    #     tempered_sampler,
    #     iterations;
    #     discard_initial = adaptation,
    #     progress=progress,
    #     verbose=false
    # )
    mc_samples = sample(
        rng,
        model,
        tempered_sampler,
        iterations;
        nadapts = adaptation,
        init_params = initial_θ_t,
        discard_initial = adaptation,
        progress=progress,
        verbose=false
    )
    stop_time = time()

    stat = map(s->s.stat, mc_samples)
    samples = map(s->s.z.θ, mc_samples)
    logposterior = map(s->s.z.ℓπ.value, mc_samples)

    @info "Resolving deterministic variables"
    chain_res = arr2nt.(samples)
    @info "Constructing chains"
    mcmcchain = DirectDetections.result2mcmcchain(system, chain_res)
    mcmcchain_with_info = MCMCChains.setinfo(
        mcmcchain,
        (;start_time, stop_time, model=system)
    )

    mean_accept = mean(getproperty.(stat, :acceptance_rate))
    num_err_frac = mean(getproperty.(stat, :numerical_error))
    mean_tree_depth = mean(getproperty.(stat, :tree_depth))
    max_tree_depth_frac = mean(getproperty.(stat, :tree_depth) .== tree_depth)

    println("Sampling report:")
    @show mean_accept
    @show num_err_frac
    @show mean_tree_depth
    @show max_tree_depth_frac

    # Report some warnings if sampling did not work well
    if num_err_frac == 1.0
        @error "Numerical errors encountered in ALL iterations. Check model and priors."
    elseif num_err_frac > 0.1
        @warn "Numerical errors encountered in more than 10% of iterations" num_err_frac
    end
    if max_tree_depth_frac > 0.1
        @warn "Maximum tree depth hit in more than 10% of iterations (reduced efficiency)" max_tree_depth_frac
    end

    return (;
        chain=mcmcchain_with_info,
        stats=stat,
        logposterior,
        adaptor
    )
end
