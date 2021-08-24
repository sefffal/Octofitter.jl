
export sample_priors
sample_priors(planet::Planet) = ourrand.(planet.priors.priors)
sample_priors(planet::Planet,N) = ourrand.(planet.priors.priors,N)

priors_fixed(planet::Planet) = typeof.(planet.priors.priors) .<: Real

function sample_priors(planet::ReparameterizedPlanet3)
    sample = NamedTuple(sample_priors(planet.planet))

    plus = sample.ω + sample.Ω
    minus = sample.ω - sample.Ω
    reparameterized = merge(
        delete(sample,  (:ω, :Ω)),
        (ωΩ⁺ = plus, ωΩ⁻=minus)
    )
    
    # Φ = sample.ω + sample.Ω + 2π*sample.τ
    # Φω⁻ = Φ - sample.ω
    # ΦΩ⁻ = Φ - sample.Ω
    # reparameterized = merge(
    #     delete(sample,  (:τ, :ω, :Ω)),
    #     (;Φ, Φω⁻, ΦΩ⁻)
    # )

    return ComponentVector(reparameterized)
end
function sample_priors(planet::ReparameterizedPlanet3, N)
    sample = NamedTuple(sample_priors(planet.planet, N))
    plus = sample.ω .+ sample.Ω
    minus = sample.ω .- sample.Ω
    reparameterized = merge(
        delete(sample,  (:ω, :Ω)),
        (ωΩ⁺ = plus, ωΩ⁻=minus)
    )

    # Φ = sample.ω .+ sample.Ω .+ 2π*sample.τ
    # Φω⁻ = Φ .- sample.ω
    # ΦΩ⁻ = Φ .- sample.Ω
    # reparameterized = merge(
    #     delete(sample,  (:τ, :ω, :Ω)),
    #     (;Φ, Φω⁻, ΦΩ⁻)
    # )
    return ComponentVector(reparameterized)
end

ourrand(n::Real) = n
ourrand(n::Real,N) = fill(n,N)
ourrand(d::Any) = rand(d)
ourrand(d::Any,N) = rand(d,N)

function sample_priors(system::System)
    sampled = ComponentVector(
        merge(NamedTuple(ourrand.(system.priors.priors)),
        # (;planets=[sample_priors(planet) for planet in system.planets])
        (;planets=[
            ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(sample_priors(planet)))]))
            for planet in system.planets
        ])
    ))
    return sampled
end
function sample_priors(system::System,N)
    sampled = ComponentVector(
        merge(NamedTuple(ourrand.(system.priors.priors,N)),
        # (;planets=[sample_priors(planet,N) for planet in system.planets])
        (;planets=[
            ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(sample_priors(planet,N)))]))
            for planet in system.planets
        ])
    ))
    return sampled
end
function priors_fixed(system)
    sampled = ComponentVector(merge(
        NamedTuple(typeof.(system.priors.priors) .<: Real),
        (;planets=[
            ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(priors_fixed(planet)))]))
            for planet in system.planets
        ])
    ))
    return sampled
end
export mean_priors
# Instead of just calling mean for the distributions, we sample and then take the mean of the data.
# This does add a little jitter, but some distributions do not directly define the mean function!
# Specifically, truncated(InverseGamma()) does not work, and this is very useful.
# mean_priors(planet::Planet) = Statistics.mean.(Statistics.rand.(planet.priors.priors,1000))
# function mean_priors(system::System)
#     priors_all = ComponentVector(;
#         NamedTuple(system.priors.priors)...,
#         planets=[planet.priors.priors for planet in system.planets]
#     )
#     return Statistics.mean.(Statistics.rand.(priors_all,1000))
# end
# mean_priors(planet::Planet) = Statistics.mean.(Statistics.rand.(planet.priors.priors,1000))
function mean_priors(system::System)
    N = 5000
    sampled = ComponentVector(
        merge(NamedTuple(mean.(rand.(system.priors.priors,N))),
        # (;planets=[mean.(sample_priors(planet,N)) for planet in system.planets])
        (;planets=[
            ComponentArray(NamedTuple([k=>mean(v) for (k,v) in pairs(NamedTuple(sample_priors(planet,N)))]))
            for planet in system.planets
        ])
    ))
    return sampled
end


function guess_starting_position(system, N=500_000)

    @info "Guessing a good starting location by sampling from priors" N
    # TODO: this shouldn't have to allocate anything, we can just loop keeping the best.
    θ0 = sample_priors(system)
    θ = sample_priors(system, N)
    ax = getaxes(θ0)
    l = length(θ0)
    A = reshape(getdata(θ), :, l)
    posts = zeros(size(A,1))
    Threads.@threads for i in eachindex(posts)
        # posts[i] = DirectDetections.ln_post(ComponentVector(view(A,i,:), ax), system)
        posts[i] = DirectDetections.ln_post(ComponentVector(view(A,i,:), ax), system)
    end
    # posts = map(eachrow(A)) do c
    #     DirectDetections.ln_post(ComponentVector(c, ax), system)
    # end
    mapv,mapi = findmax(posts)
    best = ComponentArray(reshape(getdata(θ), :, l)[mapi,:], ax)
    
    @info "Found good location" mapv a=getproperty.(best.planets, :a)

    return best
end


using Optim
function find_starting_position(system)

    initial = guess_starting_position(system, 10_000)

    ax = getaxes(initial)
    function objective(θ_dat)
        θ = ComponentArray(θ_dat, ax)
        return -DirectDetections.ln_post(θ, system)
    end
    
    result = optimize(objective, getdata(initial), GradientDescent(), Optim.Options(show_trace=true, show_every=1000, iterations=100_000), autodiff=:forward)

    display(result)

    best = ComponentArray(Optim.minimizer(result), ax)
    
    @info "Found good location" a=getproperty.(best.planets, :a)

    return best
end



function get_ωΩτ(θ_system, θ_planet)
    if haskey(θ_planet, :ωΩ⁺)
        
        # Derivation:
        # ωΩ⁺ = ω+Ω
        # ωΩ⁻ = ω-Ω
        # ωΩ⁺ + ωΩ⁻ = ω+Ω + ω-Ω = ω
        # ωΩ⁺ - ωΩ⁻ = ω+Ω - ω+Ω = Ω
        ω = θ_planet.ωΩ⁺ + θ_planet.ωΩ⁻ 
        Ω = θ_planet.ωΩ⁺ - θ_planet.ωΩ⁻
        τ = θ_planet.τ
    elseif haskey(θ_planet, :Φ)
    
        # Φ = sample.ω .+ sample.Ω .+ 2π*sample.τ
        # Φω⁻ = Φ .- sample.ω
        # ΦΩ⁻ = Φ .- sample.Ω
    
        τ2pi = θ_planet.Φ + θ_planet.Φω⁻ + θ_planet.ΦΩ⁻
        ω = θ_planet.Φ - τ2pi + θ_planet.ΦΩ⁻
        Ω = θ_planet.Φ - τ2pi + θ_planet.Φω⁻
        τ = τ2pi/2π

    else
        ω = θ_planet.ω
        τ = θ_planet.τ

        
        if hasproperty(θ_planet,:Ω)
            Ω = θ_planet.Ω
        elseif hasproperty(θ_system,:Ω)
            Ω = θ_system.Ω
        else
            error("property `Ω` not specified for planet or system")
        end
    end
    return ω, Ω, τ
end


function construct_elements(θ_system, θ_planet)
    # Handle re-parameterized models
    ω, Ω, τ = get_ωΩτ(θ_system, θ_planet)

    if hasproperty(θ_planet,:i)
        i = θ_planet.i
    elseif hasproperty(θ_system,:i)
        i = θ_system.i
    else
        error("property `i` not specified for planet or system")
    end

    return KeplerianElements((;
        θ_system.μ,
        θ_system.plx,
        i,
        Ω,
        ω,
        θ_planet.e,
        τ,
        θ_planet.a,
    ))
end
function construct_elements(θ_system, θ_planet, i)
    # Handle re-parameterized models. See above.
    if haskey(θ_planet, :ωΩ⁺)
        ω = θ_planet.ωΩ⁺[i] + θ_planet.ωΩ⁻[i]
        Ω = θ_planet.ωΩ⁺[i] - θ_planet.ωΩ⁻[i]
        τ = θ_planet.τ[i]
    elseif haskey(θ_planet, :Φ)
        τ2pi = θ_planet.Φ[i] + θ_planet.Φω⁻[i] + θ_planet.ΦΩ⁻[i]
        ω = θ_planet.Φ[i] - τ2pi + θ_planet.ΦΩ⁻[i]
        Ω = θ_planet.Φ[i] - τ2pi + θ_planet.Φω⁻[i]
        τ = τ2pi/2π
    else
        ω = θ_planet.ω[i]
        if hasproperty(θ_planet,:Ω)
            Ω = θ_planet.Ω[i]
        elseif hasproperty(θ_system,:Ω)
            Ω = θ_system.Ω[i]
        else
            error("property `Ω` not specified for planet or system")
        end
        τ = θ_planet.τ[i]
    end
    if hasproperty(θ_planet,:i)
        inc = θ_planet.i[i]
    elseif hasproperty(θ_system,:i)
        inc = θ_system.i[i]
    else
        error("property `i` not specified for planet or system")
    end

    return KeplerianElements((;
        μ=θ_system.μ[i],
        plx=θ_system.plx[i],
        i=inc,
        Ω,
        ω,
        e=θ_planet.e[i],
        τ,
        a=θ_planet.a[i],
    ))
end


function mcmc(
    system::System;
    burnin,
    numwalkers,
    numsamples_perwalker,
    thinning = 1,
    squash = true,
)
    ln_post(θ) = ln_prior(θ, system) + ln_like(θ, system)
 
    # ln_prior_system_specialized = make_ln_prior(θ, system)
    # ln_post(θ) = ln_prior_system_specialized(θ, system) + ln_like(θ, system)

    @info "Finding starting point"
    # initial_walkers = find_starting_walkers(ln_post, priors, numwalkers)
    initial_walkers = [sample_priors(system) for _ in 1:numwalkers]

    # Convert the initial walkers into static arrays for stack allocation.
    # This messy line should have no impact on the semantics of the code.
    initial_walkers_static = [
        ComponentVector{SVector{length(cv)}}(; NamedTuple(cv)...) for cv in initial_walkers
    ]

    # Run the MCMC
    thetase, _accept_ratioe = KissMCMC.emcee(
        ln_post,
        initial_walkers_static;
        nburnin = burnin * numwalkers,
        use_progress_meter = true,
        nthin = thinning,
        niter = numsamples_perwalker * numwalkers,
    )

    # Convert the output into an MCMCChains.Chain.
    # Use reinterpret to avoid re-allocating all that memory
    if squash
        thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
        reinterptted = reinterpret(reshape, eltype(first(thetase′)), thetase′)
        chains = ComponentArray(collect(eachrow(reinterptted)), getaxes(thetase′[1]))
    else
        matrix_paramxstep_per_walker = [reinterpret(reshape, eltype(first(θ)), θ) for θ in thetase]
        A = reshape(
            mapreduce(θ->reinterpret(reshape, eltype(first(θ)), θ), hcat, thetase),
            (length(thetase[1][1]), :, numwalkers,)
        )
        ax = getaxes(thetase[1][1])
        chains = ComponentArray(collect(eachslice(A,dims=1)), ax)
    end

    # return Chains(reinterptted, column_names)
    return chains
end

"""
Given a component vector and a matching component vector that is a boolean mask,
return only the values that are true
"""
function select_cv(cv, mask)
    dat = getdata(cv)
    ax = getaxes(cv)
    dat, ax
    # ComponentArray(dat[mask],ax[1][mask])
    ax[1][mask]
end


# using Zygote
# https://github.com/FluxML/Zygote.jl/issues/570
# @Zygote.adjoint (T::Type{<:SArray})(x::Number...) = T(x...), y->(nothing, y...)
function hmc(
    system::System, target_accept=0.8;
    numwalkers=1,
    burnin,
    numsamples_perwalker,
    initial_samples=100_000,
    initial_parameters=nothing
)

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(system)
    fixed = priors_fixed(system)
    D = length(initial_θ_0)

    # AdvancedHMC doesn't play well with component arrays by default, so we pass in just the underlying data
    # array and reconstruct the component array on each invocation (this get's compiled out, no perf affects)
    ax = getaxes(initial_θ_0)
    # Capture the axis into the closure for performance via the let binding.
    ℓπ = let ax=ax, system=system, initial_θ_0=initial_θ_0, fixed=fixed, notfixed = .! fixed
        function (θ)
            θ_cv = ComponentArray(θ, ax)

            # Correct fixed parameters
            θ_cv_merged = θ_cv .* notfixed .+ initial_θ_0 .* fixed

            # TODO: verify that this is still a static array
            ll = ln_post(θ_cv_merged, system)

            return ll
        end
    end

    # ℓπ_grad = let ax=ax, system=system, initial_θ_0=initial_θ_0, fixed=fixed, notfixed = .! fixed
    #     f(θ) = ln_post(θ, system)
    #     function (θ)
    #         θ_cv = ComponentArray(θ, ax)

    #         # Correct fixed parameters
    #         θ_cv_merged = θ_cv .* notfixed .+ initial_θ_0 .* fixed

    #         # TODO: verify that this is still a static array
    #         ll = ln_post(θ_cv_merged, system)

    #         ll_grad = FiniteDiff.finite_difference_gradient(f,θ_cv_merged)

    #         return ll, getdata(ll_grad)
    #     end
    # end

    chains = []
    stats = []
    # Threads.@threads
     for _ in 1:numwalkers
        # initial_θ = sample_priors(system)
        # initial_θ_cv = mean_priors(system)
        # initial_θ_cv = guess_starting_position(system,100_000)

        if isnothing(initial_parameters)
            initial_θ_cv = guess_starting_position(system,initial_samples)
        else
            initial_θ_cv = initial_parameters
        end

        # Use a static comopnent array for efficiency
        # initial_θ = ComponentVector{SVector{length(initial_θ_cv)}}(; NamedTuple(initial_θ_cv)...)
        initial_θ = initial_θ_cv

        # Define a Hamiltonian system
        metric = DenseEuclideanMetric(D)
        hamiltonian = Hamiltonian(metric, ℓπ, ForwardDiff)
        # hamiltonian = Hamiltonian(metric, ℓπ, ℓπ_grad)


        # Define a leapfrog solver, with initial step size chosen heuristically
        # if !isnothing(system.images)
        #     initial_ϵ = 0.002
        # else
            initial_ϵ = find_good_stepsize(hamiltonian, getdata(initial_θ))
            # initial_ϵ = 0.002
        # end


        integrator = Leapfrog(initial_ϵ)
        # 1.05 improves the sampling over standard leapfrog, but 4.0 is too much. It just stays stuck.
        # 1.5 seems better but seems to favour certain trajectories.
        # integrator = TemperedLeapfrog(initial_ϵ, 1.05)
        proposal = NUTS(integrator, max_depth=12) 


        # # We have change some parameters when running with image data
        # if !isnothing(system.images) && target_accept > 0.4
        #     target_accept = 0.2
        #     @info "Sampling from images, lowering target_accept to 0.2"
        # end

        adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(target_accept, integrator)) 
        # adaptor = MassMatrixAdaptor(metric)

        logger = SimpleLogger(stdout, Logging.Error)
        samples, stat = with_logger(logger) do
            sample(hamiltonian, proposal, getdata(initial_θ), numsamples_perwalker, adaptor, burnin; progress=(numwalkers==1), drop_warmup=!(adaptor isa AdvancedHMC.NoAdaptation))
        end

        sample_grid = reduce(hcat, samples);
        chain = ComponentArray(collect(eachrow(sample_grid)), ax)

        # notfixed = .! fixed
        # chain_merged = chain .* notfixed .+ initial_θ_0' .* fixed
        
        push!(chains,chain)
        push!(stats,stat)
    end
    return chains, stats
end





# using Optim, ForwardDiff

function find_starting_point(ln_post, priors)
    θ₀ = rand.(priors)
    i = 0
    while !isfinite(ln_post(θ₀))
        i += 1
        θ₀ = rand.(priors)
        if i > 1000
            error(
                "Could not find a starting point in the posterior that is finite by drawing from the priors after 1000 attempts",
            )
        end
    end
    return θ₀

    # goal(θ) = -ln_post(θ)

    # m = optimize(goal, θ₀, BFGS(), Optim.Options(show_trace=true,x_tol=-1,g_tol=1e-1,f_tol=-1); autodiff = :forward)
    # # m = optimize(goal, θ₀, Newton(), Optim.Options(show_trace=true,); autodiff = :forward)
    # # m = optimize(goal, θ₀, NelderMead(), Optim.Options(show_trace=true,); autodiff = :forward)
    # # m = optimize(goal, θ₀, Optim.SimulatedAnnealing(), Optim.Options(show_trace=true,iterations=1_000_000); autodiff = :forward)
    # # m = optimize(goal, θ₀, Optim.ParticleSwarm(), Optim.Options(show_trace=true,iterations=100_000); autodiff = :forward)

    # display(m)
    # return Optim.minimizer(m)

end



# Start walkers in a gaussian ball around the MLE, while ensuring we don't
# step outside the ranges defined by the priors
function find_starting_walkers(ln_post, priors, numwalkers)
    initial_walkers = map(1:numwalkers) do i
        return initial_position = find_starting_point(ln_post, priors)
    end
    #     initial_position
    # end
    #     # # This used to intiialize the walkers in a Gaussian ball around the MAP.
    #     # # But now we just draw starting points randomly from the priors, this isn't needed.
    #     # map(eachindex(initial_position)) do i
    #     #     p = NaN
    #     #     while !(minimum(priors[i]) < p < maximum(priors[i]))
    #     #         p = initial_position[i] + 0.01randn()*initial_position[i]
    #     #     end
    #     #     p
    #     # end
    #     # initial_position .* (1 .+ randn(length(initial_position)))
    # end
    return initial_walkers
end
