
export sample_priors
sample_priors(planet::Planet) = rand.(planet.priors.priors)
sample_priors(planet::Planet,N) = rand.(planet.priors.priors,N)

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

function sample_priors(system::System)
    sampled = ComponentVector(
        merge(NamedTuple(rand.(system.priors.priors)),
        (;planets=[sample_priors(planet) for planet in system.planets])
    ))
    return sampled
end
function sample_priors(system::System,N)
    sampled = ComponentVector(
        merge(NamedTuple(rand.(system.priors.priors,N)),
        (;planets=[sample_priors(planet,N) for planet in system.planets])
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
    N = 5
    sampled = ComponentVector(
        merge(NamedTuple(mean.(rand.(system.priors.priors,N))),
        # (;planets=[mean.(sample_priors(planet,N)) for planet in system.planets])
        (;planets=[
            ComponentArray(NamedTuple([k=>mean(v) for (k,v) in pairs(NamedTuple(sample_priors(planet,5)))]))
            for planet in system.planets
        ])
    ))
    return sampled
end



function construct_elements(θ_system, θ_planet)
    # Handle re-parameterized models
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
        Ω = θ_planet.Ω
        τ = θ_planet.τ
    end
    return KeplerianElements((;
        θ_system.μ,
        θ_system.plx,
        θ_planet.i,
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
        Ω = θ_planet.Ω[i]
        τ = θ_planet.τ[i]
    end
    return KeplerianElements((;
        μ=θ_system.μ[i],
        plx=θ_system.plx[i],
        i=θ_planet.i[i],
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
        # We can reinterpret the vector of SVectors as a matrix directly without copying!
        # This can save massive amounts of memory and time on large changes
        reinterptted =
            cat([reinterpret(reshape, eltype(first(θ)), θ) for θ in thetase]..., dims = 3)
        chains = ComponentArray(
            collect(eachslice(reinterptted, dims = 1)),
            getaxes(thetase[1][1]),
        )
    end

    # return Chains(reinterptted, column_names)
    return chains
end


function hmc(
    system::System;
    numwalkers=1,
    burnin,
    numsamples_perwalker,
)

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(system)
    D = length(initial_θ_0)

    # AdvancedHMC doesn't play well with component arrays by default, so we pass in just the underlying data
    # array and reconstruct the component array on each invocation (this get's compiled out, no perf affects)
    ax = getaxes(initial_θ_0)
    # Capture the axis into the closure for performance via the let binding.
    ℓπ = let ax=ax, system=system
        function (θ)
            θ_cv = ComponentArray(θ, ax)
            return ln_post(θ_cv, system)
            # return ln_prior(θ_cv, system)
            # return ln_like(θ_cv, system)
        end
    end

    chains = []
    stats = []
    # Threads.@threads
     for _ in 1:numwalkers
        # initial_θ = sample_priors(system)
        initial_θ = mean_priors(system)

        # Define a Hamiltonian system
        metric = DiagEuclideanMetric(D)
        # metric = DenseEuclideanMetric(D)
        hamiltonian = Hamiltonian(metric, ℓπ, ForwardDiff)

        # Define a leapfrog solver, with initial step size chosen heuristically
        initial_ϵ = find_good_stepsize(hamiltonian, getdata(initial_θ))
        # initial_ϵ = 0.02
        integrator = Leapfrog(initial_ϵ)
        # integrator = TemperedLeapfrog(initial_ϵ, 31.0)


        proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator, 10)

        adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator)) #0.08
        # adaptor = AdvancedHMC.NoAdaptation()


        logger = SimpleLogger(stdout, Logging.Error)
        samples, stat = with_logger(logger) do
            sample(hamiltonian, proposal, getdata(initial_θ), numsamples_perwalker, adaptor, burnin; progress=(numwalkers==1), drop_warmup=!(adaptor isa AdvancedHMC.NoAdaptation))
        end

        sample_grid = reduce(hcat, samples);
        chain = ComponentArray(collect(eachrow(sample_grid)), ax)
        
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
