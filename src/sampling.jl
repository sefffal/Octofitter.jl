using DiffResults
using AbstractDifferentiation
AD = AbstractDifferentiation
using LinearAlgebra

export sample_priors
sample_priors(planet::Planet) = rand.(ComponentArray(planet.priors.priors))
sample_priors(planet::Planet,N) = [sample_priors(planet) for _ in 1:N]


# function sample_priors(system::System)
#     sampled = ComponentVector(
#         merge(NamedTuple(rand.(ComponentArray(system.priors.priors))),
#         # (;planets=[sample_priors(planet) for planet in system.planets])
#         (;planets=namedtuple(collect(keys(system.planets)), [
#             ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(sample_priors(planet)))]))
#             for planet in system.planets
#         ]))
#     ))
#     return getdata(sampled)
# end


function sample_priors(system::System)
    priors_flat = map(((k,v),)->rand(v), Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
    # return rand.(priors_flat)
end


sample_priors(system::System,N) = [sample_priors(system) for _ in 1:N]





# # Instead of just calling mean for the distributions, we sample and then take the mean of the data.
# # This does add a little jitter, but some distributions do not directly define the mean function!
# # Specifically, truncated(InverseGamma()) does not work, and this is very useful.
# # mean_priors(planet::Planet) = Statistics.mean.(Statistics.rand.(planet.priors.priors,1000))
# function mean_priors(system::System)
#     priors_all = ComponentVector(;
#         NamedTuple(system.priors.priors)...,
#         planets=[planet.priors.priors for planet in system.planets]
#     )
#     # return Statistics.mean.(Statistics.rand.(priors_all,1000))
#     return Statistics.mean.(priors_all)
# end


function guess_starting_position(system, N=500_000)

    # TODO: this shouldn't have to allocate anything, we can just loop keeping the best.
    θ = sample_priors(system, N)
    arr2nt = DirectDetections.make_arr2nt(system) 

    posts = zeros(N)
    ln_prior = make_ln_prior(system)
    Threads.@threads for i in eachindex(posts)
        θ_res = arr2nt(θ[i])
        posts[i] = ln_prior(θ[i]) + ln_like(system, θ_res)
    end
    mapv,mapi = findmax(posts)
    best = θ[mapi]
    

    return best, mapv
end

# This code works but I have not found it useful. Commenting out to remove Optim.jl dependency.
# using Optim
# function optimize_starting_position(ℓπ, initial_θ_t)

#     @info "Optimizing starting location "
#     results = optimize(θ_t-> -ℓπ(θ_t), initial_θ_t, LBFGS(), autodiff = :forward, Optim.Options(iterations=100_000, time_limit=5, allow_f_increases=true))
#     # results = optimize(θ_t-> -ℓπ(θ_t), initial_θ_t, NelderMead(), Optim.Options(iterations=100_000, time_limit=5))
#     # results = optimize(θ_t-> -ℓπ(θ_t), initial_θ_t, ParticleSwarm(), Optim.Options(iterations=100_000, time_limit=5))

#     intial_objective = ℓπ(initial_θ_t)

#     minimizer = Optim.minimizer(results)

#     @info "Starting location improved" logimprovement=(ℓπ(minimizer) - intial_objective)
#     return minimizer
# end


"""
    construct_elements(θ_system, θ_planet)

Given a named tuple for of parameters from a System (θ_system) and Planet (θ_planet),
return a `VisualOrbit PlanetOrbits.jl.
"""
function construct_elements(::Type{VisualOrbit}, θ_system, θ_planet)
    return VisualOrbit((;
        θ_system.M,
        θ_system.plx,
        θ_planet.i,
        θ_planet.Ω,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    ))
end
function construct_elements(::Type{KepOrbit}, θ_system, θ_planet)
    return KepOrbit((;
        θ_system.M,
        θ_planet.i,
        θ_planet.Ω,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    ))
end
function construct_elements(::Type{ThieleInnesOrbit}, θ_system, θ_planet)
    return ThieleInnesOrbit((;
        θ_system.M,
        θ_system.plx,
        θ_planet.A,
        θ_planet.B,
        θ_planet.F,
        θ_planet.G,
        θ_planet.e,
        θ_planet.τ,
    ))
end
function construct_elements(::Type{RadialVelocityOrbit}, θ_system, θ_planet)
    return RadialVelocityOrbit((;
        θ_system.M,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    ))
end


"""
    construct_elements(chains, :b, 4)

Given a Chains object, a symbol matching the name of a planet, and an index,
construct a `VisualOrbit DirectOrbits of that planet from that
index of the chains.
"""
function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, i::Union{Integer,CartesianIndex})
    pk = string(planet_key)
    if haskey(chain, :plx) && haskey(chain, Symbol(pk*".i")) && haskey(chain, Symbol(pk*".Ω"))
        return VisualOrbit((;
            M=chain["M"][i],
            plx=chain["plx"][i],
            i=chain[pk*".i"][i],
            Ω=chain[pk*".Ω"][i],
            ω=chain[pk*".ω"][i],
            e=chain[pk*".e"][i],
            τ=chain[pk*".τ"][i],
            a=chain[pk*".a"][i],
        ))
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*".A")) && haskey(chain, Symbol(pk*".B")) && haskey(chain, Symbol(pk*".G"))&& haskey(chain, Symbol(pk*".F"))
        return ThieleInnesOrbit((;
            M=chain["M"][i],
            plx=chain["plx"][i],
            e=chain[pk*".e"][i],
            τ=chain[pk*".τ"][i],
            A=chain[pk*".A"][i],
            B=chain[pk*".B"][i],
            F=chain[pk*".F"][i],
            G=chain[pk*".G"][i],
        ))
    elseif haskey(chain, Symbol(pk*".i")) && haskey(chain, Symbol(pk*".Ω"))
        return KepOrbit((;
            M=chain["M"][i],
            i=chain[pk*".i"][i],
            Ω=chain[pk*".Ω"][i],
            ω=chain[pk*".ω"][i],
            e=chain[pk*".e"][i],
            τ=chain[pk*".τ"][i],
            a=chain[pk*".a"][i],
        ))
    elseif haskey(chain, :M) && haskey(chain, :rv)
        return RadialVelocityOrbit((;
            M=chain["M"][i],
            ω=chain[pk*".ω"][i],
            e=chain[pk*".e"][i],
            τ=chain[pk*".τ"][i],
            a=chain[pk*".a"][i],
        ))
    else
        error("Unrecognized columns")
    end
end

"""
    construct_elements(chains, :b, [4,5,10])

Given a Chains object, a symbol matching the name of a planet, and an array of indices,
construct a `VisualOrbit DirectOrbits of that planet from those indices
of the chains.
"""
function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, ii::AbstractArray{<:Union{Integer,CartesianIndex}})
    pk = string(planet_key)
    if haskey(chain, :plx) && haskey(chain, Symbol(pk*".i")) && haskey(chain, Symbol(pk*".Ω"))
        Ms=chain["M"]
        plxs=chain["plx"]
        is=chain[pk*".i"]
        Ωs=chain[pk*".Ω"]
        ωs=chain[pk*".ω"]
        es=chain[pk*".e"]
        τs=chain[pk*".τ"]
        as=chain[pk*".a"]
        return map(ii) do i
            VisualOrbit((;
                M=Ms[i],
                plx=plxs[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                τ=τs[i],
                a=as[i],
            ))
        end
    elseif haskey(chain, Symbol(pk*".i")) && haskey(chain, Symbol(pk*".Ω"))
        Ms=chain["M"]
        is=chain[pk*".i"]
        Ωs=chain[pk*".Ω"]
        ωs=chain[pk*".ω"]
        es=chain[pk*".e"]
        τs=chain[pk*".τ"]
        as=chain[pk*".a"]
        return map(ii) do i
            KepOrbit((;
                M=Ms[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                τ=τs[i],
                a=as[i],
            ))
        end
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*".A")) && haskey(chain, Symbol(pk*".B"))
        Ms=chain["M"]
        plxs=chain["plx"]
        As=chain[pk*".A"]
        Bs=chain[pk*".B"]
        Fs=chain[pk*".F"]
        Gs=chain[pk*".G"]
        es=chain[pk*".e"]
        τs=chain[pk*".τ"]
        return map(ii) do i
            ThieleInnesOrbit((;
                M=Ms[i],
                plx=plxs[i],
                e=es[i],
                τ=τs[i],
                A=As[i],
                B=Bs[i],
                F=Fs[i],
                G=Gs[i],
            ))
        end
    elseif haskey(chain, Symbol("M")) && haskey(chain, Symbol("rv"))
        Ms=chain["M"]
        ωs=chain[pk*".ω"]
        es=chain[pk*".e"]
        τs=chain[pk*".τ"]
        as=chain[pk*".a"]
        return map(ii) do i
            RadialVelocityOrbit((;
                M=Ms[i],
                ω=ωs[i],
                e=es[i],
                τ=τs[i],
                a=as[i],
            ))
        end
    else
        error("Unrecognized chain format")
    end
end
construct_elements(chain::Chains, planet_key::Union{String,Symbol}, ii::Colon) = construct_elements(chain, planet_key, 1:size(chain,1))
function construct_elements(chain, planet_key::Union{String,Symbol}, ii::AbstractArray{<:Union{Integer,CartesianIndex}})
    pk = string(planet_key)
    Ms=chain[:,"M"]
    plxs=chain[:,"plx"]
    is=chain[:,pk*".i"]
    Ωs=chain[:,pk*".Ω"]
    ωs=chain[:,pk*".ω"]
    es=chain[:,pk*".e"]
    τs=chain[:,pk*".τ"]
    as=chain[:,pk*".a"]
    return map(ii) do i
        VisualOrbit((;
            M=Ms[i],
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
construct_elements(chain::Chains, planet::Planet, args...; kwargs...) = construct_elements(chain, planet.name, args...; kwargs...) 

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
    thinning=1,
    discard_initial=adaptation,
    tree_depth=10,
    initial_samples=50_000,
    initial_parameters=nothing,
    verbosity=2,
    autodiff=:ForwardDiff,
    debug=false
)

    if ensemble != MCMCSerial()
        @warn "TODO: In-place model gradients currently not supported with MCMCThreads"
    end

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(system)
    D = length(initial_θ_0)

    ln_prior_transformed = make_ln_prior_transformed(system)
    # ln_prior = make_ln_prior(system)
    arr2nt = DirectDetections.make_arr2nt(system) 

    priors_vec = _list_priors(system)
    Bijector_invlinkvec = make_Bijector_invlinkvec(priors_vec)
    initial_θ_0_t = Bijectors.link.(priors_vec, initial_θ_0)
    arr2nt = DirectDetections.make_arr2nt(system)

    # Test out model likelihood and prior computations. This way, if they throw
    # an error, we'll see it right away instead of burried in some deep stack
    # trace from the sampler, autodiff, etc.
    ln_like(system, arr2nt(initial_θ_0))
    ln_prior_transformed(initial_θ_0_t)


    verbosity >= 1 && @info "Preparing model"
    # Capture these variables in a let binding to improve performance
    # We also set up temporary storage to reduce allocations
    # ForwardDiff is used to compute the likelihood gradients using the in-place 
    # API. This ensures type stability.
    ℓπ,∇ℓπ = let system=system,
                 ln_prior_transformed=ln_prior_transformed,
                 arr2nt=arr2nt,
                 Bijector_invlinkvec=Bijector_invlinkvec,
                 initial_θ_0_t=initial_θ_0_t
        # function ℓπcallback(θ_transformed)
        #     # Transform back from the unconstrained support to constrained support for the likelihood function
        #     θ_natural = Bijector_invlinkvec(θ_transformed)
        #     θ_structured = arr2nt(θ_natural)
        #     lp = ln_prior_transformed(θ_transformed)
        #     ll = ln_like(system, θ_structured)#*θ_structured.γ # WIP: ad-hoc tempering scheme
        #     @show ll
        #     return lp+ll
        # end
        function ℓπcallback(θ_transformed)
            # Transform back from the unconstrained support to constrained support for the likelihood function
            θ_natural = Bijector_invlinkvec(θ_transformed)
            θ_structured = arr2nt(θ_natural)
            lp = ln_prior_transformed(θ_transformed)
            ll = ln_like(system, θ_structured)#*θ_structured.γ # WIP: ad-hoc tempering scheme
            return lp+ll
        end

        if autodiff == :Enzyme
            # Enzyme mode:
            ∇ℓπcallback = let diffresult = copy(initial_θ_0_t)
                function (θ_t)
                    primal = ℓπcallback(θ_t)
                    fill!(diffresult,0)
                    # Main.Enzyme.autodiff(ℓπcallback, Main.Enzyme.Duplicated(θ_t,diffresult))
                    Main.Enzyme.autodiff(Main.Enzyme.Reverse, ℓπcallback, Main.Enzyme.Active, Main.Enzyme.Duplicated(θ_t,diffresult))
                    return primal, diffresult
                end
            end

        elseif autodiff == :FiniteDiff
                ∇ℓπcallback = let diffresult = copy(initial_θ_0_t)
                    function (θ_t)
                        primal = ℓπcallback(θ_t)
                        Main.FiniteDiff.finite_difference_gradient!(diffresult, ℓπcallback, θ_t)
                        return primal, diffresult
                    end
                end
        
        elseif autodiff == :Zygote

            # Zygote mode:
            ∇ℓπcallback = function (θ_t)
                Main.Zygote.gradient(ℓπ, θ_t)
            end

        elseif autodiff == :ForwardDiff

            # ForwardDiff mode:
            # Create temporary storage space for gradient computations
            diffresult = DiffResults.GradientResult(initial_θ_0_t)

            # Perform dynamic benchmarking to pick a ForwardDiff chunk size.
            # We're going to call this thousands of times so worth a few calls
            # to get this optimized.
            chunk_sizes = unique([1; 2:2:D; D])
            ideal_chunk_size_i = argmin(map(chunk_sizes) do chunk_size
                cfg = ForwardDiff.GradientConfig(ℓπcallback, initial_θ_0_t, ForwardDiff.Chunk{chunk_size}());
                ForwardDiff.gradient!(diffresult, ℓπcallback, initial_θ_0_t, cfg)
                t = minimum(
                    @elapsed ForwardDiff.gradient!(diffresult, ℓπcallback, initial_θ_0_t, cfg)
                    for _ in 1:10
                )

                verbosity >= 3 && @info "Timing autodiff" chunk_size t
                return t
            end)
            ideal_chunk_size =  chunk_sizes[ideal_chunk_size_i]
            verbosity >= 1 && @info "Selected auto-diff chunk size" ideal_chunk_size

            cfg = ForwardDiff.GradientConfig(ℓπcallback, initial_θ_0_t, ForwardDiff.Chunk{ideal_chunk_size}());
            ∇ℓπcallback = let cfg=cfg, diffresult=diffresult
                function (θ_transformed)
                    result = ForwardDiff.gradient!(diffresult, ℓπcallback, θ_transformed, cfg)
                    return DiffResults.value(result), DiffResults.gradient(result)
                end
            end
        else
            error("Unsupported option for autodiff: $autodiff. Valid options are :ForwardDiff (default), :Enzyme, and :Zygote.")
        end
        
        ℓπcallback, ∇ℓπcallback
    end

        
    # Display their run time. If something is egregously wrong we'll notice something
    # here in the output logs.
    ℓπ(initial_θ_0_t)
    ∇ℓπ(initial_θ_0_t) 
    verbosity >= 1 && @showtime ℓπ(initial_θ_0_t)
    verbosity >= 1 && @showtime ∇ℓπ(initial_θ_0_t)

    if debug
        return ℓπ, ∇ℓπ, initial_θ_0_t
    end

    # # Uncomment for debugging model type-stability
    # Main.InteractiveUtils.@code_warntype ℓπ(initial_θ_0_t)
    # Main.InteractiveUtils.@code_warntype ∇ℓπ(initial_θ_0_t)
    # return

    # Guess initial starting positions by drawing from priors a bunch of times
    # and picking the best one (highest likelihood).
    # Then transform that guess into our unconstrained support
    if isnothing(initial_parameters)
        verbosity >= 1 && @info "Guessing a starting location by sampling from prior" initial_samples
        initial_θ, mapv = guess_starting_position(system,initial_samples)
        verbosity > 2 && @info "Found starting location" ℓπ(θ)=mapv θ=arr2nt(initial_θ)
        # Transform from constrained support to unconstrained support
        initial_θ_t = Bijectors.link.(priors_vec, initial_θ)
    else
        initial_θ = initial_parameters
        # Transform from constrained support to unconstrained support
        initial_θ_t = Bijectors.link.(priors_vec, initial_θ)
    end

    # verbosity >= 1 && @info "Determining initial positions and metric using pathfinder"
    # # Use Pathfinder to initialize HMC.
    # # It seems to hit a PosDefException sometimes when factoring a matrix.
    # # When that happens, the next try usually succeeds.
    # start_time = time()
    # local result_pf = nothing
    # for retry in 1:5
    #     try
    #         result_pf = pathfinder(
    #             ℓπ;
    #             # ad_backend=AD.FiniteDifferencesBackend(),
    #             ad_backend=AD.ForwardDiffBackend(),
    #             init=collect(initial_θ_t),
    #             progress=true,
    #             # maxiters=5,
    #             # maxtime=5.0,
    #             # reltol=1e-4,
    #         ) 
    #         break
    #     catch ex
    #         if ex isa LinearAlgebra.PosDefException
    #             @warn "pathfinder hit a PosDefException. Retrying" exception=ex retry
    #             continue
    #         elseif ex isa InterruptException
    #             rethrow(ex)
    #         else
    #             @error "Unexpected error occured running pathfinder" exception=(ex, catch_backtrace())
    #             rethrow(ex)
    #         end
    #     end
    # end
    # if isnothing(result_pf)
    #     error("Warm up failed: pathfinder failed 5 times")
    # end
    # stop_time = time()

    # @info(
    #     "Pathfinder results",
    #     mode=arr2nt(Bijectors.invlink.(priors_vec, result_pf.fit_distribution.μ)),
    #     inv_metric=Matrix(result_pf.fit_distribution.Σ)
    # )

    # # Return a chains object with the resampled pathfinder draws
    # # Transform samples back to constrained support
    # pathfinder_samples = map(eachcol(result_pf.draws)) do θ_t
    #     Bijectors.invlink.(priors_vec, θ_t)
    # end
    # pathfinder_chain =  DirectDetections.result2mcmcchain(system, arr2nt.(pathfinder_samples))
    # pathfinder_chain_with_info = MCMCChains.setinfo(
    #     pathfinder_chain,
    #     (;
    #         start_time,
    #         stop_time,
    #         model=system,
    #         result_pf,
    #     )
    # )

    # # Start using a draw from the typical set as estimated by Pathfinder
    # initial_θ_t = result_pf.draws[:, end]

    # # Use the metric found by Pathfinder for HMC sampling
    # metric = Pathfinder.RankUpdateEuclideanMetric(result_pf.fit_distribution.Σ)
    
    # Start with found pathfinder metric then adapt a dense metric:
    # metric = DenseEuclideanMetric(Matrix(result_pf.fit_distribution.Σ))

    # Fit a dense metric from scratch
    verbosity >= 3 && @info "Creating metric"
    metric = DenseEuclideanMetric(D)

    verbosity >= 3 && @info "Creating model" 
    model = AdvancedHMC.DifferentiableDensityModel(ℓπ, ∇ℓπ)

    verbosity >= 3 && @info "Creating hamiltonian"
    hamiltonian = Hamiltonian(metric, ℓπ, ∇ℓπ)
    verbosity >= 3 && @info "Finding good stepsize"
    ϵ = find_good_stepsize(hamiltonian, initial_θ_t)
    verbosity >= 3 && @info "Found initial stepsize" ϵ 
    integrator = Leapfrog(ϵ)
    # integrator = JitteredLeapfrog(ϵ, 0.1) # 10% normal distribution on step size to help in areas of high curvature. 
    verbosity >= 3 && @info "Creating kernel"
    # κ = NUTS{MultinomialTS,GeneralisedNoUTurn}(integrator, max_depth=tree_depth)
    κ = NUTS{SliceTS,GeneralisedNoUTurn}(integrator, max_depth=tree_depth)
    
    verbosity >= 3 && @info "Creating adaptor"
    mma = MassMatrixAdaptor(metric)
    ssa = StepSizeAdaptor(target_accept, integrator)
    adaptor = StanHMCAdaptor(mma, ssa) 
    # adaptor = StepSizeAdaptor(target_accept, integrator)

    # κ = NUTS(integrator, max_depth=tree_depth) 
    verbosity >= 3 && @info "Creating sampler"
    sampler = AdvancedHMC.HMCSampler(κ, metric, adaptor)


    verbosity >= 1 && @info "Adapting sampler..."
    start_time = fill(time(), num_chains)

    # # Neat: it's possible to return a live iterator
    # # We could use this to build e.g. live plotting while the code is running
    # # once the analysis code is ported to Makie.
    # # return  AbstractMCMC.steps(
    # #     rng,
    # #     model,
    # #     sampler,
    # #     nadapts = adaptation,
    # #     init_params = initial_θ_t,
    # #     discard_initial = adaptation,
    # #     progress=progress,
    # #     verbose=false
    # # )


    last_output_time = Ref(time())
    function callback(rng, model, sampler, transition, state, iteration; kwargs...)
        if verbosity >= 1 && iteration == 1
            @info "Adaptation complete."

            # Show adapted step size and mass matrix
            if verbosity >= 3
                adapted_ss = AdvancedHMC.getϵ(adaptor)
                println("Adapated stepsize ϵ=", adapted_ss)
                # adapted_mm = AdvancedHMC.getM⁻¹(adaptor)
                # print("Adapted mass matrix M⁻¹ ")
                # display(adapted_mm)
            end
            
            @info "Sampling..."
            verbosity >= 2 && println("Progress legend: divergence iter(thread) td=tree-depth ℓπ=log-posterior-density ")
        end
        if verbosity < 2 || last_output_time[] + 1 > time()
            return
        end
        # Give different messages if the log-density is non-finite,
        # or if there was a divergent transition.
        if !isfinite(transition.z.ℓπ)
            # TODO: this never runs since any non-finite proposal is rejected during sampling.
            note = "∞" 
        elseif transition.stat.numerical_error
            note = "X"
        else
            note = " "
        end
        if transition.z.ℓπ isa AdvancedHMC.DualValue
            ℓπ = transition.z.ℓπ.value
        else
            ℓπ = transition.z.ℓπ
        end

        θ_message = ""
        if verbosity >= 3
            θ = Bijector_invlinkvec(transition.z.θ)
            θ_res = arr2nt(θ)
            # Fill the remaining width of the terminal with info
            max_width = displaysize(stdout)[2]-34
            θ_str = string(θ_res)
            θ_str_trunc = θ_str[begin:prevind(θ_str, min(end,max_width))]
            θ_message = "θ="*θ_str_trunc*"..."
        end
        
        @printf("%1s%6d(%2d) td=%2d ℓπ=%6.0f. %s\n", note, iteration, Threads.threadid(), transition.stat.tree_depth, ℓπ, θ_message)
    
        # Support for live plotting orbits as we go.
        # This code works but I found it slows down a lot as we go. Maybe it would
        # work better with Makie.
        # for p_i in keys(system.planets)
        #     kep_elements = construct_elements(θ_res, θ_res.planets[p_i])
        #     color = if transition.stat.numerical_error
        #         :red
        #     elseif iteration <= adaptation
        #         :blue
        #     else
        #         :black
        #     end
        #     Main.plot!(kep_elements,color=color, label="")
        # end
        # display(Main.Plots.current())

        last_output_time[] = time()
        return
    end

    # # For N chains, use the last N pathfinder draws as initial parameters
    # if isnothing(initial_parameters)
    #     if num_chains <= size(result_pf.draws,2)
    #         initial_parameters = [
    #             result_pf.draws[:, end-i+1]
    #             for i in 1:num_chains
    #         ]
    #     # If we have more chains to inialize than pathfinder draws, pick from them at random
    #     else
    #         initial_parameters = [
    #             result_pf.draws[:, rand(axes(result_pf.draws,2))]
    #             for _ in 1:num_chains
    #         ]
    #     end
    # else
        initial_parameters = fill(initial_θ_t, num_chains)
    # end

    mc_samples_all_chains = sample(
        rng,
        model,
        sampler,
        ensemble,
        iterations,
        num_chains;
        nadapts = adaptation,
        thinning,
        init_params=initial_parameters,
        discard_initial,
        progress=verbosity >= 1,
        callback,
        verbose=verbosity>=4
    )
    stop_time = fill(time(), num_chains)

    
    verbosity >= 1 && @info "Sampling compete. Building chains."


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
    
        verbosity >= 1 && println("""
        Sampling report for chain $i:
        mean_accept         = $mean_accept
        num_err_frac        = $num_err_frac
        mean_tree_depth     = $mean_tree_depth
        max_tree_depth_frac = $max_tree_depth_frac\
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
            states=mc_samples_all_chains,
            # pathfinder=pathfinder_chain_with_info,
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

# include("tempered-sampling.jl")
# include("zigzag.jl")


"""
Convert a vector of component arrays returned from sampling into an MCMCChains.Chains
object.
"""
function result2mcmcchain(system, chain_in)
    # `system` not currently used, but a more efficient/robust mapping in future might require it.

    # There is a specific column name convention used by MCMCChains to indicate
    # that multiple parameters form a group. Instead of planets.X.a, we adapt our to X[a] 
    # accordingly
    flattened_labels = keys(flatten_named_tuple(first(chain_in)))
    data = zeros(length(chain_in), length(flattened_labels))
    for (i, sample) in enumerate(chain_in)
        for (j, val) in enumerate(Iterators.flatten(Iterators.flatten(sample)))
            data[i,j] = val
        end
    end
    c = Chains(data, [string(l) for l in flattened_labels])
    return c
end



# Used for flattening a nested named tuple posterior sample into a flat named tuple
# suitable to be used as a Tables.jl table row.
function flatten_named_tuple(nt)
    pairs = Pair{Symbol, Float64}[]
    for key in keys(nt)
        if key != :planets
            push!(pairs, key => nt[key])
        end
    end
    for pl in keys(get(nt, :planets, (;)))
        for key in keys(nt.planets[pl])
            push!(pairs, Symbol(pl, '.', key) =>  nt.planets[pl][key])
        end
    end
    return namedtuple(pairs)

end
