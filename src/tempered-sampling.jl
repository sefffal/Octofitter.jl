using MCMCTempering



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
    sampler,
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



function MCMCTempering.compute_tempered_logdensities(
    model::DifferentiableDensityModel,
    transition,
    transition_other,
    β
)
    lp = β * model.ℓπ.ℓlikelihood(transition) * β
    # Compute for the other.
    lp_other = β * model.ℓπ.ℓlikelihood(transition_other) * β
    return lp, lp_other
end


using MCMCTempering
# Fallback when no random number generator is provided (as is usually the case)
function hmctempered(system::System, target_accept::Number=0.8, ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial(); kwargs...)
    return hmctempered(Random.default_rng(), system, target_accept, ensemble; kwargs...)
end
function hmctempered(
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
    arr2nt = Octofitter.make_arr2nt(system) 
    ln_like_generated = make_ln_like(system, arr2nt(initial_θ_0))

    priors_vec = _list_priors(system)
    Bijector_invlinkvec = make_Bijector_invlinkvec(priors_vec)
    initial_θ_0_t = Bijectors.link.(priors_vec, initial_θ_0)
    arr2nt = Octofitter.make_arr2nt(system)

    # Test out model likelihood and prior computations. This way, if they throw
    # an error, we'll see it right away instead of burried in some deep stack
    # trace from the sampler, autodiff, etc.
    ln_like_generated(system, arr2nt(initial_θ_0))
    ln_prior_transformed(initial_θ_0_t)


    verbosity >= 1 && @info "Preparing model"
    # Capture these variables in a let binding to improve performance
    # We also set up temporary storage to reduce allocations
    # ForwardDiff is used to compute the likelihood gradients using the in-place 
    # API. This ensures type stability.
    ℓπ,∇ℓπ,ℓπ_prior,ℓπ_like = let system=system,
                 ln_prior_transformed=ln_prior_transformed,
                 arr2nt=arr2nt,
                 Bijector_invlinkvec=Bijector_invlinkvec,
                 initial_θ_0_t=initial_θ_0_t,
                 ln_like=ln_like_generated

        function ℓπcallback(θ_transformed, system=system)
            # Transform back from the unconstrained support to constrained support for the likelihood function
            θ_natural = Bijector_invlinkvec(θ_transformed)
            θ_structured = arr2nt(θ_natural)
            lp = ln_prior_transformed(θ_transformed)
            ll = ln_like(system, θ_structured)
            return lp+ll
        end

        function ℓπ_prior(θ_transformed, system=system)
            return ln_prior_transformed(θ_transformed)
        end
        function ℓπ_like(θ_transformed, system=system)
            # Transform back from the unconstrained support to constrained support for the likelihood function
            θ_natural = Bijector_invlinkvec(θ_transformed)
            θ_structured = arr2nt(θ_natural)
            return ln_like(system, θ_structured)
        end

        if autodiff == :Enzyme
            # Enzyme mode:
            ∇ℓπcallback = let diffresult = copy(initial_θ_0_t), system=system, ℓπcallback=ℓπcallback
                system_tmp = deepcopy(system)
                function (θ_t)
                    likelihood = ℓπcallback(θ_t)
                    fill!(diffresult,0)
                    out= Main.Enzyme.autodiff(
                        Main.Enzyme.Reverse,
                        ℓπcallback,
                        Main.Enzyme.Active,
                        Main.Enzyme.Duplicated(θ_t,diffresult),
                        Main.Enzyme.DuplicatedNoNeed(system, system_tmp)
                    )
                    return likelihood, diffresult
                end
            end


            ## Main.Enzyme.autodiff(ℓπcallback, Main.Enzyme.Duplicated(θ_t,diffresult))
            ## Main.Enzyme.autodiff(Main.Enzyme.Reverse, ℓπcallback, Main.Enzyme.Duplicated(θ_t,diffresult))
            ## Main.Enzyme.autodiff(Main.Enzyme.Forward, ℓπcallback, Main.Enzyme.Duplicated, Main.Enzyme.Duplicated(θ_t,diffresult))


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
            ∇ℓπcallback = let cfg=cfg, diffresult=diffresult, ℓπcallback=ℓπcallback
                function (θ_transformed)
                    result = ForwardDiff.gradient!(diffresult, ℓπcallback, θ_transformed, cfg)
                    return DiffResults.value(result), DiffResults.gradient(result)
                end
            end
        else
            error("Unsupported option for autodiff: $autodiff. Valid options are :ForwardDiff (default), :Enzyme, and :Zygote.")
        end
        
        ℓπcallback, ∇ℓπcallback, ℓπ_prior, ℓπ_like
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

    verbosity >= 1 && @info "Determining initial positions and metric using pathfinder"
    # Use Pathfinder to initialize HMC.
    # It seems to hit a PosDefException sometimes when factoring a matrix.
    # When that happens, the next try usually succeeds.
    start_time = time()
    local result_pf = nothing
    for retry in 1:5
        try
            result_pf = pathfinder(
                ℓπ;
                # ad_backend=AD.FiniteDifferencesBackend(),
                ad_backend=AD.ForwardDiffBackend(),
                init=collect(initial_θ_t),
                progress=true,
                # maxiters=5,
                # maxtime=5.0,
                # reltol=1e-4,
            ) 
            break
        catch ex
            if ex isa LinearAlgebra.PosDefException
                @warn "pathfinder hit a PosDefException. Retrying" exception=ex retry
                continue
            elseif ex isa InterruptException
                rethrow(ex)
            else
                @error "Unexpected error occured running pathfinder" exception=(ex, catch_backtrace())
                rethrow(ex)
            end
        end
    end
    if isnothing(result_pf)
        error("Warm up failed: pathfinder failed 5 times")
    end
    stop_time = time()

    verbosity >= 3 && @info(
        "Pathfinder results",
        mode=arr2nt(Bijectors.invlink.(priors_vec, result_pf.fit_distribution.μ)),
        inv_metric=Matrix(result_pf.fit_distribution.Σ)
    )

    # # Return a chains object with the resampled pathfinder draws
    # # Transform samples back to constrained support
    # pathfinder_samples = map(eachcol(result_pf.draws)) do θ_t
    #     Bijectors.invlink.(priors_vec, θ_t)
    # end
    # pathfinder_chain =  Octofitter.result2mcmcchain(system, arr2nt.(pathfinder_samples))
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
    initial_θ_t = result_pf.draws[:, end]

    verbosity >= 3 && @info "Creating metric"

    # # Use the metric found by Pathfinder for HMC sampling
    # metric = Pathfinder.RankUpdateEuclideanMetric(result_pf.fit_distribution.Σ)
    
    # Start with found pathfinder metric then adapt a dense metric:
    # metric = DenseEuclideanMetric(Matrix(result_pf.fit_distribution.Σ))
    metric = DiagEuclideanMetric(diag(Matrix(result_pf.fit_distribution.Σ)))
    # metric = DenseEuclideanMetric(collect(Diagonal(Matrix(result_pf.fit_distribution.Σ))))

    # Fit a dense metric from scratch
    # metric = DenseEuclideanMetric(D)
    # metric = DiagEuclideanMetric(D)

    verbosity >= 3 && @info "Creating model" 
    # model = AdvancedHMC.DifferentiableDensityModel(ℓπ, ∇ℓπ)
    model = AdvancedHMC.DifferentiableDensityModel(
        TemperedJoint(ℓπ_prior, ℓπ_like),
        ∇ℓπ
    )

    verbosity >= 3 && @info "Creating hamiltonian"
    hamiltonian = Hamiltonian(metric, ℓπ, ∇ℓπ)
    verbosity >= 3 && @info "Finding good stepsize"
    ϵ = find_good_stepsize(hamiltonian, initial_θ_t)
    verbosity >= 3 && @info "Found initial stepsize" ϵ 
    # integrator = Leapfrog(ϵ)
    integrator = JitteredLeapfrog(ϵ, 0.1) # 10% normal distribution on step size to help in areas of high curvature. 
    verbosity >= 3 && @info "Creating kernel"
    κ = NUTS{MultinomialTS,GeneralisedNoUTurn}(integrator, max_depth=tree_depth)
    # κ = NUTS{SliceTS,GeneralisedNoUTurn}(integrator, max_depth=tree_depth)
    
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


    
    inverse_temperatures = MCMCTempering.check_inverse_temperatures([0.25, 0.5, 0.75, 1.0])
    swapstrategy = MCMCTempering.NonReversibleSwap()
    swap_every = 2
    tempered_sampler = tempered(
        sampler, inverse_temperatures, swapstrategy;
        adapt=false,
        adapt_schedule=MCMCTempering.Geometric(),
        adapt_stepsize=1,
        adapt_eta=0.66,
        adapt_target=0.234,
        swap_every=swap_every
    )


    mc_samples_all_chains = sample(
        rng,
        model,
        tempered_sampler,
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
        push!(chains, Octofitter.result2mcmcchain(system, chain_res))
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
