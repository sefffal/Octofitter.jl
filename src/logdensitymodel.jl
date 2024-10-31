
using LogDensityProblems

# Define the target distribution using the `LogDensityProblem` interface
# TODO: in future, just roll this all into the System type.
mutable struct LogDensityModel{D,Tℓπ,T∇ℓπ,TSys,TLink,TInvLink,TArr2nt,TPriSamp}
    # Dimensionality
    const D::Int
    # Auto diff backend symbol (need to keep track of this for Pigeons)
    const autodiff_backend_symbol::Symbol
    # Calculate the log-posterior density given transformed parameters
    const ℓπcallback::Tℓπ
    # Calculate the log-posterior density and gradient given transformed parameters
    const ∇ℓπcallback::T∇ℓπ
    # The underlying System object
    const system::TSys
    # Convert flat parameter vector into transformed domain
    const link::TLink
    # Convert flat transformed parameter vector back to natural domain
    const invlink::TInvLink
    # Convert a flat parameter vector into a nested named tuple structure,
    # matching the variable definitions in the System and Planet blocks
    const arr2nt::TArr2nt
    # Sample IID from the model's priors
    const sample_priors::TPriSamp
    # A set of starting points that can be sampled from to initialize a sampler, or nothing
    starting_points::Union{Nothing,Vector} 
    function LogDensityModel(system::System; autodiff=:ForwardDiff, verbosity=0, chunk_sizes=nothing)
        verbosity >= 1 && @info "Preparing model"

        sample_priors = make_prior_sampler(system)

        # Choose parameter dimensionality and initial parameter value
        initial_θ_0 = sample_priors(Random.default_rng())
        D = length(initial_θ_0)
        verbosity >= 2 && @info "Determined number of free variables" D


        ln_prior_transformed = make_ln_prior_transformed(system)
        # ln_prior = make_ln_prior(system)
        arr2nt = Octofitter.make_arr2nt(system) 
        ln_like_generated = make_ln_like(system, arr2nt(initial_θ_0))

        priors_vec = _list_priors(system)
        Bijector_invlinkvec = make_Bijector_invlinkvec(priors_vec)
        # TODO : this could be unrolled like invlink if used anywhere performance sensitive. Currently
        # it just transforms things back after all sampling is done.
        Bijector_linkvec = let priors_vec=priors_vec
            # (θ) -> Bijectors.link.(priors_vec, θ)
            function (θ)
                i = 0
                out = zeros(eltype(θ), length(θ))
                for prior in priors_vec
                    if length(prior) == 1
                        i += 1
                        out[i] = Bijectors.link(prior, θ[i])
                    else
                        param_outs = Bijectors.link(prior, θ[i+1:i+length(prior)])
                        for param_out in param_outs
                            i += 1
                            out[i] = param_out
                        end
                    end
                end
                return out
            end
        end
        initial_θ_0_t = Bijector_linkvec(initial_θ_0)
        arr2nt = Octofitter.make_arr2nt(system)

        # Test out model likelihood and prior computations. This way, if they throw
        # an error, we'll see it right away instead of burried in some deep stack
        # trace from the sampler, autodiff, etc.
        ln_like_generated(system, arr2nt(initial_θ_0))

        ln_prior_transformed(initial_θ_0,false)


        # We use let blocks to prevent type instabilities from closures
        # A function barrier would also work.
        ℓπcallback, ∇ℓπcallback = (function(
            arr2nt,
            system,
            Bijector_invlinkvec,
            ln_prior_transformed,
            ln_like_generated,
            D,
        )

            # Capture these variables in a let binding to improve performance
            # We also set up temporary storage to reduce allocations
            # ForwardDiff is used to compute the likelihood gradients using the in-place 
            # API. This ensures type stability.
            function ℓπcallback(
                θ_transformed,
                system=system,
                arr2nt=arr2nt,
                Bijector_invlinkvec=Bijector_invlinkvec,
                ln_prior_transformed=ln_prior_transformed,
                ln_like_generated=ln_like_generated;sampled=true)::eltype(θ_transformed)

                lpost = zero(eltype(θ_transformed))
                # Stop right away if we are given non-finite arguments
                if any(!isfinite, θ_transformed)
                    # @warn "non finite parameters encountered (maxlog=1)" θ_transformed maxlog=1
                    lpost += NaN
                    return lpost
                end
                # Transform back from the unconstrained support to constrained support for the likelihood function
                θ_natural = @inline Bijector_invlinkvec(θ_transformed)
                θ_structured = @inline arr2nt(θ_natural)
                lpost += @inline ln_prior_transformed(θ_natural,sampled)
                # Don't compute likelihood if we fell outside the prior bounds
                if !isfinite(lpost)
                    @warn "non finite log prior (maxlog=1)" lpost maxlog=1
                    return lpost
                end
                lpost += @inline ln_like_generated(system, θ_structured)
                # if !isfinite(lpost)
                #     # TODO: check for performance impact here
                #     # Display parameters that caused an invalid log-likelihood to be calculated
                #     # Strip off any forward diff Dual tags, as these make it impossible to see
                #     # what's going on.
                #     θ_transformed_primals = ForwardDiff.value.(θ_transformed)
                #     θ_structured = arr2nt(Bijector_invlinkvec(θ_transformed_primals))
                #     llike = ln_like_generated(system, θ_structured)
                #     @warn "Invalid log likelihood encountered. (maxlog=1)" θ=θ_structured llike θ_transformed=θ_transformed_primals  maxlog=1
                # end
                return lpost
            end

            args = (
                system,
                arr2nt,
                Bijector_invlinkvec,
                ln_prior_transformed,
                ln_like_generated,
            )
            # Test likelihood function immediately to give user a clean error
            # if it fails for some reason.
            ℓπcallback(initial_θ_0_t,args...)
            # Also Display their run time. If something is egregiously wrong we'll notice something
            # here in the output logs.
            if verbosity >= 1
                (function(ℓπcallback, θ, args)
                    @showtime ℓπcallback(θ, args...)
                end)(ℓπcallback, initial_θ_0_t, args)
            end

            if autodiff == :Enzyme
                if !isdefined(Main, :Enzyme)
                    error("To use the :Enzyme autodiff backend, load the Enzyme package first: `using Enzyme`")
                end
                # Enzyme mode:
                ∇ℓπcallback = let Enzyme=Main.Enzyme, diffresult = copy(initial_θ_0_t), system=system, system_shadow=deepcopy(system), ℓπcallback=ℓπcallback
                    oh = Enzyme.onehot(diffresult)
                    forward = function (θ_t)
                        primal, out = Enzyme.autodiff(
                            Enzyme.Forward,
                            ℓπcallback,
                            Enzyme.BatchDuplicated,
                            Enzyme.BatchDuplicated(θ_t,oh),
                            Enzyme.Const.(args)...
                            # Enzyme.Const(system),
                            # Enzyme.Const(arr2nt),
                            # Enzyme.Const(Bijector_invlinkvec),
                            # Enzyme.Const(ln_prior_transformed),#ln_prior_transformed_dup,
                            # Enzyme.Const(ln_like_generated),#ln_like_generated_dup,
                        )
                        diffresult .= tuple(out...)
                        return primal, diffresult
                    end
                    # reverse = function (θ_t)
                    #     fill!(diffresult,0)
                    #     θ_t_dup  = Enzyme.Duplicated(θ_t,diffresult)
                    #     # system_dup  = Enzyme.Duplicated(system,deepcopy(system))
                    #     # arr2nt_dup  = Enzyme.Duplicated(arr2nt, deepcopy(arr2nt))
                    #     # Bijector_invlinkvec_dup  = Enzyme.Duplicated(Bijector_invlinkvec, deepcopy(Bijector_invlinkvec))
                    #     # ln_prior_transformed_dup  = Enzyme.Duplicated(ln_prior_transformed, deepcopy(ln_prior_transformed))
                    #     # ln_like_generated_dup  = Enzyme.Duplicated(ln_like_generated, deepcopy(ln_like_generated))
                    #     out, primal = Enzyme.autodiff(
                    #         # Enzyme.Reverse,
                    #         Enzyme.ReverseWithPrimal,
                    #         (ℓπcallback),
                    #         θ_t_dup,
                    #         # Enzyme.Duplicated.(args, deepcopy.(args))...
                    #         # Enzyme.Const(system,),# system_dup, #Enzyme.Const(system,),
                    #         # Enzyme.Const(arr2nt,),# arr2nt_dup, #Enzyme.Const(arr2nt,),
                    #         Enzyme.Const(Bijector_invlinkvec,),# Bijector_invlinkvec_dup, #Enzyme.Const(Bijector_invlinkvec,),
                    #         # Enzyme.Const(ln_prior_transformed,),# ln_prior_transformed_dup, #Enzyme.Const(ln_prior_transformed,),
                    #         # Enzyme.Const(ln_like_generated,),# ln_like_generated_dup, #Enzyme.Const(ln_like_generated,),
                    #     )
                    #     return primal, diffresult
                    # end 
                    # tforward = minimum(
                    #     @elapsed forward(initial_θ_0_t)
                    #     for _ in 1:100
                    # )
                    # treverse = minimum(
                    #     @elapsed reverse(initial_θ_0_t)
                    #     for _ in 1:100
                    # )
                    # if tforward > treverse
                    #     verbosity > 2  && @info "selected reverse mode AD" tforward treverse
                    #     reverse
                    # else
                    #     verbosity > 2  && @info "selected forward mode AD" tforward treverse
                    #     forward
                    # end
                end

            elseif autodiff == :FiniteDiff
            
                ∇ℓπcallback = let diffresult = copy(initial_θ_0_t),ℓπcallback=ℓπcallback
                    function (θ_t)
                        primal = ℓπcallback(θ_t)
                        Main.FiniteDiff.finite_difference_gradient!(diffresult, ℓπcallback, θ_t)
                        return primal, diffresult
                    end
                end
            
            elseif autodiff == :Zygote

                # Zygote mode:
                ∇ℓπcallback = function (θ_t)
                    Main.Zygote.gradient(ℓπcallback, θ_t)
                end

            elseif autodiff == :ForwardDiff

                # https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/#Fixing-NaN/Inf-Issues

                # Test likelihood function gradient immediately to give user a clean error
                # if it fails for some reason.
                ForwardDiff.gradient(ℓπcallback, initial_θ_0_t)

                # ForwardDiff mode:
                # Create temporary storage space for gradient computations
                diffresults = [
                    DiffResults.GradientResult(collect(initial_θ_0_t))
                    for _ in 1:Threads.nthreads()
                ]

                # Perform dynamic benchmarking to pick a ForwardDiff chunk size.
                # We're going to call this thousands of times so worth a few calls
                # to get this optimized.
                if isnothing(chunk_sizes)
                    if D < 100
                        # If less than 20 dimensional, a single chunk with ForwardDiff
                        # is almost always optimial.
                        chunk_sizes = D
                    else
                        chunk_sizes = unique([1; D÷4; D])
                    end
                end
                ideal_chunk_size_i = argmin(map(chunk_sizes) do chunk_size
                    cfg = ForwardDiff.GradientConfig(ℓπcallback, initial_θ_0_t, ForwardDiff.Chunk{chunk_size}());
                    ForwardDiff.gradient!(diffresults[1], ℓπcallback, initial_θ_0_t, cfg, Val{false}())
                    t = minimum(
                        @elapsed ForwardDiff.gradient!(diffresults[1], ℓπcallback, initial_θ_0_t, cfg, Val{false}())
                        for _ in 1:10
                    )

                    verbosity >= 3 && @info "Tuning autodiff" chunk_size t
                    return t
                end)
                ideal_chunk_size =  chunk_sizes[ideal_chunk_size_i]
                verbosity >= 1 && @info "Selected auto-diff chunk size" ideal_chunk_size

                cfg = ForwardDiff.GradientConfig(ℓπcallback, initial_θ_0_t, ForwardDiff.Chunk{ideal_chunk_size}());
                ∇ℓπcallback = let cfg=cfg, diffresults=diffresults, ℓπcallback=ℓπcallback
                    function (θ_transformed)
                        # TODO: this is not safe for tasks that migrate across threads! Would need task-local storage instead.
                        result = ForwardDiff.gradient!(diffresults[Threads.threadid()], ℓπcallback, θ_transformed, cfg, Val{false}())
                        return DiffResults.value(result), DiffResults.gradient(result)
                        # return DiffResults.gradient(result)
                    end
                end
            else
                error("Unsupported option for autodiff: $autodiff. Valid options are :ForwardDiff (default), :Enzyme, and :Zygote.")
            end


            # Run the callback once right away. If there is a coding error in the users 
            # model, we want to surface it ASAP.
            ∇ℓπcallback(initial_θ_0_t) 
            if verbosity >= 1
                (function(∇ℓπcallback, θ)
                    @showtime ∇ℓπcallback(θ)
                end)(∇ℓπcallback, initial_θ_0_t)
            end


            ℓπcallback, ∇ℓπcallback
        end)(
            arr2nt,
            system,
            Bijector_invlinkvec,
            ln_prior_transformed,
            ln_like_generated,
            D
        )
        
        # Perform some quick diagnostic checks to warn users for performance-gtochas
        out_type_model = Core.Compiler.return_type(ℓπcallback, typeof((initial_θ_0_t,)))
        out_type_model_grad = Core.Compiler.return_type(∇ℓπcallback, typeof((initial_θ_0_t,)))
        out_type_arr2nt = Core.Compiler.return_type(arr2nt, typeof((initial_θ_0_t,)))
        out_type_prior = Core.Compiler.return_type(ln_prior_transformed, typeof((initial_θ_0,false,)))
        out_type_like = Core.Compiler.return_type(ln_like_generated, typeof((system,arr2nt(initial_θ_0),)))

        if isconcretetype(out_type_prior) &&
            isconcretetype(out_type_like) &&
            isconcretetype(out_type_arr2nt) &&
            !isconcretetype(out_type_model)
            @warn "\nThis model's log density function is not type stable, but all of its components are. \nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_model_grad)
            @warn "\nThis model's log density gradient is not type stable, but all of its components are. \nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub."
            end
        if !isconcretetype(out_type_prior)
            @warn "\nThis model's prior sampler does not appear to be type stable, which will likely hurt sampling performance.\nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_like)
            @warn "\nThis model's likelihood function does not appear to be type stable, which will likely hurt sampling performance.\nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_arr2nt)
            @warn "\nThis model specification is not type-stable, which will likely hurt sampling performance.\nCheck for global variables used within your model definition, and prepend these with `\$`.\nIf that doesn't work, you could trying running:\n`@code_warntype model.arr2nt(randn(model.D))` for a bit more information.\nFor assistance, please file an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end

        # Return fully concrete type wrapping all these functions
        new{
            D,
            typeof(ℓπcallback),
            typeof(∇ℓπcallback),
            typeof(system),
            typeof(Bijector_linkvec),
            typeof(Bijector_invlinkvec),
            typeof(arr2nt),
            typeof(sample_priors)
        }(
            D,
            autodiff,
            ℓπcallback,
            ∇ℓπcallback,
            system,
            Bijector_linkvec,
            Bijector_invlinkvec,
            arr2nt,
            sample_priors,
            nothing # no starting points set
        )
    end
end
LogDensityProblems.logdensity(p::LogDensityModel, θ) = p.ℓπcallback(θ)
LogDensityProblems.logdensity_and_gradient(p::LogDensityModel, θ) = p.∇ℓπcallback(θ)
LogDensityProblems.dimension(p::LogDensityModel{D}) where D = D
LogDensityProblems.capabilities(::Type{<:LogDensityModel}) = LogDensityProblems.LogDensityOrder{1}()

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::LogDensityModel)
    L = _count_epochs(p.system)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) and $(L) epochs with fields .ℓπcallback and .∇ℓπcallback")
end
function Base.show(io::IO, @nospecialize p::LogDensityModel)
    L = _count_epochs(p.system)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) and $(L) epochs with fields .ℓπcallback and .∇ℓπcallback")
end
