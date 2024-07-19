
using LogDensityProblems

# Define the target distribution using the `LogDensityProblem` interface
# TODO: in future, just roll this all into the System type.
struct LogDensityModel{Tℓπ,T∇ℓπ,TSys,TLink,TInvLink,TArr2nt}
    D::Int
    ℓπcallback::Tℓπ
    ∇ℓπcallback::T∇ℓπ
    system::TSys
    link::TLink
    invlink::TInvLink
    arr2nt::TArr2nt
    function LogDensityModel(system::System; autodiff=:ForwardDiff, verbosity=0, chunk_sizes=nothing)
        verbosity >= 1 && @info "Preparing model"

        # Choose parameter dimensionality and initial parameter value
        initial_θ_0 = sample_priors(system)
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
            (θ) -> Bijectors.link.(priors_vec, θ)
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
        ℓπcallback, ∇ℓπcallback = let arr2nt=arr2nt,
                                      system=system,
                                      Bijector_invlinkvec=Bijector_invlinkvec,
                                      ln_prior_transformed=ln_prior_transformed,
                                      ln_like_generated=ln_like_generated,
                                      D=D

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
                ln_like_generated=ln_like_generated;sampled=true)
                # Transform back from the unconstrained support to constrained support for the likelihood function
                θ_natural = Bijector_invlinkvec(θ_transformed)
                θ_structured = arr2nt(θ_natural)
                lprior = @inline ln_prior_transformed(θ_natural,sampled)
                # lprior = @inline ln_prior_transformed(θ_transformed)
                # CAUTION: This inline annotation is necessary for correct gradients from Enzyme. Yikes!
                llike  = @inline ln_like_generated(system, θ_structured)
                lpost = lprior+llike
                # if !isfinite(lprior)
                #     @warn "Invalid log prior encountered. This likely indicates a problem with the prior support." θ=θ_structured lprior θ_transformed maxlog=5
                # end
                # if !isfinite(llike)
                #     @warn "Invalid log likelihood encountered." θ=θ_structured llike θ_transformed maxlog=5
                # end
                return lpost
            end

            # Test likelihood function immediately to give user a clean error
            # if it fails for some reason.
            ℓπcallback(initial_θ_0_t)
            # Also Display their run time. If something is egregously wrong we'll notice something
            # here in the output logs.
            if verbosity >= 1
                (function(ℓπcallback, θ)
                    @showtime ℓπcallback(θ)
                end)(ℓπcallback, initial_θ_0_t)
            end

            if autodiff == :Enzyme
                if !isdefined(Main, :Enzyme)
                    error("To use the :Enzyme autodiff backend, load the Enzyme package first: `using Enzyme`")
                end
                # Enzyme mode:
                ∇ℓπcallback = let Enzyme=Main.Enzyme, diffresult = copy(initial_θ_0_t), system=system, system_shadow=deepcopy(system), ℓπcallback=ℓπcallback
                    # oh = Enzyme.onehot(diffresult)
                    # forward = function (θ_t)
                    #     primal, out = Enzyme.autodiff(
                    #         Enzyme.Forward,
                    #         ℓπcallback,
                    #         Enzyme.BatchDuplicated,
                    #         Enzyme.BatchDuplicated(θ_t,oh),
                    #         Enzyme.Const(system),
                    #         Enzyme.Const(arr2nt),
                    #         Enzyme.Const(Bijector_invlinkvec)
                    #         Enzyme.Const(ln_prior_transformed),#ln_prior_transformed_dup,
                    #         Enzyme.Const(ln_like_generated),#ln_like_generated_dup,
                    #     )
                    #     diffresult .= tuple(out...)
                    #     return primal, diffresult
                    # end
                    reverse = function (θ_t)
                        fill!(diffresult,0)
                        θ_t_dup  = Enzyme.Duplicated(θ_t,diffresult)
                        system_dup  = Enzyme.Duplicated(system,system_shadow)
                        arr2nt_dup  = Enzyme.Duplicated(arr2nt, deepcopy(arr2nt))
                        Bijector_invlinkvec_dup  = Enzyme.Duplicated(Bijector_invlinkvec, deepcopy(Bijector_invlinkvec))
                        ln_prior_transformed_dup  = Enzyme.Duplicated(ln_prior_transformed, deepcopy(ln_prior_transformed))
                        ln_like_generated_dup  = Enzyme.Duplicated(ln_like_generated, deepcopy(ln_like_generated))
                        out, primal = Enzyme.autodiff(
                            # Enzyme.Reverse,
                            Enzyme.ReverseWithPrimal,
                            ℓπcallback,
                            θ_t_dup,
                            system_dup, #Enzyme.Const(system,),
                            arr2nt_dup, #Enzyme.Const(arr2nt,),
                            Bijector_invlinkvec_dup, #Enzyme.Const(Bijector_invlinkvec,),
                            ln_prior_transformed_dup, #Enzyme.Const(ln_prior_transformed,),
                            ln_like_generated_dup, #Enzyme.Const(ln_like_generated,),
                        )
                        return primal, diffresult
                    end 
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
                        reverse
                    # else
                    #     verbosity > 2  && @info "selected forward mode AD" tforward treverse
                    #     forward
                    # end
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
                    Main.Zygote.gradient(ℓπcallback, θ_t)
                end

            elseif autodiff == :ForwardDiff

                # https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/#Fixing-NaN/Inf-Issues
                set_preferences!(ForwardDiff, "nansafe_mode" => true)

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
                    if D < 40
                        # If less than 20 dimensional, a single chunk with ForwardDiff
                        # is almost always optimial.
                        chunk_sizes = D
                    else
                        chunk_sizes = unique([1; 2:2:D; D])
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
        end
        

    
        # Return fully concrete type wrapping all these functions
        new{
            typeof(ℓπcallback),
            typeof(∇ℓπcallback),
            typeof(system),
            typeof(Bijector_linkvec),
            typeof(Bijector_invlinkvec),
            typeof(arr2nt),
        }(
            D,
            ℓπcallback,
            ∇ℓπcallback,
            system,
            Bijector_linkvec,
            Bijector_invlinkvec,
            arr2nt
        )
    end
end
LogDensityProblems.logdensity(p::LogDensityModel, θ) = p.ℓπcallback(θ)
LogDensityProblems.logdensity_and_gradient(p::LogDensityModel, θ) = p.∇ℓπcallback(θ) # TODO: may need to copy vector
LogDensityProblems.dimension(p::LogDensityModel) = p.D
LogDensityProblems.capabilities(::Type{<:LogDensityModel}) = LogDensityProblems.LogDensityOrder{1}()

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::LogDensityModel)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) with fields .ℓπcallback and .∇ℓπcallback")
end
function Base.show(io::IO, @nospecialize p::LogDensityModel)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) with fields .ℓπcallback and .∇ℓπcallback")
end
