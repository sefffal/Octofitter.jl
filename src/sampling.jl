using DiffResults
using AbstractDifferentiation
AD = AbstractDifferentiation
using LinearAlgebra
using Preferences

export sample_priors


sample_priors(arg::Union{Planet,System}, args...; kwargs...) = sample_priors(Random.default_rng(), arg, args...; kwargs...)
sample_priors(rng::Random.AbstractRNG, planet::Planet) = rand.(rng, ComponentArray(planet.priors.priors))
sample_priors(rng::Random.AbstractRNG, planet::Planet, N::Number) = [sample_priors(rng, planet) for _ in 1:N]

function sample_priors(rng::Random.AbstractRNG, system::System)
    priors_flat_sampled = map(((k,v),)->rand(rng, v), Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
    return priors_flat_sampled
end


sample_priors(rng::Random.AbstractRNG, system::System, N::Number) = [sample_priors(rng, system) for _ in 1:N]




function guess_starting_position(system::System, args...; kwargs...)
    return guess_starting_position(Random.default_rng(), system, args...; kwargs...)
end
function guess_starting_position(rng::Random.AbstractRNG, system::System, N=500_000)

    # TODO: this shouldn't have to allocate anything, we can just loop keeping the best.
    θ = sample_priors(rng, system, N)
    arr2nt = Octofitter.make_arr2nt(system) 
    ln_like = Octofitter.make_ln_like(system, arr2nt(sample_priors(rng, system)))

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
    if haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        return VisualOrbit((;
            M=chain["M"][i],
            plx=chain["plx"][i],
            i=chain[pk*"_i"][i],
            Ω=chain[pk*"_Ω"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            a=chain[pk*"_a"][i],
        ))
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*"_A")) && haskey(chain, Symbol(pk*"_B")) && haskey(chain, Symbol(pk*"_G"))&& haskey(chain, Symbol(pk*"_F"))
        return ThieleInnesOrbit((;
            M=chain["M"][i],
            plx=chain["plx"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            A=chain[pk*"_A"][i],
            B=chain[pk*"_B"][i],
            F=chain[pk*"_F"][i],
            G=chain[pk*"_G"][i],
        ))
    elseif haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        return KepOrbit((;
            M=chain["M"][i],
            i=chain[pk*"_i"][i],
            Ω=chain[pk*"_Ω"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            a=chain[pk*"_a"][i],
        ))
    elseif haskey(chain, :M) && haskey(chain, :rv)
        return RadialVelocityOrbit((;
            M=chain["M"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            a=chain[pk*"_a"][i],
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
    if haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        Ms=chain["M"]
        plxs=chain["plx"]
        is=chain[pk*"_i"]
        Ωs=chain[pk*"_Ω"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        τs=chain[pk*"_τ"]
        as=chain[pk*"_a"]
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
    elseif haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        Ms=chain["M"]
        is=chain[pk*"_i"]
        Ωs=chain[pk*"_Ω"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        τs=chain[pk*"_τ"]
        as=chain[pk*"_a"]
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
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*"_A")) && haskey(chain, Symbol(pk*"_B"))
        Ms=chain["M"]
        plxs=chain["plx"]
        As=chain[pk*"_A"]
        Bs=chain[pk*"_B"]
        Fs=chain[pk*"_F"]
        Gs=chain[pk*"_G"]
        es=chain[pk*"_e"]
        τs=chain[pk*"_τ"]
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
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        τs=chain[pk*"_τ"]
        as=chain[pk*"_a"]
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
    is=chain[:,pk*"_i"]
    Ωs=chain[:,pk*"_Ω"]
    ωs=chain[:,pk*"_ω"]
    es=chain[:,pk*"_e"]
    τs=chain[:,pk*"_τ"]
    as=chain[:,pk*"_a"]
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
    function LogDensityModel(system::System; autodiff=:ForwardDiff, verbosity=0)

        # Choose parameter dimensionality and initial parameter value
        initial_θ_0 = sample_priors(system)
        D = length(initial_θ_0)

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
        ln_prior_transformed(initial_θ_0_t)


        verbosity >= 1 && @info "Preparing model"
        # Capture these variables in a let binding to improve performance
        # We also set up temporary storage to reduce allocations
        # ForwardDiff is used to compute the likelihood gradients using the in-place 
        # API. This ensures type stability.
        function ℓπcallback(θ_transformed, system=system)
            # Transform back from the unconstrained support to constrained support for the likelihood function
            θ_natural = Bijector_invlinkvec(θ_transformed)
            θ_structured = arr2nt(θ_natural)
            lp = ln_prior_transformed(θ_transformed)
            ll = ln_like_generated(system, θ_structured)
            return lp+ll
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

            # https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/#Fixing-NaN/Inf-Issues
            set_preferences!(ForwardDiff, "nansafe_mode" => true)

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
                    # return DiffResults.value(result), DiffResults.gradient(result)
                    return DiffResults.gradient(result)
                end
            end
        else
            error("Unsupported option for autodiff: $autodiff. Valid options are :ForwardDiff (default), :Enzyme, and :Zygote.")
        end
        
            
        # Display their run time. If something is egregously wrong we'll notice something
        # here in the output logs.
        ℓπcallback(initial_θ_0_t)
        ∇ℓπcallback(initial_θ_0_t) 
        verbosity >= 1 && @showtime ℓπcallback(initial_θ_0_t)
        verbosity >= 1 && @showtime ∇ℓπcallback(initial_θ_0_t)
    
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
# TODO: use diff results API to avoid calculating primal twice!
LogDensityProblems.logdensity_and_gradient(p::LogDensityModel, θ) = (p.ℓπcallback(θ), p.∇ℓπcallback(θ)) # TODO: may need to copy vector
LogDensityProblems.dimension(p::LogDensityModel) = p.D
LogDensityProblems.capabilities(::Type{<:LogDensityModel}) = LogDensityProblems.LogDensityOrder{1}()

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::LogDensityModel)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) with fields .ℓπcallback and .∇ℓπcallback")
end


# Fallback when no random number generator is provided (as is usually the case)
function advancedhmc(model::LogDensityModel, target_accept::Number=0.8, ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial(); kwargs...)
    return advancedhmc(Random.default_rng(), model, target_accept, ensemble; kwargs...)
end

"""
The method signature of Octofitter.hmc is as follows:

    advancedhmc(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel
        target_accept::Number=0.8,
        ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
        num_chains=1,
        adaptation,
        iterations,
        thinning=1,
        discard_initial=adaptation,
        tree_depth=12,
        initial_samples=50_000,
        initial_parameters=nothing,
        step_size=nothing,
        verbosity=2,
    )

The only required arguments are system, adaptation, and iterations.
The two positional arguments are system, the model you wish to sample;
and target_accept, the acceptance rate that should be targeted during
windowed adaptation. During this time, the step size and mass matrix
will be adapted (see AdvancedHMC.jl for more information). The number
of steps taken during adaptation is controlled by adaptation. You can
prevent these samples from being dropped by pasing include_adaptation=false.
The total number of posterior samples produced are given by iterations.
These include the adaptation steps that may be discarded.
tree_depth controls the maximum tree depth of the sampler.
initial_parameters is an optional way to pass a starting point for the chain.
If you don't pass a default position, one will be selected by drawing
initial_samples from the priors.
The sample with the highest posterior value will be used as the starting point.
"""
function advancedhmc(
    rng::Random.AbstractRNG,
    model::LogDensityModel,
    target_accept::Number=0.8,
    ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
    num_chains=1,
    adaptation,
    iterations,
    thinning=1,
    discard_initial=adaptation,
    tree_depth=12,
    initial_samples=50_000,
    initial_parameters=nothing,
    verbosity=2,
    autodiff=:ForwardDiff,
)

    if ensemble != MCMCSerial()
        @warn "TODO: In-place model gradients currently not supported with MCMCThreads"
    end

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(rng, model.system)
    D = length(initial_θ_0)

    # ln_prior_transformed = make_ln_prior_transformed(system)
    # # ln_prior = make_ln_prior(system)
    # arr2nt = Octofitter.make_arr2nt(system) 
    # ln_like_generated = make_ln_like(system, arr2nt(initial_θ_0))

    # priors_vec = _list_priors(system)
    # Bijector_invlinkvec = make_Bijector_invlinkvec(priors_vec)
    initial_θ_0_t = model.link(initial_θ_0)
    # arr2nt = Octofitter.make_arr2nt(system)

    # # Test out model likelihood and prior computations. This way, if they throw
    # # an error, we'll see it right away instead of burried in some deep stack
    # # trace from the sampler, autodiff, etc.
    # ln_like_generated(system, arr2nt(initial_θ_0))
    # ln_prior_transformed(initial_θ_0_t)



    # # Uncomment for debugging model type-stability
    # Main.InteractiveUtils.@code_warntype ℓπ(initial_θ_0_t)
    # Main.InteractiveUtils.@code_warntype ∇ℓπ(initial_θ_0_t)
    # return

    # Guess initial starting positions by drawing from priors a bunch of times
    # and picking the best one (highest likelihood).
    # Then transform that guess into our unconstrained support
    if isnothing(initial_parameters)
        verbosity >= 1 && @info "Guessing a starting location by sampling from prior" initial_samples
        initial_θ, mapv = guess_starting_position(rng,model.system,initial_samples)
        verbosity > 2 && @info "Found starting location" θ=model.arr2nt(initial_θ)
        # Transform from constrained support to unconstrained support
        initial_θ_t = model.link(initial_θ)
    else
        initial_θ = initial_parameters
        # Transform from constrained support to unconstrained support
        initial_θ_t = model.link(initial_θ)
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

    # verbosity >= 3 && @info(
    #     "Pathfinder results",
    #     ℓπ(θ)=ℓπ(result_pf.fit_distribution.μ),
    #     mode=arr2nt(Bijectors.invlink.(priors_vec, result_pf.fit_distribution.μ)),
    #     inv_metric=Matrix(result_pf.fit_distribution.Σ)
    # )

    # # # Return a chains object with the resampled pathfinder draws
    # # # Transform samples back to constrained support
    # # pathfinder_samples = map(eachcol(result_pf.draws)) do θ_t
    # #     Bijectors.invlink.(priors_vec, θ_t)
    # # end
    # # pathfinder_chain =  Octofitter.result2mcmcchain(system, arr2nt.(pathfinder_samples))
    # # pathfinder_chain_with_info = MCMCChains.setinfo(
    # #     pathfinder_chain,
    # #     (;
    # #         start_time,
    # #         stop_time,
    # #         model=system,
    # #         result_pf,
    # #     )
    # # )

    # # Start using a draw from the typical set as estimated by Pathfinder
    # finite_pathfinder_draws = filter(row->all(isfinite, row), collect(eachrow(result_pf.draws)))
    # if length(finite_pathfinder_draws) > 0 && all(isfinite, Matrix(result_pf.fit_distribution.Σ))
        
    #     # initial_θ_t = last(finite_pathfinder_draws)

    #     verbosity >= 3 && @info "Creating metric"

    #     # # Use the metric found by Pathfinder for HMC sampling
    #     # metric = Pathfinder.RankUpdateEuclideanMetric(result_pf.fit_distribution.Σ)
        
    #     # Start with found pathfinder metric then adapt a dense metric:
    #     # metric = DenseEuclideanMetric(Matrix(result_pf.fit_distribution.Σ))
    #     # metric = DiagEuclideanMetric(diag(Matrix(result_pf.fit_distribution.Σ)))

    #     metric = DenseEuclideanMetric(collect(Diagonal(Matrix(result_pf.fit_distribution.Σ))))

    # else
    #     @warn "Pathfinder failed to provide a finite initial draw and metric. Check your model. Starting from initial guess instead."
        # Fit a dense metric from scratch
        # metric = DiagEuclideanMetric(D)
    # end

    # verbosity >= 3 && @info "Creating model" 
    # # model = AdvancedHMC.DifferentiableDensityModel(ℓπ, ∇ℓπ)
    # model = LogDensityModel(system; autodiff, verbosity)

    metric = DenseEuclideanMetric(model.D)


    verbosity >= 3 && @info "Creating hamiltonian"
    # hamiltonian = Hamiltonian(metric, model.ℓπcallback, model.∇ℓπcallback)
    hamiltonian = Hamiltonian(metric, model)
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


    verbosity >= 1 && @info "Sampling, beginning with adaptation phase..."
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
    function callback(rng, logdensitymodel, sampler, transition, state, iteration; kwargs...)
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
            θ = model.invlink(transition.z.θ)
            θ_res = model.arr2nt(θ)
            # Fill the remaining width of the terminal with info
            max_width = displaysize(stdout)[2]-34
            θ_str = string(θ_res)
            θ_str_trunc = θ_str[begin:prevind(θ_str, min(end,max_width))]
            θ_message = "θ="*θ_str_trunc*"_.."
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

    mc_samples_all_chains = AbstractMCMC.sample(
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
            θ = model.invlink(θ_t)
            return θ
        end
        chain_res = model.arr2nt.(samples)
        push!(chains, Octofitter.result2mcmcchain(chain_res))
        # push!(chains, Strapping.deconstruct(chain_res))
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
            model=model.system,
            logpost=logposts_mat,
            states=mc_samples_all_chains,
            # pathfinder=pathfinder_chain_with_info,
            # _restart=(;
            #     model,
            #     sampler,
            #     adaptor,
            #     state = last.(mc_samples_all_chains)
            # )
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
function result2mcmcchain(chain_in)

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
            push!(pairs, Symbol(pl, '_', key) =>  nt.planets[pl][key])
        end
    end
    return namedtuple(pairs)

end
