using DiffResults
using LinearAlgebra
using Preferences

export sample_priors

# TODO: consolidate all these functions. There is a lot of duplication.

sample_priors(arg::Union{Planet,System}, args...; kwargs...) = sample_priors(Random.default_rng(), arg, args...; kwargs...)
# sample_priors(rng::Random.AbstractRNG, planet::Planet) = rand.(rng, ComponentArray(planet.priors.priors))
sample_priors(rng::Random.AbstractRNG, planet::Planet, N::Number) = [sample_priors(rng, planet) for _ in 1:N]

function sample_priors(rng::Random.AbstractRNG, system::System)
    priors_flat_sampled = map(((k,v),)->rand(rng, v), Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
    return priors_flat_sampled
end

function priors_flat(system::System)
    return priors_flat = map(((k,v),)->v, Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
end


sample_priors(rng::Random.AbstractRNG, system::System, N::Number) = [sample_priors(rng, system) for _ in 1:N]


# Function to give the parameter names as a flat vector of symbols. Only returns
# active parameters (i.e.) and not any derived variables.
function list_parameter_names(system::System)
    return map(((k,v),)->k, Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
end

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

"""
    construct_elements(θ_system, θ_planet)

Given a named tuple for of parameters from a System (θ_system) and Planet (θ_planet),
return a `Visual{KepOrbit} PlanetOrbits.jl.
"""
function construct_elements(::Type{Visual{KepOrbit}}, θ_system, θ_planet)
    return Visual{KepOrbit}(;(;
        θ_system.M,
        θ_system.plx,
        θ_planet.i,
        θ_planet.Ω,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    )...)
end
function construct_elements(::Type{KepOrbit}, θ_system, θ_planet)
    return KepOrbit(;(;
        θ_system.M,
        θ_planet.i,
        θ_planet.Ω,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    )...)
end
function construct_elements(::Type{ThieleInnesOrbit}, θ_system, θ_planet)
    return ThieleInnesOrbit(;(;
        θ_system.M,
        θ_system.plx,
        θ_planet.A,
        θ_planet.B,
        θ_planet.F,
        θ_planet.G,
        θ_planet.e,
        θ_planet.τ,
    )...)
end
function construct_elements(::Type{RadialVelocityOrbit}, θ_system, θ_planet)
    return RadialVelocityOrbit(;(;
        θ_system.M,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    )...)
end
function construct_elements(::Type{CartesianOrbit}, θ_system, θ_planet)
    return CartesianOrbit(;(;
        θ_system.M,
        θ_planet.x,
        θ_planet.y,
        θ_planet.z,
        θ_planet.vx,
        θ_planet.vy,
        θ_planet.vz,
    )...)
end
function construct_elements(::Type{Visual{CartesianOrbit}}, θ_system, θ_planet)
    return Visual{CartesianOrbit}(;(;
        θ_system.M,
        θ_system.plx,
        θ_planet.x,
        θ_planet.y,
        θ_planet.z,
        θ_planet.vx,
        θ_planet.vy,
        θ_planet.vz,
    )...)
end

"""
    construct_elements(chains, :b, 4)

Given a Chains object, a symbol matching the name of a planet, and an index,
construct a `Visual{KepOrbit} DirectOrbits of that planet from that
index of the chains.
"""
function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, i::Union{Integer,CartesianIndex})
    pk = string(planet_key)
    if haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        return Visual{KepOrbit}(;(;
            M=chain["M"][i],
            plx=chain["plx"][i],
            i=chain[pk*"_i"][i],
            Ω=chain[pk*"_Ω"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            a=chain[pk*"_a"][i],
        )...)
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*"_A")) && haskey(chain, Symbol(pk*"_B")) && haskey(chain, Symbol(pk*"_G"))&& haskey(chain, Symbol(pk*"_F"))
        return ThieleInnesOrbit(;(;
            M=chain["M"][i],
            plx=chain["plx"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            A=chain[pk*"_A"][i],
            B=chain[pk*"_B"][i],
            F=chain[pk*"_F"][i],
            G=chain[pk*"_G"][i],
        )...)
    elseif haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        return KepOrbit(;(;
            M=chain["M"][i],
            i=chain[pk*"_i"][i],
            Ω=chain[pk*"_Ω"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            a=chain[pk*"_a"][i],
        )...)
    elseif haskey(chain, :M) && haskey(chain, :rv)
        return RadialVelocityOrbit(;(;
            M=chain["M"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            τ=chain[pk*"_τ"][i],
            a=chain[pk*"_a"][i],
        )...)
    else
        error("Unrecognized columns")
    end
end

"""
    construct_elements(chains, :b, [4,5,10])

Given a Chains object, a symbol matching the name of a planet, and an array of indices,
construct a `Visual{KepOrbit} DirectOrbits of that planet from those indices
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
            Visual{KepOrbit}(;(;
                M=Ms[i],
                plx=plxs[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                τ=τs[i],
                a=as[i],
            )...)
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
            KepOrbit(;(;
                M=Ms[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                τ=τs[i],
                a=as[i],
            )...)
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
            ThieleInnesOrbit(;(;
                M=Ms[i],
                plx=plxs[i],
                e=es[i],
                τ=τs[i],
                A=As[i],
                B=Bs[i],
                F=Fs[i],
                G=Gs[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_vx")) && haskey(chain, Symbol("plx"))
        M=chain["M"]
        plx=chain["plx"]

         x=chain[ pk*"_x"]
         y=chain[ pk*"_y"]
         z=chain[ pk*"_z"]
        vx=chain[pk*"_vx"]
        vy=chain[pk*"_vy"]
        vz=chain[pk*"_vz"]
        tref=chain[pk*"_tref"]
        return map(ii) do i
            Visual{CartesianOrbit}(;(;
                x = x[i],
                y = y[i],
                z = z[i],
                vx = vx[i],
                vy = vy[i],
                vz = vz[i],           
                M = M[i],
                tref = tref[i],
                plx = plx[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_vx"))
        M=chain["M"]
         x=chain[ pk*"_x"]
         y=chain[ pk*"_y"]
         z=chain[ pk*"_z"]
        vx=chain[pk*"_vx"]
        vy=chain[pk*"_vy"]
        vz=chain[pk*"_vz"]
        tref=chain[pk*"_tref"]
        return map(ii) do i
            CartesianOrbit(;(;
                x = x[i],
                y = y[i],
                z = z[i],
                vx = vx[i],
                vy = vy[i],
                vz = vz[i],           
                M = M[i],
                tref = tref[i],
            )...)
        end
    elseif haskey(chain, Symbol("M"))
        Ms=chain["M"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        τs=chain[pk*"_τ"]
        as=chain[pk*"_a"]
        return map(ii) do i
            RadialVelocityOrbit(;(;
                M=Ms[i],
                ω=ωs[i],
                e=es[i],
                τ=τs[i],
                a=as[i],
            )...)
        end
    else
        error("Unrecognized chain format")
    end
end
construct_elements(chain::Chains, planet_key::Union{String,Symbol}, ii::Colon) = construct_elements(chain, planet_key, 1:size(chain,1)*size(chain,3))
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
        Visual{KepOrbit}((;
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
    function LogDensityModel(system::System; autodiff=:Enzyme, verbosity=0, chunk_sizes=nothing)
        verbosity >= 1 && @info "Preparing model"

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
        ln_prior_transformed(initial_θ_0)


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
            function ℓπcallback(θ_transformed, system=system, arr2nt=arr2nt, Bijector_invlinkvec=Bijector_invlinkvec)
                # Transform back from the unconstrained support to constrained support for the likelihood function
                θ_natural = Bijector_invlinkvec(θ_transformed)
                θ_structured = arr2nt(θ_natural)
                lprior = @inline ln_prior_transformed(θ_natural)
                # CAUTION: This inline annotation is necessary for correct gradients from Enzyme. Yikes!
                llike  = @inline ln_like_generated(system, θ_structured)
                # llike = 0.0
                # llike += Octofitter.ln_like(system.planets.b.observations[1], θ_structured.planets.b, Octofitter.Visual{KepOrbit}(;θ_structured...,θ_structured.planets.b...))
                # llike += Octofitter.ln_like(system.planets.b.observations[2], θ_structured.planets.b, Octofitter.Visual{KepOrbit}(;θ_structured...,θ_structured.planets.b...))
                # llike += Octofitter.ln_like(system.planets.b.observations[3], θ_structured.planets.b, Octofitter.Visual{KepOrbit}(;θ_structured...,θ_structured.planets.b...))
                
                lpost = lprior+llike
                # if !isfinite(lprior)
                #     @error "Invalid log prior encountered. This likely indicates a problem with the prior support."
                #     error()
                # end
                # if !isfinite(llike)
                #     @error "Invalid log likelihood encountered. This likely indicates a problem with the prior support."
                #     error()
                # end
                return lpost
            end

            # Test likelihood function immediately to give user a clean error
            # if it fails for some reason.
            ℓπcallback(initial_θ_0_t,)
            # return (ℓπcallback,
            # initial_θ_0_t,
            # system,
            # arr2nt,
            # Bijector_invlinkvec,
            # )
            if autodiff == :Enzyme
                # Enzyme mode:
                ∇ℓπcallback = let diffresult = copy(initial_θ_0_t), system=system, ℓπcallback=ℓπcallback
                    system_tmp = deepcopy(system)
                    # oh = Enzyme.onehot(diffresult)
                    # forward = function (θ_t)
                    #     primal, out = Enzyme.autodiff(
                    #         Enzyme.Forward,
                    #         ℓπcallback,
                    #         Enzyme.Duplicated,
                    #         Enzyme.BatchDuplicated(θ_t,oh),
                    #         (system),
                    #         (arr2nt),
                    #         (Bijector_invlinkvec)
                    #     )
                    #     return primal, out
                    # end
                    reverse = function (θ_t)
                        fill!(diffresult,0)
                        primal = @inline ℓπcallback(θ_t)
                        @inline Enzyme.autodiff(
                            Enzyme.Reverse,
                            ℓπcallback,
                            Enzyme.Duplicated(θ_t,diffresult),
                            Enzyme.DuplicatedNoNeed(system,system_tmp),
                            Enzyme.Const(arr2nt),
                            Enzyme.Const(Bijector_invlinkvec)
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
                    #     reverse
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
                    Main.Zygote.gradient(ℓπ, θ_t)
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
                        result = ForwardDiff.gradient!(diffresults[Threads.threadid()], ℓπcallback, θ_transformed, cfg, Val{false}())
                        return DiffResults.value(result), DiffResults.gradient(result)
                        # return DiffResults.gradient(result)
                    end
                end
            else
                error("Unsupported option for autodiff: $autodiff. Valid options are :ForwardDiff (default), :Enzyme, and :Zygote.")
            end

            # Display their run time. If something is egregously wrong we'll notice something
            # here in the output logs.
            ℓπcallback(initial_θ_0_t)
            ∇ℓπcallback(initial_θ_0_t) 
            if verbosity >= 1
                (function(ℓπcallback, ∇ℓπcallback, θ)
                    @showtime ℓπcallback(θ)
                    @showtime ∇ℓπcallback(θ)
                end)(ℓπcallback, ∇ℓπcallback, initial_θ_0_t)
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
# TODO: use diff results API to avoid calculating primal twice!
# LogDensityProblems.logdensity_and_gradient(p::LogDensityModel, θ) = (p.ℓπcallback(θ), p.∇ℓπcallback(θ)) # TODO: may need to copy vector
LogDensityProblems.logdensity_and_gradient(p::LogDensityModel, θ) = p.∇ℓπcallback(θ) # TODO: may need to copy vector
LogDensityProblems.dimension(p::LogDensityModel) = p.D
LogDensityProblems.capabilities(::Type{<:LogDensityModel}) = LogDensityProblems.LogDensityOrder{1}()

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::LogDensityModel)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) with fields .ℓπcallback and .∇ℓπcallback")
end

"""
Test that a model returns valid probability densities for 
all input values covered by the priors.
"""
function checkmodel(@nospecialize model)

    (;system,arr2nt,link,invlink,ℓπcallback,∇ℓπcallback,) = model

    priors_kv = Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ])
    priors_flat = map(((k,v),)->v, priors_kv)
    prior_names = map(((k,v),)->k, priors_kv)

    local waserr = false

    θ0 = link(rand.(priors_flat))
    local lp = 0
    try
        lp = ℓπcallback(θ0)
    catch err
        @error "posterior density function failed on simple input" arr2nt(invlink(θ0)) exception=(err, catch_backtrace())
        waserr = true
    end
    if !isfinite(lp)
        @error "posterior density function gave non-finite answer for simple input" lp arr2nt(invlink(θ0))
        waserr = true
    end


    local gradlp = 0
    try
        primal, gradlp = ∇ℓπcallback(θ0)
    catch err
        @error "posterior density gradient failed on simple input" arr2nt(invlink(θ0)) exception=(err, catch_backtrace())
        waserr = true
    end
    if any(!isfinite, gradlp)
        @error "posterior density gradient function gave non-finite answer for simple input" lp arr2nt(invlink(θ0))
        waserr = true
    end

    # Test with extreme values near edge of the prior support
    θ0 = link(quantile.(priors_flat, 1e-5))
    if any(!isfinite, θ0)
        ii = findall(!isfinite, θ0)
        priors_with_errors = namedtuple(prior_names[ii], gradlp[ii])
        @error "lower edge of prior support included a non-finite value" priors_with_errors
    end
    try
        primal, gradlp = ∇ℓπcallback(θ0)
    catch err
        @error "posterior density gradient failed near start of prior support" arr2nt(invlink(θ0)) exception=(err, catch_backtrace()) 
        waserr = true
    end
    if any(!isfinite, gradlp)
        ii = findall(!isfinite, gradlp)
        params_with_errors = namedtuple(prior_names[ii], gradlp[ii])
        @error "posterior density function gave non-finite answer for parameters near start of prior supoprt" arr2nt(invlink(θ0))  params_with_errors
        waserr = true
    end

    # Test with extreme values near edge of the prior support
    θ0 = quantile.(priors_flat, 0.99)
    # This is a bit extreme for eccentricity
    θ0[prior_names .== :e] .= 0.999
    θ0 = link(θ0)
    if any(!isfinite, θ0)
        ii = findall(!isfinite, θ0)
        priors_with_errors = namedtuple(prior_names[ii], gradlp[ii])
        @error "upper edge of prior support included a non-finite value" priors_with_errors
    end
    try
        primal, gradlp = ∇ℓπcallback(θ0)
    catch err
        @error "posterior density gradient failed near end of prior support" arr2nt(invlink(θ0)) exception=(err, catch_backtrace()) 
        waserr = true
    end
    if any(!isfinite, gradlp)
        ii = findall(!isfinite, gradlp)
        params_with_errors = namedtuple(prior_names[ii], gradlp[ii])
        @error "posterior density function gave non-finite answer for parameters near end of prior support"  arr2nt(invlink(θ0))  params_with_errors
        waserr = true
    end

    if waserr
        error("model check failed.")
    end

    
    

    
end


# Fallback when no random number generator is provided (as is usually the case)
function advancedhmc(model::LogDensityModel, target_accept::Number=0.8; kwargs...)
    return advancedhmc(Random.default_rng(), model, target_accept; kwargs...)
end

function octofit(args...; kwargs...)
    return advancedhmc(args...; kwargs...)
end
export octofit


"""
The method signature of Octofitter.hmc is as follows:

    advancedhmc(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel
        target_accept::Number=0.8,
        ensemble::AbstractMCMC.AbstractMCMCEnsemble=MCMCSerial();
        adaptation,
        iterations,
        drop_warmup=true,
        max_depth=12,
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
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    model::LogDensityModel,
    target_accept::Number=0.8;
    adaptation,
    iterations,
    drop_warmup=true,
    max_depth=12,
    initial_samples=250_000,
    initial_parameters=nothing,
    verbosity=2,
)
    if adaptation < 1000
        @warn "At least 1000 steps of adapation are recomended for good sampling"
    end

    # Guess initial starting positions by drawing from priors a bunch of times
    # and picking the best one (highest likelihood).
    # Then transform that guess into our unconstrained support
    if isnothing(initial_parameters)
        verbosity >= 1 && @info "Guessing a starting location by sampling from prior" initial_samples
        initial_θ, mapv = guess_starting_position(rng,model.system,initial_samples)
        verbosity > 2 && @info "Found starting location" θ=stringify_nested_named_tuple(model.arr2nt(initial_θ))
        # Transform from constrained support to unconstrained support
        initial_θ_t = model.link(initial_θ)
    else
        initial_θ = initial_parameters
        # Transform from constrained support to unconstrained support
        initial_θ_t = model.link(initial_θ)
    end

    if verbosity >= 4
        @info "flat starting point" initial_θ
        @info "transformed flat starting point" initial_θ_t
    end


    # # Use Pathfinder to initialize HMC. Works but currently disabled.
    # verbosity >= 1 && @info "Determining initial positions and metric using pathfinder"
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


    # Help adaptation by starting the metric with a rough order of magnitude of the
    # variable variances along the diagonals.

    # We already sampled from the priors earlier to get the starting positon.
    # Use those variance estimates and transform them into the unconstrainted space.
    # variances_t = (model.link(initial_θ .+ sqrt.(variances)/2) .- model.link(initial_θ .- sqrt.(variances)/2)).^2
    p = priors_flat(model.system)
    variances_t = (model.link(quantile.(p, 0.85)) .- model.link(quantile.(p, 0.15))).^2
    # metric = DenseEuclideanMetric(model.D)
    metric = DenseEuclideanMetric(collect(Diagonal(variances_t)))
    if verbosity >= 3
        print("Initial mass matrix M⁻¹ from priors\n")
        display(metric.M⁻¹)
    end
    if any(v->!isfinite(v)||v==0, variances_t)
        error("failed to initialize mass matrix")
    end

    verbosity >= 3 && @info "Creating hamiltonian"
    hamiltonian = Hamiltonian(metric, model)
    verbosity >= 3 && @info "Finding good stepsize"
    ϵ = find_good_stepsize(rng, hamiltonian, initial_θ_t)
    verbosity >= 3 && @info "Found initial stepsize" ϵ 

    # Create integrator. Had best luck with JitteredLeapfrog but all perform similarily.
    integrator = Leapfrog(ϵ)
    # integrator = JitteredLeapfrog(ϵ, 0.05) # 5% normal distribution on step size to help in areas of high curvature. 
    # integrator = TemperedLeapfrog(ϵ, 1.05) 
    
    # This is an attempt at a custom integrator that scales down the step size
    # as eccentricity grows.
    # idx_ecc = only(findall(==(:e), list_parameter_names(model.system)))
    # integrator = EccentricLeapfrog1(ϵ, idx_ecc)


    verbosity >= 3 && @info "Creating kernel"
    kernel = HMCKernel(Trajectory{MultinomialTS}(integrator, GeneralisedNoUTurn(;max_depth)))

    
    verbosity >= 3 && @info "Creating adaptor"
    mma = MassMatrixAdaptor(metric)
    ssa = StepSizeAdaptor(target_accept, integrator)
    adaptor = StanHMCAdaptor(mma, ssa) 






    # If using pathfinder, take N pathfinder draws as initial parameters (disabled)
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
    # if isnothing(initial_parameters)
    #     initial_parameters = fill(initial_θ_t, num_chains)
    # else
    #     if eltype(initial_parameters) <: Number
    #         initial_parameters = fill(initial_parameters, num_chains)
    #     else
    #         # Assume they know what they're doing and are initializing multiple chains separately
    #     end
    #     initial_parameters = map(model.link, initial_parameters)
    # end

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
    num_err_frac = mean(numerical_error)
    mean_tree_depth = mean(actual_tree_depth)
    max_tree_depth_frac = mean(==(max_depth), actual_tree_depth)

    gradient_evaluations = 2sum(getproperty.(stats, :n_steps))

    verbosity >= 1 && println("""
    Sampling report for chain:
    mean_accept         = $mean_accept
    num_err_frac        = $num_err_frac
    mean_tree_depth     = $mean_tree_depth
    max_tree_depth_frac = $max_tree_depth_frac
    gradient_evaluations = $gradient_evaluations\
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
        )
    )
    return mcmcchains_with_info
end

# Helper function for displaying nested named tuples in a compact format.
function stringify_nested_named_tuple(num::Number)
    string(round(num,digits=1))*","
end
function stringify_nested_named_tuple(nt::NamedTuple)
    "(;"*join(map(keys(nt), ", ") do k
        "$k="*stringify_nested_named_tuple(nt[k])
    end)*")"
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
        for (j, val) in enumerate(Iterators.flatten(Iterators.flatten(sample)))
            data[i,j] = val
        end
    end
    c = Chains(data, [string(l) for l in flattened_labels], sectionmap)
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
