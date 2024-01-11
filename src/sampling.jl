using DiffResults
using LinearAlgebra
using Preferences
using Pathfinder
using CovarianceEstimation
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
        θ_planet.tp,
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
        θ_planet.tp,
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
        θ_planet.tp,
    )...)
end
function construct_elements(::Type{RadialVelocityOrbit}, θ_system, θ_planet)
    return RadialVelocityOrbit(;(;
        θ_system.M,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.tp,
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
            tp=chain[pk*"_tp"][i],
            a=chain[pk*"_a"][i],
        )...)
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*"_A")) && haskey(chain, Symbol(pk*"_B")) && haskey(chain, Symbol(pk*"_G"))&& haskey(chain, Symbol(pk*"_F"))
        return ThieleInnesOrbit(;(;
            M=chain["M"][i],
            plx=chain["plx"][i],
            e=chain[pk*"_e"][i],
            tp=chain[pk*"_tp"][i],
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
            tp=chain[pk*"_tp"][i],
            a=chain[pk*"_a"][i],
        )...)
    elseif haskey(chain, :M) && haskey(chain, :rv)
        return RadialVelocityOrbit(;(;
            M=chain["M"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            tp=chain[pk*"_tp"][i],
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
        tps=chain[pk*"_tp"]
        as=chain[pk*"_a"]
        return map(ii) do i
            Visual{KepOrbit}(;(;
                M=Ms[i],
                plx=plxs[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                tp=tps[i],
                a=as[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        Ms=chain["M"]
        is=chain[pk*"_i"]
        Ωs=chain[pk*"_Ω"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        tps=chain[pk*"_tp"]
        as=chain[pk*"_a"]
        return map(ii) do i
            KepOrbit(;(;
                M=Ms[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                tp=tps[i],
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
        tps=chain[pk*"_tp"]
        return map(ii) do i
            ThieleInnesOrbit(;(;
                M=Ms[i],
                plx=plxs[i],
                e=es[i],
                tp=tps[i],
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
        pkkey = pk*"_tref"
        if haskey(chain, pkkey)
            tref=chain[pkkey]
        else
            tref = fill(0.0, size(x))
        end
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
        tps=chain[pk*"_tp"]
        as=chain[pk*"_a"]
        return map(ii) do i
            RadialVelocityOrbit(;(;
                M=Ms[i],
                ω=ωs[i],
                e=es[i],
                tp=tps[i],
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
    tps=chain[:,pk*"_tp"]
    as=chain[:,pk*"_a"]
    return map(ii) do i
        Visual{KepOrbit}((;
            M=Ms[i],
            plx=plxs[i],
            i=is[i],
            Ω=Ωs[i],
            ω=ωs[i],
            e=es[i],
            tp=tps[i],
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
            ℓπcallback(initial_θ_0_t)
            # Also Display their run time. If something is egregously wrong we'll notice something
            # here in the output logs.
            if verbosity >= 1
                (function(ℓπcallback, θ)
                    @showtime ℓπcallback(θ)
                end)(ℓπcallback, initial_θ_0_t)
            end

            if autodiff == :Enzyme
                # Enzyme mode:
                ∇ℓπcallback = let diffresult = copy(initial_θ_0_t), system=system, ℓπcallback=ℓπcallback
                    oh = Enzyme.onehot(diffresult)
                    forward = function (θ_t)
                        primal, out = Enzyme.autodiff(
                            Enzyme.Forward,
                            ℓπcallback,
                            Enzyme.BatchDuplicated,
                            Enzyme.BatchDuplicated(θ_t,oh),
                            Enzyme.Const(system),
                            Enzyme.Const(arr2nt),
                            Enzyme.Const(Bijector_invlinkvec)
                        )
                        return primal, tuple(out...)
                    end
                    reverse = function (θ_t)
                        fill!(diffresult,0)
                        out, primal = Enzyme.autodiff(
                            # Enzyme.Reverse,
                            Enzyme.ReverseWithPrimal,
                            ℓπcallback,
                            Enzyme.Duplicated(θ_t,diffresult),
                            Enzyme.Const(system),
                            Enzyme.Const(arr2nt),
                            Enzyme.Const(Bijector_invlinkvec)
                        )
                        return primal, diffresult
                    end 
                    tforward = minimum(
                        @elapsed forward(initial_θ_0_t)
                        for _ in 1:100
                    )
                    treverse = minimum(
                        @elapsed reverse(initial_θ_0_t)
                        for _ in 1:100
                    )
                    if tforward > treverse
                        verbosity > 2  && @info "selected reverse mode AD" tforward treverse
                        reverse
                    else
                        verbosity > 2  && @info "selected forward mode AD" tforward treverse
                        forward
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
                        result = ForwardDiff.gradient!(diffresults[Threads.threadid()], ℓπcallback, θ_transformed, cfg, Val{false}())
                        return DiffResults.value(result), DiffResults.gradient(result)
                        # return DiffResults.gradient(result)
                    end
                end
            else
                error("Unsupported option for autodiff: $autodiff. Valid options are :ForwardDiff (default), :Enzyme, and :Zygote.")
            end


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
Base.@nospecializeinfer function advancedhmc(model::LogDensityModel, target_accept::Number=0.8; kwargs...)
    return advancedhmc(Random.default_rng(), model, target_accept; kwargs...)
end

Base.@nospecializeinfer function octofit(args...; kwargs...)
    return advancedhmc(args...; kwargs...)
end
export octofit

# Define some wrapper functions that hide type information
# so that we don't have to recompile pathfinder() with each 
# new model.
# It's worth it for the sampler, but not pathfinder.
struct CallableAny
    func::Function
end
(ca::CallableAny)(args...;kwargs...) = ca.func(args...;kwargs...)

struct LogDensityModelAny
    ldm::LogDensityModel
end
LogDensityProblems.logdensity(ldm_any::LogDensityModelAny, θ) = LogDensityProblems.logdensity(ldm_any.ldm, θ)
LogDensityProblems.logdensity_and_gradient(ldm_any::LogDensityModelAny, θ) = LogDensityProblems.logdensity_and_gradient(ldm_any.ldm, θ)
LogDensityProblems.dimension(ldm_any::LogDensityModelAny) = LogDensityProblems.dimension(ldm_any.ldm)
LogDensityProblems.capabilities(::Type{<:LogDensityModelAny}) = LogDensityProblems.LogDensityOrder{1}()



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
        initial_samples= pathfinder ? 500 : 250_000,
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
Base.@nospecializeinfer function advancedhmc(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    model::LogDensityModel,
    target_accept::Number=0.8;
    adaptation::Int=1000,
    iterations::Int=1000,
    drop_warmup::Bool=true,
    max_depth::Int=12,
    initial_parameters::Union{Nothing,Vector{Float64}}=nothing,
    verbosity::Int=2,
    pathfinder::Bool=true,
    initial_samples::Int= pathfinder ? 10_000 : 250_000,
)
    @nospecialize
    if adaptation < 1000
        @warn "At least 1000 steps of adapation are recomended for good sampling"
    end



    # Use Pathfinder to initialize HMC. Works but currently disabled.
    verbosity >= 1 && @info "Determining initial positions and metric using pathfinder"
    # It seems to hit a PosDefException sometimes when factoring a matrix.
    # When that happens, the next try usually succeeds.
    # start_time = time()




    local result_pf = nothing
    ldm_any = LogDensityModelAny(model)
    if pathfinder 
        for i in 1:8
            try
                if !isnothing(initial_parameters)
                    initial_θ = initial_parameters
                    # Transform from constrained support to unconstrained support
                    initial_θ_t = model.link(initial_θ)
                    errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
                    result_pf = with_logger(errlogger) do 
                        Pathfinder.multipathfinder(
                            ldm_any, 1000;
                            nruns=1,
                            init=collect(initial_θ_t),
                            progress=verbosity > 1,
                            maxiters=25_000,
                            # maxtime=25.0,
                            # reltol=1e-4,
                        )
                    end
                else
                    init_sampler = function(rng, x) 
                        initial_θ, mapv = guess_starting_position(rng,model.system,initial_samples)
                        initial_θ_t = model.link(initial_θ)
                        x .= initial_θ_t
                    end
                    errlogger = ConsoleLogger(stderr, verbosity >=3 ? Logging.Info : Logging.Error)
                    result_pf = with_logger(errlogger) do 
                        Pathfinder.multipathfinder(
                            ldm_any, 1000;
                            nruns=1,
                            init_sampler=CallableAny(init_sampler),
                            progress=verbosity > 1,
                            maxiters=25_000,
                            # maxtime=25.0,
                            # reltol=1e-4,
                        ) 
                    end
                end
                # Check pareto shape diagnostic to see if we have a good approximation
                # If we don't, just try again
                if result_pf.psis_result.pareto_shape > 3
                    verbosity > 3 && display(result_pf)
                    verbosity >= 4 && display(result_pf.psis_result)
                    verbosity > 2 && @warn "Restarting pathfinder" i
                    continue
                end
                S =  (cov(BiweightMidcovariance(), stack(result_pf.draws)'))
                cholesky(Symmetric(S)) # test can factor without posdef exception
                verbosity >= 3 && "Pathfinder complete"
                break
            catch ex
                result_pf = nothing
                if ex isa PosDefException
                    verbosity > 2 && @warn "Mass matrix failed to factorize. Restarting pathfinder" i
                    continue
                end
                if ex isa InterruptException
                    rethrow(ex)
                end
                @error "Unexpected error occured running pathfinder" exception=(ex, catch_backtrace())
                break
            end
        end
    end
    if isnothing(result_pf) # failed or not attempted
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

        verbosity > 1 && @warn("Falling back to initializing the diagonals with the prior interquartile ranges.")
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

    else # !isnothing(result_pf)
        # Start using a draw from the typical set as estimated by Pathfinder
        if result_pf.psis_result.pareto_shape < 3
            verbosity >= 4 && @info "PSIS result good; starting with sample from typical set"
            initial_θ_t = collect(last(eachcol(result_pf.draws))) # result_pf.fit_distribution.μ
        else
            verbosity >= 4 && @info "PSIS result bad; starting with distribution mean"
            initial_θ_t = collect(mean(result_pf.fit_distribution))
        end
        initial_θ = model.invlink(initial_θ_t)
        verbosity >= 3 && @info "Creating metric"
        # Use the metric found by Pathfinder for HMC sampling
        # metric = Pathfinder.RankUpdateEuclideanMetric(result_pf.pathfinder_results[1].fit_distribution.Σ)
        # display(metric)
        if length(result_pf.pathfinder_results) > 3
            S =  (cov(SimpleCovariance(), stack(result_pf.draws)'))
            metric = DenseEuclideanMetric(S)
        else
            # estimate it emperically by pooling draws from all pathfinder runs
            metric = DenseEuclideanMetric(collect(Matrix(result_pf.pathfinder_results[1].fit_distribution.Σ)))
        end

    end


    if verbosity >= 4
        @info "flat starting point" initial_θ
        @info "transformed flat starting point" initial_θ_t
    end
    

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



    # # Help adaptation by starting the metric with a rough order of magnitude of the
    # # variable variances along the diagonals.

    # # We already sampled from the priors earlier to get the starting positon.
    # # Use those variance estimates and transform them into the unconstrainted space.
    # # variances_t = (model.link(initial_θ .+ sqrt.(variances)/2) .- model.link(initial_θ .- sqrt.(variances)/2)).^2
    # p = priors_flat(model.system)
    # variances_t = (model.link(quantile.(p, 0.85)) .- model.link(quantile.(p, 0.15))).^2
    # # metric = DenseEuclideanMetric(model.D)
    # metric = DenseEuclideanMetric(collect(Diagonal(variances_t)))

    # if verbosity >= 3
    #     print("Initial mass matrix M⁻¹ from priors\n")
    #     display(metric.M⁻¹)
    # end
    # if any(v->!isfinite(v)||v==0, variances_t)
    #     error("failed to initialize mass matrix")
    # end

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
    # if isnothing(result_pf)
    # if metric isa Pathfinder.RankUpdateEuclideanMetric
    #     adaptor = StepSizeAdaptor(target_accept, integrator)
    # else
        mma = MassMatrixAdaptor(metric)
        ssa = StepSizeAdaptor(target_accept, integrator)
        adaptor = StanHMCAdaptor(mma, ssa) 
    # end
    # end







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
    elseif num_err_frac > 0.05
        @warn "Numerical errors encountered in more than 5% of iterations. Results are likely incorrect." num_err_frac
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
function stringify_nested_named_tuple(val::Any) # fallback
    string(val)
end
function stringify_nested_named_tuple(num::Number)
    string(round(num,digits=1))
end
function stringify_nested_named_tuple(nt::NamedTuple)
    "(;"*join(map(keys(nt)) do k
        "$k="*stringify_nested_named_tuple(nt[k])
    end, ", ")*")"
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
