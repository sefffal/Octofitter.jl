using DiffResults
using LinearAlgebra
using Preferences
using Pathfinder
using Transducers
using CovarianceEstimation
export sample_priors

sample_priors(arg::Union{Planet,System,<:LogDensityModel}, args...; kwargs...) = sample_priors(Random.default_rng(), arg, args...; kwargs...)
# Sample priors from system once
function sample_priors(rng::Random.AbstractRNG, system::System)
    priors_flat_sampled = map(((k,v),)->rand(rng, v), Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
    return priors_flat_sampled
end
# Sample priors from system many times
sample_priors(rng::Random.AbstractRNG, system::System, N::Number) = [sample_priors(rng, system) for _ in 1:N]
sample_priors(rng::Random.AbstractRNG, model::LogDensityModel, N::Number) = [model.sample_priors(rng) for _ in 1:N]


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





function guess_starting_position(model::LogDensityModel, args...; kwargs...)
    return guess_starting_position(Random.default_rng(), model, args...; kwargs...)
end
function guess_starting_position(rng::Random.AbstractRNG, model::LogDensityModel, N=500_000)

    # Seed RNGs for each thread
    rngs = [
        Xoshiro(reinterpret(UInt64,rand(rng)))
        for i in 1: Threads.nthreads()
    ]
    bestlogposts = fill(-Inf64, 1:Threads.nthreads())
    bestparams = fill(model.sample_priors(rng), 1:Threads.nthreads())

    Threads.@threads for i in 1:N
        tid = Threads.threadid()
        params = model.sample_priors(rngs[tid])
        params_t = model.link(params)
        logpost = model.ℓπcallback(params)
        if logpost > bestlogposts[tid]
            bestparams[tid] = params
            bestlogposts[tid] = logpost
        end
    end

    I_max = argmax(bestlogposts)

    return bestparams[I_max], bestlogposts[I_max]
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
function construct_elements(::Type{AbsoluteVisual{KepOrbit}}, θ_system, θ_planet)
    return AbsoluteVisual{KepOrbit}(;(;
        θ_system.M,
        θ_system.ref_epoch,
        θ_system.ra,
        θ_system.dec,
        θ_system.rv,
        θ_system.pmra,
        θ_system.pmdec,
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
function construct_elements(::Type{FixedPosition}, θ_system, θ_planet)
    return Visual{CartesianOrbit}(;(;
        θ_planet.x,
        θ_planet.y,
        z = hasproperty(θ_planet, :z) ? θ_planet.z : zero(θ_planet.x)
    )...)
end
function construct_elements(::Type{Visual{<:FixedPosition}}, θ_system, θ_planet)
    return Visual{CartesianOrbit}(;(;
        θ_system.plx,
        θ_planet.x,
        θ_planet.y,
        z = hasproperty(θ_planet, :z) ? θ_planet.z : zero(θ_planet.x)
    )...)
end

"""
    construct_elements(chains, :b, 4)

Given a Chains object, a symbol matching the name of a planet, and an index,
construct a PlanetOrbits.jl orbit object.
"""
function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, i::Union{Integer,CartesianIndex})
    pk = string(planet_key)
    if haskey(chain, :ra) && haskey(chain, :ref_epoch) && haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        o = AbsoluteVisual{KepOrbit}(;(;
            M=chain["M"][i],
            ref_epoch=chain["ref_epoch"][i],
            ra=chain["ra"][i],
            dec=chain["dec"][i],
            rv=chain["rv"][i],
            pmra=chain["pmra"][i],
            pmdec=chain["pmdec"][i],
            plx=chain["plx"][i],
            i=chain[pk*"_i"][i],
            Ω=chain[pk*"_Ω"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            tp=chain[pk*"_tp"][i],
            a=chain[pk*"_a"][i],
        )...)
        return o
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
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
    elseif haskey(chain, :ra) && haskey(chain, :ref_epoch) && haskey(chain, :plx) && haskey(chain, :M) && haskey(chain, :rv)
        return AbsoluteVisual{RadialVelocityOrbit}(;(;
            ref_epoch=chain["ref_epoch"][i],
            ra=chain["ra"][i],
            dec=chain["dec"][i],
            rv=chain["rv"][i],
            pmra=chain["pmra"][i],
            pmdec=chain["pmdec"][i],
            plx=chain["plx"][i],
            M=chain["M"][i],
            ω=chain[pk*"_ω"][i],
            e=chain[pk*"_e"][i],
            tp=chain[pk*"_tp"][i],
            a=chain[pk*"_a"][i],
        )...)
    elseif haskey(chain, :plx) && haskey(chain, :M) && haskey(chain, :rv)
        return Visual{RadialVelocityOrbit}(;(;
            pmra=chain["pmra"][i],
            pmdec=chain["pmdec"][i],
            plx=chain["plx"][i],
            M=chain["M"][i],
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
    if haskey(chain, :ra) && haskey(chain, :ref_epoch) && haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
        Ms=chain["M"]
        ref_epochs=chain["ref_epoch"]
        ras=chain["ra"]
        decs=chain["dec"]
        rvs=chain["rv"]
        pmras=chain["pmra"]
        pmdecs=chain["pmdec"]
        plxs=chain["plx"]
        plxs=chain["plx"]
        is=chain[pk*"_i"]
        Ωs=chain[pk*"_Ω"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        tps=chain[pk*"_tp"]
        as=chain[pk*"_a"]
        return map(ii) do i
            AbsoluteVisual{KepOrbit}(;(;
                M=Ms[i],
                ref_epoch=ref_epochs[i],
                ra=ras[i],
                dec=decs[i],
                rv=rvs[i],
                pmra=pmras[i],
                pmdec=pmdecs[i],
                plx=plxs[i],
                i=is[i],
                Ω=Ωs[i],
                ω=ωs[i],
                e=es[i],
                tp=tps[i],
                a=as[i],
            )...)
        end
    elseif haskey(chain, :plx) && haskey(chain, Symbol(pk*"_i")) && haskey(chain, Symbol(pk*"_Ω"))
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
    elseif haskey(chain, Symbol(pk*"_x")) && haskey(chain,  Symbol(pk*"_y")) && haskey(chain, Symbol("plx"))
        plx=chain["plx"]
         x=chain[ pk*"_x"]
         y=chain[ pk*"_y"]
        pkkey = pk*"_z"
        if haskey(chain, pkkey)
            z=chain[pkkey]
        else
            z = fill(0.0, size(x))
        end
        return map(ii) do i
            Visual{FixedPosition}(;(;
                x = x[i],
                y = y[i],
                z = z[i],
                plx = plx[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_sep")) && haskey(chain,  Symbol(pk*"_pa")) && haskey(chain, Symbol("plx"))
        plx=chain["plx"]
        sep=chain[ pk*"_sep"]
        pa=chain[ pk*"_pa"]
        pkkey = pk*"_z"
        if haskey(chain, pkkey)
            z=chain[pkkey]
        else
            z = fill(0.0, size(sep))
        end
        return map(ii) do i
            Visual{FixedPosition}(;(;
                sep = sep[i],
                pa = pa[i],
                z = z[i],
                plx = plx[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_ra")) && haskey(chain,  Symbol(pk*"_dec")) && haskey(chain, Symbol("plx"))
        plx=chain["plx"]
        ra=chain[ pk*"_ra"]
        dec=chain[ pk*"_dec"]
        pkkey = pk*"_z"
        if haskey(chain, pkkey)
            z=chain[pkkey]
        else
            z = fill(0.0, size(ra))
        end
        return map(ii) do i
            Visual{FixedPosition}(;(;
                ra = ra[i],
                dec = dec[i],
                z = z[i],
                plx = plx[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_x")) && haskey(chain,  Symbol(pk*"_y"))
         x=chain[ pk*"_x"]
         y=chain[ pk*"_y"]
        pkkey = pk*"_z"
        if haskey(chain, pkkey)
            z=chain[pkkey]
        else
            z = fill(0.0, size(x))
        end
        return map(ii) do i
            FixedPosition(;(;
                x = x[i],
                y = y[i],
                z = z[i],
            )...)
        end
    elseif haskey(chain, Symbol(pk*"_sep")) && haskey(chain,  Symbol(pk*"_pa"))
        sep=chain[ pk*"_sep"]
        pa=chain[ pk*"_pa"]
       pkkey = pk*"_z"
       if haskey(chain, pkkey)
           z=chain[pkkey]
       else
           z = fill(0.0, size(sep))
       end
       return map(ii) do i
           FixedPosition(;(;
               sep = sep[i],
               pa = pa[i],
               z = z[i],
           )...)
       end

    elseif haskey(chain, :M) && haskey(chain, :ra) && haskey(chain, :ref_epoch) && haskey(chain, :plx)
        Ms=chain["M"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        tps=chain[pk*"_tp"]
        as=chain[pk*"_a"]
        plx=chain[:plx]
        ref_epochs=chain["ref_epoch"]
        ras=chain["ra"]
        decs=chain["dec"]
        rvs=chain["rv"]
        pmras=chain["pmra"]
        pmdecs=chain["pmdec"]
        return map(ii) do i
            AbsoluteVisual{RadialVelocityOrbit}(;(;
                M=Ms[i],
                ω=ωs[i],
                e=es[i],
                tp=tps[i],
                a=as[i],
                plx=plx[i],
                ref_epoch=ref_epochs[i],
                ra=ras[i],
                dec=decs[i],
                rv=rvs[i],
                pmra=pmras[i],
                pmdec=pmdecs[i],
            )...)
        end

    elseif haskey(chain, :M) && haskey(chain, :plx)
        Ms=chain["M"]
        ωs=chain[pk*"_ω"]
        es=chain[pk*"_e"]
        tps=chain[pk*"_tp"]
        as=chain[pk*"_a"]
        plx=chain[:plx]
        return map(ii) do i
            Visual{RadialVelocityOrbit}(;(;
                M=Ms[i],
                ω=ωs[i],
                e=es[i],
                tp=tps[i],
                a=as[i],
                plx=plx[i],
            )...)
        end
    elseif haskey(chain, :M)
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
# function construct_elements(chain, planet_key::Union{String,Symbol}, ii::AbstractArray{<:Union{Integer,CartesianIndex}})
#     pk = string(planet_key)
#     Ms=chain[:,"M"]
#     plxs=chain[:,"plx"]
#     is=chain[:,pk*"_i"]
#     Ωs=chain[:,pk*"_Ω"]
#     ωs=chain[:,pk*"_ω"]
#     es=chain[:,pk*"_e"]
#     tps=chain[:,pk*"_tp"]
#     as=chain[:,pk*"_a"]
#     return map(ii) do i
#         Visual{KepOrbit}((;
#             M=Ms[i],
#             plx=plxs[i],
#             i=is[i],
#             Ω=Ωs[i],
#             ω=ωs[i],
#             e=es[i],
#             tp=tps[i],
#             a=as[i],
#         ))
#     end
# end
construct_elements(chain::Chains, planet::Planet, args...; kwargs...) = construct_elements(chain, planet.name, args...; kwargs...) 


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
        initial_samples= pathfinder ? 500 : 250_000,  # deprecated
        initial_parameters=nothing, # deprecated
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
    verbosity::Int=2,

    # deprecated options: see model.starting_points and Octofitter.default_initializer!
    pathfinder::Bool=true,
    initial_parameters::Union{Nothing,Vector{Float64}}=nothing,
    initial_samples::Int= pathfinder ? 10_000 : 250_000,
)
    @nospecialize
    if adaptation < 1000
        @warn "At least 1000 steps of adapation are recomended for good sampling"
    end

    if !pathfinder
        error("pathfinder option is deprecated. Use model.starting_points and/or Octofitter.default_initializer!.")
    end
    if !isnothing(initial_parameters)
        error("initial_parameters option is deprecated. Use model.starting_points and/or Octofitter.default_initializer!.")
    end
    if initial_samples != 10_000
        error("initial_samples option is deprecated. Use model.starting_points and/or Octofitter.default_initializer!.")
    end
    
    # inialize if not already done or set by user
    get_starting_point!!(model)

    local metric
    try
        # This can fail, triggering an exception
        S =  (cov(SimpleCovariance(), stack(model.starting_points)'))
        metric = DenseEuclideanMetric(S)
    catch err
        display(err)

        verbosity > 1 && @warn("Falling back to initializing the diagonals with the prior interquartile ranges.")
        # We already sampled from the priors earlier to get the starting positon.
        # Use those variance estimates and transform them into the unconstrainted space.
        # variances_t = (model.link(initial_θ .+ sqrt.(variances)/2) .- model.link(initial_θ .- sqrt.(variances)/2)).^2
        p = _list_priors(model.system)
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
    end

    initial_θ_t = rand(rng, model.starting_points)
    initial_θ = model.invlink(initial_θ_t)

    if verbosity >= 4
        @info "flat starting point" initial_θ
        @info "transformed flat starting point" initial_θ_t
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
    ratio_divergent_transitions = mean(numerical_error)
    mean_tree_depth = mean(actual_tree_depth)
    max_tree_depth_frac = mean(==(max_depth), actual_tree_depth)

    gradient_evaluations = 2sum(getproperty.(stats, :n_steps))

    verbosity >= 1 && println("""
    Sampling report for chain:
    mean_accept                 = $mean_accept
    ratio_divergent_transitions = $ratio_divergent_transitions
    mean_tree_depth             = $mean_tree_depth
    max_tree_depth_frac         = $max_tree_depth_frac
    gradient_evaluations        = $gradient_evaluations\
    """)

    # Report some warnings if sampling did not work well
    if ratio_divergent_transitions == 1.0
        @error "Numerical errors encountered in ALL iterations. Check model and priors."
    elseif ratio_divergent_transitions > 0.05
        @warn "Numerical errors encountered in more than 5% of iterations. Results are likely incorrect. Your model might have high curvature, and could be improved. Otherwise, increasing the target acceptance rate (second argument to `octofit`) might help" ratio_divergent_transitions
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
            samples_transformed=mc_samples
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

"""
    mcmcchain2result(model, chain_in,)

Does the opposite of result2mcmcchain: given a model and a chain, return a vector of named tuples.
"""
function mcmcchain2result(model, chain,)

    # Quickly construct a named tuple template
    θ = sample_priors(model.system)
    nt = model.arr2nt(θ)

    planetkeys = string.(keys(model.system.planets))
        # for pk in planetkeys
            # if startswith(k, pk*"_")

    # These are the columns we expect in the Chains object
    flattened_labels = keys(flatten_named_tuple(nt))
    # These are the labels corresponding to the flattened named tuple without the *planet_key* prepended
    return broadcast(1:size(chain,1),1:size(chain,3)') do i,j
        # Take existing NT and recurse through it. replace elements
        nt_sys = Dict(
            k => el 
            for (k,el) in zip(flattened_labels, Array(chain[i,:,j]))
            if !any(map(pk->startswith(string(k),pk*"_"), planetkeys))
                # this search operation could be majorly spread up by computing a set
                # of valid keys *once*
        )
        nt_planets = map(collect(planetkeys)) do pk
            return Symbol(pk)=>namedtuple(Dict(
                replace(string(k), r"^"*string(pk)*"_" =>"") => el 
                for (k,el) in zip(flattened_labels, Array(chain[i,:,j]))
                if startswith(string(k),pk*"_")
                    # this search operation could be majorly spread up by computing a set
                    # of valid keys *once*
            ))
        end
        return merge(namedtuple(nt_sys), (;planets=namedtuple(nt_planets)))
        nt_sys
    end
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

include("octoquick.jl")