
struct Astrometry{N,T<:Number}
    # observations::NamedTuple{(:epoch, :ra, :dec, :σ_ra, :σ_dec), NTuple{5, T}}
    epoch::SVector{N,T}
    ra::SVector{N,T}
    dec::SVector{N,T}
    σ_ra::SVector{N,T}
    σ_dec::SVector{N,T}
end
export Astrometry

function Astrometry(observations::NamedTuple...)
    return Astrometry(
        SVector(getproperty.(observations, :epoch)),
        SVector(getproperty.(observations, :ra)),
        SVector(getproperty.(observations, :dec)),
        SVector(getproperty.(observations, :σ_ra)),
        SVector(getproperty.(observations, :σ_dec)),
    )
end
Base.show(io::IO, ::MIME"text/plain", astrom::Astrometry) = print(
    io, """
        Astrometry[$(length(astrom.epoch))]
        epoch   \tra\tdec\tσ_ra\tσ_dec
        ────────────────────────────────────────────────
        $(join(["$(round(astrom.epoch[i],digits=2))\t$(round(astrom.ra[i],digits=2))\t$(round(astrom.dec[i],digits=2))\t$(round(astrom.σ_ra[i],digits=2))\t$(round(astrom.σ_dec[i],digits=2))\t"
             for i in eachindex(astrom.epoch)],"\n"))
        ────────────────────────────────────────────────
        """)
@recipe function f(astrom::Astrometry) where T
    xerror := astrom.σ_ra
    yerror := astrom.σ_dec
    xflip --> true
    return astrom.ra, astrom.dec
end

struct ProperMotionAnom{N,T<:Number}
    # observations::NamedTuple{(:epoch, :ra, :dec, :σ_ra, :σ_dec), NTuple{5, T}}
    ra_epoch::SVector{N,T}
    dec_epoch::SVector{N,T}
    pm_ra::SVector{N,T}
    pm_dec::SVector{N,T}
    σ_pm_ra::SVector{N,T}
    σ_pm_dec::SVector{N,T}
end
export ProperMotionAnom 
function ProperMotionAnom(observations::NamedTuple...)
    T = promote_type(typeof.(values(first(observations)))...)
    return ProperMotionAnom{length(observations),T}(
        SVector(getproperty.(observations, :ra_epoch)),
        SVector(getproperty.(observations, :dec_epoch)),
        SVector(getproperty.(observations, :pm_ra)),
        SVector(getproperty.(observations, :pm_dec)),
        SVector(getproperty.(observations, :σ_pm_ra)),
        SVector(getproperty.(observations, :σ_pm_dec)),
    )
end

Base.show(io::IO, ::MIME"text/plain", pma::ProperMotionAnom) = print(
    io, """
        ProperMotionAnom[$(length(pma.ra_epoch))]
        ra_epoch   \tdec_epoch   \tpm_ra\tpm_dec\tσ_pm_ra\tσ_pm_dec
        ──────────────────────────────────────────────────────────────────
        $(join(["$(round(pma.ra_epoch[i],digits=2)) \t$(round(pma.dec_epoch[i],digits=2))   \t$(round(pma.pm_ra[i],digits=2))\t$(round(pma.pm_dec[i],digits=2))\t$(round(pma.σ_pm_ra[i],digits=2))\t$(round(pma.σ_pm_dec[i],digits=2))\t"
             for i in eachindex(pma.ra_epoch)],"\n"))
        ──────────────────────────────────────────────────────────────────
        """)



struct Images{TImg,TCont}
    band::Vector{Symbol}
    image::Vector{TImg}
    epoch::Vector{Float64}
    platescale::Vector{Float64}
    contrast::Vector{TCont}
end
export Images
function Images(observations::NamedTuple...)
    observations_prepared = map(observations) do obs
        if !hasproperty(obs, :contrast)
            @info "Measuring contrast from image"
            contrast = contrast_interp(obs.image)
            obs = merge(obs, (;contrast))
        end
        return obs
    end
    return Images(
        collect(getproperty.(observations_prepared, :band)),
        collect(getproperty.(observations_prepared, :image)),
        collect(getproperty.(observations_prepared, :epoch)),
        collect(getproperty.(observations_prepared, :platescale)),
        collect(getproperty.(observations_prepared, :contrast)),
    )
end
Base.show(io::IO, ::MIME"text/plain", is::Images) = print(
    io, """
        Images[$(length(is.images))]
        epoch\tband\tplatescale
        ───────────────────────────
        $(join(["$(round(is.epoch[i],digits=2))\t$(is.band[i])\t$(round(is.platescale[i],digits=2))" for i in eachindex(is.epoch)],"\n"))
        ───────────────────────────
        """)


struct Priors{T}
    priors::T
end
export Priors 
# function Priors(;priors...)
#     ln_prior = make_ln_prior(priors)

#     # Compile and test result
#     𝓁prior = ln_prior(rand.(priors_cv))
#     if !isfinite(𝓁prior)
#         error("Test of ln_prior calculation returned $𝓁prior")
#     end
#     return Priors{typeof(priors)}(
#         priors,
#         ln_prior,
#     )
# end
Priors(;priors...) = Priors{typeof(priors)}(priors)

function Base.show(io::IO, mime::MIME"text/plain", priors::Priors)
    println(io, "Priors:")
    for (k,prior) in zip(keys(priors.priors), priors.priors)
        print(io, "\t- ", k, ":\t")
        show(io, mime, prior)
        print(io, "\n")
    end
end

struct Deterministic{T}
    variables::T
end
const Derived = Deterministic
export Deterministic, Derived
function Deterministic(;variables...)
    # Basically just wrap a named tuple
    return Deterministic(
        NamedTuple(variables),
    )
end
function Base.show(io::IO, mime::MIME"text/plain", det::Deterministic)
    println(io, "Deterministic:")
    for k in keys(det.variables)
        print(io, "\t- ", k, "\n")
    end
end

abstract type AbstractPlanet end
struct Planet{TD<:Union{Deterministic,Nothing},TP<:Union{Priors,Nothing},TA<:Union{Astrometry,Nothing}} <: AbstractPlanet
    deterministic::TD
    priors::TP
    astrometry::TA
    name::Symbol
end
export Planet
Planet(det::Deterministic,priors::Priors,astrometry::Union{Astrometry,Nothing}=nothing; name) = Planet(det,priors, astrometry, name)
Planet(priors::Priors,astrometry::Union{Astrometry,Nothing}=nothing; name) = Planet(nothing,priors, astrometry, name)
function Base.show(io::IO, mime::MIME"text/plain", p::AbstractPlanet)
    print(io, "Planet $(p.name)\n")
    show(io, mime, p.deterministic)
    show(io, mime, p.priors)
    if hasproperty(p, :astrometry) && !isnothing(p.astrometry)
        show(io, mime, p.astrometry)
    end
    print(io, "\n")
    println(io)
end

struct ReparameterizedPlanet{TPlanet} <: AbstractPlanet
    planet::TPlanet
    priors::Priors
    name::Symbol
end

astrometry(planet::Planet) = planet.astrometry
astrometry(planet::ReparameterizedPlanet) = planet.planet.astrometry


function reparameterize(planet::Planet)
    function reparameterized_ln_prior(θ_cv)
        # τ2pi = θ_cv.Φ + θ_cv.Φω⁻ + θ_cv.ΦΩ⁻
        return planet.priors.ln_prior(merge(
            NamedTuple(θ_cv), (;
                ω = θ_cv.ωΩ⁺ + θ_cv.ωΩ⁻,
                Ω = θ_cv.ωΩ⁺ - θ_cv.ωΩ⁻
            )
             
            

            # NamedTuple(θ_cv), (;
            #     ω = θ_cv.Φ - τ2pi + θ_cv.ΦΩ⁻,
            #     Ω = θ_cv.Φ - τ2pi + θ_cv.Φω⁻,
            #     τ = τ2pi/2π,
            # )

        ))
    end
    priors = Priors(planet.priors.priors, reparameterized_ln_prior)
    return ReparameterizedPlanet(planet, priors, planet.name)
end
export reparameterize

struct System{TDet<:Union{Deterministic,Nothing}, TPriors<:Priors,TPMA<:Union{ProperMotionAnom,Nothing}, TImages<:Union{Nothing,Images},TModels,TPlanet}
    deterministic::TDet
    priors::TPriors
    propermotionanom::TPMA
    images::TImages
    models::TModels
    planets::TPlanet
    name::Symbol
    function System(
        system_det::Union{Deterministic,Nothing},
        system_priors::Union{Priors,Nothing},
        propermotionanom::Union{ProperMotionAnom,Nothing},
        images::Union{Images,Nothing},
        planets::AbstractPlanet...;
        models=nothing,
        name
    )
        planets_nt = namedtuple(
            getproperty.(planets, :name),
            planets
        )
        return new{typeof(system_det), typeof(system_priors), typeof(propermotionanom), typeof(images), typeof(models),typeof(planets_nt)}(
            system_det, system_priors, propermotionanom, images, models, planets_nt, name
        )
    end
end
export System
# Argument standardization:
Trest = Union{ProperMotionAnom,Images,Planet}
System(priors::Priors, args...; kwargs...) =
    System(nothing, priors, args...,; kwargs...)
System(priors::Priors, planets::Planet...; kwargs...) =
    System(nothing, priors, nothing, nothing, planets...,; kwargs...)
System(det::Deterministic, priors::Priors, planets::Planet...; kwargs...) =
    System(det, priors, nothing, nothing, planets...; kwargs...)


Trest = Union{Images,Planet}
System(det::Union{Deterministic,Nothing}, priors::Union{Priors,Nothing}, propermotionanom::ProperMotionAnom, planets::Planet...; kwargs...) =
    System(det, priors, propermotionanom, nothing, planets...; kwargs...)

System(det::Union{Deterministic,Nothing}, priors::Union{Priors,Nothing}, images::Images, planets::Planet...; kwargs...) =
    System(det, priors, nothing, images, planets...; kwargs...)

System(det::Union{Deterministic,Nothing}, priors::Union{Priors,Nothing}, images::Images, propermotionanom::ProperMotionAnom, planets::Planet...; kwargs...) =
    System(det, priors, propermotionanom, images, planets...; kwargs...)

# System(det::Union{Deterministic,Nothing}, priors::Union{Priors,Nothing}, propermotionanom::Union{ProperMotionAnom,Nothing}, images::Union{Images,Nothing}, planets::AbstractPlanet...; kwargs...)    = (println("D"); true) &&
#     System(det, priors, propermotionanom, images, planets...; kwargs...)



function Base.show(io::IO, mime::MIME"text/plain", sys::System)
    print(io, "System model $(sys.name)\n")
    show(io, mime, sys.priors)
    for planet in sys.planets
        show(io, mime, planet)
    end
    if !isnothing(sys.propermotionanom)
        show(io, mime, sys.propermotionanom)
    end
    if !isnothing(sys.images)
        show(io, mime, sys.images)
    end
    if !isnothing(sys.models)
        show(io, mime, keys(sys.models))
    end
    println(io)
end





# Copied from ModellingToolkit.

export @named 

macro named(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        esc(_named_idxs(name, idxs, :($(gensym()) -> $call)))
    else
        check_name(name)
        esc(:($name = $(_named(name, call))))
    end
end

macro named(name::Symbol, idxs, call)
    esc(_named_idxs(name, idxs, call))
end

function _named(name, call, runtime=false)
    has_kw = false
    call isa Expr || throw(Meta.ParseError("The rhs must be an Expr. Got $call."))
    if length(call.args) >= 2 && call.args[2] isa Expr
        # canonicalize to use `:parameters`
        if call.args[2].head === :kw
            call.args[2] = Expr(:parameters, Expr(:kw, call.args[2].args...))
            has_kw = true
        elseif call.args[2].head === :parameters
            has_kw = true
        end
    end

    if !has_kw
        param = Expr(:parameters)
        if length(call.args) == 1
            push!(call.args, param)
        else
            insert!(call.args, 2, param)
        end
    end

    kws = call.args[2].args

    if !any(kw->(kw isa Symbol ? kw : kw.args[1]) == :name, kws) # don't overwrite `name` kwarg
        pushfirst!(kws, Expr(:kw, :name, runtime ? name : Meta.quot(name)))
    end
    call
end

function _named_idxs(name::Symbol, idxs, call)
    if call.head !== :->
        throw(ArgumentError("Not an anonymous function"))
    end
    if !isa(call.args[1], Symbol)
        throw(ArgumentError("not a single-argument anonymous function"))
    end
    sym, ex = call.args
    ex = Base.Cartesian.poplinenum(ex)
    ex = _named(:(Symbol($(Meta.quot(name)), :_, $sym)), ex, true)
    ex = Base.Cartesian.poplinenum(ex)
    :($name = $map($sym->$ex, $idxs))
end

check_name(name) = name isa Symbol || throw(Meta.ParseError("The lhs must be a symbol (a) or a ref (a[1:10]). Got $name."))





function split_assign(expr)
    if !(expr isa Expr && expr.head === :(=) && expr.args[2].head === :call)
        throw(ArgumentError("expression should be of the form `sys = foo(a, b)`"))
    end
    name, call = expr.args
end