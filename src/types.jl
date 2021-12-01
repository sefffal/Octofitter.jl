
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
    if any(<(1), abs.(getproperty.(observations, :ra))) ||
        any(<(1), abs.(getproperty.(observations, :dec)))
        @warn "Confirm that astrometry is entered in milliarcseconds (very small values detected)"
    end
    return Astrometry(
        SVector(getproperty.(observations, :epoch)),
        SVector(getproperty.(observations, :ra)),
        SVector(getproperty.(observations, :dec)),
        SVector(getproperty.(observations, :σ_ra)),
        SVector(getproperty.(observations, :σ_dec)),
    )
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

"""
    Images(...)

A block of images of a system. Pass a vector of named tuples with the following fields:
* band
* image 
* epoch 
* platescale 
* contrast 

For example:
```julia
Images([
    (; epoch=1234.0, band=:J, image=readfits("abc.fits"), platescale=19.4)
])
```
Contrast can be a function that returns the 1 sigma contrast of the image from a separation in mas to the same units as the image file.
Or, simply leave it out and it will be calculated for you.
Epoch is in MJD.
Band is a symbol which matches the one used in the planet's `Priors()` block.
Platescale is in mas/px.
"""
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

"""
    Priors(key=dist, ...)

A group of zero or more priors passed by keyword arguments.
The priors must be univariate distributions from the Distributions.jl 
package.
"""
struct Priors{T}
    priors::T
end
export Priors
# Basically just wrap a named tuple
Priors(;priors...) = Priors{typeof(priors)}(priors)

function Base.show(io::IO, mime::MIME"text/plain", priors::Priors)
    println(io, "Priors:")
    for (k,prior) in zip(keys(priors.priors), values(priors.priors))
        print(io, "\t- ", k, ":\t")
        show(io, mime, prior)
        print(io, "\n")
    end
end

"""
    Derived(key=func, ...)

A group of zero or more functions that resolve to a parameter for the model.
Derived parameters must be functions accepting one argument if applied to a System,
or two arguments if applied to a Planet.
Must always return the same value for the same input (be a pure function) and 
support autodifferentiation.
"""
struct Derived{T}
    variables::T
end
const Deterministic = Derived
export Derived, Deterministic
function Derived(;variables...)
    # Basically just wrap a named tuple
    return Derived(
        NamedTuple(variables),
    )
end
function Base.show(io::IO, mime::MIME"text/plain", det::Derived)
    println(io, "Derived:")
    for k in keys(det.variables)
        print(io, "\t- ", k, "\n")
    end
end

abstract type AbstractPlanet end
"""
    Planet([derived,] priors, [astrometry,], name=:symbol)

A planet (or substellar companion) part of a model.
Must be constructed with a block of priors, and optionally
additional derived parameters and/or astrometry.
`name` must be a symbol, e.g. `:b`.
"""
struct Planet{TD<:Union{Derived,Nothing},TP<:Union{Priors,Nothing},TA<:Union{Astrometry,Nothing}} <: AbstractPlanet
    deterministic::TD
    priors::TP
    astrometry::TA
    name::Symbol
end
export Planet
Planet(det::Derived,priors::Priors,astrometry::Union{Astrometry,Nothing}=nothing; name) = Planet(det,priors, astrometry, name)
Planet(priors::Priors,astrometry::Union{Astrometry,Nothing}=nothing; name) = Planet(nothing,priors, astrometry, name)


astrometry(planet::Planet) = planet.astrometry


"""
    System([derived,] priors, [images,] [propermotionanom,] planets..., name=:symbol)

Construct a model of a system.
Must be constructed with a block of priors, and optionally
additional derived parameters.
You may provide `ProperMotionAnom()` and/or `Images()` of the system.
Finally, planet models are listed last.
`name` must be a symbol e.g. `:HD56441`.
"""
struct System{TDet<:Union{Derived,Nothing}, TPriors<:Priors,TPMA<:Union{ProperMotionAnom,Nothing}, TImages<:Union{Nothing,Images},TModels,TPlanet}
    deterministic::TDet
    priors::TPriors
    propermotionanom::TPMA
    images::TImages
    models::TModels
    planets::TPlanet
    name::Symbol
    function System(
        system_det::Union{Derived,Nothing},
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

# Argument standardization / method cascade.
# Allows users to pass arguments to System in any convenient order.
Trest = Union{ProperMotionAnom,Images,Planet}
System(priors::Priors, args...; kwargs...) =
    System(nothing, priors, args...,; kwargs...)
System(priors::Priors, planets::Planet...; kwargs...) =
    System(nothing, priors, nothing, nothing, planets...,; kwargs...)
System(det::Derived, priors::Priors, planets::Planet...; kwargs...) =
    System(det, priors, nothing, nothing, planets...; kwargs...)
Trest = Union{Images,Planet}
System(det::Union{Derived,Nothing}, priors::Union{Priors,Nothing}, propermotionanom::ProperMotionAnom, planets::Planet...; kwargs...) =
    System(det, priors, propermotionanom, nothing, planets...; kwargs...)
System(det::Union{Derived,Nothing}, priors::Union{Priors,Nothing}, images::Images, planets::Planet...; kwargs...) =
    System(det, priors, nothing, images, planets...; kwargs...)
System(det::Union{Derived,Nothing}, priors::Union{Priors,Nothing}, images::Images, propermotionanom::ProperMotionAnom, planets::Planet...; kwargs...) =
    System(det, priors, propermotionanom, images, planets...; kwargs...)

#### Show methods

Base.show(io::IO, ::MIME"text/plain", astrom::Astrometry) = print(
    io, """
        Astrometry[$(length(astrom.epoch))]
        epoch   \tra\tdec\tσ_ra\tσ_dec
        ────────────────────────────────────────────────
        $(join(["$(round(astrom.epoch[i],digits=2))\t$(round(astrom.ra[i],digits=2))\t$(round(astrom.dec[i],digits=2))\t$(round(astrom.σ_ra[i],digits=2))\t$(round(astrom.σ_dec[i],digits=2))\t"
             for i in eachindex(astrom.epoch)],"\n"))
        ────────────────────────────────────────────────
        """)
Base.show(io::IO, ::MIME"text/plain", pma::ProperMotionAnom) = print(
    io, """
        ProperMotionAnom[$(length(pma.ra_epoch))]
        ra_epoch   \tdec_epoch   \tpm_ra\tpm_dec\tσ_pm_ra\tσ_pm_dec
        ──────────────────────────────────────────────────────────────────
        $(join(["$(round(pma.ra_epoch[i],digits=2)) \t$(round(pma.dec_epoch[i],digits=2))   \t$(round(pma.pm_ra[i],digits=2))\t$(round(pma.pm_dec[i],digits=2))\t$(round(pma.σ_pm_ra[i],digits=2))\t$(round(pma.σ_pm_dec[i],digits=2))\t"
             for i in eachindex(pma.ra_epoch)],"\n"))
        ──────────────────────────────────────────────────────────────────
        """)
Base.show(io::IO, ::MIME"text/plain", is::Images) = print(
    io, """
        Images[$(length(is.image))]
        epoch\tband\tplatescale
        ───────────────────────────
        $(join(["$(round(is.epoch[i],digits=2))\t$(is.band[i])\t$(round(is.platescale[i],digits=2))" for i in eachindex(is.epoch)],"\n"))
        ───────────────────────────
        """)
function Base.show(io::IO, mime::MIME"text/plain", p::AbstractPlanet)
    print(io, "Planet $(p.name)\n")
    if !isnothing(p.deterministic)
        show(io, mime, p.deterministic)
    end
    show(io, mime, p.priors)
    if hasproperty(p, :astrometry) && !isnothing(p.astrometry)
        show(io, mime, p.astrometry)
    end
    print(io, "\n")
    println(io)
end
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





# function reparameterize(planet::Planet)
#     function reparameterized_ln_prior(θ_cv)
#         # τ2pi = θ_cv.Φ + θ_cv.Φω⁻ + θ_cv.ΦΩ⁻
#         return planet.priors.ln_prior(merge(
#             NamedTuple(θ_cv), (;
#                 ω = θ_cv.ωΩ⁺ + θ_cv.ωΩ⁻,
#                 Ω = θ_cv.ωΩ⁺ - θ_cv.ωΩ⁻
#             )
             
            

#             # NamedTuple(θ_cv), (;
#             #     ω = θ_cv.Φ - τ2pi + θ_cv.ΦΩ⁻,
#             #     Ω = θ_cv.Φ - τ2pi + θ_cv.Φω⁻,
#             #     τ = τ2pi/2π,
#             # )

#         ))
#     end
#     priors = Priors(planet.priors.priors, reparameterized_ln_prior)
#     return ReparameterizedPlanet(planet, priors, planet.name)
# end
# export reparameterize