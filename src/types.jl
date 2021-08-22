
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
        $(typeof(astrom))
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
        $(typeof(pma))
        ra_epoch   \tdec_epoch   \tpm_ra\tpm_dec\tσ_pm_ra\tσ_pm_dec
        ──────────────────────────────────────────────────────────────────
        $(join(["$(round(pma.ra_epoch[i],digits=2)) \t$(round(pma.dec_epoch[i],digits=2))   \t$(round(pma.pm_ra[i],digits=2))\t$(round(pma.pm_dec[i],digits=2))\t$(round(pma.σ_pm_ra[i],digits=2))\t$(round(pma.σ_pm_dec[i],digits=2))\t"
             for i in eachindex(pma.ra_epoch)],"\n"))
        ──────────────────────────────────────────────────────────────────
        """)



struct Images
    band::Vector{Symbol}
    image::Vector
    epoch::Vector{Float64}
    platescale::Vector{Float64}
    contrast::Vector
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
        $(typeof(is))
        epoch\tband\tplatescale
        ───────────────────────────
        $(join(["$(round(is.epoch[i],digits=2))\t$(is.band[i])\t$(round(is.platescale[i],digits=2))" for i in eachindex(is.epoch)],"\n"))
        ───────────────────────────
        """)


struct Priors{N}
    priors::ComponentVector#TODO
    ln_prior::Function
end
export Priors 
function Priors(;priors...)
    priors_cv = ComponentVector(priors)
    ln_prior = make_ln_prior(priors_cv)

    # Compile and test result
    𝓁prior = ln_prior(mean.(priors_cv))
    if !isfinite(𝓁prior)
        error("Test of ln_prior calculation returned $𝓁prior")
    end
    return Priors{length(priors)}(
        priors_cv,
        ln_prior
    )
end

function Base.show(io::IO, mime::MIME"text/plain", priors::Priors)
    print(io, "Priors:")
    for (k,prior) in zip(keys(priors.priors), priors.priors)
        print(io, "\t- ", k, ":\t")
        show(io, mime, prior)
        print(io, "\n")
    end
end

abstract type AbstractPlanet{T} end
struct Planet{T} <: AbstractPlanet{T}
    priors::Priors
    astrometry::T
end
export Planet
Planet(priors::Priors) = Planet(priors, nothing)
function Base.show(io::IO, mime::MIME"text/plain", p::AbstractPlanet{T}) where T
    print(io, typeof(p), " model")
    if T == Nothing
        print(io, " with no associated astrometry")
    else
        print(io, " with associated astrometry")
    end
    print(io, "\n")
    show(io, mime, p.priors)
    println(io)
end

struct ReparameterizedPlanet3{T} <: AbstractPlanet{T}
    planet::Planet{T}
    priors::Priors
end
astrometry(planet::Planet) = planet.astrometry
astrometry(planet::ReparameterizedPlanet3) = planet.planet.astrometry


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
    priors = Priors{length(planet.priors.priors)}(planet.priors.priors, reparameterized_ln_prior)
    return ReparameterizedPlanet3(planet, priors)
end
export reparameterize

struct System{TPMA<:Union{ProperMotionAnom,Nothing}, TImages<:Union{Nothing,Images},TPlanet}
    priors::Priors
    propermotionanom::TPMA
    images::TImages
    planets::TPlanet
    function System(system_priors, propermotionanom, images, planets::AbstractPlanet...)
        if length(planets) > 1 && !all(==(keys(first(planets).priors.priors)), [keys(planet.priors.priors) for planet in planets])
            error("All planets in the system model must have priors for the same properties defined")
        end
        return new{ typeof(propermotionanom), typeof(images), typeof(planets)}(
            system_priors, propermotionanom, images, planets
        )
    end
end
export System
System(system_priors::Priors, images::Images, propermotionanom::ProperMotionAnom, planets::AbstractPlanet...) = 
    System(
        system_priors, propermotionanom, images, planets...
    )
System(system_priors::Priors, propermotionanom::ProperMotionAnom, planets::AbstractPlanet...) =
    System(
        system_priors, propermotionanom, nothing, planets...
    )
System(system_priors::Priors, images::Images, planets::AbstractPlanet...) = 
    System(
        system_priors, nothing, images, planets...
    )
System(system_priors::Priors, planets::AbstractPlanet...) = 
    System(system_priors, nothing, nothing, planets...)

function Base.show(io::IO, mime::MIME"text/plain", sys::System)
    # print(io, "System model with $(length(sys.planets)) planets:\n")
    print(io, "System model:\n")
    show(io, mime, sys.priors)
    print(io, "with $(length(sys.planets)) planets:\n")
    for planet in sys.planets
        print(io, "- ")
        show(io, mime, planet)
    end
    if !isnothing(sys.propermotionanom)
        print(io, "- associated proper motion anomaly:")
        show(io, mime, sys.propermotionanom)
    end
    if !isnothing(sys.images)
        print(io, "- associated images\n")
        show(io, mime, sys.images)
    end

    # if T == Nothing
    #     print(io, " with no associated astrometry")
    # else
    #     print(io, " with associated astrometry")
    # end
    println(io)
end









