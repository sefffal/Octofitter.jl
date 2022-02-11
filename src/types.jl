
struct Astrometry{N,T<:Number}
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


struct Photometry{N,T<:Number}
    band::SVector{N,Symbol}
    phot::SVector{N,T}
    σ_phot::SVector{N,T}
end
export Photometry


function Photometry(observations::NamedTuple...)
    return Photometry(
        SVector(getproperty.(observations, :band)),
        SVector(getproperty.(observations, :phot)),
        SVector(getproperty.(observations, :σ_phot)),
    )
end


struct ProperMotionAnom{N,T<:Number}
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
Images(
    (; epoch=1234.0, band=:J, image=readfits("abc.fits"), platescale=19.4)
)
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

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize priors::Priors)
    println(io, "Priors:")
    for (k,prior) in zip(keys(priors.priors), values(priors.priors))
        print(io, "  ", k, "\t")
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
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize det::Derived)
    print(io, "Derived:\n  ")
    for k in keys(det.variables)
        print(io, k, ", ")
    end
    print(io, "\n")
end

abstract type AbstractPlanet end
"""
    Planet([derived,] priors, [astrometry,], name=:symbol)

A planet (or substellar companion) part of a model.
Must be constructed with a block of priors, and optionally
additional derived parameters and/or astrometry.
`name` must be a symbol, e.g. `:b`.
"""
struct Planet{TD<:Union{Derived,Nothing},TP<:Union{Priors,Nothing},TA<:Union{Astrometry,Nothing},TPhot<:Union{Photometry,Nothing}} <: AbstractPlanet
    deterministic::TD
    priors::TP
    astrometry::TA
    photometry::TPhot
    name::Symbol
end
export Planet
# There has got to be a better way...
Planet(det::Derived,priors::Priors,astrometry::Union{Astrometry,Nothing}=nothing, photometry::Union{Photometry,Nothing}=nothing; name) = Planet(det,priors, astrometry, photometry, name)
Planet(priors::Priors,astrometry::Union{Astrometry,Nothing}=nothing, photometry::Union{Photometry,Nothing}=nothing; name) = Planet(nothing,priors, astrometry,photometry, name)
Planet(det::Derived,priors::Priors, photometry::Photometry, astrometry::Union{Astrometry,Nothing}=nothing; name) = Planet(det,priors, astrometry, photometry, name)
Planet(priors::Priors, photometry::Photometry, astrometry::Union{Astrometry,Nothing}=nothing; name) = Planet(nothing,priors, astrometry,photometry, name)
Planet(priors::Priors, det::Derived, args...; name) = Planet(det, priors, args...; name)


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
struct System{TDet<:Union{Derived,Nothing}, TPriors<:Priors, TPMA<:Union{ProperMotionAnom,Nothing}, TImages<:Union{Nothing,Images},TPlanet}
    deterministic::TDet
    priors::TPriors
    propermotionanom::TPMA
    images::TImages
    planets::TPlanet
    name::Symbol
    function System(
        system_det::Union{Derived,Nothing},
        system_priors::Union{Priors,Nothing},
        propermotionanom::Union{ProperMotionAnom,Nothing},
        images::Union{Images,Nothing},
        planets::AbstractPlanet...;
        name
    )
        if isempty(planets)
            planets_nt = (;)
        else
            planets_nt = namedtuple(
                getproperty.(planets, :name),
                planets
            )
        end
        return new{typeof(system_det), typeof(system_priors), typeof(propermotionanom), typeof(images), typeof(planets_nt)}(
            system_det, system_priors, propermotionanom, images, planets_nt, name
        )
    end
end
export System

# Argument standardization / method cascade.
# Allows users to pass arguments to System in any convenient order.
System(planets::Planet...; kwargs...) = System(Priors(), planets...; kwargs...)
System(priors::Priors, args...; kwargs...) =
    System(nothing, priors, args...,; kwargs...)
System(priors::Priors, planets::Planet...; kwargs...) =
    System(nothing, priors, nothing, nothing, planets...,; kwargs...)
System(priors::Priors, det::Derived, args...; kwargs...) = System(det, priors, args...; kwargs...)
System(det::Derived, priors::Priors, planets::Planet...; kwargs...) =
    System(det, priors, nothing, nothing, planets...; kwargs...)
System(det::Union{Derived,Nothing}, priors::Union{Priors,Nothing}, propermotionanom::ProperMotionAnom, planets::Planet...; kwargs...) =
    System(det, priors, propermotionanom, nothing, planets...; kwargs...)
System(det::Union{Derived,Nothing}, priors::Union{Priors,Nothing}, images::Images, planets::Planet...; kwargs...) =
    System(det, priors, nothing, images, planets...; kwargs...)
System(det::Union{Derived,Nothing}, priors::Union{Priors,Nothing}, images::Images, propermotionanom::ProperMotionAnom, planets::Planet...; kwargs...) =
    System(det, priors, propermotionanom, images, planets...; kwargs...)

#### Show methods
Base.show(io::IO, ::MIME"text/plain", @nospecialize phot::Photometry) = print(
    io, """
        Astrometry[$(length(phot.band))]
        band\tphot\tσ_phot
        ──────────────────────
        $(join(["$(phot.band[i])\t$(phot.phot[i])\t$(phot.σ_phot[i])"
             for i in eachindex(phot.band)],"\n"))
        ──────────────────────
        """)
Base.show(io::IO, ::MIME"text/plain", @nospecialize astrom::Astrometry) = print(
    io, """
        Astrometry[$(length(astrom.epoch))]
        epoch   \tra\tdec\tσ_ra\tσ_dec
        ────────────────────────────────────────────────
        $(join(["$(round(astrom.epoch[i],digits=2))\t$(round(astrom.ra[i],digits=2))\t$(round(astrom.dec[i],digits=2))\t$(round(astrom.σ_ra[i],digits=2))\t$(round(astrom.σ_dec[i],digits=2))\t"
             for i in eachindex(astrom.epoch)],"\n"))
        ────────────────────────────────────────────────
        """)
Base.show(io::IO, ::MIME"text/plain", @nospecialize pma::ProperMotionAnom) = print(
    io, """
        ProperMotionAnom[$(length(pma.ra_epoch))]
        ra_epoch   \tdec_epoch   \tpm_ra\tpm_dec\tσ_pm_ra\tσ_pm_dec
        ──────────────────────────────────────────────────────────────────
        $(join(["$(round(pma.ra_epoch[i],digits=2)) \t$(round(pma.dec_epoch[i],digits=2))   \t$(round(pma.pm_ra[i],digits=2))\t$(round(pma.pm_dec[i],digits=2))\t$(round(pma.σ_pm_ra[i],digits=2))\t$(round(pma.σ_pm_dec[i],digits=2))\t"
             for i in eachindex(pma.ra_epoch)],"\n"))
        ──────────────────────────────────────────────────────────────────
        """)
Base.show(io::IO, ::MIME"text/plain", @nospecialize is::Images) = print(
    io, """
        Images[$(length(is.image))]
        epoch\tband\tplatescale
        ───────────────────────────
        $(join(["$(round(is.epoch[i],digits=2))\t$(is.band[i])\t$(round(is.platescale[i],digits=2))" for i in eachindex(is.epoch)],"\n"))
        ───────────────────────────
        """)
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::AbstractPlanet)
    print(io, "Planet $(p.name)\n")
    # if !isnothing(p.deterministic)
    #     show(io, mime, p.deterministic)
    # end
    # show(io, mime, p.priors)
    # if hasproperty(p, :astrometry) && !isnothing(p.astrometry)
    #     show(io, mime, p.astrometry)
    # end
    # print(io, "\n")
    # println(io)
end
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize sys::System)
    print(io, "System model $(sys.name)\n")
    # show(io, mime, sys.priors)
    # for planet in sys.planets
    #     show(io, mime, planet)
    # end
    # if !isnothing(sys.propermotionanom)
    #     show(io, mime, sys.propermotionanom)
    # end
    # if !isnothing(sys.images)
    #     show(io, mime, sys.images)
    # end
    # println(io)
end


"""
    _list_priors(system::System)

Given a System, return a flat vector of prior distributions for all parameters in the System and 
any Planets, according to the same order as `make_arr2nt`.
"""
function _list_priors(system::System)

    priors_vec = []
    # System priors
    for prior_distribution in values(system.priors.priors)
        push!(priors_vec,prior_distribution)
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            push!(priors_vec, prior_distribution)
        end
    end

    # narrow the type
    return map(identity, priors_vec)
end


"""
    make_arr2nt(system::System)

Returns a function that maps an array of parameter values (e.g. sampled from priors,
sampled from posterior, dual numbers, etc.) to a nested named tuple structure suitable
for passing into `ln_like`.
In the process, this evaluates all Deterministic variables defined in the model.

Example:
```julia
julia> arr2nt = DirectDetections.make_arr2nt(system);
julia> params = DirectDetections.sample_priors(system)
9-element Vector{Float64}:
  1.581678216418196
 29.019567489624023
  1.805551918844831
  4.527607604709967
 -2.4432825071288837
 10.165695048003027
 -2.8667829383306063
  0.0006589349005787781
 -1.2600092725864576
julia> arr2nt(params)
(M = 1.581678216418196, plx = 29.019567489624023, planets = (B = (τ = 1.805551918844831, ω = 4.527607604709967, i = -2.4432825071288837, Ω = 10.165695048003027, loge = -2.8667829383306063, loga = 0.0006589349005787781, logm = -1.2600092725864576, e = 0.001358992505357737, a = 1.0015184052910453, mass = 0.054952914078000334),))
```
"""
function make_arr2nt(system::System)

    # Roll flat vector of variables e.g. drawn from priors or posterior
    # into a nested NamedTuple structure for e.g. ln_like
    i = 0
    body_sys_priors = Expr[]

    # This function uses meta-programming to unroll all the code at compile time.
    # This is a big performance win, since it avoids looping over all the functions
    # etc. In fact, most user defined e.g. Derived callbacks can get inlined right here.

    # Priors
    for key in keys(system.priors.priors)
        i += 1
        ex = :(
            $key = arr[$i]
        )
        push!(body_sys_priors,ex)
    end

    # Deterministic variables for the system
    body_sys_determ = Expr[]
    if !isnothing(system.deterministic)
        for (key,func) in zip(keys(system.deterministic.variables), values(system.deterministic.variables))
            ex = :(
                $key = $func(sys)
            )
            push!(body_sys_determ,ex)
        end
    end

    # Planets: priors & deterministic variables
    body_planets = Expr[]
    for planet in system.planets
        
        # Priors
        body_planet_priors = Expr[]
        for key in keys(planet.priors.priors)
            i += 1
            ex = :(
                $key = arr[$i]
            )
            push!(body_planet_priors,ex)
        end

        if isnothing(planet.deterministic)
            ex = :(
                $(planet.name) = (;
                    $(body_planet_priors...)
                )
            )
        # Resolve deterministic vars.
        else
            body_planet_determ = Expr[]
            for (key,func) in zip(keys(planet.deterministic.variables), values(planet.deterministic.variables))
                ex = :(
                    $key = $func(sys_res, (;$(body_planet_priors...)))
                )
                push!(body_planet_determ,ex)
            end

            ex = :(
                $(planet.name) = (;
                    $(body_planet_priors...),
                    $(body_planet_determ...)
                )
            )
        end
        push!(body_planets,ex)
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    func = @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array. Got ", length(arr))
        end
        # Expand system variables from priors
        sys = (;$(body_sys_priors...))
        # Resolve deterministic system variables
        sys_res = (;
            sys...,
            $(body_sys_determ...)
        )
        # Get resolved planets
        pln = (;$(body_planets...))
        # Merge planets into resolved system
        sys_res_pln = (;sys_res..., planets=pln)
        return sys_res_pln
    end))

    return func
end
