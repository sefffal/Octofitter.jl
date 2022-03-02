abstract type AbstractObs end
TypedTables.Table(obs::AbstractObs) = obs.table

const astrom_cols = (:epoch, :ra, :dec, :σ_ra, :σ_dec)
struct Astrometry{TTable<:Table} <: AbstractObs
    table::TTable
    function Astrometry(observations...)
        table = Table(observations...)
        if !issubset(astrom_cols, Tables.columnnames(table))
            error("Expected columns $astrom_cols")
        end
        return new{typeof(table)}(table)
    end
end
Astrometry(observations::NamedTuple...) = Astrometry(observations)
export Astrometry


const phot_cols = (:band, :phot, :σ_phot)
struct Photometry{TTable<:Table} <: AbstractObs
    table::TTable
    function Photometry(observations...)
        table = Table(observations...)
        if !issubset(phot_cols, Tables.columnnames(table))
            error("Expected columns $phot_cols")
        end
        return new{typeof(table)}(table)
    end
end
Photometry(observations::NamedTuple...) = Photometry(observations)
export Photometry



const pma_cols = (:ra_epoch, :dec_epoch, :dt, :pm_ra, :pm_dec, :σ_pm_ra, :σ_pm_dec)
struct ProperMotionAnom{TTable<:Table} <: AbstractObs
    table::TTable
    function ProperMotionAnom(observations...)
        table = Table(observations...)
        if !issubset(pma_cols, Tables.columnnames(table))
            error("Expected columns $pma_cols")
        end
        return new{typeof(table)}(table)
    end
end
ProperMotionAnom(observations::NamedTuple...) = ProperMotionAnom(observations)
export ProperMotionAnom

struct ProperMotionAnomHGCA{TTable<:Table} <: AbstractObs
    table::TTable
    function ProperMotionAnomHGCA(observations...)
        table = Table(observations...)
        # if !issubset(pma_cols, Tables.columnnames(table))
        #     error("Expected columns $pma_cols")
        # end
        return new{typeof(table)}(table)
    end
end
ProperMotionAnomHGCA(observations::NamedTuple...) = ProperMotionAnomHGCA(observations)
export ProperMotionAnomHGCA

const images_cols = (:band, :image, :epoch, :platescale, :contrast, )

"""
    Images(...)

A block of images of a system. Pass a vector of named tuples with the following fields:
$images_cols

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
struct Images{TTable<:Table} <: AbstractObs
    table::TTable
    function Images(observations...)
        table = Table(observations...)
        # Fallback to calculating contrast automatically
        if !in(:contrast, columnnames(table))
            @info "Measuring contrast from image"
            contrast = contrast_interp(obs.image)
            table = Table(table, contrast=contrast)
        end
        if !issubset(images_cols, columnnames(table))
            error("Expected columns $images_cols")
        end
        return new{typeof(table)}(table)
    end
end
Images(observations::NamedTuple...) = Images(observations)
export Images


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
export Derived
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

"""
A list of variables drawn from prior distributions or functions thereof.

Example:
```julia
    Variables(
        a = Uniform(1, 100),
        e = Uniform(0, 1),
        i = Sine(),
        τ = Uniform(0, 1),
        mass = Uniform(0, 100),
        Ωpω = Uniform(0, 2π),
        Ωmω = Normal(0, π/2),
        Ω = (sys, pl) -> (pl.Ωpω + pl.Ωmω)/2,
        ω = (sys, pl) -> (pl.Ωpω - pl.Ωmω)/2,
    )
```
"""
function Variables(; kwargs...)
    isdist(v) = typeof(v) <: Distribution
    priors = filter(isdist ∘ last, kwargs)
    derived = filter(!(isdist ∘ last), kwargs)
    # Convert anything <: Number to a function that returns that number
    map!(values(derived)) do v
        if typeof(v) <: Number
            Returns(v)
        else
            # Assume it's callable. We check if it's
            # <: Function but this excludes callable objects.
            v
        end
    end
    return Priors(priors), Derived(derived)
end
export Variables


"""
    Planet([derived,] priors, [astrometry,], name=:symbol)

A planet (or substellar companion) part of a model.
Must be constructed with a block of priors, and optionally
additional derived parameters and/or astrometry.
`name` must be a symbol, e.g. `:b`.
"""
struct Planet{TP<:Priors,TD<:Union{Derived,Nothing},TObs<:NTuple{N,<:AbstractObs} where N}
    priors::TP
    deterministic::TD
    observations::TObs
    name::Symbol
end
export Planet
Planet((priors,det)::Tuple{Priors,Derived}, args...; kwargs...) = Planet(priors, det, args...; kwargs...)
Planet(priors::Priors, obs::AbstractObs...; name) = Planet(priors, nothing, obs, name)
Planet(priors::Priors, det::Derived, obs::AbstractObs...; name) = Planet(priors, det, obs, name)


"""
    System([derived,] priors, [images,] [propermotionanom,] planets..., name=:symbol)

Construct a model of a system.
Must be constructed with a block of priors, and optionally
additional derived parameters.
You may provide `ProperMotionAnom()` and/or `Images()` of the system.
Finally, planet models are listed last.
`name` must be a symbol e.g. `:HD56441`.
"""
struct System{TPriors<:Priors, TDet<:Union{Derived,Nothing},TObs<:NTuple{N,AbstractObs} where N, TPlanet}
    priors::TPriors
    deterministic::TDet
    observations::TObs
    planets::TPlanet
    name::Symbol
    function System(
        system_priors::Priors,
        system_det::Union{Derived,Nothing},
        observations::NTuple{N,AbstractObs} where N,
        planets::NTuple{M,Planet} where M;
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
        return new{typeof(system_priors), typeof(system_det), typeof(observations), typeof(planets_nt)}(
            system_priors, system_det, observations, planets_nt, name
        )
    end
end
export System

# Argument standardization / method cascade.
# Allows users to pass arguments to System in any convenient order.
System((priors,det)::Tuple{Priors,Derived}, args...; kwargs...) = System(priors, det, args...; kwargs...)
System(planets::Planet...; kwargs...) = System(Priors(), nothing, planets...; kwargs...)
System(priors::Priors, args::Union{AbstractObs,Planet}...; kwargs...) = System(priors, nothing, args...; kwargs...)
System(priors::Priors, det::Union{Derived,Nothing}, args::Union{AbstractObs,Planet}...; kwargs...) = System(priors, det, group_obs_planets(args)...; kwargs...)

function group_obs_planets(args)
    observations = filter(o->typeof(o) <: AbstractObs, args)
    planets = filter(p->p isa Planet, args)
    return observations, planets
end

## Helpers for accessing the first observation of a certain type in a planet or system
function astrometry(planet::Planet)
    for obs in planet.observations
        if obs isa Astrometry
            return obs
        end
    end
    return nothing
end
export astrometry
function propermotionanom(system::System)
    for obs in system.observations
        if obs isa ProperMotionAnom || obs isa ProperMotionAnomHGCA
            return obs
        end
    end
    return nothing
end
function images(system::System)
    for obs in system.observations
        if obs isa Images
            return obs
        end
    end
    return nothing
end


#### Show methods
for ObsType in (Astrometry, Photometry, ProperMotionAnom, ProperMotionAnomHGCA, Images)
    @eval function Base.show(io::IO, mime::MIME"text/plain", @nospecialize obs::($ObsType))
        print(io, "$($ObsType) ")
        Base.show(io::IO, mime, Table(obs))
    end
end


function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::Planet)
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
