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

const rv_cols = (:epoch, :rv, :σ_rv)
struct RadialVelocity{TTable<:Table} <: AbstractObs
    table::TTable
    function RadialVelocity(observations...)
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Ecpected columns $rv_cols")
        end
        return new{typeof(table)}(table)
    end
end
RadialVelocity(observations::NamedTuple...) = RadialVelocity(observations)
export RadialVelocity

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
struct Priors
    priors::Dict{Symbol,Distribution}
end
export Priors
# Basically just wrap a named tuple
Priors(;priors...) = Priors(Dict(priors))

function Base.show(io::IO, mime::MIME"text/plain", priors::Priors)
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
struct Derived
    variables::Dict{Symbol,Base.Callable}
end
export Derived
Derived(;variables...) = Derived(Dict(variables))
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
additional derived parameters and/or sastrometry.
`name` must be a symbol, e.g. `:b`.
"""
struct Planet{TElem<:AbstractOrbit, TP<:Priors,TD<:Union{Derived,Nothing},TObs<:NTuple{N,<:AbstractObs} where N}
    priors::TP
    derived::TD
    observations::TObs
    name::Symbol
    Planet{O}(priors::Priors, det::Derived, obs::NTuple{N,<:AbstractObs}, name::Symbol) where {O<:AbstractOrbit,N} = new{O,typeof(priors),typeof(det),typeof(obs)}(priors,det,obs,name)
end
export Planet
# Planet(O::Type{<:AbstractOrbit}, (priors,det)::Tuple{Priors,Derived}, args...; kwargs...) = Planet(O, priors, det, args...; kwargs...)
# Planet(O::Type{<:AbstractOrbit}, priors::Priors, obs::AbstractObs...; name) = Planet(O, priors, nothing, obs, name)
# Planet(O::Type{<:AbstractOrbit}, priors::Priors, det::Derived, obs::AbstractObs...; name) = Planet(O, priors, det, obs, name)

Planet{O}((priors,det)::Tuple{Priors,Derived}, args...; kwargs...) where {O<:AbstractOrbit} = Planet{O}(priors, det, args...; kwargs...)
Planet{O}(priors::Priors, obs::AbstractObs...; name) where {O<:AbstractOrbit} = Planet{O}(priors, nothing, obs, name)
Planet{O}(priors::Priors, det::Derived, obs::AbstractObs...; name) where {O<:AbstractOrbit} = Planet{O}(priors, det, obs, name)


"""
    orbittype(planet::Planet)

Returns the type of orbital elements that should be used to represent this planet.

Example:
```julia
julia> d = Planet(KeplerianElements, Variables(), name=:test)
Planet test
julia> DirectDetections.orbittype(d)
KeplerianElements
```
"""
orbittype(::Planet{TElem}) where TElem = TElem

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
    derived::TDet
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
for ObsType in (Astrometry, Photometry, ProperMotionAnom, ProperMotionAnomHGCA, RadialVelocity, Images)
    @eval function Base.show(io::IO, mime::MIME"text/plain", @nospecialize obs::($ObsType))
        print(io, "$($ObsType) ")
        Base.show(io::IO, mime, Table(obs))
    end
end


function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::Planet)
    println(io, "Planet $(p.name)")
    # if !isnothing(p.derived)
    #     show(io, mime, p.derived)
    # end
    # show(io, mime, p.priors)
    # if hasproperty(p, :astrometry) && !isnothing(p.astrometry)
    #     show(io, mime, p.astrometry)
    # end
    # print(io, "\n")
    # println(io)
end
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize sys::System)
    println(io, "System model $(sys.name)")
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
    if !isnothing(system.derived)
        for (key,func) in zip(keys(system.derived.variables), values(system.derived.variables))
            ex = :(
                $key = $func(sys)
            )
            push!(body_sys_determ,ex)
        end
    end

    # Planets: priors & derived variables
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

        if isnothing(planet.derived)
            ex = :(
                $(planet.name) = (;
                    $(body_planet_priors...)
                )
            )
        # Resolve derived vars.
        else
            body_planet_determ = Expr[]
            for (key,func) in zip(keys(planet.derived.variables), values(planet.derived.variables))
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
        # Resolve derived system variables
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





# This is a straight forward implementation that unfortunately is not type stable.
# This is because we are looping over a heterogeneous container
# function make_ln_prior(priors)
#     return function ln_prior(params)
#         lp = zero(first(params))
#         for i in eachindex(params)
#             pd = priors[i]
#             param = params[i]
#             lp += logpdf(pd, param)
#         end
#         return lp 
#     end
# end

function make_ln_prior(system::System)

    # This function uses meta-programming to unroll all the code at compile time.
    # This is a big performance win, since it avoids looping over all the different
    # types of distributions that might be specified as priors.
    # Otherwise we would have to loop through an abstract vector and do runtime dispatch!
    # This way all the code gets inlined into a single tight numberical function in most cases.

    i = 0
    prior_evaluations = Expr[]

    # System priors
    for prior_distribution in values(system.priors.priors)
        i += 1
        ex = :(
            lp += $logpdf($prior_distribution, arr[$i])
        )
        push!(prior_evaluations,ex)
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            i += 1
            # Work around for Beta distributions.
            # Outside of the range [0,1] logpdf returns -Inf.
            # This works fine, but AutoDiff outside this range causes a DomainError.
            if typeof(prior_distribution) <: Beta
                ex = :(
                    lp += 0 <= arr[$i] < 1 ? $logpdf($prior_distribution, arr[$i]) : -Inf
                )
            else
                ex = :(
                    lp += $logpdf($prior_distribution, arr[$i])
                )
            end
            push!(prior_evaluations,ex)
        end
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        lp = zero(first(arr))
        # Add contributions from planet priors
        @inbounds begin
           $(prior_evaluations...) 
        end
        return lp
    end))
end

# Same as above, but assumes the input to the log prior was sampled
# using transformed distributions from Bijectors.jl
# Uses logpdf_with_trans() instead of logpdf to make the necessary corrections.
function make_ln_prior_transformed(system::System)

    i = 0
    prior_evaluations = Expr[]

    # System priors
    for prior_distribution in values(system.priors.priors)
        i += 1
        ex = :(
            lp += $logpdf_with_trans($prior_distribution, arr[$i], true)
        )
        push!(prior_evaluations,ex)
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            i += 1
            ex = :(
                lp += $logpdf_with_trans($prior_distribution, arr[$i], true)
            )
            push!(prior_evaluations,ex)
        end
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        lp = zero(first(arr))
        # Add unrolled prior evaluations
        @inbounds begin
           $(prior_evaluations...) 
        end
        return lp
    end))
end


# # Replaces `θ = Bijectors.invlink.(priors_vec, θ_t)` with a type stable
# # unrolled version.
# function make_Bijector_invlinkvec(priors_vec)

#     i = 0
#     parameter_transformations = Expr[]

#     # System priors
#     for prior_distribution in priors_vec
#         i += 1
#         ex = :(
#             theta_out[$i] = $(Bijectors.invlink)($prior_distribution, arr[$i])
#         )
#         push!(parameter_transformations, ex)
#     end

#     # Here is the function we return.
#     # It maps an array of parameters into our nested named tuple structure
#     # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
#     # The RuntimeGeneratedFunctions package avoids these in all cases.
#     return @RuntimeGeneratedFunction(:(function (arr)
#         l = $i
#         theta_out = @MVector zeros(eltype(arr), l)
#         # theta_out = zeros(eltype(arr), l)
#         @boundscheck if length(arr) != l
#             error("Expected exactly $l elements in array (got $(length(arr)))")
#         end
#         # Add unrolled parameter transformations to fill theta_out
#         @inbounds begin
#            $(parameter_transformations...) 
#         end
#         return theta_out
#     end))
# end


# Replaces `θ = Bijectors.invlink.(priors_vec, θ_t)` with a type stable
# unrolled version.
function make_Bijector_invlinkvec(priors_vec)

    i = 0
    parameter_transformations = Expr[]

    # System priors
    for prior_distribution in priors_vec
        i += 1
        ex = :(
            $(Bijectors.invlink)($prior_distribution, arr[$i])
        )
        push!(parameter_transformations, ex)
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        # theta_out = zeros(eltype(arr), l)
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        # Add unrolled parameter transformations to fill theta_out
        @inbounds begin
            # theta_out = SVector{l,eltype(arr)}(
            # theta_out = MVector{l,eltype(arr)}(
            theta_out = tuple(
                $(parameter_transformations...) 
            )
        end
        return theta_out
    end))
end
