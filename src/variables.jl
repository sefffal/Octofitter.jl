abstract type AbstractLikelihood end
TypedTables.Table(like::AbstractLikelihood) = like.table

"""
    ln_like(likeobj::AbstractLikelihood,  θ_planet, orbit, orbit_solutions, i_orbsol_start)

Compute the natural log of the likelihood of the data given by `likeobj` being generated
given the model parameters `θ_planet`. 

The remaining parameters are pre-calculated cached values generated `θ_planet`.
"""
function ln_like end

"""
    likeobj_from_epoch_subset(likeobj::AbstractLikelihood,  observation_indices)

Given a likelihood object that wraps an observation table, construct a new one only containing
observations specified by the index or indices in `observation_indices`.

This allows sub-setting data for various statistical checks.
"""
function likeobj_from_epoch_subset end



"""
    Priors(key=dist, ...)

A group of zero or more priors passed by keyword arguments.
The priors must be univariate distributions from the Distributions.jl 
package.
"""
struct Priors
    priors::OrderedDict{Symbol,Distribution}
end
# Basically just wrap a named tuple
Priors(;priors...) = Priors(OrderedDict(priors))

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
    variables::OrderedDict{Symbol,Base.Callable}
end
Derived(;variables...) = Derived(OrderedDict(variables))
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
        M = 1.0,
        i = Sine(),
        τ = UniformCircular(1.0),
        mass = Uniform(0, 100),
        Ωpω = Uniform(0, 2π),
        Ωmω = Normal(0, π/2),
        Ω = (sys, pl) -> (pl.Ωpω + pl.Ωmω)/2,
        ω = (sys, pl) -> (pl.Ωpω - pl.Ωmω)/2,

    )
```
"""
function Variables(; kwargs...)
    # Start by expaning parameterizations.
    # Users can make anything they want using functions, but we
    # have a nice mechanism for pre-canned parameterizations
    # that can introduce auxillary variables without much fuss.
    kwargs_expanded = Pair[]
    extra_obs_likelihoods = []
    for (k,v) in pairs(kwargs)
        expanded, extra_obs_likelihood = expandparam(k,v)
        # Each input can expand to any number of inputs.
        for (k,v) in pairs(expanded)
            push!(kwargs_expanded, k => v)
        end
        push!(extra_obs_likelihoods, extra_obs_likelihood)
    end
    kwargs_dict = OrderedDict(kwargs_expanded)
    # Now divide the list into priors (random variables) and derived (functions of random variables)
    isdist(v) = typeof(v) <: Distribution
    priors = filter(isdist ∘ last, kwargs_dict)
    derived = filter(!(isdist ∘ last), kwargs_dict)
    observation_likelihoods = filter(!isnothing, extra_obs_likelihoods)
    
    return (Priors(priors), Derived(derived), observation_likelihoods...)
end



abstract type Parameterization end
"""
    UniformCircular(domain=2π)

Creates a variable parameterized on a continuous, periodic domain.
Creates two normally diastributed random variables :Ωx and :Ωy 
and maps them to a circular domain using a derived variable 
according to `atan(:Ωy, :Ωx)`.
This may be more efficient than using e.g. `Ω = Uniform(-pi, pi)`
because the sampler can wrap around freely.

Example:
```julia
Variables(;
    a = Uniform(1, 10),
    e = Uniform(0, 1),
    UniformCircular(:Ω)...,
)
```
"""
struct UniformCircular <: Parameterization
    domain::Float64
end
UniformCircular() = UniformCircular(2π)
export UniformCircular
expandparam(var, n::Number) = OrderedDict(var => Returns(n)), nothing
expandparam(var, f::Base.Callable) = OrderedDict(var => f), nothing
expandparam(var, d::Distribution) = OrderedDict(var => d), nothing


function expandparam(var, p::UniformCircular)

    varx = Symbol("$(var)y")
    vary = Symbol("$(var)x")

    callback_inner = @RuntimeGeneratedFunction(:(
        # This parameterization needs to work for either a planet
        # or a system as a whole.
        function (body)
            atan(body.$vary, body.$varx)/2pi*$(p.domain)
        end
    ))
    callback(sys,pl) = callback_inner(pl)
    callback(sys) = callback_inner(sys)

    # We need to create a "prior" on the length of the unit vector so that it doesn't get pinched at (0,0)
    # This has no observable effect on the results of the model as whole, it just improves sampling.
    prior = UnitLengthPrior(varx, vary)

    return OrderedDict(
        varx => Normal(0,1),
        vary => Normal(0,1),
        var  => callback,
    ), prior
end

# This likelihood can be attached to put a prior on the two variables repsenting a uniform circular distribtion.
# It makes sure they keep away from the origin for better sampling.
struct UnitLengthPrior{X,Y} <: AbstractLikelihood where {X,Y}
    UnitLengthPrior(xsymbol, ysymbol) = new{xsymbol, ysymbol}()
end
_isprior(::UnitLengthPrior) = true
function likeobj_from_epoch_subset(obs::UnitLengthPrior{X,Y}, obs_inds) where {X,Y}
    return UnitLengthPrior(X,Y)
end
TypedTables.Table(like::UnitLengthPrior) = nothing

function ln_like(::UnitLengthPrior{X,Y}, θ_system::NamedTuple, _args...) where {X,Y}
    x = getproperty(θ_system, X)
    y = getproperty(θ_system, Y)
    vector_length = sqrt(x^2 + y^2)
    return logpdf(LogNormal(log(1.0), 0.1), vector_length);
end
function ln_like(::UnitLengthPrior{X,Y}, θ_system::NamedTuple, θ_planet::NamedTuple, _args...) where {X,Y}
    x = getproperty(θ_planet, X)
    y = getproperty(θ_planet, Y)
    vector_length = sqrt(x^2 + y^2)
    return logpdf(LogNormal(log(1.0), 0.1), vector_length);
end
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize like::UnitLengthPrior{X,Y}) where {X,Y}
    T = typeof(like)
    println(io, "$(T): √($X^2+$Y^2) ~ LogNormal(log(1), 0.02)")
end
generate_from_params(like::UnitLengthPrior, θ_planet, orbit) = like



"""
    Planet([derived,] priors, [astrometry,], name=:symbol)

A planet (or substellar companion) part of a model.
Must be constructed with a block of priors, and optionally
additional derived parameters and/or sastrometry.
`name` must be a symbol, e.g. `:b`.
"""
struct Planet{TElem<:AbstractOrbit, TP<:Priors,TD<:Union{Derived,Nothing},TObs<:Tuple}
    priors::TP
    derived::TD
    observations::TObs
    name::Symbol
    Planet{O}(priors::Priors, det::Derived, like::Tuple; name::Symbol) where {O<:AbstractOrbit} = new{O,typeof(priors),typeof(det),typeof(like)}(priors,det,like,name)
end
export Planet
Planet{O}((priors,det)::Tuple{Priors,Derived}, args...; kwargs...) where {O<:AbstractOrbit} = Planet{O}(priors, det, args...; kwargs...)
Planet{O}(priors::Priors, like::AbstractLikelihood...; name) where {O<:AbstractOrbit} = Planet{O}(priors, nothing, like; name)
Planet{O}(priors::Priors, det::Derived, like::AbstractLikelihood...; name) where {O<:AbstractOrbit} = Planet{O}(priors, det, like; name)


"""
    orbittype(planet::Planet)

Returns the type of orbital elements that should be used to represent this planet.

Example:
```julia
julia> d = Planet(VisualElements, Variables(), name=:test)
Planet test
julia> Octofitter.orbittype(d)
VisualElements
```
"""
orbittype(::Planet{TElem}) where TElem = TElem

"""
    System([derived,] priors, [images,] [propermotionanom,] planets..., name=:symbol)

Construct a model of a system.
Must be constructed with a block of priors, and optionally
additional derived parameters.
You may provide `ProperMotionAnomLikelihood()` and/or `Images()` of the system.
Finally, planet models are listed last.
`name` must be a symbol e.g. `:HD56441`.
"""
struct System{TPriors<:Priors, TDet<:Union{Derived,Nothing},TObs<:NTuple{N,AbstractLikelihood} where N, TPlanet}
    priors::TPriors
    derived::TDet
    observations::TObs
    planets::TPlanet
    name::Symbol
    function System(
        system_priors::Priors,
        system_det::Union{Derived,Nothing},
        observations::NTuple{N,AbstractLikelihood} where N,
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
System(priors::Priors, args::Union{AbstractLikelihood,Planet}...; kwargs...) = System(priors, nothing, args...; kwargs...)
System(priors::Priors, det::Union{Derived,Nothing}, args::Union{AbstractLikelihood,Planet}...; kwargs...) = System(priors, det, group_obs_planets(args)...; kwargs...)

_isprior(::AbstractLikelihood) = false

function _count_likeobj(system::System)::Int
    likeobj_count = 0
    for obj in [system; system.planets...]
        for obs in obj.observations
            if !_isprior(obs)
                likeobj_count += 1
            end
        end
    end
    return likeobj_count
end

function _count_epochs(system::System)::Int
    observation_count = 0
    for obj in [system; system.planets...]
        for obs in obj.observations
            if hasproperty(obs, :table)
                observation_count += Tables.rowcount(obs.table)
            else
                observation_count += 1
            end
        end
    end
    return observation_count
end

# Function to give the parameter names as a flat vector of symbols. Only returns
# active parameters (i.e.) and not any derived variables.
function list_parameter_names(system::System)
    return map(((k,v),)->k, Iterators.flatten([
        system.priors.priors,
        [planet.priors.priors for planet in system.planets]...
    ]))
end


function group_obs_planets(args)
    observations = filter(o->typeof(o) <: AbstractLikelihood, args)
    planets = filter(p->p isa Planet, args)
    return observations, planets
end

#### Show methods
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize like::AbstractLikelihood)
    ObsType = typeof(like)
    liketype_name = string(ObsType)
    if occursin("{",liketype_name)
        liketype_name = liketype_name[1:findfirst(==('{'),collect(liketype_name))-1]
    end
    print(io, "$(liketype_name) ")
    if hasproperty(like, :table)
        Base.show(io::IO, mime, Table(like))
    else
        println(io)
    end
end


function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::Planet)
    println(io, "Planet $(p.name)")
    if !isnothing(p.derived)
        show(io, mime, p.derived)
    end
    show(io, mime, p.priors)
    for like in p.observations
        show(io, mime, like)
    end
    print(io, "\n")
    println(io)
end
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize sys::System)
    println(io, "System model $(sys.name)")
    if !isnothing(sys.derived)
        show(io, mime, sys.derived)
    end
    show(io, mime, sys.priors)
    for planet in sys.planets
        show(io, mime, planet)
    end
    for like in sys.observations
        show(io, mime, like)
    end
    print(io, "\n")
    println(io)
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
julia> arr2nt = Octofitter.make_arr2nt(system);
julia> params = Octofitter.sample_priors(system)
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
        if length(system.priors.priors[key]) > 1
            ex_is = []
            # Handle vector-valued distributions
            for _ in 1:length(system.priors.priors[key])
                i += 1
                ex_i = :(arr[$i])
                push!(ex_is, ex_i)
            end
            ex  = :(
                $key = ($(ex_is...),)
            )
        else
            i += 1
            ex  = :(
                $key = arr[$i]
            )
        end
        push!(body_sys_priors,ex)
    end

    # Deterministic variables for the system
    body_sys_determ = Expr[]
    if isnothing(system.derived)
        push!(body_sys_determ,:(sys = sys0))
    else
        for (j,(key,func)) in enumerate(zip(keys(system.derived.variables), values(system.derived.variables)))
            ex = :(
                $(Symbol("sys$j")) = (; $(Symbol("sys$(j-1)"))..., $key = $func($(Symbol("sys$(j-1)"))))
            )
            push!(body_sys_determ,ex)
        end
        l = length(keys(system.derived.variables))
        push!(body_sys_determ,:(sys = $(Symbol("sys$l"))))
    end

    # Planets: priors & derived variables
    body_planets = Expr[]
    for planet in system.planets
        
        # Priors
        body_planet_priors = Expr[]
        for key in keys(planet.priors.priors)
            if length(planet.priors.priors[key]) > 1
                ex_is = []
                # Handle vector-valued distributions
                for _ in 1:length(planet.priors.priors[key])
                    i += 1
                    ex_i = :(arr[$i])
                    push!(ex_is, ex_i)
                end
                ex  = :(
                    $key = ($(ex_is...),)
                )
            else
                i += 1
                ex  = :(
                    $key = arr[$i]
                )
            end


            push!(body_planet_priors,ex)
        end
        j = 0
        body_planet_determ = Expr[]
        if !isnothing(planet.derived)
            # Resolve derived vars.
            for (key,func) in zip(keys(planet.derived.variables), values(planet.derived.variables))
                ex = :(
                    $(Symbol("planet$(j+1)")) = (; $(Symbol("planet$j"))..., $key = $func(sys, $(Symbol("planet$j"))))
                )
                push!(body_planet_determ,ex)
                j += 1
            end

            ex = :(
                $(planet.name) = (;
                    $(body_planet_priors...),
                    $(body_planet_determ...)
                )
            )
        end
        ex = :($(planet.name) = begin
            planet0 = (;$(body_planet_priors...));
            $(body_planet_determ...);
            $(Symbol("planet$(j)"))
        end)
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
        sys0 = (;$(body_sys_priors...))
        # Resolve derived system variables
        $(body_sys_determ...)
        # Get resolved planets
        pln = (;$(body_planets...))
        # Merge planets into resolved system
        sys_res_pln = (;sys..., planets=pln)
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
            lp += $logpdf($prior_distribution, arr[$i]);
            if !isfinite(lp)
                println("invalid prior value encountered for prior: Distributions.logpdf(", $prior_distribution, ", ", arr[$i], ")")
                error()
            end
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
                    lp += 0 <= arr[$i] < 1 ? $logpdf($prior_distribution, arr[$i]) : -Inf;
                    if !isfinite(lp)
                        println("invalid prior value encountered for prior: Distributions.logpdf(", $prior_distribution, ", ", arr[$i], ")")
                        error()
                    end
                )
            else
                ex = :(
                    lp += $logpdf($prior_distribution, arr[$i]);
                    if !isfinite(lp)
                        println("invalid prior value encountered for prior: Distributions.logpdf(", $prior_distribution, ", ", arr[$i], ")")
                    end
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
# Note! The array of values passed in must already be transformed back to the correct domain.
function make_ln_prior_transformed(system::System)

    i = 0
    prior_evaluations = Expr[]

    # System priors
    for prior_distribution in values(system.priors.priors)
        # prior_unconstrained = Bijectors.transformed(prior_distribution)
        if length(prior_distribution) > 1
            samples = []
            for _ in 1:length(prior_distribution)
                i += 1
                samples = [samples; :(arr[$i])]
            end
            samples = :(SVector($(samples...),))
        else
            i += 1
            samples = :(arr[$i])
        end
        ex = :(
            p = $logpdf_with_trans($prior_distribution, $samples, sampled);
            # Try and "heal" out of bounds values.
            # Since we are sampling from the unconstrained space they only happen due to insufficient numerical 
            # precision. 
            if !isfinite(p) && $(eltype(prior_distribution) <: AbstractFloat)
                # println("invalid prior value encountered for prior: Bijectors.logpdf_with_trans(", $prior_distribution, ", ", arr[$i], ", $sampled)=", p)
                if sign(p) > 1
                    return prevfloat(typemax($(eltype(prior_distribution))))
                else
                    return nextfloat(typemin($(eltype(prior_distribution))))
                end
            end;
            lp += p
        )
        push!(prior_evaluations,ex)
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            # prior_distribution = Bijectors.transformed(prior_distribution)
            if length(prior_distribution) > 1
                samples = []
                for _ in 1:length(prior_distribution)
                    i += 1
                    samples = [samples; :(arr[$i])]
                end
                samples = :(SVector($(samples...),))
            else
                i += 1
                samples = :(arr[$i])
            end
            ex = :(
                p = $logpdf_with_trans($prior_distribution, $samples, sampled);
                # Try and "heal" out of bounds values.
                # Since we are sampling from the unconstrained space they only happen due to insufficient numerical 
                # precision. 
                if !isfinite(p) && $(eltype(prior_distribution) <: AbstractFloat)
                    # println("invalid prior value encountered for prior: Bijectors.logpdf_with_trans(", $prior_distribution, ", ", arr[$i], ", $sampled)=", p)
                    if sign(p) > 1
                        return prevfloat(typemax($(eltype(prior_distribution))))
                    else
                        return nextfloat(typemin($(eltype(prior_distribution))))
                    end
                end;
                lp += p
                )
            push!(prior_evaluations,ex)
        end
    end
    if isempty(prior_evaluations)
        error("Model includes no free variables")
    end
    # # System priors
    # for prior_distribution in values(system.priors.priors)
    #     i += 1
    #     ex = :(
    #         lp += $logpdf_with_trans($prior_distribution, arr[$i], sampled)
    #     )
    #     push!(prior_evaluations,ex)
    # end

    # # Planet priors
    # for planet in system.planets
    #     # for prior_distribution in values(planet.priors.priors)
    #     for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
    #         i += 1
    #         ex = :(
    #             lp += $logpdf_with_trans($prior_distribution, arr[$i], sampled)
    #         )
    #         push!(prior_evaluations,ex)
    #     end
    # end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr,sampled)
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

# Generate a callback function to efficiently sample from a system's prior distributions.
function make_prior_sampler(system::System)

    # This function uses meta-programming to unroll all the code at compile time.
    # This is a big performance win, since it avoids looping over all the different
    # types of distributions that might be specified as priors.
    # Otherwise we would have to loop through an abstract vector and do runtime dispatch!
    # This way all the code gets inlined into a single tight numberical function in most cases.
    prior_sample_expressions = Expr[]




    # System priors
    for prior_distribution in values(system.priors.priors)
        # Performance: Instead of splatting, loop through according to the
        # statically known distribution length.
        push!(prior_sample_expressions, :(sample = $rand(rng, $prior_distribution)))
        for i in 1:length(prior_distribution)
            push!(prior_sample_expressions, :(prior_samples = (prior_samples..., sample[$i])))
        end
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for prior_distribution in values(planet.priors.priors)
            # Performance: Instead of splatting, loop through according to the
            # statically known distribution length.
            push!(prior_sample_expressions, :(sample = $rand(rng, $prior_distribution)))
            for i in 1:length(prior_distribution)
                push!(prior_sample_expressions, :(prior_samples = (prior_samples..., sample[$i])))
            end
        end
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (rng)
        prior_samples = ()
        @inbounds begin
           $(prior_sample_expressions...) 
        end
        return prior_samples
    end))
end


# Replaces `θ = Bijectors.invlink.(priors_vec, θ_t)` with a type stable
# unrolled version.
function make_Bijector_invlinkvec(priors_vec)

    i = 0
    parameter_transformations = Expr[]

    # System priors
    for prior_distribution in priors_vec
        if length(prior_distribution) > 1
            # Vector-valued distribution
            array_access_ex = []
            for _ in 1:length(prior_distribution)
                i += 1
                push!(array_access_ex, :(arr[$i]))
            end
            ex = :(
                $(Bijectors.invlink)($prior_distribution, ($(array_access_ex...),))...
            )
        else
            i += 1
            ex = :(
                $(Bijectors.invlink)($prior_distribution, arr[$i])
            )
        end
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
            theta_out = tuple(
                $(parameter_transformations...) 
            )
        end
        return theta_out
    end))
end




"""
Sample system parameters from prior distributions.
"""
function drawfrompriors(system::System)
    θ = Octofitter.sample_priors(system)
    arr2nt = Octofitter.make_arr2nt(system)
    θnt = arr2nt(θ)
    return θnt
end
export drawfrompriors


# Helper function for determining the umber type to use in a likelihood function. 
# It loops through all system and planet variables in a nested named tuple structure
# and promotes them
_system_number_type(T::NamedTuple) = _system_number_type(typeof(T))
@generated function _system_number_type(::Type{NamedTuple{Keys,Vals}}) where {Keys, Vals}
    T = Bool # The narrowest number type
    for (K,V) in zip(Keys,fieldtypes(Vals))
        if V <: Number
            T = promote_type(T, V)
        elseif V <: NamedTuple
            T = promote_type(T, _system_number_type(V))
        elseif V <: Tuple
            T = promote_type(T, eltype(V))
        end
    end
    return float(T) # Now promote to an appropriate floating point type
end
_planet_orbit_type(::Planet{OrbitType}) where {OrbitType} = OrbitType