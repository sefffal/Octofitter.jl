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
    likelihoodname(likeobj::AbstractLikelihood)

Return the name for a likelihood object. 
Most likelihood objects have a `name` field, but some specialized
types may override this function to provide their name differently.
"""
likelihoodname(likeobj::AbstractLikelihood) = likeobj.name
export likelihoodname



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

A group of zero or more expressions that resolve to a parameter for the model.
Captures support access to constant variables.
Must always return the same value for the same input (be a pure expression) and 
support autodifferentiation.
"""
struct Derived
    variables::OrderedDict{Symbol}
    captured_names::Tuple
    captured_vals::Tuple
end

function Derived(;_captured=nothing, variables...)
    captured_names = isnothing(_captured) ? () : _captured[1]
    captured_vals = isnothing(_captured) ? () : _captured[2]
    return Derived(OrderedDict(variables), captured_names, captured_vals)
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
    display(kwargs_dict)
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
Creates two normally distributed random variables :Ωx and :Ωy 
and maps them to a circular domain using a derived variable 
according to `atan(:Ωy, :Ωx)`.
This may be more efficient than using e.g. `Ω = Uniform(-pi, pi)`
because the sampler can wrap around freely.

Example:
```julia
@variables begin
    a ~ Uniform(1, 10)
    e ~ Uniform(0, 1)
    Ω ~ UniformCircular()
end
```
"""
struct UniformCircular <: Parameterization
    domain::Float64
end
UniformCircular() = UniformCircular(2π)
export UniformCircular

# We need to create a "prior" on the length of the unit vector so that it doesn't get pinched at (0,0)
struct UnitLengthPrior{X,Y} <: AbstractLikelihood
    varx::Symbol
    vary::Symbol
end
likelihoodname(::UnitLengthPrior{X,Y}) where {X,Y} = "unitlengthprior_$(X)_$(Y)"

# Generic expandparam for regular distributions/values
expandparam(var, d::Distribution) = (priors=[(var, d)], derived=[], likelihoods=[])
expandparam(var, n::Number) = (priors=[], derived=[(var, n)], likelihoods=[])
expandparam(var, f::Expr) = (priors=[], derived=[(var, f)], likelihoods=[])

# Expand UniformCircular into priors and derived variables
function expandparam(var::Symbol, p::UniformCircular)
    varx = Symbol("$(var)x")
    vary = Symbol("$(var)y")
    
    # The derived expression to compute the angle
    derived_expr = :(atan($vary, $varx) / 2π * $(p.domain))
    
    # Create the unit length prior
    unit_length_prior = UnitLengthPrior{varx,vary}(varx, vary)
    
    return (
        priors = [
            (varx, Normal(0, 1)),
            (vary, Normal(0, 1))
        ],
        derived = [
            (var, derived_expr)
        ],
        likelihoods = [unit_length_prior]
    )
end

_isprior(::UnitLengthPrior) = true
instrument_name(::UnitLengthPrior{X,Y}) where {X,Y} = "unitlengthprior-$X-$Y" 
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
    θ = merge(θ_system, θ_planet)
    x = getproperty(θ, X)
    y = getproperty(θ, Y)
    vector_length = sqrt(x^2 + y^2)
    return logpdf(LogNormal(log(1.0), 0.1), vector_length);
end
function Base.show(io::IO, mime::MIME"text/plain", @nospecialize like::UnitLengthPrior{X,Y}) where {X,Y}
    T = typeof(like)
    println(io, "$(T): √($X^2+$Y^2) ~ LogNormal(log(1), 0.02)")
end
generate_from_params(like::UnitLengthPrior, θ_planet, orbit) = like


# User-defined likelihood for expressions that should follow a distribution
struct UserLikelihood{TSym_LHS, TSym_RHS} <: AbstractLikelihood
    priors::Priors
    derived::Derived
    name::String
end
# Constructor to embed the symbol as a type parameter
UserLikelihood(sym_lhs::Symbol, sym_rhs::Symbol, name::String) =  UserLikelihood{sym_lhs, sym_rhs}(Priors(), Derived(), name)
UserLikelihood(sym_lhs::Symbol, sym_rhs::Symbol, name::Symbol) =  UserLikelihood{sym_lhs, sym_rhs}(Priors(), Derived(), String(name))

# Required AbstractLikelihood interface methods
likelihoodname(like::UserLikelihood) = like.name
_isprior(::UserLikelihood) = true
likeobj_from_epoch_subset(like::UserLikelihood, obs_inds) = like
TypedTables.Table(::UserLikelihood) = nothing
generate_from_params(like::UserLikelihood, θ_planet, orbit) = like

# System-level likelihood
function ln_like(user_like::UserLikelihood{TSym_LHS, TSym_RHS}, θ_system::NamedTuple, _args...) where {TSym_LHS, TSym_RHS}
    lhs = getproperty(θ_system, TSym_LHS)
    rhs = getproperty(θ_system, TSym_RHS)
    if rhs isa NTuple{N,<:Number} where N
        rhs = SVector(rhs)
    end
    if lhs isa Distribution
        return logpdf(lhs, rhs)
    elseif rhs isa Distribution
        return logpdf(rhs, lhs)
    else
        error("neither the left nor right hand side of the `~` expression evaluated to a distribution")
    end
end

# Planet-level likelihood
function ln_like(user_like::UserLikelihood{TSym_LHS, TSym_RHS}, θ_system::NamedTuple, θ_planet::NamedTuple, θ_obs::NamedTuple, _args...) where {TSym_LHS, TSym_RHS}
    θ = merge(θ_planet, θ_obs)
    lhs = getproperty(θ, TSym_LHS)
    rhs = getproperty(θ, TSym_RHS)
    if rhs isa NTuple{N,<:Number} where N
        rhs = SVector(rhs)
    end
    return logpdf(lhs, rhs)
end

# Show method
function Base.show(io::IO, mime::MIME"text/plain", like::UserLikelihood{TSym_LHS, TSym_RHS}) where {TSym_LHS, TSym_RHS}
    println(io, "UserLikelihood: $TSym_LHS ~ $TSym_RHS")
end

# We need a blank likelihood type to hold variables when constructing
# e.g. prior only models.
# TODO: In future if/when we get RHS ~ working in our models, we can probably 
# use those objects for this purpose too.
struct BlankLikelihood <: AbstractLikelihood
    priors::Priors
    derived::Derived
    name::String
    function BlankLikelihood(
        variables::Tuple{Priors,Derived}=(Priors(),Derived()),
        name=""
    )
        (priors,derived)=variables
        return new(priors,derived,name)
    end
end
function Octofitter.ln_like(::BlankLikelihood, θ_system, args...)
    T = Octofitter._system_number_type(θ_system)
    return zero(T)
end

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
end
export Planet
function Planet(;
    name::Union{Symbol,AbstractString},
    basis::Type,
    variables::Tuple,
    likelihoods=()
)
    (priors,derived,additional_likelihoods...)=variables
    name = Symbol(name)
    # Type asserts
    priors::Priors
    derived::Derived
    for l in additional_likelihoods
        l::AbstractLikelihood
    end
    likes = (likelihoods..., additional_likelihoods...)
    
    # Check for duplicate observation/likelihood names on this planet
    like_names = String[]
    for like in likes
        like_name = likelihoodname(like)
        if like_name in like_names
            error("Planet $name: Duplicate observation/likelihood name '$like_name'. Each observation/likelihood attached to a planet must have a unique name.")
        end
        push!(like_names, like_name)
    end
    
    # Check for duplicate variables in Prior and Derived blocks
    prior_keys = Set(keys(priors.priors))
    if !isnothing(derived)
        derived_keys = Set(keys(derived.variables))
        overlap = intersect(prior_keys, derived_keys)
        if !isempty(overlap)
            error("Planet $name: Variables $(collect(overlap)) are defined as both a prior (~) and derived variable (=). Each variable must be defined only once.")
        end
    end
    
    return Planet{
        basis,
        typeof(priors),typeof(derived),typeof(likes)
    }(priors, derived, likes, name)
end


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
end
export System
function System(;
    name::Union{Symbol,AbstractString},
    variables::Tuple,
    companions=(),
    likelihoods=()
)
    (priors,derived,additional_likelihoods...)=variables
    name = Symbol(name)
    # Type asserts
    priors::Priors
    derived::Derived
    for l in additional_likelihoods
        l::AbstractLikelihood
    end
    for p in companions
        p::Planet
    end
    likes = (likelihoods..., additional_likelihoods...)
    
    # Check for duplicate observation/likelihood names at system level
    like_names = String[]
    for like in likes
        like_name = likelihoodname(like)
        if like_name in like_names
            error("System $name: Duplicate observation/likelihood name '$like_name'. Each observation/likelihood attached to a system must have a unique name.")
        end
        push!(like_names, like_name)
    end
    
    # Check for duplicate variables in Prior and Derived blocks
    prior_keys = Set(keys(priors.priors))
    if !isnothing(derived)
        derived_keys = Set(keys(derived.variables))
        overlap = intersect(prior_keys, derived_keys)
        if !isempty(overlap)
            error("System $name: Variables $(collect(overlap)) are defined as both a prior (~) and derived variable (=). Each variable must be defined only once.")
        end
    end
    
    if isempty(companions)
        planets_nt = (;)
    else
        planets_nt = namedtuple(
            getproperty.(companions, :name),
            companions
        )
    end
    return System{
        typeof(priors),typeof(derived),typeof(likes),typeof(planets_nt)
    }(priors, derived, likes, planets_nt, name)
end

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
        push!(priors_vec, prior_distribution)
    end

    # System observation priors
    for obs in system.observations
        if !hasproperty(obs, :priors)
            continue
        end
        for prior_distribution in values(obs.priors.priors)
            push!(priors_vec, prior_distribution)
        end
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            push!(priors_vec, prior_distribution)
        end
        
        # Planet observation priors
        for obs in planet.observations
            if !hasproperty(obs, :priors)
                continue
            end
            for prior_distribution in values(obs.priors.priors)
                push!(priors_vec, prior_distribution)
            end
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
        for (j,(key,expr)) in enumerate(zip(keys(system.derived.variables), values(system.derived.variables)))
            # Create let bindings with the actual captured values
            captured_bindings = [:($(system.derived.captured_names[i]) = $(system.derived.captured_vals[i])) 
                               for i in 1:length(system.derived.captured_names)]
            
            # Get all variable names available up to this point
            prior_keys = collect(keys(system.priors.priors))
            derived_keys_so_far = j == 1 ? Symbol[] : collect(keys(system.derived.variables))[1:j-1]
            all_available_keys = vcat(prior_keys, derived_keys_so_far)
            
            ex = :(
                $(Symbol("sys$j")) = let _prev=$(Symbol("sys$(j-1)"))
                    (; _prev..., $key = let $(captured_bindings...)
                        # Make previous variables available
                        $([:($(k) = _prev.$k) for k in all_available_keys]...)
                        result = $expr
                        result
                    end)
                end
            )
            push!(body_sys_determ,ex)
        end
        l = length(keys(system.derived.variables))
        push!(body_sys_determ,:(sys = $(Symbol("sys$l"))))
    end

    # System observations
    body_observations = Expr[]
    for obs in system.observations
        
        # Priors
        body_obs_priors = Expr[]
        if hasproperty(obs, :priors)
            for key in keys(obs.priors.priors)
                if length(obs.priors.priors[key]) > 1
                    ex_is = []
                    # Handle vector-valued distributions
                    for _ in 1:length(obs.priors.priors[key])
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
                push!(body_obs_priors,ex)
            end
        end
        
        j = 0
        body_obs_determ = Expr[]
        if hasproperty(obs, :priors) && !isnothing(obs.derived)
            # Build let bindings with actual captured values
            captured_bindings = [:($(obs.derived.captured_names[i]) = $(obs.derived.captured_vals[i])) 
                            for i in 1:length(obs.derived.captured_names)]
            
            for (key,expr) in zip(keys(obs.derived.variables), values(obs.derived.variables))
                # Get all observation variable names available up to this point
                obs_prior_keys = hasproperty(obs, :priors) ? collect(keys(obs.priors.priors)) : Symbol[]
                obs_derived_keys_so_far = j == 0 ? Symbol[] : collect(keys(obs.derived.variables))[1:j]
                all_obs_keys = vcat(obs_prior_keys, obs_derived_keys_so_far)
                
                ex = :(
                    $(Symbol("obs$(j+1)")) = let system=sys, _prev=$(Symbol("obs$j"))
                        (; _prev..., $key = let $(captured_bindings...)
                            # Make previous observation variables available
                            $([:($(k) = _prev.$k) for k in all_obs_keys]...)
                            result = $expr
                            if result isa Distributions.Distribution
                                error("System observation derived variable '$($(Meta.quot(key)))' evaluated to a Distribution object ($result). Did you mean to sample from it using '~' instead of '='?")
                            end
                            result
                        end)
                    end
                )
                push!(body_obs_determ,ex)
                j += 1
            end
        end
        
        if isempty(body_obs_priors) && isempty(body_obs_determ)
            continue
        end
        name = normalizename(likelihoodname(obs))
        ex = :($name = begin
            obs0 = (;$(body_obs_priors...));
            $(body_obs_determ...);
            (;$(Symbol("obs$(j)"))...)
        end)
        push!(body_observations,ex)
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
            # Build let bindings with actual captured values
            captured_bindings = [:($(planet.derived.captured_names[i]) = $(planet.derived.captured_vals[i])) 
                            for i in 1:length(planet.derived.captured_names)]
            
            for (key,expr) in zip(keys(planet.derived.variables), values(planet.derived.variables))
                # Get all variable names available up to this point
                prior_keys = collect(keys(planet.priors.priors))
                derived_keys_so_far = j == 0 ? Symbol[] : collect(keys(planet.derived.variables))[1:j]
                all_available_keys = vcat(prior_keys, derived_keys_so_far)
                
                ex = :(
                    $(Symbol("planet$(j+1)")) = let system=sys, _prev=$(Symbol("planet$j"))
                        (; _prev..., $key = let $(captured_bindings...)
                            # Make previous planet variables available
                            $([:($(k) = _prev.$k) for k in all_available_keys]...)
                            result = $expr
                            if result isa Distributions.Distribution
                                error("Planet derived variable '$($(Meta.quot(key)))' evaluated to a Distribution object ($result). Did you mean to sample from it using '~' instead of '='?")
                            end
                            result
                        end)
                    end
                )
                push!(body_planet_determ,ex)
                j += 1
            end
        end

        # Planet observations
        planet_observations = Expr[]
        for obs in planet.observations
            
            # Priors
            planet_obs_priors = Expr[]
            if hasproperty(obs, :priors)
                for key in keys(obs.priors.priors)
                    if length(obs.priors.priors[key]) > 1
                        ex_is = []
                        # Handle vector-valued distributions
                        for _ in 1:length(obs.priors.priors[key])
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
                    push!(planet_obs_priors,ex)
                end
            end
            
            k = 0
            planet_obs_determ = Expr[]
            # For planet observations with derived variables:
            if hasproperty(obs, :derived) && !isnothing(obs.derived)
                # Build let bindings with actual captured values
                captured_bindings = [:($(obs.derived.captured_names[i]) = $(obs.derived.captured_vals[i])) 
                                for i in 1:length(obs.derived.captured_names)]
                
                for (key,expr) in zip(keys(obs.derived.variables), values(obs.derived.variables))
                    # Get all planet and observation variable names available
                    planet_prior_keys = collect(keys(planet.priors.priors))
                    planet_derived_keys = !isnothing(planet.derived) ? collect(keys(planet.derived.variables)) : Symbol[]
                    obs_prior_keys = hasproperty(obs, :priors) ? collect(keys(obs.priors.priors)) : Symbol[]
                    obs_derived_keys_so_far = k == 0 ? Symbol[] : collect(keys(obs.derived.variables))[1:k]
                    
                    # Include both planet variables and observation variables
                    ex = :(
                        $(Symbol("obs$(k+1)")) = let planet=$(Symbol("planet$(j)")), _prev=$(Symbol("obs$k"))
                            (; _prev..., $key = let $(captured_bindings...)
                                # Make planet and previous observation variables available
                                $([:($(pk) = planet.$pk) for pk in vcat(planet_prior_keys, planet_derived_keys)]...)
                                $([:($(ok) = _prev.$ok) for ok in vcat(obs_prior_keys, obs_derived_keys_so_far)]...)
                                result = $expr
                                if result isa Distributions.Distribution
                                    error("Planet observation derived variable '$($(Meta.quot(key)))' evaluated to a Distribution object ($result). Did you mean to sample from it using '~' instead of '='?")
                                end
                                result
                            end)
                        end
                    )
                    push!(planet_obs_determ,ex)
                    k += 1
                end
            end

            if isempty(planet_obs_priors) && isempty(planet_obs_determ)
                continue
            end
            
            name = normalizename(likelihoodname(obs))
            ex = :($name = begin
                obs0 = (;$(planet_obs_priors...));
                $(planet_obs_determ...);
                (;$(Symbol("obs$(k)"))...)
            end)
            push!(planet_observations,ex)
        end

        ex = :($(planet.name) = begin
            planet0 = (;$(body_planet_priors...));
            $(body_planet_determ...);
            (;(;$(Symbol("planet$(j)"))...)..., observations=(;$(planet_observations...)))
        end)
        push!(body_planets,ex)
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problems"
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
        # Get resolved observations
        # obs = $(body_observations...)
        # Get resolved planets
        pln = (;$(body_planets...))
        # Merge planets into resolved system
        sys_res_pln = (;sys..., observations=(;$(body_observations...)), planets=pln)
        return sys_res_pln
    end))

    return func
end

# From CSV.jl:
const RESERVED = Set(["local", "global", "export", "let",
    "for", "struct", "while", "const", "continue", "import",
    "function", "if", "else", "try", "begin", "break", "catch",
    "return", "using", "baremodule", "macro", "finally",
    "module", "elseif", "end", "quote", "do"])
function normalizename(name::String)::Symbol
    uname = strip(Base.Unicode.normalize(name))
    id = Base.isidentifier(uname) ? uname : map(c->Base.is_id_char(c) ? c : '_', uname)
    cleansed = string((isempty(id) || !Base.is_id_start_char(id[1]) || id in RESERVED) ? "_" : "", id)
    return Symbol(replace(cleansed, r"(_)\1+"=>"_"))
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

    # System observation priors
    for obs in system.observations
        if !hasproperty(obs, :priors)
            continue
        end
        for prior_distribution in values(obs.priors.priors)
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

        # Planet observation priors
        for obs in planet.observations
            for prior_distribution in values(obs.priors.priors)
                i += 1
                # Apply same Beta distribution workaround if needed
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
        # Add contributions from all priors
        # @inbounds begin
           $(prior_evaluations...) 
        # end
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

    # System observation priors
    for obs in system.observations
        if !hasproperty(obs, :priors)
            continue
        end
        for prior_distribution in values(obs.priors.priors)
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

        # Planet observation priors
        for obs in planet.observations
            if !hasproperty(obs, :priors)
                continue
            end
            for prior_distribution in values(obs.priors.priors)
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
    end
    
    if isempty(prior_evaluations)
        error("Model includes no free variables")
    end

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
    for obs in system.observations
        if !hasproperty(obs, :priors)
            continue
        end
        for prior_distribution in values(obs.priors.priors)
            # Performance: Instead of splatting, loop through according to the
            # statically known distribution length.
            push!(prior_sample_expressions, :(sample = $rand(rng, $prior_distribution)))
            for i in 1:length(prior_distribution)
                push!(prior_sample_expressions, :(prior_samples = (prior_samples..., sample[$i])))
            end
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
        for obs in planet.observations
            if !hasproperty(obs, :priors)
                continue
            end
            for prior_distribution in values(obs.priors.priors)
                # Performance: Instead of splatting, loop through according to the
                # statically known distribution length.
                push!(prior_sample_expressions, :(sample = $rand(rng, $prior_distribution)))
                for i in 1:length(prior_distribution)
                    push!(prior_sample_expressions, :(prior_samples = (prior_samples..., sample[$i])))
                end
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