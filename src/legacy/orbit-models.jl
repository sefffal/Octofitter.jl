#=

All the code defining the different Keplerian orbit types, their parameterizations
and calculations are in PlanetOrbits.jl.

This file contains an additional subtype of PlanetOrbits.AbstractOrbit that 
has only a fixed position. This is a useful parameterization when you only
have one epoch of data.

It hooks into all the machinery of PlanetOrbits.jl and just very simply returns
the saved position.

It purposely doesn't define any methods for velocity or accleration, only positions.
=#
using PlanetOrbits

struct FixedPosition{T<:Number} <: AbstractOrbit{T}
    # Cartesian XYX positions in AU
    x::T
    y::T
    z::T
    function FixedPosition(x,y,z)
        x,y,z=promote(x,y,z)
        return new{typeof(x)}(x,y,z)
    end
end
FixedPosition(;x,y,z, kwargs...) = FixedPosition(x,y,z)

PlanetOrbits.period(::FixedPosition) = Inf
PlanetOrbits.meanmotion(::FixedPosition) = 0.0
PlanetOrbits.eccentricity(::FixedPosition) = 0.0
PlanetOrbits.totalmass(::FixedPosition) = 0.0
PlanetOrbits.semimajoraxis(::FixedPosition) = 0.0
PlanetOrbits.periastron(::FixedPosition) = 0.0


struct ObitSolutionFixedPosition{T<:Number,TEl} <: AbstractOrbitSolution
    elem::TEl
    EA::T
    t::T
end
PlanetOrbits.meananom(::ObitSolutionFixedPosition{T}) where T = zero(T)
PlanetOrbits.trueanom(::ObitSolutionFixedPosition{T}) where T = zero(T)
function PlanetOrbits.orbitsolve(o::FixedPosition, t, method::PlanetOrbits.AbstractSolver=PlanetOrbits.Auto())
    return ObitSolutionFixedPosition(o, promote(0, t)...)
end
function PlanetOrbits.orbitsolve_ν(o::FixedPosition, ν, EA=0, t=0)
    return ObitSolutionFixedPosition(o, promote(EA,  t)...)
end
function PlanetOrbits.orbitsolve_eccanom(o::FixedPosition, EA)
    return ObitSolutionFixedPosition(o, promote(EA, 0)...)
end
function PlanetOrbits.orbitsolve_meananom(o::FixedPosition, MA)
    return ObitSolutionFixedPosition(o, promote(MA, 0)...)
end


PlanetOrbits.posx(o::ObitSolutionFixedPosition{T,<:FixedPosition}) where T = o.elem.x
PlanetOrbits.posy(o::ObitSolutionFixedPosition{T,<:FixedPosition}) where T = o.elem.y
PlanetOrbits.posz(o::ObitSolutionFixedPosition{T,<:FixedPosition}) where T = o.elem.z

## Visual Wrapper
function Visual{OrbitType}(;
    plx=100,
    x=nothing,
    y=nothing,
    z=0,
    ra=nothing,
    dec=nothing,
    sep=nothing,
    pa=nothing,
    kwargs...
) where {OrbitType <: FixedPosition}
    dist = 1000/plx * PlanetOrbits.pc2au
    if !isnothing(x) && !isnothing(y)
        T = promote_type(typeof(x),typeof(y),typeof(z))
        return PlanetOrbits.VisualOrbit{T,FixedPosition{T}}(FixedPosition(x,y,z),plx,dist)
    end
    if !isnothing(ra) && !isnothing(dec)
        x = ra/plx
        y = dec/plx
        T = promote_type(typeof(x),typeof(y),typeof(z))
        return PlanetOrbits.VisualOrbit{T,FixedPosition{T}}(FixedPosition(x,y,z),plx,dist)
    end
    if !isnothing(sep) && !isnothing(pa)
        x = sep*sin(pa)/plx
        y = sep*cos(pa)/plx
        T = promote_type(typeof(x),typeof(y),typeof(z))
        return PlanetOrbits.VisualOrbit{T,FixedPosition{T}}(FixedPosition(x,y,z),plx,dist)
    end
    error("Invalid keyword parameters. Pass either x & y, ra & dec, or sep and pa.")
end
function PlanetOrbits.orbitsolve(o::Visual{<:FixedPosition}, t, method::PlanetOrbits.AbstractSolver=PlanetOrbits.Auto())
    sol = orbitsolve(o.parent, t)
    return PlanetOrbits.OrbitSolutionVisual(o, sol)
end
function PlanetOrbits.orbitsolve_ν(o::PlanetOrbits.VisualOrbit{<:FixedPosition}, _::Number)
    sol = orbitsolve_ν(o.parent, 0.0)
    return PlanetOrbits.OrbitSolutionVisual(o, sol)
end
function PlanetOrbits.orbitsolve_eccanom(o::Visual{<:FixedPosition}, _::Number)
    sol = orbitsolve_eccanom(o.parent, 0.0)
    return PlanetOrbits.OrbitSolutionVisual(o, sol)
end
function PlanetOrbits.orbitsolve_meananom(o::Visual{<:FixedPosition}, _::Number)
    sol = orbitsolve_meananom(o.parent, 0.0)
    return PlanetOrbits.OrbitSolutionVisual(o, sol)
end

PlanetOrbits.posx(o::ObitSolutionFixedPosition{T,O}) where T where O <: Visual{<:FixedPosition} = o.elem.parent.x
PlanetOrbits.posy(o::ObitSolutionFixedPosition{T,O}) where T where O <: Visual{<:FixedPosition} = o.elem.parent.y
PlanetOrbits.posz(o::ObitSolutionFixedPosition{T,O}) where T where O <: Visual{<:FixedPosition} = o.elem.parent.z
