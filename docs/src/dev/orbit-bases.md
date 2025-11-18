# Orbit Bases and Parameterizations

This document explains the different orbit bases (parameterizations) available in Octofitter via the [PlanetOrbits.jl](https://sefffal.github.io/PlanetOrbits.jl/) package, how they're implemented, and when to use each one.

## Overview

Octofitter re-exports `PlanetOrbits.jl`, which provides several orbit parameterizations optimized for different use cases. The choice of orbit basis affects:

- **Sampling efficiency**: Some parameterizations avoid degeneracies at low eccentricity/inclination
- **Physical interpretation**: Different bases emphasize different observable quantities
- **Computational performance**: Some calculations are faster in certain bases

## Available Orbit Types

### 1. Visual{KepOrbit} - Campbell Elements

**Most common parameterization for visual astrometry.**

#### Parameters

```julia
@planet b Visual{KepOrbit} begin
    # Orbital elements
    a ~ LogUniform(1, 100)           # Semi-major axis [AU]
    e ~ Uniform(0.0, 0.5)            # Eccentricity
    i ~ Sine()                        # Inclination [radians]
    ω ~ UniformCircular()            # Argument of periastron [radians]
    Ω ~ UniformCircular()            # Longitude of ascending node [radians]

    # Time reference
    θ ~ UniformCircular()            # Mean anomaly at reference epoch
    tp = θ_at_epoch_to_tperi(θ, ref_epoch, ...)  # Periastron time

    # Or alternatively:
    # τ ~ Uniform(0, 1)               # Fraction of orbit period
    # tp = θ_epoch_to_tperi_years(τ, ref_epoch, ...)
end
```

#### Coordinate System

The `Visual` wrapper transforms from 3D Cartesian coordinates to on-sky projected positions:

```julia
# Under the hood:
elem = KepOrbit(M=1.0, plx=50.0, a=10.0, e=0.2, i=π/4, ω=π/6, Ω=π/3, tp=50000.0)
visual_elem = Visual(elem, plx=50.0)

# Solving returns projected positions:
sol = orbitsolve(visual_elem, epoch)
ra_offset = raoff(sol)   # milliarcseconds
dec_offset = decoff(sol) # milliarcseconds
```

#### Degeneracies

**Warning**: Campbell elements have degeneracies:

- As `e → 0`: `ω` becomes undefined (circular orbit has no periastron)
- As `i → 0`: `Ω` becomes undefined (face-on orbit has no ascending node)
- As `i → π`: `Ω` and `ω` trade roles

These can make sampling difficult for low-eccentricity or low-inclination systems. Consider using `ThieleInnesOrbit` instead.

### 2. ThieleInnesOrbit - Thiele-Innes Elements

**Recommended for low-eccentricity orbits and faster sampling.**

#### Parameters

```julia
@planet b ThieleInnesOrbit begin
    # Thiele-Innes elements [milliarcseconds]
    A ~ Normal(0, 1000)
    B ~ Normal(0, 1000)
    F ~ Normal(0, 1000)
    G ~ Normal(0, 1000)

    e ~ Uniform(0.0, 0.5)            # Eccentricity

    # Time reference (same as Campbell)
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(θ, ref_epoch; system.plx, M, e, A, B, F, G)
end
```

#### Advantages

1. **No degeneracies**: A, B, F, G are well-defined at all eccentricities and inclinations
2. **Faster sampling**: Linear relationship to observables makes HMC/NUTS more efficient
3. **Natural for astrometry**: Directly relates to on-sky motion

#### Relationship to Campbell Elements

The Thiele-Innes constants relate to Campbell elements via:

```
A = a (cos(ω)cos(Ω) - sin(ω)sin(Ω)cos(i))
B = a (cos(ω)sin(Ω) + sin(ω)cos(Ω)cos(i))
F = a (-sin(ω)cos(Ω) - cos(ω)sin(Ω)cos(i))
G = a (-sin(ω)sin(Ω) + cos(ω)cos(Ω)cos(i))
```

where `a` is projected semi-major axis in milliarcseconds.

#### Converting Back to Campbell

After fitting with Thiele-Innes, you can convert to Campbell elements:

```julia
# Get Thiele-Innes orbit objects
orbits_ti = Octofitter.construct_elements(model, results, :b, :)

# Convert to Campbell
orbits_campbell = Visual{KepOrbit}.(orbits_ti)

# Extract traditional elements
a = semimajoraxis.(orbits_campbell)
e = eccentricity.(orbits_campbell)
i = inclination.(orbits_campbell)
ω = argumentofperiastron.(orbits_campbell)
Ω = longitudeofascendingnode.(orbits_campbell)
```

See the [Thiele-Innes tutorial](../thiele-innes.md) for a complete example.

### 3. RadialVelocityOrbit - RV-Only Parameterization

**Optimized for radial velocity data.**

#### Parameters

```julia
@planet b RadialVelocityOrbit begin
    K ~ LogUniform(0.1, 100)         # RV semi-amplitude [m/s]
    e ~ Uniform(0.0, 0.5)            # Eccentricity
    ω ~ UniformCircular()            # Argument of periastron [radians]

    # Time reference
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(θ, ref_epoch, ...)
end
```

#### Key Differences

- Uses `K` (RV semi-amplitude) instead of `a` (semi-major axis)
- No inclination `i` or position angle `Ω` (degenerate with `K` for RV-only)
- Cannot compute projected positions (no `raoff`, `decoff`)

#### Relationship to Visual Orbit

When both RV and astrometry are available:

```
K = (2π / P) * (a * sin(i)) / sqrt(1 - e²)
```

where `a * sin(i)` is the projected semi-major axis in AU.

#### Usage

```julia
using OctofitterRadialVelocity

rvlike = StarAbsoluteRVLikelihood(
    rv_data,
    name="instrument",
    variables=@variables begin
        offset ~ Normal(0, 100)      # Systemic velocity offset [m/s]
        jitter ~ LogUniform(0.1, 10) # Instrumental jitter [m/s]
    end
)

@planet b RadialVelocityOrbit begin
    K ~ LogUniform(1, 100)
    e ~ Uniform(0, 0.5)
    ω ~ UniformCircular()
    P ~ LogUniform(1, 1000)          # Period [days]
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(θ, ref_epoch, P=P, e=e)
end
```

See the [RV tutorials](../rv-1.md) for complete examples.

### 4. KepOrbit - Raw Keplerian Elements

**Low-level interface without projection.**

#### Parameters

Same as `Visual{KepOrbit}`, but without the visual projection wrapper. Returns 3D Cartesian positions rather than projected on-sky offsets.

```julia
@planet b KepOrbit begin
    # System mass required
    M ~ truncated(Normal(1.2, 0.1), lower=0.1)  # System mass [M☉]

    # Standard Keplerian elements
    a ~ LogUniform(1, 100)           # Semi-major axis [AU]
    e ~ Uniform(0.0, 0.5)            # Eccentricity
    i ~ Sine()                        # Inclination [radians]
    ω ~ UniformCircular()            # Argument of periastron [radians]
    Ω ~ UniformCircular()            # Longitude of ascending node [radians]
    tp ~ UniformCircular()           # Time of periastron [MJD]
end
```

#### When to Use

- **Transit fitting**: Need 3D positions for transit geometry
- **Direct imaging**: Working with absolute positions rather than relative
- **Special calculations**: Need access to 3D coordinates

#### Coordinate Access

```julia
elem = KepOrbit(M=1.0, a=10.0, e=0.2, ...)
sol = orbitsolve(elem, epoch)

# 3D Cartesian coordinates [AU]
x = posx(sol)
y = posy(sol)
z = posz(sol)

# 3D Cartesian velocities [AU/day]
vx = velx(sol)
vy = vely(sol)
vz = velz(sol)
```

### 5. CartesianOrbit - State Vectors

**Cartesian position and velocity at a reference epoch.**

#### Parameters

```julia
@planet b CartesianOrbit begin
    # Position at reference epoch [AU]
    x ~ Normal(0, 10)
    y ~ Normal(0, 10)
    z ~ Normal(0, 10)

    # Velocity at reference epoch [AU/day]
    vx ~ Normal(0, 0.1)
    vy ~ Normal(0, 0.1)
    vz ~ Normal(0, 0.1)

    # System mass for propagation
    M ~ truncated(Normal(1.2, 0.1), lower=0.1)
end
```

#### When to Use

- **Astrometric acceleration**: Fitting proper motion anomaly with insufficient data for full orbit
- **Very long period orbits**: When you only see a small arc
- **Integration with N-body**: Initial conditions for dynamical simulations

### 6. FixedPosition - Single Epoch

**Special case for single-epoch imaging.**

**Location**: [`src/orbit-models.jl`](../../src/orbit-models.jl)

#### Parameters

```julia
@planet b Visual{FixedPosition} begin
    # Direct position specification [AU]
    x ~ Normal(0, 10)
    y ~ Normal(0, 10)
    z = 0  # Usually fixed at zero

    # Or use on-sky coordinates [mas]
    # ra ~ Normal(100, 10)
    # dec ~ Normal(50, 10)

    # Or use separation and position angle
    # sep ~ LogUniform(10, 100)  # mas
    # pa ~ UniformCircular()      # radians
end
```

#### Implementation Detail

`FixedPosition` is defined in Octofitter (not PlanetOrbits) and implements the `AbstractOrbit` interface:

```julia
struct FixedPosition{T<:Number} <: AbstractOrbit{T}
    x::T
    y::T
    z::T
end

# "Solving" just returns the fixed position
orbitsolve(o::FixedPosition, t) = ObitSolutionFixedPosition(o, 0, t)

# Position queries return the stored values
posx(sol::ObitSolutionFixedPosition) = sol.elem.x
posy(sol::ObitSolutionFixedPosition) = sol.elem.y
posz(sol::ObitSolutionFixedPosition) = sol.elem.z
```

#### When to Use

- **Single epoch data**: Only one observation available
- **Extremely wide orbits**: Period >> observation baseline, treat as fixed
- **Position priors**: Constrain planet position from other sources

## Orbit Type Implementation Details

### Type Hierarchy

```julia
# From PlanetOrbits.jl
abstract type AbstractOrbit{T<:Number} end

# Specific orbit types
struct KepOrbit{T} <: AbstractOrbit{T}
    M::T    # Total system mass [M☉]
    plx::T  # Parallax [mas]
    a::T    # Semi-major axis [AU]
    e::T    # Eccentricity
    i::T    # Inclination [rad]
    ω::T    # Argument of periastron [rad]
    Ω::T    # Longitude of ascending node [rad]
    tp::T   # Time of periastron passage [MJD]
end

struct ThieleInnesOrbit{T} <: AbstractOrbit{T}
    M::T    # Total system mass [M☉]
    plx::T  # Parallax [mas]
    A::T    # Thiele-Innes constant [mas]
    B::T    # Thiele-Innes constant [mas]
    F::T    # Thiele-Innes constant [mas]
    G::T    # Thiele-Innes constant [mas]
    e::T    # Eccentricity
    tp::T   # Time of periastron passage [MJD]
end

# Wrapper for projected coordinates
struct VisualOrbit{T,O<:AbstractOrbit{T}} <: AbstractOrbit{T}
    parent::O
    plx::T
    dist::T
end

const Visual{O} = VisualOrbit{<:Any,O}
```

### How Octofitter Uses Orbit Types

#### 1. Type Parameter in Planet

```julia
struct Planet{TElem<:AbstractOrbit,TTable,TDerived}
    name::Symbol
    basis::Type{TElem}              # e.g., Visual{KepOrbit}
    variables::Tuple{Priors,TDerived,UserLikelihoods}
    observations::TTable
end
```

#### 2. Code Generation in `make_ln_like`

**Location**: [`src/likelihoods/system.jl:108-109`](../../src/likelihoods/system.jl)

```julia
# Extract orbit type from planet's type parameter
OrbitType = _planet_orbit_type(system.planets[i])

# Generate construction code
planet_construction = quote
    # Merges system parameters (M, plx) with planet parameters (a, e, ...)
    $planet_key = $(OrbitType)(;merge(θ_system, θ_system.planets[$i])...)
end
```

At runtime, this becomes (for example):

```julia
planet_b = Visual{KepOrbit}(;
    M=θ_system.M,
    plx=θ_system.plx,
    a=θ_system.planets.b.a,
    e=θ_system.planets.b.e,
    i=θ_system.planets.b.i,
    ω=θ_system.planets.b.ω,
    Ω=θ_system.planets.b.Ω,
    tp=θ_system.planets.b.tp
)
```

#### 3. Orbit Solving Interface

All orbit types implement a common solving interface:

```julia
# Main solving function
sol = orbitsolve(orbit::AbstractOrbit, epoch::Number)

# Returns AbstractOrbitSolution with methods:
meananom(sol)    # Mean anomaly
trueanom(sol)    # True anomaly
eccanom(sol)     # Eccentric anomaly

# Position (type-dependent)
raoff(sol)       # RA offset [mas] (Visual only)
decoff(sol)      # Dec offset [mas] (Visual only)
posx(sol)        # X position [AU] (Cartesian)
posy(sol)        # Y position [AU] (Cartesian)
posz(sol)        # Z position [AU] (Cartesian)

# Velocity (not all types)
radvel(sol)      # Radial velocity [m/s]
pmra(sol)        # Proper motion in RA [mas/yr]
pmdec(sol)       # Proper motion in Dec [mas/yr]
```

### Performance Considerations

#### Type Stability

All orbit construction and solving is type-stable:

```julia
# Compiler knows exact return type
planet_b::Visual{KepOrbit{Float64}}
sol::OrbitSolutionVisual{Float64, KepOrbit{Float64}}
ra::Float64 = raoff(sol)
```

This enables:
- LLVM optimization
- Efficient memory layout
- Fast automatic differentiation

#### Pre-compilation

Orbit types are typically instantiated once during `make_ln_like` generation:

```julia
# At model construction time:
θ_example = sample_priors(rng)
ln_like = make_ln_like(system, θ_example)

# The above triggers compilation of:
# - Orbit construction for each planet's specific type
# - orbitsolve for those specific types
# - All likelihood calculations with concrete types
```

Once compiled, repeated likelihood evaluations reuse the specialized code.

## Choosing an Orbit Basis

### Decision Tree

```
Do you have only RV data?
├─ YES → Use RadialVelocityOrbit
└─ NO → Do you have only one epoch of imaging?
    ├─ YES → Use Visual{FixedPosition}
    └─ NO → Is eccentricity likely < 0.3?
        ├─ YES → Use ThieleInnesOrbit (faster sampling)
        └─ NO → Is this for transits or special 3D calculations?
            ├─ YES → Use KepOrbit
            └─ NO → Use Visual{KepOrbit} (default)
```

### Performance Comparison

Based on typical multi-planet systems:

| Orbit Type | Relative Speed | Sampling Efficiency | Use Case |
|------------|---------------|---------------------|----------|
| `Visual{KepOrbit}` | 1.0x (baseline) | 1.0x | General astrometry |
| `ThieleInnesOrbit` | 1.0x | **2-5x better** | Low-e orbits |
| `RadialVelocityOrbit` | **0.8x** | 1.0x | RV-only |
| `KepOrbit` | 1.0x | 1.0x | Transits |
| `CartesianOrbit` | 1.1x | 0.5x | Wide orbits, PMA |
| `FixedPosition` | **0.5x** | N/A | Single epoch |

*Speed: Likelihood evaluation time. Efficiency: Effective sample size per iteration.*

## Example: Multi-Basis System

You can mix orbit bases in the same system:

```julia
# Inner planet: well-constrained, use Campbell
@planet b Visual{KepOrbit} begin
    a ~ LogUniform(1, 5)
    e ~ Uniform(0.0, 0.9)
    # ... other elements
end

# Middle planet: low eccentricity, use Thiele-Innes
@planet c ThieleInnesOrbit begin
    A ~ Normal(0, 1000)
    B ~ Normal(0, 1000)
    F ~ Normal(0, 1000)
    G ~ Normal(0, 1000)
    e ~ Uniform(0.0, 0.2)
    # ... time reference
end

# Outer planet: only RV detection
@planet d RadialVelocityOrbit begin
    K ~ LogUniform(1, 50)
    P ~ LogUniform(1000, 10000)
    e ~ Uniform(0.0, 0.5)
    # ... time reference
end

system = System(
    companions=[b, c, d],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)
```

Each planet uses the orbit type best suited to its available data and orbital characteristics.

## See Also

- **User Documentation**:
  - [Thiele-Innes Tutorial](../thiele-innes.md) - Complete workflow with conversion
  - [RV Fitting](../rv-1.md) - Using RadialVelocityOrbit
  - [Quick Start](../quick-start.md) - Basic Visual{KepOrbit} usage

- **Developer Documentation**:
  - [Architecture Overview](architecture.md) - How orbit types fit into code generation
  - [Epoch Tables and Kepler Solving](epoch-tables-kepler.md) - How orbits are solved

- **External**:
  - [PlanetOrbits.jl Documentation](https://sefffal.github.io/PlanetOrbits.jl/) - Full API reference
