# Kepler Solver and Orbit Types

The heart of this package is being able to take a set of Keplerian elements and output relative positions, velocities, etc.
For this, we use [PlanetOrbits.jl](https://github.com/sefffal/PlanetOrbits.jl) which adopts the same conventions as Orbitize!.

For full documentation on orbit types, coordinate conventions, and available functions, see the PlanetOrbits.jl documentation:
- [Introduction](https://sefffal.github.io/PlanetOrbits.jl/dev/introduction/) - Getting started with orbit creation and solving
- [Coordinate Conventions](https://sefffal.github.io/PlanetOrbits.jl/dev/conventions/) - Coordinate system and orbital element definitions
- [API Reference](https://sefffal.github.io/PlanetOrbits.jl/dev/api/) - Complete function reference

## Kepler Solver

The Kepler solver used to go from mean anomaly to eccentric anomaly is a tweaked version copied from [AstroLib.jl](http://juliaastro.github.io/AstroLib.jl/stable/ref/#AstroLib.kepler_solver).

From AstroLib.jl:

> Many different numerical methods exist to solve Kepler's equation. This function implements the algorithm proposed in Markley (1995) Celestial Mechanics and Dynamical Astronomy, 63, 101 ([DOI:10.1007/BF00691917](http://dx.doi.org/10.1007/BF00691917)). This method is not iterative, requires only four transcendental function evaluations, and has been proved to be fast and efficient over the entire range of elliptic motion 0≤e≤10.

On my laptop, this solves for a single eccentric anomaly in just 47 ns.
Since it is implemented in pure Julia, there is no overhead from calling into a C or Cython compiled function and no need for vectorization.

## Choosing an Orbit Type in Octofitter

When defining a `Planet` in Octofitter, specify the orbit type via the `basis` parameter:

| Basis | Use Case |
|-------|----------|
| `Visual{KepOrbit}` | Relative astrometry only |
| `AbsoluteVisual{KepOrbit}` | Astrometry + proper motion anomaly (HGCA, Hipparcos) |
| `ThieleInnesOrbit` | Low-eccentricity or face-on orbits with astrometry |
| `RadialVelocityOrbit` | RV-only fitting |

For generating synthetic data using PlanetOrbits functions, see the [Data Simulation](@ref data-simulation) tutorial.
