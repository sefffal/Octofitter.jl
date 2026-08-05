# Kepler Solver and Orbit Types

The heart of this package is being able to take a set of Keplerian elements and output relative positions, velocities, etc.
For this, we use [PlanetOrbits.jl](https://github.com/sefffal/PlanetOrbits.jl) which adopts the same conventions as Orbitize!.

For a quick summary of the coordinate system and orbital-element conventions ($a$, $e$, $i$, $\omega$, $\Omega$, $t_p$, ...), see [What conventions does Octofitter use for orbital elements?](@ref) in the FAQ.

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

## Choosing a Parameterization

There is no `basis=` keyword in Octofitter any more, and no orbit *types* to choose
between. Two independent choices replace it.

**The frame** is chosen by which frame variables the `System` block defines:

| System block defines | Frame | Observables available |
|---|---|---|
| nothing | none | physical units only (AU, AU/yr, m/s) |
| `plx` | parallax | angular observables in mas — the usual choice for relative astrometry |
| `plx, ra, dec, pmra, pmdec, rv, ref_epoch` | absolute | rigorously propagated space motion — needed for Gaia/Hipparcos absolute astrometry |

A *partial* absolute frame (some of those variables but not all) is an error rather than
a silent downgrade. For an RV-only fit, omit `plx` and fix `i = π/2`, `Ω = 0` — radial
velocities constrain only m·sin i.

**The parameterization** is chosen by which orbital elements a `Body` block declares.
Supply exactly one alternative from each group; supplying two, or none, is a mechanical
error from the constructor rather than a silently ignored keyword.

| Group | Alternatives |
|---|---|
| size | `a` [AU] or `P` [**days**] |
| shape | (`e`, `ω`) or (`secosω`, `sesinω`) or (`ecosω`, `esinω`) |
| phase | `tp` or `M0` + `epoch` or `θ` + `epoch` |
| orientation | `i` [rad], `Ω` [rad] |
| joint | `x, y, z, vx, vy, vz` + `epoch` replaces every group above |

`secosω` = √e·cosω and `ecosω` = e·cosω sample the eccentricity disc rather than the
half-plane, which removes the ω degeneracy as e → 0. `θ` is the planet's sky-plane
position angle at `epoch` and is usually far better constrained by imaging than `tp`.

!!! warning "`P` is in days"
    `P` matches `period(sys)` so the two round-trip. If you think in years, multiply:
    `P = P_years * 365.25`.

!!! note "`τ` is gone"
    v1's `τ ~ UniformCircular(1.0)` needed hidden period and reference-epoch state and has
    no clean meaning under N-body integration. Use `tp`, or `M0` + `epoch`, or
    `θ` + `epoch`.

A Thiele-Innes fit is written as a derived line rather than an orbit type — see
[Fit with a Thiele-Innes Basis](@ref).

## Choosing a Propagator

By default, orbits are superposed Keplerians (`PlanetOrbits.KeplerianApprox()`), which is
what every version of Octofitter did. You can instead integrate the bodies' mutual
gravity, without changing anything else in the model:

```julia
System(...; method=PlanetOrbits.AHL21(h=40.0, t0=57388.0))
```

Note that the elements mean different things under the two propagators — constant versus
osculating at `t0` — so a chain fitted with one is not element-for-element comparable
with the other.

For generating synthetic data using PlanetOrbits functions, see the [Data Simulation](@ref data-simulation) tutorial.
