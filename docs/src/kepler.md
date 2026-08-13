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


### Orbitize!'s `tau`, written out

There is no `τ` element, because a bare `τ` carries hidden period and reference-epoch
state and has no clean meaning under N-body integration. But orbitize!'s
`tau` — the fraction of a period elapsed since a reference epoch — is two ordinary
lines, and it keeps its usual `Uniform(0, 1)` prior:

```julia
b = Body(name="b", about=A, variables=@variables begin
    P ~ LogUniform(365.0, 365_000.0)      # days
    τ ~ Uniform(0, 1)                     # orbitize!'s tau
    tp = τ * P + 58849                    # orbitize!'s default tau_ref_epoch, MJD
    # …e, i, ω, Ω as usual
end)
```

`58849` is orbitize!'s default `tau_ref_epoch`; use whatever your file records, and note
that [`Octofitter.loadhdf5`](@ref) reads it out of the file's attributes for you when
importing an orbitize! chain (see [Compatibility with Orbitize!](@ref compat-orbitize)). `τ` stays in the
chain as a sampled column, so posteriors remain directly comparable with orbitize!'s.

Sampling `τ` with [`UniformCircular`](@ref)`(1.0)` instead of `Uniform(0, 1)` avoids the
hard boundary at the wrap point, at the cost of one extra dimension.

A Thiele-Innes fit is written as a derived line rather than an orbit type — see
[Fit with a Thiele-Innes Basis](@ref).

## Choosing a Propagator

By default, orbits are superposed Keplerians (`PlanetOrbits.KeplerianApprox()`). You can
instead integrate the bodies' mutual gravity, without changing anything else in the
model:

```julia
System(...; method=PlanetOrbits.AHL21(h=40.0, t0=57388.0))
```

`h` is the integrator timestep in days and `t0` is the epoch at which the declared
elements are the *osculating* ones. Pick `h` well below the shortest period in the
system, and `t0` near the middle of your data. The elements mean different things under
the two propagators — constant versus osculating at `t0` — so a chain fitted with one is
not element-for-element comparable with the other. This is the machinery a transit-timing
fit needs; see the
[N-body integration](https://sefffal.github.io/PlanetOrbits.jl/dev/nbody/) page in the
PlanetOrbits manual for the details, and expect a worked TTV tutorial here later.

## Precision opt-outs

Two corrections are applied on every solve and can each be turned off independently.
They gate different physics, so they are separate switches rather than one "fast mode":

| keyword | what it does | when it matters |
|---|---|---|
| `observing_geometry=false` | skips the observing-geometry pass | every term scales with the system's angular extent ρ, so it is negligible for a compact system and is not for a wide one |
| `barycentric_lighttime=false` | skips the barycentric light-travel solve | a whole-system *timing* correction that scales with proximity and proper motion, not with ρ |

They are keywords on `orbitsolve`:

```julia
sols = orbitsolve(posys, epochs; observing_geometry=false, barycentric_lighttime=false)
```

Turning both off recovers the older, cheaper sky path exactly, and is what the
performance comparisons in the migration guide are measured against. Read
[Precision opt-outs](https://sefffal.github.io/PlanetOrbits.jl/dev/precision/) before
setting either in a fit you intend to publish: these change *results*, not just runtime,
and a system whose frame declares no `plx` skips the angular corrections anyway.

For generating synthetic data using PlanetOrbits functions, see the [Data Simulation](@ref data-simulation) tutorial.
