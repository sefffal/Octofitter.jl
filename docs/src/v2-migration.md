# Migrating to Octofitter v2

Octofitter v2 is built on PlanetOrbits v2. The change is not cosmetic: a model
is now a **flat list of bodies plus a flat list of observations**, and the
orbital physics lives entirely in PlanetOrbits. This page is the map from the
old surface to the new one, and the reasoning where it is not obvious.

There is no compatibility shim. Much of the v1 surface (`basis=`, `M`,
per-companion observations) has no meaning in v2, so a silent shim would have
to guess — and guessing wrong changes numbers rather than raising an error.

## The shape of a model

```julia
using Octofitter, Distributions

A = Body(name="A", variables=@variables begin
    mass ~ truncated(Normal(1.0, 0.1), lower=0.2)   # M⊙
end)

b = Body(name="b", about=A, variables=@variables begin
    mass_jup ~ LogUniform(0.1, 100)
    mass = mass_jup * mjup                          # M⊙
    a ~ LogUniform(1, 100)                          # AU
    e ~ Uniform(0, 0.99)
    ω ~ Uniform(0, 2pi)
    i ~ Sine()
    Ω ~ Uniform(0, 2pi)
    tp ~ Uniform(50_000, 60_000)                    # MJD
end)

astrom = RelAstromObs(data; target=b, ref=A, name="GPI")

sys = System(name="HD 12345", bodies=[A, b], observations=[astrom],
    variables=@variables begin
        plx ~ Uniform(1, 100)                       # mas
    end)

model = Octofitter.LogDensityModel(sys)
chain = octofit(model)
```

Three things moved, and each removes a class of mistake.

**The host star is a body.** It has a mass like everything else, so a row's
gravitating mass is the total mass of the bodies it binds — computed, not
declared. `M`, `M_pri`/`M_sec`, `M = M_pri + M_sec*mjup2msol`,
`mass = system.M_sec` and the per-planet `M` bookkeeping all disappear.

!!! warning "Results will shift where the bookkeeping was subtly wrong"
    Automatic per-row masses reproduce what careful users intended, and
    silently fix models where the hand computation was off. The
    `fit-coplanar` tutorial, for instance, gave its inner planet
    `M = M_pri + M_b` where a Jacobi chain wants `M_pri + M_c`.

**The topology is declared, never inferred.** `about=A` (astrocentric) and
`about=(A, b)` (Jacobi) are both spellable and are *materially different
models* under the Keplerian approximation — the rows are the approximation.
v1 picked between them at runtime by comparing semi-major axes, which broke
for crossing orbits, equal `a`, and anything that is not a planet orbiting a
star. Moons, circumbinary planets and 2+2 quadruples now work.

**Observations name their references.** Every observation says what it
observes and what it is measured against, so it needs to reconstruct nothing:

```julia
RelAstromObs(data;      target=b, ref=A)                  # relative astrometry
RelAstromObs(data;      target=c, ref=b)                  # companion vs companion
RadialVelocityObs(rvs;  target=A, ref=Barycentre)         # stellar reflex
RadialVelocityObs(rvs;  target=b, ref=A)                  # relative RV
GaiaDR4AstromObs(scans; target=Photocentre, ref=Barycentre)
```

`raoff(sol, b, A)` is the whole astrometric model, for any hierarchy and
either propagator. The epicyclic superposition v1 wrote out by hand in six
likelihoods is gone.

## Migration table

| v1 | v2 |
|---|---|
| `Planet(name=…, basis=Visual{KepOrbit}, …)` | `Body(name=…, about=…, …)` + `plx` in the system block |
| `basis=AbsoluteVisual{KepOrbit}` | define `plx, ra, dec, pmra, pmdec, rv, ref_epoch` in the system block |
| `basis=RadialVelocityOrbit` | omit `plx`; set `i = π/2`, `Ω = 0` |
| `basis=ThieleInnesOrbit` | derive the elements in the body block (below) |
| `basis=CartesianOrbit`, `FixedPosition` | give `x, y, z, vx, vy, vz, epoch` in the body block |
| `companions=[b, c]` | `bodies=[A, b, c]` — the star is a node |
| `M = 1.2` (system) | `mass ~ …` on each body |
| `mass ~ …` in Mjup | `mass = mass_jup * mjup` — masses are M⊙ throughout |
| `a = cbrt(M*P^2)` | `P ~ …` directly (**days**, matching `period`) |
| `tp = θ_at_epoch_to_tperi(θ, ep; M, e, a, i, ω, Ω)` | `θ ~ …` and `epoch = ep` |
| `τ ~ UniformCircular(1.0)` | `tp`, `M0`+`epoch`, or `θ`+`epoch` |
| `Planet(…, observations=[astrom])` | `System(…, observations=[astrom])` with `target=`/`ref=` |
| `fluxratio` on the observation | `flux` / `flux_<band>` on the body |
| `construct_elements(chain, :b, i)` | `construct_system(model, chain, i)` |
| chain column `b_a` | chain column `b_a` (unchanged) |

Angles are radians and epochs MJD, as before. Element keywords are the Unicode
ones (`ω`, `Ω`) with no ASCII aliases.

## Parametrizations are constructor groups

Supply exactly one alternative from each group; the wrong combination is a
mechanical error rather than a silently ignored keyword.

| group | alternatives |
|---|---|
| size | `a` [AU] or `P` [**days**] |
| shape | (`e`, `ω`) or (`secosω`, `sesinω`) or (`ecosω`, `esinω`) |
| phase | `tp` or `M0`+`epoch` or `θ`+`epoch` |
| orientation | `i`, `Ω` |
| joint | `x,y,z,vx,vy,vz`+`epoch` replaces every group above |

Any other name in a body block is an ordinary local, so intermediates
(`mass_jup`, `logmass`, …) work as they always did. That is also how
Thiele-Innes fits are written now — a derived line rather than an orbit type:

```julia
b = Body(name="b", about=A, variables=@variables begin
    mass = 5mjup
    TI_A ~ Normal(0, 100)      # mas
    TI_B ~ Normal(0, 100)
    TI_F ~ Normal(0, 100)
    TI_G ~ Normal(0, 100)
    ti = PlanetOrbits.ThieleInnes(A=TI_A, B=TI_B, F=TI_F, G=TI_G, plx=system.plx)
    a = ti.a; i = ti.i; ω = ti.ω; Ω = ti.Ω
    e ~ Uniform(0, 0.9)
    tp ~ Uniform(50_000, 60_000)
end)
```

The v2 inverse has no quadrant fixups and returns the `Ω ∈ [0, π)` branch of
the genuine ±180° node degeneracy; v1's `ThieleInnesOrbit` carried two
documented π errors for `Ω ≥ π` and `ω + Ω > 2π`.

## Evaluation order

> Bodies look up and outward; the system block looks down and inward, after
> the fact; siblings never see each other.

1. System block, lines that do not mention a body → visible as `system.*`.
2. Every body's block, in declaration order.
3. System block, lines that *do* mention a body — **deferred automatically**.
4. Observation blocks.
5. Build the `PlanetOrbits.System`, solve once, evaluate the likelihoods.

Shared parameters are *hoisted* to the system block; that is the correct
expression of exact coupling, not a workaround:

```julia
sys = System(name="coplanar", bodies=[A, b, c], variables=@variables begin
    plx ~ Uniform(1, 100)
    i_shared ~ Sine()          # hoisted: one inclination, two users
    mut_inc = b.Ω - c.Ω        # deferred: mentions a body, so it runs last
    mut_inc ~ Normal(0, 0.1)   # …and can then be constrained
end)
```

Deferral needs no new syntax and is detected by a static walk over the stored
expressions. A body reading a deferred system variable is a genuine cycle and
is reported with both sides named.

Sibling references inside a *body* block (`i = b.i + di`) are not supported.
Everything expressible that way is expressible by hoisting.

## Choosing the frame and the propagator

Which frame variables the system block defines chooses the frame. There is no
`basis=`:

| system block defines | frame | observables |
|---|---|---|
| nothing | none | physical units only (AU, AU/yr, m/s) |
| `plx` | parallax | angular observables in mas |
| `plx, ra, dec, pmra, pmdec, rv, ref_epoch` | absolute | rigorous 3D space-motion propagation |

A partial absolute frame is an error, not a silent downgrade.

The propagator, and the two precision opt-outs, are `System` keywords:

```julia
System(…; method=PlanetOrbits.AHL21(h=40.0, t0=57388.0),
          observing_geometry=true,        # per-body LTT, depth scale, LOS projection
          barycentric_lighttime=true)     # whole-system light-travel timing
```

`AHL21` integrates every body's mutual gravity instead of superposing
Keplerians. Nothing else in the model changes. Note that the elements mean
different things under the two propagators — constant versus osculating at
`t0` — so a chain fitted with one is not element-for-element comparable with
the other.

## Performance

Same 220-transit Gaia DR4-like model, same data, same parameters, same
machine (`examples/gaia-dr4-v2.jl`):

| | v1 ℓπ | v2 ℓπ | v1 ∇ℓπ | v2 ∇ℓπ |
|---|---|---|---|---|
| 1 companion (D=13) | 45.6 µs | **10.7 µs** (4.3×) | 574 µs | **75 µs** (7.6×) |
| 2 companions (D=20) | 87.3 µs | **16.9 µs** (5.2×) | 581 µs | **173 µs** (3.4×) |

measured with `observing_geometry=false`, which selects exactly v1's sky
geometry — the log likelihoods then agree to ~1e-15 relative, which is the
gate `test/v2/v1-regression.jl` holds. The v2 default adds per-body
light-travel retardation, per-body depth scaling and line-of-sight projection,
and costs 13.9 / 22.8 µs and 195 / 509 µs; it is more accurate, not
equivalent. Both are allocation-free.

## What is not ported yet

`src/legacy/` holds the v1 code that has not moved to the new surface, kept
unmodified so each port is a diff rather than a rewrite. Its tests are beside
it in `test/legacy/`.

- Likelihoods: HGCA and its linear-fit variant, Hipparcos IAD, G23H, images,
  interferometry, likelihood maps, photometry, and the observable / planet
  order / non-crossing priors.
- Analysis and IO: cross-validation, SBC, completeness, the NSS catalogue
  reader, orbitize!/HDF5 IO, and the plotting-facing helpers in `analysis.jl`.
- The v1 package extensions (the old Makie recipe set, PairPlots, Pigeons,
  Dynesty), parked in `src/legacy/ext/`. Every v1 plot recipe reaches into
  `system.planets`; they are not declared in `Project.toml` while unported.

Plotting itself is **not** on that list any more: `ext/OctofitterMakieExt.jl`
is the new v2 plotting layer, built on the backend-agnostic API in
`src/plotting-api.jl` (`octoplot`, `PosteriorSeries`, `ObservableQuery`,
`plotchannels`/`Octofitter.residuals`, `timeseriespanel!`, `skypanel!`).
v1 plot names not yet reproduced (`hgcaplot`, `pmaplot`, `rvpostplot`'s
phase-folding, `octocorner`, …) return with their observation types' ports.

Porting a likelihood is mostly deletion: replace the per-companion
superposition loop with one `raoff(sol, target, ref)`, take the references as
constructor keywords, and drop the `PlanetObservationContext` /
`SystemObservationContext` distinction for the single `ObsContext`.

Until PlanetOrbits v2 is registered, Octofitter declares no `[compat]` bound
on it and both packages must be `Pkg.develop`ed together.
