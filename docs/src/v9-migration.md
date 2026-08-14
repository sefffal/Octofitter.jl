# [Migrating to Octofitter v9](@id v9-migration)

Octofitter v9 is built on PlanetOrbits v2. The change is not cosmetic: a model
is now a **flat list of bodies plus a flat list of observations**, and the
orbital physics lives entirely in PlanetOrbits. This page is the map from the
old surface to the new one, and the reasoning where it is not obvious.

There is no compatibility shim. Much of the v8 surface (`basis=`, `M`,
per-companion observations) has no meaning in v9, so a silent shim would have
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
v8 picked between them at runtime by comparing semi-major axes, which broke
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
either propagator. The epicyclic superposition v8 wrote out by hand in six
likelihoods is gone.

## Migration table

| v8 | v9 |
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
| `fluxratio` / `flux` on the observation | `flux` / `flux_<band>` on the body |
| `G23HObs(gaia_id=…)` | `G23HObs(gaia_id=…, host=A, companions=(b, c))` |
| `construct_elements(chain, :b, i)` | `construct_system(model, chain, i)` |
| `θ_at_epoch_to_tperi(θ, ep; …)` | *(removed)* — declare `θ` and `epoch` as elements |
| chain column `b_a` | chain column `b_a` (unchanged) |
| chain column `b_GPI_jitter` | chain column `GPI_jitter` |
| chain column `b_SPHERE_flux` | chain column `b_flux_H` |

Angles are radians and epochs MJD, as before. Element keywords are the Unicode
ones (`ω`, `Ω`) with no ASCII aliases.

### Retired names

This is a breaking major release, so old names are **not** aliased to their
replacements — a name that resolves to something with different arguments fails
one line later with a worse message, and holding a name hostage stops it being
reused. Rename them:

| v8 name | v9 |
|---|---|
| `Planet` | [`Body`](@ref), and the observation moves to `System`'s `observations=` |
| `θ_at_epoch_to_tperi` | declare `θ` and `epoch` as orbital elements |
| `PlanetRelAstromObs` | [`RelAstromObs`](@ref)`(tab; target, ref)` |
| `StarAbsoluteRVObs` | [`RadialVelocityObs`](@ref)`(tab; target=A, ref=Barycentre)` |
| `PlanetRelativeRVObs` | [`RadialVelocityObs`](@ref)`(tab; target=b, ref=A)` |
| `MarginalizedStarAbsoluteRVObs` | [`MarginalizedRVObs`](@ref)`(tab; target, ref)` |
| `masspostplot` | [`octocorner`](@ref), or `hist(vec(chain[:b_mass]))` |
| `PhotometryLikelihood` | [`PhotometryObs`](@ref) |
| `PlanetOrderPrior` | [`OrbitOrderPrior`](@ref), taking `Body` nodes or `Symbol`s |
| `ObsPriorAstromONeil2019` | [`ObsPriorONeil2019`](@ref) |
| `InterferometryLikelihood` | [`InterferometryObs`](@ref) |
| `GRAVITYWideKPLikelihood` | [`GRAVITYWideKPObs`](@ref) |
| `HGCAInstantaneousObs`, `GaiaCatalogFitObs` | [`HGCAObs`](@ref) / [`G23HObs`](@ref) |

Four of them are still *defined*, purely so the error message names the
replacement rather than being a bare `UndefVarError`: `Planet`,
`θ_at_epoch_to_tperi`, and the HGCA pair. They are the ones an old script hits
first.

### Two spellings that catch nearly everyone

**`Barycentre` and `Photocentre` are declarations, not queries.** They are
singleton *specs*, valid in `target=`, `ref=` and [`ObservableQuery`](@ref).
They are not points, because their weights come from the bodies' masses or
fluxes — so resolving one needs a solved system:

```julia
radvel(sol, :A, Barycentre)                                # error, with this advice
radvel(sol, :A, Octofitter.resolveref(posys, Barycentre))  # correct
radvel(sol, :A, PlanetOrbits.barycentre(posys))            # equivalent
```

**`$` interpolation works on `=` lines and is rejected on `~` lines.** A derived
expression is quoted and evaluated later, inside the model, where your local
variables are gone — hence the `$`. A prior's right-hand side is evaluated where
you wrote it, so it already sees them:

```julia
@variables begin
    ra    = $(cat.ra)               # derived: interpolate
    pmra  ~ Normal(cat.pmra, 10)    # prior: no `$`, `cat` is already in scope
end
```

Both mistakes used to surface as errors pointing into Octofitter's internals
rather than at the line you wrote.

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

The v9 inverse has no quadrant fixups and returns the `Ω ∈ [0, π)` branch of
the genuine ±180° node degeneracy; v8's `ThieleInnesOrbit` carried two
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
          observing_geometry=:auto,       # per-body LTT, depth scale, LOS projection
          barycentric_lighttime=:auto,    # whole-system light-travel timing
          einstein_rv=:on)                # spectroscopic vs kinematic RV
```

The two opt-outs are **tri-state** and default to `:auto`, which measures
whether each correction moves *your* predictions by anything comparable to
your uncertainties and decides once, at build, printing what it decided.
`:on`/`:off` — and `true`/`false` — still set them by hand. See
[Letting Octofitter decide: `:auto`](@ref corrections) for how the decision is
made and how to re-check it after sampling.

!!! warning "Two radial-velocity changes to know about"
    `radvel` is now the *spectroscopic* velocity: it includes the Einstein
    (second-order Doppler plus gravitational-redshift) term. And an absolute
    RV series now models the perspective-acceleration drift by default —
    `secular_acceleration=:model`. If your pipeline already removed secular
    acceleration, you must say so with `secular_acceleration=:data_corrected`,
    or it will be counted twice. Both are covered, with magnitudes, in
    [How Octofitter Computes Orbits](@ref orbit-computation).

!!! note "Catalog priors describe the primary, not the barycentre"
    The frame's `ra/dec/plx/pmra/pmdec/rv` are the **system barycentre's**.
    A catalog measured the *primary body* (or, for an unresolved pair, the
    photocentre). For a planet host the difference is negligible. For a
    stellar binary it is the full reflex velocity — potentially mas/yr — so a
    tight `pmra ~ Normal(catalog…)` on the frame applies a measurement of body
    A to the barycentre and can actively bias the orbit.

    There is no anchoring option yet. Until there is: for stellar binaries,
    either widen those priors, or fit the catalog's underlying data (DR4
    epochs, Hipparcos IAD) rather than its 5-parameter summary. Note the
    summary is an epoch-average with the orbit aliased into it whichever way
    you parametrize, so a tight catalog prior on an accelerating source is
    doubly suspect.

`AHL21` integrates every body's mutual gravity instead of superposing
Keplerians. Nothing else in the model changes. Note that the elements mean
different things under the two propagators — constant versus osculating at
`t0` — so a chain fitted with one is not element-for-element comparable with
the other.

## Performance

Same 220-transit Gaia DR4-like model, same data, same parameters, same
machine (`examples/gaia-dr4-v2.jl`):

| | v8 ℓπ | v9 ℓπ | v8 ∇ℓπ | v9 ∇ℓπ |
|---|---|---|---|---|
| 1 companion (D=13) | 45.6 µs | **10.7 µs** (4.3×) | 574 µs | **75 µs** (7.6×) |
| 2 companions (D=20) | 87.3 µs | **16.9 µs** (5.2×) | 581 µs | **173 µs** (3.4×) |

measured with `observing_geometry=false`, which selects exactly v8's sky
geometry — the log likelihoods then agree to ~1e-15 relative, which is the
gate `test/v2/v1-regression.jl` holds. The v9 default adds per-body
light-travel retardation, per-body depth scaling and line-of-sight projection,
and costs 13.9 / 22.8 µs and 195 / 509 µs; it is more accurate, not
equivalent. Both are allocation-free.

## G23H — what changed in the port

`G23HObs` is ported (`src/likelihoods/g23h.jl`). Four differences are
user-visible.

**Membership is declared.** The constructor takes `host=` and
`companions=(…)`:

```julia
absastrom = G23HObs(; gaia_id=…, host=A, companions=(b, c, d), …)
```

`companions=` fixes the order the flux-ratio vectors are indexed in — v8
indexed them by position in `system.planets`, which was implicit and easy to
get wrong. The names are validated against the system's bodies at model-build
time.

**Flux ratios come from the bodies, with a per-draw override.** By default
G23H reads each companion's `flux_G` (Gaia G) and `flux_Hp` (Hipparcos Hp)
body variables, like every other v9 observation — give the host `flux_G = 1.0`
and each companion its contrast ratio. A model that declares neither leaves the
companions dark (ratio 0), which is the right answer for a dark companion but
*is* a change from v8, where the vector was mandatory.

The observation-level `fluxratio` / `fluxratio_hip` variables still exist as an
override tier, and they are the only tier that can read deferred system
variables — so a sampled resolved-flag latent can gate a companion out of the
Gaia photocentre without the blending state round-tripping through a body's
flux. When supplied they must be length-`length(companions)` containers; a bare
scalar is accepted only when there is exactly one companion.

**Multiple luminous companions give different numbers, on purpose.** v8
normalized each companion's photocentre term by its own `(1 + f_k)` and
superposed them. The photocentre of several luminous bodies is the
flux-weighted mean of their *apparent positions*, which normalizes by
`(1 + Σ f_k)` — the non-superposition point BINARYS makes. v9 builds one
`WeightedPoint` per draw and asks for `raoff(sol, photocentre, barycentre)`,
which is exact for any number of companions and any hierarchy. With a single
luminous companion the two are algebraically identical, and v9 reproduces v8
to roundoff (gate: `test/v2/g23h.jl`). With two, on a representative draw,
they differ by 0.30 mas of sky path — 1.3% of the signal, and 0.10 mas/yr on
the fitted DR3 proper motion, about 5σ of the catalog uncertainty. There is
no legacy-compatibility mode.

**Absolute-frame propagation goes through the frame observables.** v8's
`propagate_astrom` carried a differential-light-travel term its own source
flagged as double counting the proper motion. It is deleted: the
absolute-frame branch reads `frame_ra`/`frame_dec`/`frame_pmra`/`frame_pmdec`
at the catalog reference epochs, which are part of the observation's declared
epochs and so are solved with everything else.

Two smaller notes. The Hipparcos grating response (the BINARYS atan2
"Hippacentre", its resolution taper and σ inflation) is *not* a photocentre
and stays G23H-internal — it is an instrument response, and lives with the
observation that owns the instrument. And `ueva_mode=:none` is now accepted
alongside `:RUWE`/`:EAN`, for stars whose `sig_AL`/`sig_att_radec`/`sig_cal`
calibration is absent: it drops the `:ueva_dr3` datum and the UEVA-driven DR3
covariance deflation and leaves every other channel untouched.

!!! warning "The Gaia RV channel is not differentiable"
    `Distributions.NoncentralChisq` cannot evaluate its log-density at a
    `ForwardDiff.Dual`, so with the `:rv_dr3` channel active the whole log
    density is `-Inf` under any gradient-based sampler, while the primal path
    stays finite. This is inherited from v8, where the failure was silent; it
    now warns once. Use a derivative-free sampler (the production driver uses
    Pigeons), or `include_rv=false`.

## Observation types, one by one

Every v8 observation type has now crossed to the new surface (or been
deliberately retired). The recurring pattern is the same in each case: the
observation moves from a `Planet` to the `System`, names its `target` and `ref`,
and reads brightness from the bodies rather than from a positionally-indexed
vector of its own.

### Relative astrometry

```julia
PlanetRelAstromObs(dat, name="GPI")                       # v8, attached to a planet
RelAstromObs(dat; target=b, ref=A, name="GPI")            # v9
```

`ref` is a real reference: `ref=Barycentre` and `ref=b` (one companion measured
against another) are now expressible. The observation variables `jitter`,
`platescale` and `northangle` are unchanged.

### Radial velocity

`StarAbsoluteRVObs` and `PlanetRelativeRVObs` are gone; there is one type,
distinguished by its references, and it lives in **core Octofitter** (no
`using OctofitterRadialVelocity` needed):

```julia
RadialVelocityObs(dat; target=A, ref=Barycentre, name="HARPS")  # stellar reflex
RadialVelocityObs(dat; target=b, ref=A,          name="CRIRES") # relative RV
```

!!! warning "`offset` and `jitter` are no longer injected for you"
    v8's `StarAbsoluteRVObs` silently added `offset ~ Uniform(-1000, 1000)` and
    `jitter ~ LogUniform(0.001, 100)` when constructed with no `variables=`
    block. v9 never invents a prior. **A v8 RV model copied across without those
    two lines fits with no zero point and no jitter** — it will not error, it
    will change your answers. This is the single most likely way a v8 RV fit
    changes on migration.

`gaussian_process=` and `trend_function=` are now available on relative RV as
well as absolute. With a Gaussian process, epochs must be strictly increasing:
duplicated epochs are a construction-time error instead of a silent `-Inf`.

`MarginalizedStarAbsoluteRVObs` is renamed [`MarginalizedRVObs`](@ref)
(OctofitterRadialVelocity), now *requires* `target=`, and errors if you declare
an `offset` — the zero point is integrated out analytically, so declaring one
was always a mistake.

The RV-specific Makie extension is gone: the single-draw RV summary lives in core now,
built from [`octoplot`](@ref)'s own generic panels. It is also **renamed** — it draws one
posterior draw, not the posterior, so `rvpostplot` is now [`rvplot`](@ref) and
`rvpostplot_animated` is [`rvplot_animated`](@ref). The old names still work and forward
to the new ones with a deprecation warning.

`octoplot` no longer puts several RV instruments on one panel: each gets its own, with
its data uncalibrated and each draw's model curve carrying that draw's own offset and
trend. `rvplot` is where they share an axis, which its single draw is what makes
legitimate. See [`sharepanel`](@ref).

### Photometry

The flux parameter moved from the observation to the body:

```julia
# v8
PhotometryObs(dat, name="NIRC2", variables=@variables begin flux ~ Uniform(0, 10) end)
# v9
PhotometryObs(dat; target=b, band=:H, name="NIRC2")   # reads b's `flux_H`
```

Consequence worth planning for: two instruments observing the same body in the
same band now share **one** parameter, where v8 gave each its own. The
per-point arithmetic is bit-for-bit unchanged. `target` must be a single body;
`epoch` columns are ignored (photometry forces no orbit solves).

### Images and likelihood maps

```julia
# v8: one ImageObs per planet, attached to that planet
ImageObs(dat, name="SPHERE", variables=@variables begin flux ~ Normal(3.8, 0.5) end)
# v9: ONE observation naming every source in the image
ImageObs(dat; targets=(b, c), ref=A, band=:H, name="SPHERE")
```

with `flux_H` declared on `b` and on `c`. A `flux` variable on the observation
is now a hard error naming the fix. v8's per-planet spelling counted the same
image's background once per companion.

`LogLikelihoodMapObs(dat; target=b, ref=A, name="GRAVITY")` is the same edit
with a single `target`: a precomputed logL surface is one-companion by
construction.

Both types are recentred on `ref`, which need not be the star. Multi-planet
image fits are **not** bit-identical to v8, because the epicyclic superposition
is gone; single-planet fits are (pinned to 1e-12 in the test suite).

### Interferometry

```julia
InterferometryObs(dat; targets=(A, b), ref=A, band=:K, name="NIRISS-AMI")
```

`targets` names every source in the visibility sum, **including the host**,
which is now an ordinary source with its own `flux_K`. Giving the host
`flux_K = 1.0` and the companions their old v8 contrast ratios reproduces v8's
likelihood bit-for-bit. Two gotchas:

* `platescale` is now a **divisor**, matching relative astrometry and images.
  v8's interferometry likelihood multiplied by it. `platescale = 1` is
  unaffected; any other fitted value from v8 must be inverted.
* Fibre coupling changed: each source's throughput is evaluated at its own
  offset from `fiber_pointing` (default `Photocentre(band)`), where v8 applied
  the host-to-photocentre throughput to the companion and left the host at 1.0.
  Numbers move. `fiber_pointing=A` is "fibre on the host".

`GRAVITYWideKPObs` is no longer a type but a preset function returning an
`InterferometryObs`, so `obs isa GRAVITYWideKPObs` no longer compiles.

### Absolute astrometry

[`G23HObs`](@ref) is the joint Gaia + Hipparcos model, and is now the way HGCA
data is fit: [`HGCAObs`](@ref) is a **helper function** returning a `G23HObs`
restricted to the six HGCA channels, with `ueva_mode=:none` and
`include_rv=false`.

```julia
pma = HGCAObs(; gaia_id=…, host=A, companions=(b,), ref=Barycentre)
pma isa G23HObs    # true
```

!!! warning "HGCA fits are not bit-identical to v8"
    Those six channels are now modelled by G23H's epoch-selection model, its
    five-parameter refits and its exact flux-weighted photocentre, rather than by
    the HGCA's own cross-calibration.

`HGCAInstantaneousObs` and `GaiaCatalogFitObs` are **gone** — calling them
raises an error naming the replacement, because the modelling code behind them
was subsumed rather than ported. There is no `N_ave` equivalent; the
instantaneous approximation is what G23H replaced by modelling the actual scan
epochs.

[`HipparcosIADObs`](@ref) is a standalone Hipparcos fit:

```julia
hip = HipparcosIADObs(; hip_id=1475, host=A, companions=(b,), ref=Barycentre)
```

By default `iad_Δplx = 0`, so the system's own `plx` sets the abscissa's
parallax signature — that is what makes a Hipparcos-only fit *measure* the
parallax. Write `iad_Δplx ~ Uniform(-10, 10)` to make it a nuisance instead.
Results are not bit-identical to v8, which compared a compensated absolute-frame
position against data reconstructed with the tangent-plane model; v9 uses the
tangent-plane forward model on both halves.

The `iad_Δ*` block is a general five-parameter frame offset, for astrometry
taken on a frame other than Gaia's. There are deliberately **no** per-source
offsets for Gaia channels: several `G23HObs` in one system share one frame, and
that shared frame is what constrains a wide orbit.

### Prior terms

| v8 | v9 |
|---|---|
| `ObsPriorAstromONeil2019(obs)` | [`ObsPriorONeil2019`](@ref)`(obs)` |
| `PlanetOrderPrior(b, c)` | [`OrbitOrderPrior`](@ref)`(b, c)`, taking `Body` nodes or `Symbol`s |
| `NonCrossingPrior()` | same, plus an optional `bodies=` |
| `LimitClosestApproachAUPrior(…)` | same, plus an optional `bodies=` |
| `HillStabilityPrior()` | same, plus an optional `bodies=` |

All of them go in the **System**'s `observations=` list.

`ObsPriorONeil2019` now needs to know which orbit it applies to. It defaults to
the wrapped observation's `target`, which is right for relative astrometry and
relative RV; a stellar-reflex `RadialVelocityObs(…; target=A, ref=Barycentre)`
must pass `orbit=b` explicitly. `orbit=(b, c)` sums over several orbits. It also
now covers RV data, not just astrometry.

!!! warning "`HillStabilityPrior` is not bit-compatible with v8"
    v8 overwrote the outer planet's parameters with the inner planet's, so the
    outer companion's mass never entered the criterion, and `M★` was a per-planet
    total mass that v9 has no equivalent of. The port implements the intended
    Gladman criterion. `NonCrossingPrior` and `LimitClosestApproachAUPrior` are
    numerically identical to v8.

## Analysis, IO, and other behaviour changes

* **`Octofitter.pointwise_like` results change, and v8's were wrong.** v8 passed
  the *complement* of the wanted rows to `likeobj_from_epoch_subset`, so each v8
  "pointwise" column actually held the likelihood of all the data except that
  point. Any published PSIS-LOO number produced with v8 is affected. Columns for
  prior-shaped terms are also no longer emitted, so the columns now sum exactly
  to the model log-likelihood minus those terms — which is what PSIS-LOO wants.
* **SBC data are now noisy by default.** `Octofitter.calibrationhmc` defaults to
  `add_noise=true`; v8 simulated noiseless data, so its rank histograms were not
  the SBC statistic. Pass `add_noise=false` to reproduce the old behaviour.
* **Completeness `inject` overrides nest under `bodies=`**, not `planets=`
  (which raises an explicit error), and masses are in **M⊙**.
* **Flux-table interpolators take solar masses by default.**
  `sonora_photometry_interpolator`, `sonora_cooling_interpolator` and
  `Octofitter.bhac15_mass_age_interpolator` now default to `mass_unit=:Msol`;
  pass `mass_unit=:Mjup` to keep v8 call sites verbatim. This is silent — a v8
  call passing `12.0` now means 12 M⊙, lands off-grid, and returns `NaN` rather
  than erroring. They return **absolute magnitudes**, so the
  `flux_<band> = 10^(-0.4 * mag)` step is mandatory.
* **`Octofitter.savehdf5(fname, model, chain)` takes three positional
  arguments**, and now writes the `parameter_labels` attribute its own reader
  wants. `loadhdf5` gained `host=` / `bodynames=`, imports companion masses in
  M⊙, no longer synthesises a `<planet>_M` column, and finally works with
  `numchains > 1`.
* **`Octofitter.loadchain(fname; model)`** is the recommended spelling: without
  the model, a chain saved against a different model (or by v8) loads silently
  and yields `missing` downstream. `Octofitter.checkchain(model, chain)` does the
  check on its own. v8 chains still load, with a warning about the
  `<planet>_<obs>_<var>` → `<obs>_<var>` rename.
* **`Octofitter.Whereistheplanet_astrom` now returns its vector** of
  `RelAstromObs` — in v8 the `return` was commented out, so the documented
  destructuring silently unpacked `nothing`. It gained `target=`/`ref=`.
* **`gaia_plx(; gaia_id)`** reads the Gaia DR3 catalog directly instead of the
  retired HGCA data dependency. The numbers agree for any source in both.
* **`getplotdata` and `calibrationplots` are removed.** Both were dead code in v8
  that could not have run.
* **`octofit_pigeons` needs `using Pigeons`.** Its methods live in a package
  extension, so without that import the function exists with no methods and a
  call raises a `MethodError`. This is unchanged from v8 in principle, but worth
  restating because the extension was briefly unported during the v9 work.
* **`GaiaDR4AstromObs` no longer takes `gaia_id=` and no longer carries
  `obs.gaia_sol`.** The observation models a sky path and does not need the
  catalog row. Read it with [`gaia_dr3_solution`](@ref)`(; gaia_id)`, which is
  the same query [`gaia_plx`](@ref) uses and caches to
  `_gaia_dr3_final/source-<id>.csv`.
* **`GaiaDR4AstromObs` ingests `scan_pos_angle` in degrees, not radians.** The
  Gaia archive publishes the DR4 scan angle in degrees (the VOTABLE declares
  `unit="deg"`), so a table read straight out of the archive now needs no unit
  conversion; the conversion to radians happens once, internally, at
  construction. **Delete the `deg2rad.` you applied in v8** — nothing errors if
  you leave it in, you just get a wrong posterior from scan angles compressed
  into a ±3.14° wedge. `scan_pos_angle` is the only column whose unit changed:
  epochs are still MJD, `centroid_pos_al` and `centroid_pos_error_al` still mas,
  `parallax_factor_al` still dimensionless. `obs.table.scan_pos_angle` keeps
  your degrees verbatim; the radian sin/cos the likelihood projects with are
  precomputed as `obs.sinψ`/`obs.cosψ`.
* **`generate_from_params(::RelAstromObs, …)` now carries the `:cor` column
  through and draws correlated noise.** v8 (and the first v9 port) dropped the
  column and drew independent `randn()` per component, so a replicate of
  correlated astrometry was drawn from the wrong distribution *and* was fitted
  with a diagonal covariance. Simulated data, posterior-predictive checks, SBC
  and completeness on correlated astrometry all change.
* **`generate_from_params` on an `ObsPriorONeil2019` returns the wrapper**, not
  the bare inner likelihood. Unwrapping produced a different model — no Jacobian,
  and `obspri_<name>` renamed to `<name>` — so a chain fitted to the replicate
  did not line up with one fitted to the original.
* **`initialize!` / `startingpoints!` locate parameters by layout.** They used to
  find an override's index by searching the flat parameter vector for a matching
  *value*, which picks the wrong slot when two variables draw the same number and
  cannot address one element of a vector-valued prior at all. Overrides also nest
  under `bodies=`; `planets=` raises an explicit error.
* **`Octofitter.savehdf5` derives `tp` when the model did not sample it.** The
  recommended `θ` + `epoch` phase spelling produces no `b_tp` column, and `tp` is
  determined by the elements, so the export rebuilds it per draw rather than
  refusing.
* **A GOST forecast with no scans is an error, not a cache entry.** A transient
  header-only response used to be written to `GOST-<ra>-<dec>-<baseline>.csv` and
  then reused forever, failing far downstream as a `DimensionMismatch`.

## What is not ported yet

`src/legacy/` holds the v8 code that has not moved to the new surface, kept
unmodified so each port is a diff rather than a rewrite.

- **The Dynesty extension**, and the v8 Makie recipe set parked in
  `src/legacy/ext/`. (The **Pigeons** extension *is* ported —
  [`octofit_pigeons`](@ref) works as it did in v8 once you `using Pigeons` — but
  its MPI/cluster path has not been re-verified on a real cluster.)
- **`octoplot`'s `mark_epochs_mjd`.** Every shipped observation type now declares
  `plotchannels`, so `octoplot` overlays all of their data, and the named v8
  figures are back: [`rvplot`](@ref)/[`rvplot_animated`](@ref),
  [`dotplot`](@ref), [`gaiastarplot`](@ref), [`gaiatimeplot`](@ref),
  [`skytrackplot`](@ref) and [`hipparcosplot`](@ref). `hgcaplot`, `pmaplot`,
  `absastromplot`, `astromplot`, `physorbplot` and `rvtimeplot` are not separate
  functions any more — they are `octoplot`'s generic sky, time-series and
  phase-folded panels, which now carry the data those recipes drew. `masspostplot`
  is deliberately dropped. The one v8 view with no v9 equivalent is `hgcaplot`'s
  and `pmaplot`'s proper-motion *vector* panels (μα⋆ against μδ with 1σ crosses);
  the same information is in the `pmra`/`pmdec` panels and their whitened residual
  strips.
- **`gp_predict` for the AbstractGPs backend**, so GP cross-validation works with
  the Celerite backend only (a pre-existing v8 hole, now a clear error).
- **The Gaia NSS helpers** — `initialize_from_nss!`, `query_nss`,
  `nss_to_starting_point` and `nss_to_model_chain`, for seeding a fit from a
  published Non-Single Star solution. Parked in `src/legacy/nss.jl`. Convert an
  NSS solution to orbital elements yourself and pass them to [`initialize!`](@ref)
  or `startingpoints!` in the meantime.

Plotting in general is **not** on that list: `ext/OctofitterMakieExt.jl` is the
new v9 plotting layer, built on the backend-agnostic API in
`src/plotting-api.jl` ([`octoplot`](@ref), [`PosteriorSeries`](@ref),
[`ObservableQuery`](@ref), [`plotchannels`](@ref)/`Octofitter.residuals`,
[`timeseriespanel!`](@ref), [`skypanel!`](@ref)), and [`octocorner`](@ref) works.

Until PlanetOrbits v2 is registered, Octofitter declares no `[compat]` bound
on it and both packages must be `Pkg.develop`ed together.
