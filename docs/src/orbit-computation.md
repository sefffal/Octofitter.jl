# [How Octofitter Computes Orbits](@id orbit-computation)

Every prediction Octofitter makes follows the same path. It takes the orbital
elements of the current sample, works out where each body is at each epoch, and
then turns those positions into the quantity your instrument actually
measured — a separation and position angle, a radial velocity, a Gaia
along-scan abscissa.

This page is about the two dials on that calculation:

* **how the bodies are moved** — superposed Keplerian orbits by default, or a
  full N-body integration;
* **how carefully positions are turned into observables** — a handful of small
  corrections that matter a great deal for some data and are invisible in the
  rest.

The two are really the same question asked twice: *how precise does this
calculation need to be for my data?* A 3 mas imaging point and a Gaia transit
at 30 µas are not asking for the same arithmetic, and neither are a 10 m/s and
a 10 cm/s radial velocity.

!!! tip "If you just want to fit an orbit"
    You do not have to set anything on this page. The defaults are the accurate
    ones, and Octofitter measures which corrections your data can actually feel
    and quietly declines the ones it cannot. Come back when you want to know
    what it decided, when you are working at microarcsecond or cm/s precision,
    or when your planets are close enough to tug on each other.

## The baseline: a Keplerian orbit

By default each row of the hierarchy — each `Body` with an `about=`, each
`Orbit` — follows its own Keplerian orbit, and the orbits are superposed. That
is exact for two bodies, and it is the classic approximation for hierarchical
systems where every companion is much lighter than what it orbits. The
machinery for this lives in
[PlanetOrbits.jl](https://github.com/sefffal/PlanetOrbits.jl), which uses the
same conventions as orbitize!.

For a quick summary of the coordinate system and orbital-element conventions
($a$, $e$, $i$, $\omega$, $\Omega$, $t_p$, ...), see
[What conventions does Octofitter use for orbital elements?](@ref) in the FAQ.
Fuller references live in the PlanetOrbits manual:

* [Introduction](https://sefffal.github.io/PlanetOrbits.jl/dev/introduction/) — creating and solving orbits
* [Coordinate Conventions](https://sefffal.github.io/PlanetOrbits.jl/dev/conventions/) — coordinate system and element definitions
* [API Reference](https://sefffal.github.io/PlanetOrbits.jl/dev/api/) — complete function reference

### Solving Kepler's equation

Going from mean anomaly to eccentric anomaly means solving Kepler's equation,
which has no closed-form solution, so essentially all of the arithmetic in a
Keplerian model goes here. The default solver is a tweaked version of the one
in [AstroLib.jl](http://juliaastro.github.io/AstroLib.jl/stable/ref/#AstroLib.kepler_solver):

> Many different numerical methods exist to solve Kepler's equation. This
> function implements the algorithm proposed in Markley (1995) Celestial
> Mechanics and Dynamical Astronomy, 63, 101
> ([DOI:10.1007/BF00691917](http://dx.doi.org/10.1007/BF00691917)). This method
> is not iterative, requires only four transcendental function evaluations, and
> has been proved to be fast and efficient over the entire range of elliptic
> motion 0≤e≤10.

That comes to about 47 ns for a single eccentric anomaly on a laptop. It is pure
Julia, so there is no overhead from calling into C or Cython, and no need to
vectorize by hand — PlanetOrbits already batches the solve across epochs and
uses SIMD kernels where it can.

Other solvers (including arbitrary-precision ones, and the hyperbolic branch
for unbound orbits) are selectable through the propagator; see
[Kepler solvers](https://sefffal.github.io/PlanetOrbits.jl/dev/kepler/) in the
PlanetOrbits manual. Most people never need to change it.

## [Choosing a parameterization](@id parameterization)

There is no `basis=` keyword in Octofitter any more, and no orbit *types* to
choose between. Two independent choices replace it.

**The frame** is chosen by which frame variables the `System` block defines:

| System block defines | Frame | Observables available |
|---|---|---|
| nothing | none | physical units only (AU, AU/yr, m/s) |
| `plx` | parallax | angular observables in mas — the usual choice for relative astrometry |
| `plx, ra, dec, pmra, pmdec, rv, ref_epoch` | absolute | rigorously propagated space motion — needed for Gaia/Hipparcos absolute astrometry |

A *partial* absolute frame (some of those variables but not all) is an error
rather than a silent downgrade. For an RV-only fit, omit `plx` and fix
`i = π/2`, `Ω = 0` — radial velocities constrain only m·sin i.

**The parameterization** is chosen by which orbital elements a `Body` block
declares. Supply exactly one alternative from each group; supplying two, or
none, is a mechanical error from the constructor rather than a silently ignored
keyword.

| Group | Alternatives |
|---|---|
| size | `a` [AU] or `P` [**days**] |
| shape | (`e`, `ω`) or (`secosω`, `sesinω`) or (`ecosω`, `esinω`) |
| phase | `tp` or `M0` + `epoch` or `θ` + `epoch` |
| orientation | `i` [rad], `Ω` [rad] |
| joint | `x, y, z, vx, vy, vz` + `epoch` replaces every group above |

`secosω` = √e·cosω and `ecosω` = e·cosω sample the eccentricity disc rather
than the half-plane, which removes the ω degeneracy as e → 0. `θ` is the
planet's sky-plane position angle at `epoch` and is usually far better
constrained by imaging than `tp`.

### Orbitize!'s `tau`, written out

There is no `τ` element, because a bare `τ` carries hidden period and
reference-epoch state and has no clean meaning under N-body integration. But
orbitize!'s `tau` — the fraction of a period elapsed since a reference epoch —
is two ordinary lines, and it keeps its usual `Uniform(0, 1)` prior:

```julia
b = Body(name="b", about=A, variables=@variables begin
    P ~ LogUniform(365.0, 365_000.0)      # days
    τ ~ Uniform(0, 1)                     # orbitize!'s tau
    tp = τ * P + 58849                    # orbitize!'s default tau_ref_epoch, MJD
    # …e, i, ω, Ω as usual
end)
```

`58849` is orbitize!'s default `tau_ref_epoch`; use whatever your file records,
and note that [`Octofitter.loadhdf5`](@ref) reads it out of the file's
attributes for you when importing an orbitize! chain (see
[Compatibility with Orbitize!](@ref compat-orbitize)). `τ` stays in the chain as
a sampled column, so posteriors remain directly comparable with orbitize!'s.

Sampling `τ` with [`UniformCircular`](@ref)`(1.0)` instead of `Uniform(0, 1)`
avoids the hard boundary at the wrap point, at the cost of one extra dimension.

A Thiele-Innes fit is written as a derived line rather than an orbit type — see
[Fit with a Thiele-Innes Basis](@ref).

## Why a fit might need more than a plain Kepler orbit

A Keplerian orbit gives you the bodies' positions and velocities relative to
their common barycentre. Comparing those to data means asking a second
question: what does an observer at the solar system barycentre, at this epoch,
actually see? Several small things happen between the two, and they are small
in different ways.

Some intuition first, in plain terms.

**The system is not standing still, and it is not seen from one fixed
direction.** A nearby high-proper-motion star drifts across the sky over your
baseline, so the direction you are projecting the orbit onto slowly rotates. The
same drift swings the line of sight, which turns some of the star's transverse
motion into apparent radial velocity — the *perspective* or *secular*
acceleration. Both effects are tiny for a distant, slow star and quite large for
a nearby fast one.

**Light takes time to arrive, and not the same time for everything.** The
system's distance changes as it moves, so light reaching you tonight left at a
slightly different time than a simple `t_em = t_obs` would suggest. And within a
wide system, a companion on the far side of its orbit is genuinely farther away
than its host, so its light is retarded a little more.

**Clocks in the system do not run at your rate.** A body moving fast, or sitting
deep in its companion's gravitational potential, is redshifted relative to you —
the transverse Doppler shift plus gravitational redshift, together the
*Einstein* term. It shows up in spectroscopic radial velocities.

**Bodies can pull on each other.** Superposed Keplerian orbits cannot express
transit-timing variations, resonant interactions, or secular exchange of
eccentricity between planets. That is the N-body question, at the end of this
page.

### When does each of these matter?

Rough guidance, to be read as "look into this", not as a rule:

| your data | typical σ | what tends to matter |
|---|---|---|
| imaging / relative astrometry | 1–10 mas | nothing here — a plain Kepler orbit is plenty |
| GRAVITY-class relative astrometry | 10–50 µas | observing geometry, for wide or high-proper-motion systems |
| Gaia DR4 / Hipparcos absolute astrometry | 30–100 µas per transit | per point, usually far below σ — but the epochs pile up, and a bias common to all of them grows as √n rather than averaging down (see below) |
| radial velocity | m/s | secular acceleration for a high-proper-motion host; the Einstein term for eccentric stellar or brown-dwarf companions, and for relative RVs |
| radial velocity | cm/s | the above, plus the line-of-sight projection term for wide, high-proper-motion pairs |
| transit timing | seconds | barycentric light travel — and N-body, if the planets interact |

The one that surprises people is the second-to-last row: for a *relative*
radial velocity of a wide, fast-moving pair, one of the "sky geometry"
corrections lands squarely in m/s. More on that below.

## The corrections, one at a time

Here is the whole set, and how each is controlled:

| correction | what it is | how it is controlled |
|---|---|---|
| observing geometry (four terms) | how barycentric positions are *viewed* from the solar system barycentre | `observing_geometry` on `System` — `:auto` by default |
| barycentric light travel | when the light you are seeing was emitted | `barycentric_lighttime` on `System` — `:auto` by default |
| the Einstein term in `radvel` | transverse Doppler + gravitational redshift | always modelled; `einstein_rv=:off` exists only as a counterfactual |
| perspective (secular) acceleration | the line of sight swinging as the star moves | `secular_acceleration=` on each absolute RV series |

The first two are PlanetOrbits' precision opt-outs, and its
[Precision opt-outs](https://sefffal.github.io/PlanetOrbits.jl/dev/precision/)
page carries the full derivations and tables. The magnitudes below are the
short version.

### Observing geometry

Four separate corrections separate a barycentric state from what an observer at
the solar system barycentre sees: rotating into the viewing direction at each
epoch, per-body (differential) light-travel time, a per-body rather than
per-system AU→mas scale, and the line-of-sight projection of velocity. One flag
gates all four.

| term | size | observable |
|---|---|---|
| viewing-direction rotation | 4.85e-3 · ρ[mas] · μ[″/yr] · T[yr] µas | position |
| differential (per-body) light-travel time | 0.099 · ρ[mas] · √(M/a[AU]) µas | position |
| depth scaling | 4.85e-6 · ρ[mas]² µas | position |
| line-of-sight projection | 0.023 · μ[″/yr] · a[AU] m/s | **radial velocity** |

The first three scale with ρ, the angular extent the system actually spans on
the sky. For *absolute* astrometry ρ is the photocentre reflex, not the relative
orbit: a Jupiter analog at 10 pc gives ρ = 0.475 mas, so those three terms come
out at 0.005 / 0.021 / 1e-6 µas — four orders of magnitude below a 30–100 µas
per-epoch precision. For a resolved system observed at GRAVITY+ or CRIRES+
precision, ρ is the full separation, and the same formulas give tens of µas.

!!! note "The fourth term is the one that catches people out"
    It does not scale with ρ, and it is not in µas. Two bodies sit at different
    sky positions, so they are seen along slightly different directions, and a
    velocity common to both — above all the barycentre's transverse space
    velocity — projects differently onto each. That contributes about
    `0.023 · μ[″/yr] · a[AU]` m/s to a **relative** radial velocity: 24 cm/s for
    a Barnard-like host with a companion at 1 AU, and around **239 m/s** at
    1000 AU.

    So if you are fitting relative radial velocities of a high-proper-motion
    system — a wide pair especially — this row is the one to size, not the µas
    table above. No reduction pipeline can have removed it for you, because it
    depends on the barycentre frame and each body's direction, both of which are
    being fitted.

### Barycentric light travel time

The system as a whole moves in three dimensions, so its distance changes, so
light observed at `t_obs` was emitted at a slightly different `t_em`.
`barycentric_lighttime` decides whether that is solved for or `t_em = t_obs` is
used.

The raw shift is large — a tenth of a day out to well over a day — but most of
it is unobservable: its constant part is degenerate with `tp` and its linear
part with the period. What is left is the curvature:

| system | raw shift | after removing the degenerate linear part |
|---|---|---|
| 125 pc, slow | 0.061 d | 0.01 s |
| 41 pc, 0.1″/yr | 0.30 d | 0.04 s |
| 8 pc, 1.7″/yr | 0.73 d | 2.0 s |
| 1.8 pc, 10.4″/yr (Barnard-like) | 1.34 d | 15.9 s |

Translate that into your own observable before worrying about it. For
astrometry a timing error δt displaces a companion by roughly `ρ · 2π · δt/P` —
even the Barnard-like extreme is about 0.002 µas at ρ = 1 mas and P = 1 yr. For
transit timing, δt *is* the observable and seconds are within reach of current
facilities.

Turning this off does not stop the frame being propagated: `frame_pmra`,
`frame_pmdec` and `frame_rv` are still evaluated at every epoch, perspective
acceleration included. That curvature is the absolute-astrometry signal, and
deleting it would change your model rather than merely cheapen it.

### The Einstein term in `radvel`

`radvel` returns the **spectroscopic** radial velocity — what a spectrograph
reports. That is the kinematic projection plus the Einstein term: the
second-order Doppler (time dilation) shift and the gravitational redshift,
differenced between the two references you named. (`velz` is the purely
kinematic quantity, for dynamics.)

There is no provenance keyword for this one, because there is nothing a pipeline
could have asserted about it. The part that varies depends on the sampled orbit
(e, the masses, r(t)), so no reduction can have removed it, and the constant part
is absorbed by your fitted `offset` whether it is modelled or not.

The sizes, briefly. For a **stellar reflex** RV, the varying part is sub-cm/s
with a planetary companion — but reaches the m/s level, phase-locked to the
orbit, for a stellar or brown-dwarf companion on an eccentric orbit. For a
**relative** RV the roles swap, since the companion sits deep in the star's
potential: a Jupiter at 1 AU with e = 0.4 varies by 5.6 m/s over its orbit. That
last case is the one most likely to matter, and the one v1 results are most
likely to shift on. The full tables are on PlanetOrbits'
[Precision opt-outs](https://sefffal.github.io/PlanetOrbits.jl/dev/precision/)
page, and the correction report prints the number for your own series.

`System` does carry `einstein_rv=:off`, which predicts every radial velocity
from the kinematic `velz` instead. It is a counterfactual — a way to measure
what the term does to a posterior — not a statement about your data, which is
why it sits on the system beside the other solve settings rather than on an
observation. It is deliberately *not* tri-state: `:auto` would mean "let the
data decide whether to model it", and the data cannot, since no reduction can
have removed the part that varies with the orbit. Leaving it alone is the right
default.

### Perspective acceleration, and who already removed what

This one is different from the others, because it is the only place where
Octofitter genuinely cannot know the answer and has to ask you.

A radial velocity or an astrometric position is never the raw thing a telescope
recorded. Somewhere between the detector and your table, a pipeline removed the
Earth's motion, and possibly more. Octofitter models whatever the pipeline did
*not* remove. Getting that boundary wrong double-counts a term or drops one, and
in both cases quietly.

The organising principle is *what a term depends on*:

| term | who handles it | what you declare |
|---|---|---|
| Rømer delay to the SSB, aberration, Earth-motion Einstein term, annual `v_t·ϖ` | **the pipeline**, always assumed done | nothing |
| perspective (secular) acceleration | **either** | `secular_acceleration=` per RV series |
| source-side observing geometry, barycentric light travel | **the model** | `observing_geometry=`, `barycentric_lighttime=` (both `:auto`) |
| the `radvel` Einstein term | **the model**, unconditionally | nothing — no pipeline can have removed it |

Earth-side terms depend only on the Earth ephemeris and the catalog direction —
no fit parameters at all. A pipeline can therefore remove them once, perfectly,
and every pipeline does, so Octofitter models none of them: the epochs you supply
are barycentric and the velocities are already referred to the solar system
barycentre. Source-side terms depend on the parameters being sampled, so nobody
could have removed them for you.

Perspective acceleration sits right on that boundary, and that is why it needs a
declaration. It depends only on the source's own kinematics (μ, ϖ, rv) — so a
pipeline *can* remove it, using frozen catalog values, and many do. But those
are exactly the quantities Octofitter samples, especially in a joint fit with
absolute astrometry. Only you know which happened:

```julia
RadialVelocityObs(tab; target=A, ref=Barycentre, name="HARPS",
    secular_acceleration = :model,          # the default
    variables = @variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
    end)
```

* `:model` — Octofitter adds the propagated drift to the prediction. Its
  constant part washes into `offset`; the drift and curvature are the payload —
  roughly **4.5 m/s/yr** for a Barnard-class star, **3 cm/s/yr** for a generic
  30 km/s star at 30 pc.
* `:data_corrected` — your pipeline already subtracted it, so the model adds
  nothing.

`:model` is the default for two reasons: its signature — a smooth trend with a
little curvature — is degenerate with a long-period companion, and it is what
makes a joint RV + absolute-astrometry fit self-consistent, since the same
sampled (μ, ϖ) then predict both data types.

The term is definitionally zero unless the system has a full absolute frame (all
of `ra`, `dec`, `plx`, `pmra`, `pmdec`, `rv`, `ref_epoch`). With `plx` alone,
`:model` costs nothing and does nothing; the correction report below is what
tells you a high-proper-motion target would like the fuller frame.

!!! warning "A stitched series of unknown provenance"
    The declaration is per observation, because provenance is per instrument. If
    you genuinely do not know whether a series was corrected, a reasonable
    fallback is `:data_corrected` plus a fitted per-series `trend_function` to
    absorb the difference — or re-reduce the data. Averaging the two options is
    the one thing that does not work.

!!! warning "Wide-pair boundness tests"
    Two *absolute* series, one per star, each carry their own declaration — and
    there is a trap here. A pipeline that removed secular acceleration per star
    using **that star's own catalog solution** used a proper motion contaminated
    by orbital motion, if the pair is bound. `:data_corrected` then bakes part of
    the signal you are testing for into the subtraction.

    The clean workflow is RVs with *no* secular-acceleration removal, plus
    `:model` on both series: the shared sampled frame and the per-body
    line-of-sight projection then model everything self-consistently.

#### Why a relative RV series rejects the keyword

A relative (body-vs-body) series rejects `secular_acceleration` at
construction, and the error says more than "it cancels" — because that is only
half true:

* What *does* cancel is the **common-mode** barycentre drift, which both bodies
  share. That is the only part this keyword governs, so there is nothing left to
  declare.
* What does **not** cancel is the **differential** term — the line-of-sight
  projection from the observing-geometry table above, `≈ 0.023 · μ["/yr] ·
  a[AU]` m/s. It is already modelled, from the sampled frame velocity and each
  body's fitted direction, which is also why no pipeline can have removed it.

So for a high-proper-motion relative-RV target, the thing to check is not this
keyword but that `observing_geometry` is on. `:auto` keeps it on for data like
that, because the change in prediction is far larger than σ.

## [Letting Octofitter decide: `:auto`](@id corrections)

```julia
System(name="HD 12345", bodies=[A, b], observations=[astrom, rvs],
       observing_geometry=:auto,        # the default
       barycentric_lighttime=:auto,     # the default
       variables=sysvars)               # your @variables block, kept last
```

Both flags also take `:on` and `:off` (and `true`/`false`) if you would rather
set them yourself. `:auto`, the default, measures whether they matter to *your*
data and decides.

### How `:auto` decides

Not by a rule of thumb. The ρ-scaling formulas above were derived for
star-plus-planet topologies, and v2's whole point is hierarchies — moons, 2+2
quadruples, N-body — where "the" ρ is not well defined. So `:auto` mechanizes
the honest recipe instead: solve both ways and compare against your own
uncertainties.

1. Draw ~300 parameter sets from the model priors, with a recorded seed. Draws
   that fail to build are skipped; if too few succeed — a third of the requested
   draws, capped at 100, so 100 of the default 300 — the corrections stay on and
   you get a warning.
2. Solve each draw with the flag on and off, and ask **each observation** for the
   change in *its own* model predictions — through the same `simulate` machinery
   the likelihood uses, so transit timing, astrometry and RV each automatically
   test the quantity they fit.
3. Per observation, the impact is the 99th percentile over draws of
   `max_epochs |Δ prediction| / σ`, where σ is that observation's own tightest
   uncertainty.
4. The flag resolves **off** only if every observation's *accumulated bias* stays
   under **0.1σ**. Needs are OR'd: µas astrometry keeps `observing_geometry` on
   regardless of what the RVs say.
5. A model with no observations — or none whose type can report predictions —
   resolves both **on**, without taking a draw.

Cost is a few hundred solves at ~10 µs each: negligible at build, and skipped
entirely when the answer is already determined.

### Why the threshold is on accumulated bias, not per point

The obvious criterion — "the correction must be some large factor below my error
bar" — turns out to be the wrong shape, because these corrections are
**common-mode**: they push every epoch the same way. A random per-point error of
`b·σ` averages down as `1/√n`; a *systematic* one does not, and shifts a fitted
parameter by roughly `b·√n·σ`. So a 0.001σ per-point bias is invisible in ten
points and a 0.3σ shift in 10⁵ Gaia transits, and no single flat margin can be
right for both.

`:auto` therefore tests `impact × √n` against a single 0.1σ limit — the "you
would not see this in a corner plot" line. The effective per-point margin comes
out around 100× for a hundred-point dataset and 300× for a DR4-scale one,
automatically, rather than from a number somebody picked. The build log prints
both numbers so you can check the arithmetic.

The decision is made **once, at build**, and what the model stores is the
resolved `Bool`. It is deliberately not re-decided per draw: that would make the
likelihood and its gradient discontinuous across parameter space, which HMC pays
for in divergences.

Three knobs on `System` tune the test if you need them — `correction_draws`,
`correction_seed` and `correction_threshold` — and `verbosity=0` silences the
log.

!!! note "An observation type that cannot answer"
    `correction_impact` has a conservative default — "cannot say", which resolves
    the flag **on**. So `:auto` can only ever make a model cheaper than the old
    always-on default, never less accurate.

    Types declare their capability with `has_correction_impact`, which is checked
    *before* any draws are taken: if nothing in your model can answer, both flags
    resolve `:on` immediately and the log names the types responsible. Today
    `RelAstromObs`, `RadialVelocityObs`, `GaiaDR4AstromObs`, `HipparcosIADObs`
    and `PhotometryObs` answer (and an observable prior inherits the answer of
    the likelihood it wraps); G23H, images and interferometry do not yet.

### What it tells you

One line per decision, at `verbosity >= 1`:

```
[HD 12345] observing_geometry = false (auto): worst accumulated bias 2.9e-4σ,
           345× inside the 0.1σ limit — over 300 prior draws (seed 0x...)
[HD 12345] "HARPS": Einstein term, peak-to-peak [m/s] ≈ 0.0027 (0.0009× its tightest σ = 3.0)
```

The second kind of line is an *advisory*: a magnitude worth knowing, with no
decision attached. Every RV series gets one for the Einstein term, and one for
the secular-acceleration drift when that is being modelled.

There is also a [`CorrectionReport`](@ref Octofitter.CorrectionReport) on the
built system (`model.system.corrections`), and it travels into the chain's
metadata. A posterior is not really separable from the corrections that produced
it, and two models under comparison should be able to show their verdicts side
by side.

`chain.info` carries three entries: `corrections` (the report), and
`observing_geometry` / `barycentric_lighttime` (the resolved `Bool`s). Only the
last two survive `savechain`/`loadchain` — the report does not fit in a FITS
header card — so a chain loaded from disk tells you *what* was used but not the
impact numbers behind it. Keep the script that built the model if you want those.

### What it looks like on real data

A one-companion Hipparcos IAD model of HIP 1475 — 145 real transits, ϖ = 279
mas, with genuinely broad priors (`a ~ LogUniform(0.5, 50)`,
`e ~ Uniform(0, 0.9)`, isotropic orientation):

```
observing_geometry    = false (auto): worst accumulated bias 0.0883σ,
                        1.13× inside the 0.1σ limit — over 300 prior draws
barycentric_lighttime = false (auto): changes no prediction at all
```

That is worth reading twice. The observing-geometry verdict is `off` by 13% —
a *marginal* call on a nearby star with a wide-companion prior, not a foregone
one, which is what a threshold in about the right place looks like. The
light-travel verdict is exactly zero because this model declares only `plx`:
without a full absolute frame there is no 3D motion for that correction to act
on. Resolving both cost 70 ms on top of a 60 ms build.

### Re-checking after sampling

```julia
recheck_corrections(model, chain)
```

This re-runs the identical test on draws from the **target** chain. (Target only:
tempered chains sit closer to the prior, and would re-trigger exactly the
broad-prior inclusions the posterior has ruled out.)

* **De-escalation — a hint, and the common case.** A flag resolved *on* because
  broad priors covered parameter space the posterior never visits, and it changed
  nothing. Set it `:off` next time and sample faster.
* **Escalation — a warning, and expected to be rare.** A flag resolved *off*
  under the priors would resolve *on* under the posterior: it concentrated
  somewhere the correction matters, the results may be biased, and the fit is
  worth re-running with the flag on.

### Setting them by hand

`:on` and `:off` (or `true`/`false`) skip the measurement and are recorded in
the report as `(user)` — the build log still prints them, because "the user
asked for it" is part of the answer to "what physics produced this posterior?".

Turning both off recovers the older, cheaper sky path exactly, and that is what
the performance comparisons in the migration guide are measured against. Do read
PlanetOrbits'
[Precision opt-outs](https://sefffal.github.io/PlanetOrbits.jl/dev/precision/)
page before setting either by hand in a fit you intend to publish: these change
*results*, not just runtime. A system whose frame declares no `plx` skips the
angular corrections anyway.

One level down, in PlanetOrbits itself, the same two settings are keywords on
`orbitsolve` rather than tri-state fields — PlanetOrbits cannot see how precise
your data is, so there is nothing there for `:auto` to measure:

```julia
sols = orbitsolve(posys, epochs; observing_geometry=false, barycentric_lighttime=false)
```

## When to reach for N-body instead

Superposed Keplerians cannot express anything that arises from the bodies
actually pulling on each other. If your model has planets close enough to
interact — transit-timing variations, mean-motion resonances, a moon coupled to
its planet and star, secular exchange of eccentricity — integrate their mutual
gravity instead. It is one keyword, and nothing else in the model changes:

```julia
System(name="HD 12345", bodies=[A, b, c], observations=[ttvs],
       method=PlanetOrbits.AHL21(h=40.0, t0=57388.0),
       variables=sysvars)               # your @variables block, kept last
```

`h` is the integrator timestep in days and `t0` is the epoch at which the
declared elements are the *osculating* ones. Guidance for the timestep is
`h ≲ P_min/20` for the shortest period in the system; past that the error grows
as `h²` and a warning is emitted. Put `t0` near the middle of your data — it
defaults to the frame's `ref_epoch` when the system has a full absolute frame,
and must be given explicitly otherwise.

Two things to keep in mind. Cost is real: the integrator does one map evaluation
per `h` of timespan covered, so a compact system with a short inner period is the
expensive case, while wide imaging systems with periods of years need
proportionally few steps. And the elements mean different things under the two
propagators — constants of the motion versus osculating at `t0` — so a chain
fitted with one is **not** element-for-element comparable with the other.
Re-evaluate rather than translate.

This is the machinery a transit-timing fit needs. See
[N-body integration](https://sefffal.github.io/PlanetOrbits.jl/dev/nbody/) in the
PlanetOrbits manual for the details (and please cite Agol, Hernandez & Langford
2021 when you publish with it); a worked TTV tutorial is planned here.

## Under the hood: where the observer's position enters

PlanetOrbits can compute an apparent position from an explicit observer position
(the observer-aware [`raoff`](https://sefffal.github.io/PlanetOrbits.jl/dev/api/)
and `decoff`), which is what turns on the annual–orbital (Kopeikin) coupling and
exact parallax factors. Ephemerides stay Octofitter's concern — PlanetOrbits
takes positions, not spacecraft.

The invariant every absolute-astrometry observation type obeys: it either

1. consumes **SSB observables plus explicit parallax factors** (from the catalog
   data itself — Gaia's `parallax_factor_al`, Hipparcos' IAD column), or
2. consumes **observer-aware observables and no factors**,

never both, and which one it does is literal code in the type, reviewable per
type.

Today every type is mode (i), and that is the right choice for them: the coupling
the factors omit is `≈ 4.85 · z[AU]/d[pc]²` µas per AU of observer displacement,
which is sub-µas for any Gaia DR4 or Hipparcos target. Migrating a type to mode
(ii) is a separate change, gated on validating the exact-geometry projection
against the existing factor math at first order — a check PlanetOrbits' own suite
now carries. Nothing here is a silent fallback: a likelihood that needs observer
positions and has no source for them errors at construction rather than
degrading at solve time.

## Coming from v1

* If you followed v1's advice to *un-remove* secular acceleration from your RV
  data, the new default `:model` gives you identical behaviour — that is what
  the advice was working around.
* If your data **is** pipeline-corrected, you now need to say so:
  `secular_acceleration=:data_corrected`. v1 modelled the term; v2's first
  drafts did not; the default now does again.
* `radvel` gained the Einstein term. See the sizes above and PlanetOrbits'
  migration guide — an ordinary reflex-RV planet fit will not notice, a
  relative-RV fit of an eccentric companion may.
* `observing_geometry=true` / `false` still work, and are recorded as `(user)`.

For generating synthetic data using PlanetOrbits functions, see the
[Data Simulation](@ref data-simulation) tutorial.
