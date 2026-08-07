# Corrections and data provenance

A radial velocity or an astrometric position is never the raw thing a
telescope recorded. Somewhere between the detector and your table, a pipeline
removed the Earth's motion, and possibly more. Octofitter models whatever the
pipeline did *not* remove — and getting that boundary wrong double-counts a
term or omits one, in both cases silently.

This page is the map of that boundary: who owns what, what you have to declare,
and what Octofitter decides for you.

## Who owns what

| term | owner | you declare |
|---|---|---|
| Rømer delay to the SSB, aberration, Earth-motion Einstein term, annual `v_t·ϖ` | **the pipeline**, always assumed done | nothing |
| perspective (secular) acceleration | **either** | `secular_acceleration=` per RV series |
| source-side observing geometry, barycentric light-travel | **the model** | `observing_geometry=`, `barycentric_lighttime=` (default `:auto`) |
| the `radvel` Einstein term | **the model**, unconditionally | nothing — it cannot be pipeline-removed |

The organising principle is *what the term depends on*.

**Earth-side terms depend only on the Earth ephemeris and the catalog
direction** — zero fit parameters. A pipeline can therefore remove them once,
perfectly, and every pipeline does. Octofitter assumes that and models none of
them: the epochs you supply are barycentric, and the velocities are already
referred to the solar-system barycentre. Modelling them here would
double-count.

**Source-side terms depend on the parameters being sampled.** Nobody can have
removed them for you, because they are not knowable until the orbit is.

**Perspective acceleration sits on the boundary, and that is why it needs a
declaration.** It depends only on the source's own kinematics (μ, ϖ, rv) — so a
pipeline *can* remove it, using frozen catalog values, and many do. But those
are exactly the quantities Octofitter samples, especially in a joint fit with
absolute astrometry. Only you know which happened.

## Perspective acceleration: `secular_acceleration`

```julia
RadialVelocityObs(tab; target=A, ref=Barycentre, name="HARPS",
    secular_acceleration = :model,          # the default
    variables = @variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
    end)
```

- `:model` — Octofitter adds the propagated drift, `frame_rv(sol) − rv`, to the
  prediction. Its constant part washes into `offset`; the drift and curvature
  are the payload — roughly **4.5 m/s/yr** for a Barnard-class star, **3 cm/s/yr**
  for a generic 30 km/s star at 30 pc.
- `:data_corrected` — your pipeline already subtracted it. The model adds
  nothing, which reproduces Octofitter v2's earlier behaviour bit-for-bit.

`:model` is the default because its signature — a smooth trend with slight
curvature — is degenerate with a long-period companion, and because it makes a
joint RV + absolute-astrometry fit self-consistent: the same sampled (μ, ϖ)
predict both data types.

The term is definitionally zero unless the system has a full absolute frame
(all of `ra`, `dec`, `plx`, `pmra`, `pmdec`, `rv`, `ref_epoch`). With `plx`
alone, `:model` costs nothing and does nothing; the correction report below is
what tells you a high-proper-motion target wants the fuller frame.

!!! warning "Stitched series of unknown provenance"
    The declaration is per observation because provenance is per instrument.
    If you genuinely do not know whether a series was corrected, declare
    `:data_corrected` and add a fitted per-series `trend_function` to absorb
    the difference — or re-reduce. Do not average the two.

!!! warning "Wide-pair boundness tests"
    Two *absolute* series, one per star, each carry their own declaration —
    and there is a trap. A pipeline that removed secular acceleration per star
    using **that star's own catalog solution** used a proper motion
    contaminated by orbital motion, if the pair is bound. `:data_corrected`
    then bakes part of the signal you are testing for into the subtraction.

    The clean workflow is RVs with *no* secular-acceleration removal, plus
    `:model` on both series: the shared sampled frame and the per-body
    line-of-sight projection then model everything self-consistently.

### It is rejected for relative RV, on purpose

A relative (body-vs-body) series rejects the keyword at construction, and the
error explains why in more detail than "it cancels" — because that is only
half true:

- What *does* cancel is the **common-mode** barycentre drift, which both bodies
  share. That is the only part this keyword governs, so there is nothing left
  to declare.
- What does **not** cancel is the **differential** term: two bodies sit at
  different sky positions, so the barycentre's transverse space velocity
  projects differently onto each — `≈ 0.023 · μ["/yr] · a[AU]` m/s, i.e. 24 cm/s
  for a Barnard-like host with a 1 AU companion and **~239 m/s** at 1000 AU. It
  is already modelled, as the fourth correction of PlanetOrbits' observing
  geometry, computed from the sampled frame velocity and each body's fitted
  direction. That is also why no pipeline can have removed it.

So for a high-proper-motion relative-RV target, the thing to check is not this
keyword but that `observing_geometry` is on. `:auto` keeps it on for data like
that, because Δprediction ≫ σ by a wide margin.

## The Einstein term: not declarable

`radvel` returns the spectroscopic radial velocity — the kinematic projection
plus the second-order Doppler and gravitational-redshift difference between its
two references. There is no provenance keyword for it, because there is nothing
a pipeline could have asserted: the varying part depends on the sampled orbit
(e, the masses, r(t)), and the constant part is absorbed by `offset` either way.

Sizes, in brief (PlanetOrbits' "Precision opt-outs" page has the full tables):
sub-cm/s of *variation* for a stellar reflex with a planetary companion, but
several m/s for the relative RV of a close-in eccentric one. The correction
report prints the number for your own series.

`System` does carry `einstein_rv=:off`, which predicts every radial velocity
from the kinematic `velz` instead. It is a counterfactual, for measuring what
the term does to a posterior — not a declaration about your data, which is why
it sits on the system beside the other solve settings rather than on an
observation beside `secular_acceleration`. It is deliberately **not**
tri-state: `:auto` would mean "let the data decide whether to model it", and
the data cannot, since no reduction can have removed the part that varies with
the orbit. Leave it alone otherwise.

## `observing_geometry` and `barycentric_lighttime`: `:auto`

```julia
System(name="HD 12345", bodies=[A, b], observations=[astrom, rvs],
       observing_geometry=:auto,        # the default
       barycentric_lighttime=:auto,     # the default
       variables=@variables begin … end)
```

These are PlanetOrbits' two precision opt-outs, which gate genuinely different
corrections (see its [Precision opt-outs](https://sefffal.github.io/PlanetOrbits.jl/dev/precision/)
page). `:on` and `:off` — or `true`/`false` — set them yourself. `:auto`, the
default, measures whether they matter to *your* data and decides.

### How `:auto` decides

Not by a rule of thumb. The ρ-scaling formulas were derived for star-plus-planet
topologies, and v2's whole point is hierarchies — moons, 2+2 quadruples, N-body
— where "the" ρ is not well defined. So `:auto` mechanizes the manual's own
recipe instead: solve both ways and compare against your own uncertainties.

1. Draw ~300 parameter sets from the model priors, with a recorded seed. Draws
   that fail to build are skipped; if fewer than 100 succeed, the corrections
   stay on and you get a warning.
2. Solve each draw with the flag on and off, and ask **each observation** for
   the change in *its own* model predictions — the same `simulate` machinery
   the likelihood uses, so transit timing, astrometry and RV each automatically
   test the quantity they fit.
3. Per observation, the impact is the 99th percentile over draws of
   `max_epochs |Δ prediction| / σ`.
4. The flag resolves **off** only if every observation's *accumulated bias*
   stays under **0.1σ**. Needs are OR'd: µas astrometry forces
   `observing_geometry` on regardless of what the RVs say.
5. A model with no observations — or none whose type can report predictions —
   resolves both **on**, without taking a draw.

Cost is a few hundred solves at ~10 µs each: negligible at build, and skipped
entirely when the answer is already determined.

### Why the threshold is on accumulated bias, not per point

The obvious criterion — "the correction must be some large factor below your
error bar" — is the wrong shape, because these corrections are **common-mode**:
they push every epoch the same way. A random per-point error of `b·σ` averages
down as `1/√n`; a *systematic* one does not, and shifts a fitted parameter by
roughly `b·√n·σ`. So a 0.001σ per-point bias is invisible in ten points and a
0.3σ shift in 10⁵ Gaia transits, and one flat margin cannot be right for both.

`:auto` therefore tests `impact × √n` against a single 0.1σ limit — the "you
would not see this in a corner plot" line. The effective per-point margin comes
out around 100× for a hundred-point dataset and 300× for a DR4-scale one,
automatically, rather than from a number somebody picked. The build log prints
both numbers so you can check the arithmetic.

The decision is made **once, at build**, and what the model stores is the
resolved `Bool`. It is deliberately not re-decided per draw: that would make
the likelihood and its gradient discontinuous across parameter space, which HMC
pays for in divergences.

!!! note "An observation type that cannot answer"
    `correction_impact` has a conservative default — "cannot say", which
    resolves the flag **on**. So `:auto` can only ever make a model cheaper
    than the old always-on default, never less accurate.

    Types declare their capability with `has_correction_impact`, which is
    checked *before* any draws are taken: if nothing in your model can answer,
    both flags resolve `:on` immediately and the log names the types
    responsible. Today `RelAstromObs`, `RadialVelocityObs`,
    `GaiaDR4AstromObs`, `HipparcosIADObs` and `PhotometryObs` answer; G23H,
    images and interferometry do not yet.

### What it tells you

One line per decision, at `verbosity >= 1`:

```
[HD 12345] observing_geometry = false (auto): worst accumulated bias 2.9e-4σ,
           345× inside the 0.1σ limit — over 300 prior draws (seed 0x...)
[HD 12345] "HARPS": Einstein term, peak-to-peak [m/s] ≈ 0.0027 (0.0009× its tightest σ = 3.0)
```

plus a [`CorrectionReport`](@ref Octofitter.CorrectionReport) on the built
system (`model.system.corrections`) that travels into the chain's metadata.
A posterior is not separable from the corrections that produced it, and two
models under comparison show their verdicts side by side.

`chain.info` carries three entries: `corrections` (the report), and
`observing_geometry` / `barycentric_lighttime` (the resolved `Bool`s). Only the
last two survive `savechain`/`loadchain` — a FITS header card cannot hold the
report — so a chain loaded from disk tells you *what* was used but not the
impact numbers behind it. Keep the script that built the model if you need
those.

### What it looks like on real data

A one-companion Hipparcos IAD model of HIP 1475 — 145 real transits, ϖ = 279
mas, with genuinely broad priors (`a ~ LogUniform(0.5, 50)`, `e ~ Uniform(0,
0.9)`, isotropic orientation):

```
observing_geometry    = false (auto): worst accumulated bias 0.0883σ,
                        1.13× inside the 0.1σ limit — over 300 prior draws
barycentric_lighttime = false (auto): changes no prediction at all
```

Worth reading twice. The observing-geometry verdict is `off` by 13% — that is
a *marginal* call on a nearby star with a wide-companion prior, not a
foregone one, which is what a threshold in the right place looks like. The
light-travel verdict is exactly zero because this model declares only `plx`:
without a full absolute frame there is no 3D motion for the correction to act
on. Resolving both cost 70 ms on top of a 60 ms build.

### Re-checking after sampling

```julia
recheck_corrections(model, chain)
```

Re-runs the identical test on draws from the **target** chain. (Target only:
tempered chains sit closer to the prior and would re-trigger exactly the
broad-prior inclusions the posterior has ruled out.)

- **De-escalation — a hint, and the common case.** A flag resolved *on* because
  broad priors covered parameter space the posterior never visits, and it
  changed nothing. Set it `:off` next time and sample faster.
- **Escalation — a warning, and expected to be rare.** A flag resolved *off*
  under the priors would resolve *on* under the posterior: it concentrated
  somewhere the correction matters, the results may be biased, and the fit
  should be re-run with the flag on.

## Where the observer goes: one layer, never two

PlanetOrbits can compute an apparent position from an explicit observer
position (the observer-aware [`raoff`](https://sefffal.github.io/PlanetOrbits.jl/dev/api/)
and `decoff`), which is what turns on the annual–orbital (Kopeikin) coupling
and exact parallax factors. Ephemerides stay Octofitter's concern —
PlanetOrbits takes positions, not spacecraft.

The invariant every absolute-astrometry observation type obeys: it either

1. consumes **SSB observables plus explicit parallax factors** (from the
   catalog data itself — Gaia's `parallax_factor_al`, Hipparcos' IAD column),
   or
2. consumes **observer-aware observables and no factors**,

never both, and which one is literal code in the type, reviewable per type.

Today every type is mode (i), and that is correct for them: the coupling the
factors omit is `≈ 4.85 · z[AU]/d[pc]²` µas per AU of observer displacement,
which is sub-µas for any Gaia DR4 or Hipparcos target. Migrating a type to
mode (ii) is a separate change, gated on validating the exact-geometry
projection against the existing factor math at first order — a check
PlanetOrbits' own suite now carries. Nothing here is a silent fallback: a
likelihood that needs observer positions and has no source for them errors at
construction rather than degrading at solve time.

## Migrating from v1

- If you followed v1's advice to *un-remove* secular acceleration from your RV
  data, the new default `:model` gives you identical behaviour — that is what
  the advice was working around.
- If your data **is** pipeline-corrected, you must now say so:
  `secular_acceleration=:data_corrected`. v1 modelled the term; v2's first
  drafts did not; the default now does again.
- `radvel` gained the Einstein term. See the sizes above and PlanetOrbits'
  migration guide — an ordinary reflex-RV planet fit will not notice, a
  relative-RV fit of an eccentric companion may.
- `observing_geometry=true` / `false` still work, and are recorded as `(user)`.
- The Einstein switch moved: it is `System(…; einstein_rv=:off)`, not a
  keyword on `RadialVelocityObs`.
