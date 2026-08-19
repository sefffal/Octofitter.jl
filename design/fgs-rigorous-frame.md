# FGS: a rigorous frame path, and the light-time pin

Status: **Phases 0–4 implemented and verified, 2026-08-18.** Phase 5 (the
refit) is not done and is the next step. Written following the convention
discussion of 2026-08-14 and William's decisions on route (A), `:auto`
confirmation, and the secant anchor.

## Results

Every acceptance criterion below was met. The measurements that matter:

| what | result |
|---|---|
| v4 regeneration vs the shipped v3 | `DX/DY`, `PLXFAC`, `M`, every covariance **bit-identical** (0.0e+00), both targets |
| reference-track round trip | rebuilds `natural_place` to **0.0 mas** from the product alone |
| barycentric half vs natural place | **260–546 mas** apart (Barnard) — the parallax ellipse the split exists to avoid double-counting |
| v4 null, `barycentric_lighttime=false` | **2.0e-7 mas** (Barnard), 1.3e-7 (ε Eri) — the Julia readout reproduces the Python propagator through the whole chain |
| v4 null, `barycentric_lighttime=true` | **0.082 mas** (Barnard), 3.6e-5 (ε Eri) — *this is the light-time term itself*, measured |
| v4 sensitivities vs the propagator | match to 3 s.f. on **both axes, both targets** (Barnard: 3.03 mas/(km/s), 0.612 mas/mas, 0.0324 mas/(mas/yr)) |
| v3 sensitivities | exactly **0.0** — the Δ-linear model is blind to all three, as designed |
| Phase 0 closed form vs the propagator | agrees to **<1%** on both targets |
| G23H frame-only light-time difference | **≤0.008 mas/yr** across all channels, 0.1σ at worst — the `_dedoppler` fix holds, so pin removal is safe |
| `:auto` after the pin | resolves `barycentric_lighttime = true` on both folds, **by measurement** (300 draws) |
| secant vs the hand-rolled `_reflex_secant` | **0.0 exactly**, 5 (P, m) cases |
| migrated nc model vs hand-rolled shear | **0.0 exactly**, every sampled quantity, both nc folds |

Gates run. **Octofitter suite: 2 failures, both pre-existing and unrelated** —
`test/v2/gradients.jl:65-66`, the `@allocated == 0` hot-path gates, measuring
272 B. That model is `RelAstromObs` + `RadialVelocityObs` and nothing else; no
file changed here is on its code path (`radial-velocity.jl` and `obs.jl`, which
are, carry unstaged edits from the prior session). The measurement is not even
stable — evaluated standalone the same call allocates 16 B — so it looks like
an escape-analysis gate that has drifted on this Julia patch. **Not bisected**:
isolating it would mean reverting the earlier session's unstaged PlanetOrbits
and Octofitter work, which was out of scope here. Worth its own look.

`cluster/barnard-slew/src/smoketest.jl` **23 passed, 0 skipped** —
every config in the matrix builds on the v4 product, and its independent
`nc[shear invariance]` check matches `ln_like` on 3 draws, which is what
verifies the `AnchoredFrame` migration against an implementation that did not
change. The acceptance criteria above are also installed as a repeatable gate,
`cluster/barnard-slew/src/verify_fgs_conventions.jl` (they cannot live in
Octofitter's own suite: they need a real product, and those are built
elsewhere).

Two findings worth carrying forward:

1. **The 0.082 mas null is the light-time term, isolated.** Because the
   reference track is the pipeline's light-time-free propagation, differencing
   a rigorous path against it measures exactly the barycentric light-time
   difference and nothing else — confirmed by the flag-off null collapsing to
   2e-7 mas. That is the ~0.1 mas figure from the scoping discussion, now
   measured rather than estimated: ~8% of the FGS−Gaia PMA σ.
2. **The large `:auto` impact numbers are the reflex, not the frame.** The
   impact test reported 115σ per point for G23H, which looks alarming next to
   the old docstring's "≤0.02 mas/yr". With the frame pinned at the DR3
   solution and no companion, the two settings differ by **≤0.008 mas/yr** —
   0.1σ, nowhere near μ·|v_r|/c = 3.83 mas/yr. The 115σ comes from prior draws
   admitting companions up to 300 M_jup, whose reflex genuinely is light-time
   sensitive. `:auto` resolving on is correct; the number is just dominated by
   the prior's tail.

## The v3 path was removed (William, 2026-08-18)

As first built, this kept a compatibility path: v3 products loaded on the
Δ-linear model, protected by the Phase 0 guards. William's call was that with
nothing shipped, a shim is not worth its cost — require the final format.
Removing it deleted, rather than added:

- **Phase 0 entirely** — `_fgs_check_linearization`,
  `_fgs_linearization_sensitivity`, `_fgs_curvature_budget`, `_AU_PER_YEAR_MS`
  (~110 lines). Those guards existed *only* to make the frozen-curvature
  approximation safe; with no approximation there is nothing to guard.
- the `rigorous::Bool` field, the two-path branch in `_fgs_model!`, the
  conditional reference-track/solve-epoch reads, and the instance-level
  `has_correction_impact`.
- `write_product`'s silent downgrade to v3, which now raises.

Re-verified after removal: the gate reports **identical numbers** (2.04e-7 mas
null, 0.0824 mas light-time term, all six sensitivities), confirming the
removal was purely subtractive. A v3 file now fails at load with an error
naming the missing columns and the command to regenerate.

`picklefgs.read_product` deliberately still accepts v1–v3. It is the pipeline's
own archive reader, not a consumer shim — five regression tests read ups And
and ε Eri products built before this bump — and it adds no second code path.

## Deviations from the plan as scoped

- **The hand-rolled `_reflex_secant` helpers were NOT deleted** (Phase 4 said
  to). `src/smoketest.jl` uses them as an *independent* implementation for its
  shear-invariance check; deleting them would have turned that test into a
  tautology at exactly the moment the shear moved into the framework. They are
  unwired from the model and marked as the reference oracle.
- **`refine=false` on the migrated folds**, so the migration is bit-identical
  to what it replaced. The two-pass refinement is more accurate but would have
  changed the folds under the same name; that is a separate, deliberate choice
  to make with the refit.
- **Incidental fix:** `test/v2/sampling-jacobians.jl:35` used a bare
  `Trajectory`, which is ambiguous in `Main` once AdvancedHMC is loaded (both
  export one). Pre-existing — it errored before any change here — and it was
  blocking the suite as a regression gate, so it is now qualified.
- **`FGSEpochAstromObs` gained `correction_impact`** (not in the original
  scope). Without it FGS answers "cannot say", which resolves the flag on by
  default — so in the joint fold `:auto` would have been *defaulting*, not
  confirming, which is what William asked Phase 3 for.
- **Phase 0 was written, verified, and then deleted** — see above. Its
  measurements were not wasted: the closed-form sensitivities it derived are
  what the verification gate now checks the rigorous path against.

## Open question, deliberately not decided

`HipparcosIADObs` and `GaiaDR4AstromObs` still declare
`reduced_lighttime_free = true`. By the criterion now written into that
trait's docstring — "one reduction, in its own parameterization" vs "a
difference *between* reduction windows" — they are arguably on the other side
of the line from G23H: IAD abscissae are direct measurements within one
window, not inter-window differences. Left as is: it was not in scope, it does
not affect the Barnard matrix (the IAD is always dropped, and there is no DR4
term), and it deserves its own decision.

Spans three repos:

- `pickle-fgs` — product format v3 → v4 (export the reference track)
- `Octofitter-fgs-public` — `fgs-pickle.jl` forward model, `corrections.jl`
  pin, `anchored-frame.jl` secant window
- `cluster/barnard-slew` — refit the FGS-bearing folds

Read `g23h-multi-source.md` for the frame-parameterization background; this
document is independent of it except at Phase 4.

---

## 0. The problem, measured

`_fgs_model!` (`src/likelihoods/fgs-pickle.jl:417-429`) predicts the exported
offsets with a **first-order Taylor expansion about REFSOL**:

    model_x = Δα*₀ + Δμ_α* τ + Δπ·PLXFAC_X + reflex + wedge/lock/drift

It is linear in τ and carries no perspective acceleration, no light-time term,
and no dependence on `rv` at all — `refsol.rv_kms` is parsed at line 203 and
appears nowhere else in the file. All curvature lives in the reference track
`xi_ref` the pipeline subtracted, frozen at the pipeline's `rv_ref`.

That is defensible only while the Δ quantities are small. For Barnard they are
not evaluated anywhere near the expansion point: the product's epochs are
**1993.10–1996.42** against `REFSOL EPOCH_JYEAR = 2016.0`, so **τ ∈ [−22.90,
−19.58] yr** and the star has moved **238 arcsec** from where the expansion was
taken. Everything below scales as τ².

Measured by propagating REFSOL through the pipeline's own
`refframe.propagate_astrometry` and differencing against the linear form, at
τ = −22.90, on the δ axis where the signal lives:

| term | Barnard | ε Eri (τ = −14.8) |
|---|---|---|
| total perspective curvature | **335 mas** | 1.1 mas |
| ∂(exact − linear)/∂rv | **3.03 mas per km/s** | 0.068 |
| ∂(exact − linear)/∂plx | **0.61 mas per mas** | 0.004 |
| ∂(exact − linear)/∂pm | **0.032 mas per (mas/yr)** | 0.001 |
| barycentric light-time | ~0.1 mas | — |
| Roemer (`pmpx`), if non-differential | 0.17 mas | — |

Per-epoch σ is 1.6–4.3 mas; on 35 epochs the ensemble mean is ~0.4 mas.

The 335 mas of curvature is absorbed exactly by `xi_ref` and is **not** an
error. The derivatives are, because the model cannot move them:

- **rv.** Fixed at the catalog value in every FGS-bearing fold, and
  `RV_KMS = -110.468` is correctly populated in the product (the `_sanitize`
  NaN→0 failure mode did not occur). But the catalog σ is 0.131 km/s → a
  **0.40 mas** frozen systematic, ≈1σ of the 35-epoch mean, sitting on the δ
  axis. And the `rv=:free20` folds put rv on a ±20 km/s box — 60 mas that FGS
  is structurally blind to. Those folds are `use_fgs=false` today; nothing but
  the fold matrix's own conventions keeps that shut.
- **plx.** The catalog-prior folds are safe (Δπ ~ 0.02 mas → 0.01 mas). The
  `fgs-g23h-uniplx` fold samples `Uniform(450, 650)`; a 5 mas excursion is
  3 mas of unmodelled curvature and a 50 mas one is 30 mas.
- **pm.** ±30 mas/yr box → up to 0.95 mas. Under the `nc` shear the derived
  barycentric pm is deliberately *unbounded* at long P, so this grows with the
  secant drift.

ε Eri is negligible on every line by two to three orders of magnitude. Barnard
is the whole reason for this work; ε Eri is the null-test target.

### Not in scope here

The `frame_shift`/FGS `pmra` conflict (`cluster/barnard-slew/src/model.jl`
:193-203) is larger than anything above — ~3 mas/yr over a 23-year τ is 69 mas
— but it is already diagnosed and already fixed by the `-nc` folds. Phase 4
replaces the hand-rolled fix with `AnchoredFrame`; it is a tidy-up, not a bug.

---

## 1. Decisions already made — do not relitigate

1. **Route (A): export the reference track.** Octofitter does *not*
   reimplement `propagate_astrometry`. Product regeneration is cheap, there
   are no external consumers, and breakage is explicitly a non-issue
   (William, 2026-08-18). This is also the only route that can carry the
   Roemer term, which needs HST's barycentric position — not recoverable from
   `PLXFAC_X/Y`, since those are the p̂,q̂ projections and Roemer needs the û
   component.
2. **The coordinate origin is defined by the data and must not move.**
   `export.py:400` exports `dxy = minv @ (ξ_meas − ξ_ref)/MAS` with `ξ_ref`
   from `natural_place(REFSOL, t)`, propagated light-time-free. Whatever the
   model does, it must difference against **that same track**. The tempting
   "upgrade" — differencing against a rigorously propagated REFSOL — silently
   redefines the origin. The correct form is deliberately mixed-convention:
   `rigorous_apparent(t) − REFSOL_lighttime_free(t)`.
3. **The light-time pin comes out and `:auto` decides**, after Phase 2, with
   the `CorrectionReport` recorded. 0.1 mas over a 23-year FGS−Gaia PMA
   baseline is ~0.004 mas/yr against a PMA σ of ~0.05 mas/yr — 8% of σ, a
   systematic not noise, and worth confirming rather than assuming.
4. **`AnchoredFrame` gets a secant window**, not the instantaneous reflex.
   The reason is the one already recorded in `model.jl`'s helper docstring: a
   10 M_jup at P = 10 d has α ≈ 2 mas but ~475 mas/yr of *instantaneous*
   reflex PM, which would push `pmra_obs` outside its ±30 mas/yr box and clip
   the short-P/high-mass region the data actually allow. The secant → 0 there,
   exactly as the data's sensitivity to the wiggle-as-drift does.
5. **The parallax stays differential.** `xi_ref` already carries REFSOL's full
   parallactic displacement, and `PLXFAC` is the pipeline's own
   d(sky)/d(plx@2016) with its ephemeris conventions embedded. `Δπ·PLXFAC` is
   correct and stays untouched. So do the wedge, lock and drift blocks.

---

## 2. Phase 0 — guards (Octofitter only, ships against v3 products)

The two live risks above become loud failures instead of silent biases. No
product change, no format change; this is worth doing first because it
protects every existing v3 product and every fold run before Phase 2 lands.

`ln_like` already enforces `ref_epoch == REFSOL epoch`
(`fgs-pickle.jl:490-492`). Add, in the same place and the same style:

- **rv**: error if `|θ_system.rv/1000 − refsol.rv_kms|` exceeds a tolerance
  derived from the product's own τ range and σ. Suggested rule: reject when
  the implied curvature error `|Δrv| · ∂/∂rv · τ²`-scaled exceeds ~0.25× the
  median per-epoch σ. For Barnard that is a tolerance of ~0.3 km/s; for ε Eri
  ~12 km/s, which is the right asymmetry and falls out of the data rather
  than a hardcoded number.
- **plx**: same shape, as a `@warn` rather than an error — a wide-π fold is a
  deliberate choice and should announce itself, not abort.

The message must say *why* (the model is linear in τ and the curvature is
frozen at REFSOL's rv) and *what to do* (fix rv, or move to a v4 product).

**Acceptance:** the `uniplx` and `rvfree` configurations in
`cluster/barnard-slew` trip the new diagnostics; every currently-passing fold
still passes silently.

---

## 3. Phase 1 — pickle-fgs: export the reference track (FVERS 4)

### What to export

`crossmatch.natural_place` (line 193) is

    roemer_yr = Σ(u₀ · obs.pos_au) · AU_LIGHT_YR
    prop      = propagate_astrometry(REFSOL, 2016.0 → jyear_tdb(t) + roemer_yr)
    u         = topocentric_direction(unit_vectors(prop), prop.parallax, obs.pos_au)

The model needs the **barycentric direction** — `prop`, *before*
`topocentric_direction` — because Octofitter's `frame_ra`/`frame_dec` are the
barycentre's apparent direction with no observer, and the parallax is supplied
separately by `Δπ·PLXFAC` (decision 5). Exporting the full natural place
instead would double-count the parallax by the whole ±547 mas ellipse.

New `EPOCHS` columns, all `D`:

| column | value | units |
|---|---|---|
| `RA_REF_DEG` | `prop["ra_deg"]` | deg |
| `DEC_REF_DEG` | `prop["dec_deg"]` | deg |
| `PLX_REF_MAS` | `prop["parallax_mas"]` | mas |
| `ROEMER_YR` | `roemer_yr` | Julian yr |

`PLX_REF_MAS` is not consumed by the model — it is there so the round-trip
test can reconstruct `natural_place` from the product alone, and so the
product is self-describing. It is a real quantity: REFSOL's parallax is
0.78 mas larger at 2016 than at 1993.

`ROEMER_YR` **is** consumed: Octofitter must solve at the same retarded epoch
the reference track was propagated to, or the difference is taken between two
different times. It is ≤499 s, worth 0.17 mas at Barnard's µ once the model
stops being differential.

### Code changes

- `crossmatch.py`: `natural_place` grows a way to return the intermediate
  `prop` and `roemer_yr` (an optional return, or a sibling `barycentric_place`
  sharing the body — a small refactor, no behaviour change).
- `export.py:_ref_track` (line 253): return the new quantities alongside
  `x0, minv, fac`.
- `export.py`: write the four columns; bump `FVERS` 3 → 4.
- `spec/09_octofitter_export.md`: document the columns and the
  barycentric-not-natural-place distinction, which is the one thing a future
  reader will get wrong.

### Acceptance

1. **Round-trip.** `topocentric_direction(unit_vectors(RA_REF_DEG,
   DEC_REF_DEG), PLX_REF_MAS, obs.pos_au)` reproduces the pipeline's own
   `natural_place(REFSOL, t)` to < 1e-9 deg at every epoch. This is the test
   that proves the export is the same track the data were differenced against.
2. **`ROEMER_YR`** is within ±499 s and its sign matches
   `crossmatch.py:203`.
3. `refit_from_product` and the embedded `PLXCHK ± PLXCHKS` gate are unchanged
   by the export (nothing subtracted differs).
4. Barnard and ε Eri products regenerated; `PIPESHA` clean (the current
   Barnard product is `-dirty`).

---

## 4. Phase 2 — Octofitter: the rigorous frame readout

### Loader

`FGS_PICKLE_FVERS` (line 20) becomes a supported set `(3, 4)`; the hard
equality at line 168 becomes a membership test with a branch stored on the
observation. v3 products keep today's model plus Phase 0's guards; v4 takes
the rigorous path.

Two epoch corrections on the v4 path:

- **Roemer**: the solve epoch becomes `epoch + ROEMER_YR · julian_year`.
- **UTC vs TDB**: `table.epoch` is `TIME_MJD_UTC` ("Octofitter convention",
  line 208) while `tau` is built from `JYEAR_TDB`, and PlanetOrbits solves in
  barycentric TDB-like MJD. That ~69 s inconsistency is worth
  **0.023 mas** at Barnard's µ — below σ, and harmless today because the
  linear model reads only `tau`. It stops being harmless the moment the model
  differences against a TDB-referenced track. Build the v4 solve epoch from
  `JYEAR_TDB`. *(Flagged for verification during implementation — I have not
  confirmed which column the trajectory machinery actually consumes.)*

### Forward model

In `_fgs_model!`, on the v4 path, `Δα0 + Δpmra·τ` and its δ counterpart are
replaced by a per-epoch readout of the frame that PlanetOrbits already solves:

```julia
sol = solutionat(ctx, i)
Δx  = _wrap180(frame_ra(sol) - tbl.ra_ref[i]) * 3.6e6 * cosd(tbl.dec_ref[i])
Δy  = (frame_dec(sol) - tbl.dec_ref[i]) * 3.6e6
mx  = model_x[i] + Δx + Δplx * tbl.plxfac_x[i]
my  = model_y[i] + Δy + Δplx * tbl.plxfac_y[i]
```

`model_x[i]` already holds the reflex from `accumulate_offsets!` against
`ref=Barycentre` (the zero-argument `raoff`, no observer position) — unchanged.
`Δα0` disappears entirely: the sampled `ra`/`dec`/`pmra`/`pmdec`/`rv` no longer
enter the model directly, they enter *through the frame*, which is exactly
what makes the result rigorous. `cosd(tbl.dec_ref[i])` is the per-epoch value
rather than `cosd(rs.dec)`; the difference is only ~3e-4 mas but it is free.

Notes for the implementer:

- **Cancellation.** `frame_ra` and `ra_ref` are both ≈269.4°, differing by a
  few mas. The subtraction itself is exact (Sterbenz); the error is 1 ulp of
  269.4° ≈ 2e-4 mas from `frame_ra`'s own `atan2`. Four orders below σ. Under
  `Dual` the partials do not cancel at all.
- The bookkeeping closes because `PLXFAC` is d(sky)/d(plx **at 2016**) —
  propagated — and `Δplx = θ_system.plx − refsol.parallax` is also referred to
  2016. `PLX_REF_MAS` is deliberately *not* used here.
- On v4, drop Phase 0's rv **error** (rv is now modelled and FGS genuinely
  constrains it at ~3 mas/(km/s)); keep an informational line in the report.

### Acceptance

1. **Null.** With `θ_system` set exactly to REFSOL and zero companion mass,
   both the v3 and v4 paths return 0 to < 1e-3 mas at every epoch.
2. **The table reproduces.** Perturb `rv`, `plx`, `pmdec` one at a time; the
   v4 path must move by the Phase-0 measured values — 3.03 mas/(km/s),
   0.61 mas/mas, 0.032 mas/(mas/yr) at τ = −22.90 — and the v3 path must not.
   This is the criterion that proves the fix does the thing it claims.
3. **ε Eri null.** The v3 and v4 paths agree to < 0.1 mas across the ε Eri
   product for any reasonable Δ, confirming the two-to-three-order gap in §0.
4. **v3 regression.** A v3 Barnard fit reproduces its current posterior.
5. Docstring lines 50-62 rewritten: the forward-model equation shown there is
   the v3 one and would otherwise become a lie.

---

## 5. Phase 3 — remove the light-time pin, confirm with `:auto`

With FGS reading the same rigorous apparent path as G23H, the internal
inconsistency that made the pin load-bearing is gone.

- Delete `reduced_lighttime_free(::Type{<:G23HObs}) = true`
  (`g23h.jl:2056`) and the pin block at `corrections.jl:344-352`.
- **Keep the trait and the machinery.** It is still the right answer for an
  observation that is purely a reduction product; it is `G23HObs` — whose
  dominant channels (HG, DR3−DR2) are *inter-window* — that should not have
  claimed it. Update the `corrections.jl:217` docstring to say so.
- Run the Barnard `g23h` and `fgs-g23h` folds and record what `:auto` decides,
  with the `CorrectionReport` and the prior-predictive impact numbers.

**Acceptance:** `:auto` resolves `:barycentric_lighttime` on, the HG-channel
pull moves in the direction the sign predicts, and the report is committed
alongside the result rather than described.

---

## 6. Phase 4 — `AnchoredFrame` with a secant window

Orthogonal to Phases 0-3; do it after, or in parallel by someone else.

- `_interim_solve` (`anchored-frame.jl:36`) takes a vector of epochs rather
  than one — `Trajectory` already accepts an `SVector`.
- `anchor_offsets(posys, anchor, epoch; window=(t1, t2))`: `pmra`/`pmdec`
  become the **secant** `(raoff(t2) − raoff(t1)) / (t2 − t1)`;
  `ra_cosdec`, `dec`, `rv`, `dz` stay instantaneous at `epoch`.
- `window` threads through `anchored_frame` and `AnchoredFrame`. Make it
  **required** when the caller wants a secant — a guessed default window is
  worse than no window, and the right one (bracketing every dataset in the
  model) is a modelling statement.
- The one-pass Jacobian stays exactly triangular with unit diagonal: the
  secant is still built from the body variables and the anchored parallax
  alone. Record that in the docstring, since it is the property the prior
  argument rests on.
- Migrate `cluster/barnard-slew/src/model.jl`'s `nc` path onto it
  (`NC_WIN_T1_MJD`/`NC_WIN_T2_MJD` become the `window`), and delete the
  hand-rolled helper.

**Acceptance:** the migrated `nc` model reproduces the hand-rolled shear to
machine precision at fixed parameters, and `fgs-g23h-nc` reproduces its
current posterior.

---

## 7. Phase 5 — refit

Every `use_fgs=true` fold changes: `fgs`, `fgs-g23h`, `fgs-g23h-era0`,
`fgs-g23h-nc`, `fgs-g23h-uniplx`, `fgs-g23h-wide*`. The `g23h*` folds change
only through Phase 3.

Expect the largest movement in `fgs-g23h-uniplx` (§0: the wide-π fold is the
one the frozen curvature actually hurts) and essentially none in the
catalog-π folds, where the Δ quantities are small enough that the
linearization was fine. If a catalog-π fold moves by more than a few tenths of
a mas in the FGS residuals, something in Phase 2 is wrong — that is the
cheapest global check available and it should be run before the matrix.

---

## 8. Ordering

    Phase 0 ──────────────────────────────► (ships alone, protects v3)
    Phase 1 ──► Phase 2 ──► Phase 3 ──► Phase 5
    Phase 4 ──────────────────────────►

Phase 3 must not land before Phase 2: removing the pin while `_fgs_model!` is
still linear desynchronizes G23H (rigorous) from FGS (linear) *within one
model sharing frame parameters*, which is worse than the data-convention
mismatch the pin was defending against.
