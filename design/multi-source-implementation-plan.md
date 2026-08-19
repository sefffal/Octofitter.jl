# Multi-source absolute astrometry — implementation plan

Status: **agreed design, ready for implementation.** Written 2026-08-12,
following the problem statement in `g23h-multi-source.md` (referenced below by
its section numbers, e.g. "§4") and a design discussion with William. This
document is the hand-off to the implementing agent: it records the decisions,
the full change list in phases, and the acceptance criteria per item.

Read `g23h-multi-source.md` first. It is the *why*; this is the *what*.

Line references were verified against the working tree on 2026-08-12:

- `src/model/codegen.jl:15-17` — evaluation order: system non-deferred →
  body blocks → system deferred lines (which see `system.*` and body
  namespaces).
- `src/model/codegen.jl:540` — the generated builder returns the
  `PlanetOrbits.System`; `codegen.jl:632` — `evaluate` solves it once per
  sample.
- `src/macros.jl:202` — `LL += expr` inside `@variables` compiles to a
  `DirectLLObs`: an arbitrary additive log-likelihood term. This is the
  potential-term mechanism used throughout this plan.
- `src/macros.jl:58-107` — `expr ~ Distribution` compiles to a
  `UserLikelihood`: a logpdf constraint on a derived expression.
- `PlanetOrbits.jl/src/system.jl:90-146` — the frame is a single field
  `frame::F` on `System`; masses, rows, specs, `Ainv`, fluxes are all
  frame-independent. `_convert_frame` / `_frame_scalar_type`
  (`system.jl:295-307`) handle scalar-type promotion.
- `PlanetOrbits.jl/src/frame.jl:15-87` — the `NoFrame` / `Parallax` /
  `AbsoluteFrame` hierarchy. Body–barycentre kinematics are well defined
  without an `AbsoluteFrame`; this is what makes the interim solve possible.
- `src/likelihoods/g23h.jl` — `table.kind` names channels (line 127);
  `_g23h_channel_kinds` (223) is the canonical 11-component ordering; the
  hardcoded-width mask is at 1776 (`ntuple(i -> ... , Val(11))`); the moments
  assembly is `_g23h_catalog_moments` (2033); the TAP fetch path is
  `_query_gaia_dr3` (470).

---

## 0. Decisions already made — do not relitigate

1. **The degeneracy fix is a parameterization-layer change, not a likelihood
   change.** The likelihood layer predicts from the barycentric frame, purely
   physically. All degeneracy-breaking (the job `frame_shift` does today,
   per-observation, incorrectly for multi-source — §4) moves into explicit
   sampled-observable → physical maps declared in the model.
2. **Anchored frame parameterization**: sample the anchor source's observed
   catalog quantities (α, δ, plx, μα*, μδ, and optionally systemic RV, at the
   catalog reference epoch); derive the barycentric `AbsoluteFrame` inputs by
   subtracting the model's own motion of the anchor body relative to the
   system barycentre at that epoch. The joint map is triangular with unit
   diagonal ⇒ |Jacobian| = 1: priors may be placed on the sampled anchor
   quantities directly, or on the derived barycentric values via
   `expr ~ Dist`, with no correction term.
3. **Per-source absolute channels, not pairwise relative objects.** Gaia
   cross-source covariances are ~block-diagonal, so each Gaia source gets its
   own observation object ("this object is Gaia 1234; that one is Gaia
   1235"), and each source's own 5×5 pos/plx/pm correlations live inside its
   own (widened) coupled block.
4. **The coupled block becomes kind-driven, not index-driven.** The existing
   `table.kind` column is the source of truth; the hardcoded 11-wide index
   arithmetic goes away. (William, explicitly: use `:kind` directly.)
5. **Wide stellar companions are parameterized in Cartesian relative phase
   space** (Δx, Δy, Δz, Δvx, Δvy, Δvz at a reference epoch), with elements
   derived. Prior choice is orthogonal, expressed via `LL +=` potential
   terms (see Phase 2.3 for the closed forms).
6. **No sharing of per-source nuisances** (`transit_priorities`, `σ_AL`,
   `σ_att`, `σ_calib`) across sources (§9). The calibration terms are
   magnitude-dependent (factors of 3–6 between ups And A and B) and the
   transit selections are made by different mechanisms per source.
7. **Doc examples must drop catalog-derived priors on plx, pmra, pmdec, and
   RV** (and narrow ra/dec boxes) wherever the corresponding quantity becomes
   likelihood data — the prior and the channel would double-count the same
   catalog measurement.
8. **The barycentric parameterization remains supported and remains the
   default** until the Phase 6 measurement says otherwise. Anchoring is a
   modeling choice; for two well-measured sources plus the `mass_single_msms`
   mass ratio the barycentric form may well sample fine (§12).

---

## Phase 1 — PlanetOrbits.jl

### 1.1 `reframe(sys, frame)`

New function returning a `System` identical to `sys` but with the frame field
replaced. Scalar type must be re-promoted (`_frame_scalar_type` of the new
frame joins the promotion; `_convert_frame` applied). Masses, rows, specs,
`Ainv`, fluxes are reused as-is — they are frame-independent by construction.

Acceptance:
- `orbitsolve` on `reframe(sys, fr)` agrees with a freshly constructed
  `System(...; frame kwargs...)` to bit-identical results (same inputs, same
  epochs), for all three frame levels in both directions
  (Parallax → AbsoluteFrame and back).
- ForwardDiff through `reframe` w.r.t. frame quantities works and matches
  finite differences.

### 1.2 Cartesian-state orbit constructor

Construct an orbit spec from a relative state at an epoch:
`(x, y, z, vx, vy, vz)` [AU, AU/yr] at `t_ref` [MJD], `about=` the usual
reference grammar. Closed-form osculating-element extraction; must be AD-safe
(no branches that break dual numbers on the generic path) and must fail
loudly (not NaN) on unbound/near-parabolic states — a clear error naming the
specific energy is fine for now.

Sign/axis convention must match the package's (x=RA/East, y=Dec/North,
z=LOS/+receding) triad convention (`frame.jl:94-99`).

Acceptance:
- Round trip: elements → solve at `t_ref` → extract state → constructor →
  elements, to tolerance, across a broad randomized grid of bound orbits
  (e in [0, 0.95], all angles, both z signs).
- Solve test: an orbit built from a state reproduces that state via
  `orbitsolve` at `t_ref`.
- AD through the constructor matches finite differences.

Out of scope (recorded so nobody "helpfully" adds them): Lambert two-position
constructor; any analytic marginalization of frame parameters.

---

## Phase 2 — Octofitter: interim solve, anchored frame, prior helpers

### 2.1 Codegen: expose a frame-free interim solve in deferred scope

The enabling fact: the anchored map needs the anchor body's motion about the
system barycentre at one epoch, and that is well defined at `Parallax` level
— the interim solve never needs the frame it is helping to construct.

Change to the generated builder (`codegen.jl`, around the system-construction
at 540): after body blocks and before the deferred system lines, construct a
**Parallax-level** `PlanetOrbits.System` from the already-resolved bodies and
whatever non-deferred parallax variable the model designates (the sampled
anchor plx in the anchored pattern). Expose it in deferred scope under a
documented name (suggestion: `system_interim`), then construct the final
frame from `system.ra/dec/plx/pmra/pmdec/rv` — which may now be
deferred-derived — via `reframe` (1.1), and proceed to observation blocks.

Rules, to be documented and enforced with clear errors:
- Deferred lines evaluate in declaration order; lines feeding the frame must
  not reference observation variables (observations evaluate last).
- `system_interim` is only available in deferred system lines.
- Building it must be skipped entirely (zero cost) for models that never
  reference it — key off whether the symbol appears in any deferred
  expression.

Consistency is the point of doing this in codegen rather than documenting a
user-space pattern: the interim system is built from the *same* body
declarations as the final one, so the two cannot drift apart.

### 2.2 `AnchoredFrame` helper

Sugar over 2.1. Given an anchor body, a catalog reference epoch, and which
quantities are anchored (position / pm / plx / rv, individually selectable —
mixing anchors per quantity is legal), expand into:

- sampled variables for the anchor's observed quantities at the catalog
  epoch (e.g. `pmra_A_dr3`), with whatever priors the user declares;
- deferred derived `ra/dec/plx/pmra/pmdec/rv` subtracting the anchor body's
  barycentre-relative position/velocity at `t_ref`, read from
  `system_interim` (converted AU/yr → mas/yr with the sampled plx).

Document the unit-Jacobian property (decision 2) and that a physical prior on
a *derived* barycentric quantity is one `expr ~ Dist` line.

Acceptance:
- **Equivalence test (single source):** an anchored-frame model reproduces
  the `frame_shift=true` behavior it replaces. Reuse the §2.2 verification
  method from the problem doc: likelihood maxima agree to ~1e-7 relative,
  maximizer displaced by exactly the model's own `Δpm_dr3`.
- **Multi-source correctness test:** with two `G23HObs` on one frame and the
  anchored parameterization, the predicted DR3 proper motions of the two
  sources *differ* by the model's relative velocity (the §4 bug is
  structurally impossible), and a synthetic-data injection (see 6.2) recovers
  an injected relative proper motion.
- ForwardDiff through the whole anchored path.

### 2.3 Prior/Jacobian helper library (uses `LL +=`)

Sampling coordinates and prior coordinates are decoupled by potential terms.
Ship, with docstrings deriving each:

- `logjac_cartesian_to_campbell(a, e, Mtot)` returning
  `-(3/2)·log(Mtot) - log(e) - (1/2)·log(a)` (plus constant), i.e. minus the
  log of the phase-space→Campbell Jacobian. Derivation to include in the
  docstring: Delaunay variables are canonical, so
  d³x d³v = dL dG dH dM dω dΩ with L = √(GMa), G = L√(1−e²), H = G cos i;
  transforming (L,G,H)→(a,e,cos i) gives
  **p(a, e, cos i, M, ω, Ω) ∝ M^{3/2} · e · √a**, flat in cos i and the
  angles.
- Documentation of the default: **flat-Cartesian sampling with no potential
  term is the Jeans/Ambartsumian phase-space prior** — thermal eccentricity
  p(e) ∝ e, isotropic orientation (flat cos i), p(a) ∝ √a, time-uniform
  phase — *and* an M^{3/2} tilt multiplying whatever prior the user declared
  on total mass. The mass tilt is the surprise; call it out in bold in the
  docs.
- Uniform-in-Campbell: sample Cartesian, add
  `LL += logjac_cartesian_to_campbell(...)`, then declare the desired
  Campbell priors as `expr ~ Dist` lines on the derived elements.
- O'Neil-style observable-based priors: the epoch-Jacobian potential,
  evaluated at the data epochs. Check whether v9 already carries a port of
  v8's `ObsPriors` before writing anything new; if it does, provide the
  glue and a doc example, not a reimplementation.

Acceptance:
- Prior-only sampling (no data) of a flat-Cartesian wide companion at fixed
  Mtot: KS tests confirm p(e) ∝ e, flat cos i, p(a) ∝ √a within the box
  image.
- Same with the `logjac_cartesian_to_campbell` term added: flat in
  (a, e, cos i) within the support.

---

## Phase 3 — g23h.jl kind-driven refactor (behavior-identical)

Replace the hardcoded 11-wide machinery with assembly driven by `table.kind`:

- `_g23h_channel_kinds` (223) remains the canonical ordering of the coupled
  block, but every consumer — the mask at 1776, `_g23h_catalog_moments`
  (2033), `_g23h_restrict`, `likeobj_from_epoch_subset` (968),
  `generate_from_params` — derives widths, masks, and index maps from the
  kind list generically (Val-parameterized width or a generated function),
  never from literal integers.
- This phase adds **no channels**. It exists so Phase 4 is a table edit, not
  surgery, and so channel addition is never again a cross-cutting edit to the
  most numerically delicate code in the file (§6).

Acceptance (gate before Phase 4 may start):
- Exact log-likelihood agreement with pre-refactor code (tight tolerance;
  bit-identical where achievable) on: ups And A (full channel set including
  `:iad_hip` and `:rv_dr3`), ups And B (Gaia-only, no Hipparcos), a
  Hipparcos-epoch-subset rebuild, and a `channels=`-restricted subset.
- Gradient agreement via ForwardDiff on the same cases.
- The UEVA deflation, `rho_dr2_dr3` cross-covariance, BINARYS σ-inflation,
  and ε·|Δpm| epistemic term all survive with values unchanged (they are the
  §6 hazard; test them by name).

---

## Phase 4 — new channels and schema

### 4.1 Catalog schema (versioned artifact — bump to G23H-v1.1)

Add per source, for DR3 and DR2: the full remaining 5×5 correlation columns
(pos↔pm, pos↔plx, plx↔pm — ten correlations per catalog release), and mean
`radial_velocity` / `radial_velocity_error` if any row lacks them. Fetch via
the existing `_query_gaia_dr3` TAP path (g23h.jl:470). The loader must accept
v1.0 (position channels then unavailable — clear error if requested) and
v1.1. The feather is shared with other projects: additive columns only, no
renames.

### 4.2 New channel kinds

Naming: the existing `:ra_dr3` etc. are **proper-motion** channels; position
channels need distinct names. Suggestion: `:pos_ra_dr3`, `:pos_dec_dr3`,
`:pos_ra_dr2`, `:pos_dec_dr2`, `:plx_dr3`, `:plx_dr2`, `:meanrv_dr3`. They
join the coupled block (Phase 3 makes this a table edit) with full
within-source correlations from 4.1.

Model predictions: apparent position of the observation's host body at the
catalog reference epoch (position channels); instantaneous source parallax,
i.e. barycentre distance plus the body's line-of-sight offset (plx channels);
line-of-sight velocity of the host including orbital motion at the DR3 RV
mean epoch (mean-RV channel).

**Position and plx channels ship default-off** (opt-in via `channels=`) until
4.3 is validated — an RUWE-7 source's position error double-counts the very
companions being modeled (§8), and turning the channels on before the
deflation exists would bias exactly the interesting systems.

### 4.3 Position analogue of the UEVA deflation (§8)

Extend the `deflation_factor_dr3` treatment to the position (and plx) block.
The derivation — how much of the astrometric excess noise that inflated the
five-parameter solution's position uncertainty is explained by the modeled
companions — must be written up in the PR and **needs William's sign-off
before the channels default on**. Implement the mechanism now; treat the
statistical model as provisional.

### 4.4 Mean-RV channel zero point

Per-source RV zero-point nuisance (gravitational redshift + convective
blueshift; several hundred m/s between an F8V and an M dwarf), with an
informative prior, entering the mean-RV channel prediction. Note in the docs
that the frame propagation is insensitive to zero-point errors at this scale
(perspective terms scale as v_r·μ·plx), so the nuisance belongs to the
channel, not the frame.

### 4.5 Docs: remove double-counting priors

Update every doc example that sets catalog-derived priors on `plx`, `pmra`,
`pmdec`, or RV (and narrow ra/dec boxes) to remove them wherever the
corresponding channel is active, with a short explanation of why. (Decision
7; William's explicit list.)

---

## Phase 5 — multi-source hygiene

- **5.1 `frame_shift` deprecation path (§4).** Keep the kwarg. Error —
  not warn — when `frame_shift=true` and more than one `G23HObs` shares the
  system frame, with a message pointing at `AnchoredFrame`. Single-source
  default behavior unchanged for now; full removal is a later release, after
  Phase 6.
- **5.2 Fluxratio tier-2 (§10).** The system-level `fluxratio` fallthrough in
  `_g23h_fluxratios` (1021) must either become per-observation-keyed or raise
  an error that names the observation and its companion count when two
  observations with different companion counts both reach tier 2. No silent
  behavior.
- **5.3 No-Hipparcos audit (§11).** Tests: `_g23h_catalog_moments` with
  `dist_hip === nothing` (mask/`n_components` bookkeeping);
  `generate_from_params` for a Gaia-only source; and make the perspective-
  acceleration guard (1769) explicit on finiteness of
  `nonlinear_dpmra`/`nonlinear_dpmdec` rather than coinciding with `has_hip`.
- **5.4 §9 decision recorded in docs.** Per-source nuisances stay per-source
  (decision 6). Document the silent index-misalignment hazard of any future
  sharing scheme (two GOST queries need not return equal row counts; index
  *i* then denotes different transits — a silent failure).

---

## Phase 6 — measurement (decides defaults; do not skip)

- **6.1 The §12 four-cell.** τ_int on the frame block (`pmra`, `pmdec`) and
  on companion masses, on the `rv-g23h` fold of the ups And model:
  {barycentric, anchored} × {A only, A + B}. This is the number that decides
  whether anchored becomes the recommended default in the docs, and whether
  `frame_shift` can be removed outright.
- **6.2 Synthetic wide-orbit recovery.** ups And B is a structural test, not
  a sensitivity test (§13). Use `generate_from_params` to inject a wide pair
  with a genuinely detectable relative proper motion and (with Phase 4
  channels on) a position offset; verify the posterior recovers the injected
  relative velocity and separation. This is the end-to-end test that the
  whole stack — anchored frame, kind-driven block, new channels, Cartesian
  wide orbit — actually constrains what it claims to.

---

## Suggested order and gating

Phases 1 and 3 are independent and can proceed in parallel. Phase 2 needs
1.1; Phase 4 needs 3 (gate: behavior-identical acceptance) and 4.2 benefits
from 2.2 existing for its tests. Phase 5 is independent of 4 except 5.1's
docs pointer. Phase 6 last.

Per-phase, the acceptance criteria above are the definition of done. Where a
criterion says "bit-identical" and that proves unachievable, document why and
the achieved tolerance rather than silently loosening it.
