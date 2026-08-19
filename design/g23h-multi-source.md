# `G23HObs` on multiple sources — goals, fixes applied, and open problems

Status: **problem statement only.** Written 2026-08-12 during the ups And
A(bcd) + B(Bb) session.

This document deliberately **does not propose solutions.** Every section below
states what is true of the code today and what follows from it. The tradeoffs
between the candidate designs are the thing to think through carefully, and
pre-committing to a fix in the same breath as the diagnosis tends to foreclose
that. Where an obvious remedy exists it is left unstated on purpose.

Line numbers refer to `src/likelihoods/g23h.jl` at the time of writing.

---

## 1. The goal

Support hierarchical **stellar** systems: two or more Gaia catalog sources
belonging to one physical system, modelled as one Octofitter `System` sharing
one absolute reference frame, each carrying its own `G23HObs`.

The driving case, and the one every number below is drawn from:

| | ups And A | ups And B |
|---|---|---|
| Gaia DR3 | 348020482735930112 | 348020242217448576 |
| G | 3.97 | 12.50 |
| Hipparcos | HIP 7513 | **none** |
| RUWE (DR3) | 7.25 | 1.52 |
| matched transits (DR3) | 22 | 33 |
| `mass_single_msms` | 1.290 ± 0.016 M⊙ | 0.2256 ± 0.0015 M⊙ |
| companions to model | b, c, d | Bb (candidate) |

Separation 55.62″ = 750 au projected at 13.478 pc; relative proper motion
2.082 mas/yr.

The intended payoff is recorded in the code already, at `skypath.jl:237-253`:

> What is emphatically *not* wanted is a frame offset per Gaia source. Several
> `G23HObs` in one system share one frame, and that shared frame is what binds
> the system together — the wide pair's relative astrometry constrains the wide
> orbit for free. Inventing per-source offsets would dilute exactly the
> constraint the joint fit exists to exploit. So: offsets for other frames,
> never between sources on the same one.

Sections 4-13 are largely an account of the distance between that paragraph and
what the code currently delivers.

### What already works

The v9 reference grammar is structurally ready. `G23HObs` takes `host=`,
`companions=` and `ref=` (constructor at line 411), reports them through
`refspecs` (line 202), and `resolveref` constant-folds each spec against the
solved system. Codegen evaluates observation blocks last, after the system
block and all body blocks (`model/codegen.jl:15-18`), so an observation's
variables can read `system.*` including deferred system variables. Nothing in
the plumbing prevents several `G23HObs` in one system.

---

## 2. Fixes applied in this session

### 2.1 `_g23h_channel_table` NaN guard — a v8→v9 regression

`G23HObs` could not be constructed for **any** Gaia-only source.

The Hipparcos and Hipparcos-Gaia rows are built unconditionally and dropped
only at the end of `_g23h_channel_table` (`has_hip || append!(drop, 1:5)`). For
a source with no Hipparcos entry the four `epoch_*_hip` / `epoch_*_hg` catalog
fields are `NaN`, and `years2mjd(NaN)` throws `InexactError: Int64(NaN)` from
`Dates.Date` long before that drop is reached.

v8 guarded each conversion individually (`isnothing(hip_like) ? NaN :
years2mjd(catalog.epoch_ra_hip)`, v8 g23h.jl:543-546). Folding all eleven into
a single `_e(y) = years2mjd(y)` helper during the v9 port dropped the guard.

Fixed at line 715 by keying on the value rather than on `has_hip`, which also
covers a row that has a `hip_id` but non-finite Hipparcos-Gaia epochs. ups And
B now constructs: 34-transit pool, channels `[:ra_dr2, :dec_dr2, :ra_dr32,
:dec_dr32, :ra_dr3, :dec_dr3, :ueva_dr3]`.

### 2.2 `frame_shift` keyword

New `G23HObs(...; frame_shift::Bool=true)` (line 430), stored on the struct,
threaded through both struct-rebuild sites (`likeobj_from_epoch_subset` and the
Hipparcos-subset path), and applied at line 1658. Default `true` preserves
existing behaviour exactly. See §4 for what it is and why it needed a switch.

Verified: the two settings give different likelihoods, and the likelihood
surface over `(pmra, pmdec)` is the same surface translated — maxima agree to
1.3 × 10⁻⁷ relative, maximizer displaced by (+0.557, +0.011) mas/yr,
i.e. by the model's own `Δpmra_dr3`.

---

## 3. Problem: no channel carries absolute position, parallax, or mean RV

`_g23h_channel_kinds` (line 223) is eleven components: five proper-motion pairs
(`hip`, `hg`, `dr2`, `dr32`, `dr3`) plus `ueva_dr3`. `:iad_hip` and `:rv_dr3`
sit outside the coupled block. **Every one of them is a velocity or a scatter
statistic.**

Absolute position does enter the absolute-frame branch (lines 1528-1552), but
only as a difference across a baseline:

```julia
Δα_hg_prop = (α_dr3₀ - α_h₀) * 60 * 60 * 1000 * cosd((δ_dr3₀ + δ_h₀) / 2)
pmra_hg = (Δα_dr3_p - Δα_h + Δα_hg_prop) / Δt_hg_ra * julian_year
```

which reconstructs a proper motion (correctly, including the changing cos δ and
the perspective term). A source's large constant offset from the reference
point lands in `Δα_dr3_p` and cancels in `Δα_dr3_p − Δα_dr2_p`.

Consequences:

- For a wide pair the **separation is unconstrained**; only the relative
  velocity is. The `skypath.jl` note's "relative astrometry constrains the wide
  orbit for free" is true of velocity and aspirational for position.
- The catalog `parallax` / `parallax_error` are never used as data, although
  they are present for both sources and two independent parallaxes of one
  physical system are a real consistency constraint.
- The catalog mean `radial_velocity` is never used as data. Note the existing
  `:rv_dr3` channel is **not** the mean RV — it is the Chance et al. (2022) RV
  *variability* statistic, a sample variance compared through a
  `NoncentralChisq` (line 1822 ff.). For a wide pair the mean-RV difference
  (A − B = 0.429 km/s here) is the only handle on the line-of-sight relative
  velocity, the one phase-space direction neither relative astrometry nor
  relative proper motion can reach.
- Of the six relative phase-space dimensions of a wide pair, the current
  likelihood constrains two (tangential velocity) and leaves four to priors.

---

## 4. Problem: the frame shift is per-observation but redefines a shared parameter

At the end of `_g23h_simulate!` (line 1658, comment from the original author):

```julia
# Shift the whole reference frame so the model proper motion refers to the
# primary rather than the barycentre. This vastly improves sampling
# efficiency, and is why Δpmra_dr3 is subtracted from *every* channel
# including DR3's own.
shift = @SVector [Δpmra_dr3, Δpmdec_dr3]
```

The DR3 channel is `μ_dr3 = pmra_dr3₀ + Δpmra_dr3` where `pmra_dr3₀` is the
frame `pmra` (lines 1516-1518). Subtracting the shift makes the model
prediction **exactly `pmra`**.

For one source this is a pure reparameterization: a constant subtracted from
every channel is absorbed by the frame proper motion, which stops meaning
"barycentre proper motion" and starts meaning "the primary's proper motion as
DR3 measured it" — the direction the data pin hardest rather than the one
degenerate with the reflex. The HG channel keeps its acceleration content
because the shift is common-mode across epochs.

For several `G23HObs` on one frame it is a correctness bug. Each observation
subtracts its *own* `Δpmra_dr3` while redefining the *same* system-level
`pmra`, so every source's DR3 channel predicts one and the same number. For ups
And the catalog proper motions differ by (0.545, 2.009) mas/yr against B's DR3
uncertainties of (0.042, 0.037) — **13σ and 55σ**. The relative proper motion
is not merely uninformative under `frame_shift=true`; it is removed from the
likelihood, and it is the entire wide-orbit signal.

Two further observations about the shift, independent of the multi-source case:

- **It is a reparameterization only because the frame proper motion carries a
  wide uniform prior.** A model that put an informative prior on `pmra`/`pmdec`
  would have its posterior silently changed by `frame_shift=true`, even with a
  single source. Nothing checks for this.
- `frame_shift=false` is a switch, not a resolution. There is currently no way
  to have both the sampling benefit and multi-source correctness, and the cost
  of turning it off is unmeasured (see §12).

---

## 5. Problem: `Barycentre` is system-global, so adding a body silently redefines existing observations

`Barycentre` with no members resolves to the mass-weighted barycentre of the
**whole** system (`model/refs.jl:120`). Adding a wide stellar companion
therefore moves the reference point of every `ref=Barycentre` observation
already in the model:

| system | f_B | barycentre offset from A |
|---|---|---|
| A + b,c,d vs B | 0.149 | **8.26″** at (+4.28″, −7.06″) |
| A + b,c,d vs B + Bb | 0.183 | **10.17″** |

For the ups And model as it stands this breaks two things at once: the system
`ra`/`dec` are pinned to A's catalog position (or, in the FGS folds, to a 60 mas
box about the FGS reference solution), and both are the *primary's* position,
not the barycentre's, once B exists.

The interaction with §3 is the part worth dwelling on: **this is currently
invisible precisely because no channel constrains absolute position.** The two
problems mask each other. Anything that fixes §3 will expose §5 in every
existing multi-body model, not only in new ones.

---

## 6. Problem: the coupled covariance block is fixed-width by construction

`_g23h_catalog_moments` (line 2033) assembles an `SVector{11}` of catalog means,
an `SVector{11}` of model values, and an 11×11 covariance, with
`mask = ntuple(i -> _g23h_channel_kinds[i] ∈ kinds, Val(11))` and hardcoded
index arithmetic. The same width is assumed in `_g23h_channel_kinds` (223),
`_g23h_restrict`, `likeobj_from_epoch_subset` (968) and
`generate_from_params`.

Consequence: adding any channel is a cross-cutting edit to the most numerically
delicate code in the file — the block that also carries the UEVA-driven DR3
deflation (`deflation_factor_dr3`), the DR2↔DR3 cross-covariance via
`rho_dr2_dr3`, the BINARYS σ-inflation of Σ_h and the ε·|Δpm| epistemic term.
The cost of §3 is not conceptual, it is concentrated here.

---

## 7. Problem: the catalog schema lacks the correlations a position channel would need

Present in `G23H-v1.0.feather`: `ra_error_central_dr3`,
`dec_error_central_dr3`, `ra_dec_corr_central_dr3` (and the DR2 equivalents),
`rho_dr2_dr3`, `parallax_error`, `radial_velocity_error`.

Absent: every correlation between position, parallax and proper motion — the
rest of the DR3 5×5. Those quantities were jointly fitted by AGIS and are
correlated; a position channel added from the columns on hand would be treated
as independent of the proper-motion channels it was fitted with.

(`_query_gaia_dr3` exists at line 470 as a TAP fetch path for missing columns,
so the data are reachable — but this is a schema question, not a lookup
question, and the catalog is a versioned artifact shared with other projects.)

---

## 8. Problem: per-source astrometric quality is modelled for proper motion only

The UEVA channel exists to ask how much of Gaia's uncertainty inflation the
companion model explains, and to deflate the DR3 covariance accordingly
(`deflation_factor_dr3`, applied as `Σ_dr3 * d^2` and fed back into the DR3−DR2
cross term). There is no analogue for position.

ups And A has RUWE 7.25 and `ra_error_central_dr3` = 0.253 mas — inflated by a
Gaia five-parameter fit that did not know about b, c and d. B has RUWE 1.52 and
0.037 mas. If position channels were added, A's position uncertainty would
double-count the companion signal in exactly the way the proper-motion channels
are explicitly corrected not to, and the two sources would be affected very
unequally (a factor ~7 in σ, and A is the one with the planets).

---

## 9. Problem: nuisance dimensionality scales with source count, and sharing is not obviously correct

Each `G23HObs` samples `transit_priorities` — one dimension per forecast
transit (`MvNormal(zeros(len_epochs), I)`, line 895) — plus `σ_AL`, `σ_att`,
`σ_calib`, plus `u_dup_dr2` where the DR2 duplicate count is marginalized.
Pools here: A 40 transits (cached GOST), B 34. Under `freeze_epochs=true` these
are frozen constants and the cost vanishes; under sampling it is ~40 dimensions
per source.

The two sources are 55.62″ apart — far inside Gaia's 0.7° field — so their
forecast pools are near-identical and the temptation to share
`transit_priorities` is strong. **It is not clearly correct, and the evidence
cuts against it for this pair specifically:**

- The genuinely common losses (dead time, decontamination, the OBMT gap lists)
  are *already* removed by `dr2_ok_mask` / `dr3_ok_mask` before priorities are
  applied. What `transit_priorities` marginalizes is the *residual* loss:
  gating, window class, saturation, per-source AGIS outlier rejection.
- Those are source-specific, and this pair is 8.5 magnitudes apart. A has
  **22** matched transits; B has **33**. The brighter star has *fewer* usable
  transits, because bright-star handling discards them.
- Sharing one priorities vector forces the smaller selection to be the top-k
  subset of the larger — the `partialsortperm` construction in
  `_g23h_default_variables` gives maximal epoch overlap by design, and says so
  in its own comment. That asserts A's 22 nest inside B's 33, when the two sets
  are selected by different mechanisms.

Separately, a hazard of *any* sharing scheme: `transit_priorities` indexes each
observation's own pool, built by `_g23h_forecast_pool` from a per-source GOST
query. Two queries at different positions need not return the same row count,
and if they differ by one row, index *i* denotes a different transit in each.
That failure is silent, not loud.

And regardless of priorities: `σ_AL`, `σ_att`, `σ_calib` are
magnitude-dependent calibration terms and must not be shared. Catalog values
here: σ_AL 0.018 (A) vs 0.054 (B); σ_cal 1.231 (A) vs 0.208 (B), a factor of 6.

---

## 10. Problem: the system-level flux-ratio fallthrough is keyed by name, not by source

`_g23h_fluxratios` (line 1021) resolves companion flux ratios in three tiers:
`θ_obs.fluxratio`, then `θ_system.fluxratio`, then the bodies' own
`flux_<band>` variables. Tier 2 is keyed by the bare name.

With two `G23HObs` in one system, a system-level `fluxratio` vector is read by
**both**, each validating it against its own companion count via
`_g23h_asratios(v, Val(N), T)`. A model where A has three companions and B has
one cannot use tier 2 at all: whichever observation is checked second raises
"received 3 flux ratios but the observation declares 1 companion".

This is a latent trap rather than a live bug — the ups And model uses tier 3 —
but tier 2 exists specifically for backwards compatibility with models written
against the old default-variables block, so it is exactly the path an older
single-source model would be on when someone adds a second source.

---

## 11. Problem: the no-Hipparcos path is unaudited beyond construction

§2.1 got a Gaia-only source to *build*. It did not audit the rest of the path.
Known unknowns:

- `_g23h_catalog_moments` with `dist_hip === nothing` takes the
  `(@SVector[0.0,0.0], @SMatrix zeros(2,2))` branch; the resulting zero rows are
  masked out, but the interaction with the `mask`/`n_components` bookkeeping is
  untested for this case.
- `generate_from_params` has its own `_g23h_isabs(...) && has_hip` branch that
  has not been exercised for a Gaia-only source.
- The HGCA perspective-acceleration correction (line 1769) is applied under
  `absolute && !isnothing(dist_hip)`. B's `nonlinear_dpmra`/`nonlinear_dpmdec`
  are `NaN` in the catalog, and B is skipped because it has no Hipparcos — so
  the guard happens to coincide with where those columns are finite in the two
  rows checked. Whether that is an invariant of the catalog, or a coincidence
  that would bite on a source with Hipparcos but null nonlinear terms, is
  **unverified**.
- In a wide pair the two sources have different channel sets *and* different
  time baselines (A: 1991-2016; B: 2014-2017). Anything that implicitly assumes
  a common baseline across observations is worth re-reading with that in mind.

---

## 12. Problem: the sampling cost of `frame_shift=false` is unmeasured

The shift was introduced because the frame proper motion is degenerate with the
reflex, and the author reports the reparameterization "vastly improves sampling
efficiency". Turning it off returns the sampler to that degeneracy.

The counter-argument for a two-source fit is that `pmra + reflex_A` and
`pmra + reflex_B` must satisfy two catalog proper motions simultaneously, with
the reflex ratio fixed by a mass ratio that `mass_single_msms` constrains to
0.7% — so the degeneracy that motivated the shift is broken by the same second
source that forces the shift off. **This is an argument, not a measurement.**
No run has been done either way, and τ_int on `pmra`/`pmdec` in the first
multi-source fit is the number that settles it.

There is a related unknown: a wide companion's orbit over a 25-year baseline is
close to linear motion, so B's orbital velocity is itself nearly degenerate
with the frame proper motion. How much of that degeneracy the mass-ratio
constraint actually removes is not known analytically here.

---

## 13. Problem: the ups And B case exercises none of the above at a level that would validate it

Worth recording so it is not over-claimed later. The astrometric perturbation
of A by B is small:

| | value |
|---|---|
| apparent acceleration of A | 1.17 × 10⁻³ mas/yr² |
| Δμ over Hip→DR3 (24.75 yr) | 0.029 mas/yr |
| vs σ(μ_HG) = 0.0155 | 1.9σ |
| vs σ(μ_Hip) = 0.389 | 0.07σ |
| vs observed \|μ_DR3 − μ_HG\| = 0.838 | 3.5% |

All upper bounds: projected separation is a lower bound on true separation and
g ∝ 1/r², and the full acceleration vector is used rather than its tangential
part. The pair's relative velocity is only ~10% of circular at 750 au, so it is
near apoapsis or much wider in 3D — both of which reduce g further.

To produce the observed proper-motion anomaly at 750 au would take
M_B ≈ 6.5 M⊙, about 29× the actual mass. So ups And B is a good *structural*
test of multi-source `G23HObs` — two sources, one frame, asymmetric Hipparcos
coverage, an 8.5-magnitude contrast — and a poor *sensitivity* test. It will
not validate that the machinery recovers a wide orbit, because the signal it
carries is at the 1.9σ level in one channel. A system with a genuinely
detectable wide-orbit acceleration would be needed for that.

---

## 14. Summary of the problem set

| # | Problem | Blocks multi-source? |
|---|---|---|
| 3 | No position / parallax / mean-RV channel | Yes — separation unconstrained |
| 4 | Per-observation frame shift redefines a shared parameter | Was blocking; now switchable |
| 5 | `Barycentre` is global, silently redefined by a new body | Yes, once §3 is fixed |
| 6 | Covariance block fixed at 11 components | Cost of fixing §3 |
| 7 | Catalog lacks pos↔pm, pos↔plx correlations | Cost of fixing §3 |
| 8 | No position analogue of the UEVA deflation | Cost of fixing §3 |
| 9 | Nuisance dimensionality per source; sharing not clearly valid | Efficiency, plus a silent-misalignment hazard |
| 10 | System-level `fluxratio` fallthrough keyed by name only | Latent trap for migrated models |
| 11 | No-Hipparcos path unaudited past construction | Unknown |
| 12 | Cost of `frame_shift=false` unmeasured | Unknown |
| 13 | ups And B is a structural but not a sensitivity test | Validation gap |

§3 is the root: §5, §6, §7 and §8 are all either masked by it or are its cost.
