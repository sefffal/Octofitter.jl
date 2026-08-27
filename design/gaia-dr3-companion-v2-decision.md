# `GaiaDR3CompanionObs` on v2 — decision: do not port

**Status:** recommendation to @sefffal. Nothing is merged; PR #104 is untouched.
**Branch:** `v2-gaia-dr3-companion` (worktree only, not pushed).
**Evidence:** `test/v2/g23h-companion.jl` (10 gates, offline, all green).

## The ask

On PR #140 you wrote:

> I think we should revive the branch `feature/gaia-dr3-companion` so we can
> attach a full set of Gaia observations to a particular companion without the
> overhead of `G23HObs` etc.

`feature/gaia-dr3-companion` is PR #104 (`7a10e2f`, June 2026, +774 lines):
a planet-level `GaiaDR3CompanionObs` carrying the published DR3 five-parameter
solution of a resolved bound companion, plus a docs page and unit/integration
tests.

## What PR #104 actually was

Not epoch astrometry. It consumed **one archive `gaia_source` row** — position,
parallax, proper motion and their 5×5 correlations — and compared it against:

- the companion's absolute position at the DR3 reference epoch, from
  `sol.compensated.ra2/dec2` plus `raoff(sol) + raoff(sol, M_companion)`;
- a **three-point central finite difference** (±180 d) for the proper motion;
- optionally the DR3 `radial_velocity`.

It required an `AbsoluteVisual{...}` basis, attached to a `Planet`, and inflated
**all five** published uncertainties by 1.37.

Two things about it do not survive contact with v9:

1. **The finite-difference proper motion is not what a catalog publishes.** A
   Gaia five-parameter solution is a least-squares fit to along-scan abscissae
   over the release baseline, not the instantaneous derivative at the reference
   epoch. The two agree only in the limit where the source's path is linear
   across the baseline — which is precisely the limit in which the observation
   carries no orbital information. Reproducing the reduction (project on scan,
   refit with the same design matrix) is what `skypath.jl`'s
   `prepare_A_5param` / `fit_5param_pinv` exist for, and what `G23HObs`
   already does for its `:ra_dr3` / `:dec_dr3` channels. A ported
   `GaiaDR3CompanionObs` would either keep the finite difference (and be
   wrong where it matters) or duplicate that math.

2. **The blanket 1.37 inflation is wrong on the parallax.** `_g23h_check_5x5`
   verifies against the G23H catalog on real sources that position and
   proper-motion errors carry exactly 1.370000 while the **parallax carries
   none**. PR #104 inflates the parallax too.

## The v9 recipe, and that it works

v9 removed the type's reason to exist. `host=` already names "the body this
catalog source is centred on", and the docstring's own "Several sources in one
system" section is this feature. A resolved secondary's own DR3 row is:

```julia
A = Body(name="A", variables=@variables begin mass ~ …; flux_G = 1.0 end)
B = Body(name="B", about=:A, variables=@variables begin mass ~ …; a ~ …; … end)

obs_B = G23HObs(;
    host = :B, companions = (:A,), ref = Barycentre,   # ← see Gap 1 and Gap 2
    gaia_id = <B's source_id>, catalog = <B's row>,
    channels = (:ra_dr3, :dec_dr3), ueva_mode = :none,
    include_rv = false, include_iad = false, frame_shift = false)

sys = System(name="pair", bodies=[A, B], observations=[obs_B],
    variables=@variables begin
        plx = …; ra = …; dec = …; pmra = …; pmdec = …; rv = 0.0
        ref_epoch = 57388.5
        fluxratio = [0.0]        # ← resolved: A contributes no light to B's source
    end)
```

Add a second `G23HObs(host=:A, companions=(:B,), …)` and an `AnchoredFrame(A; …)`
and you have the two-source wide-pair model, which is the case the anchored
frame was built for.

**Verified** (`test/v2/g23h-companion.jl`, offline against the `fixtures/`
catalog row with `hip_id = NaN` so the source looks like a faint secondary):

- the likelihood responds monotonically to the companion's separation;
- data generated at truth with `generate_from_params(…; add_noise=false)` are
  recovered: the profile likelihood peaks exactly at the injected `a = 750 au`
  over `[500, 1100]`, and at the injected host mass `1.0 M☉` over `[0.6, 1.7]`,
  with the truth point above every perturbed one.

So the feature works today, with no new code. But three things are wrong with
the ergonomics, and the first is severe.

## Gap 1 — `companions=()` silently drops the entire orbital signal

`G23HObs(host=:B, companions=())` is the obvious spelling for "this source is
just body B". It constructs, evaluates, and returns a likelihood that is
**bit-identical for every orbit**: `a = 10`, `750` and `3000 au` give the same
`ln_like` to the last digit. No error, no warning.

The cause is `_g23h_simulate!` (`src/likelihoods/g23h.jl:~1714`):

```julia
active = map(c -> !iszero(PlanetOrbits._primal(masses[c.idx])), comps)
any_active = any(active)
…
else
    # "Every per-transit perturbation is exactly zero, and a five-parameter
    #  fit to zero data returns zero parameters and zero χ²."
```

That comment is true when `host` is the body the barycentre coincides with once
every companion is inactive — the configuration the fast path was written for.
It is false when `host` is itself a companion, where
`raoff(sol, host, barycentre) ≠ 0` with no companions declared at all. `active`
is keyed on `companions=`, so it never sees that.

This is a silent-wrongness trap on the exact spelling a user following your PR
#140 comment would reach for first. Pinned as a test so a later fix fails
loudly rather than quietly.

**Not fixed here** — it is a change to the numerical fast path of a
release-gating likelihood, so it is yours to call. The narrow fix is to fold
"is the host distinct from the reference under all-inactive companions" into
`any_active`.

## Gap 2 — a resolved source has no clean spelling

The workaround for Gap 1 is `companions=(A,)` — declaring the *primary* as a
possible blend into the secondary's source, purely so that `active` is
non-empty. That works, but the blending weight now matters and the default is
wrong:

- **`fluxratio = [0.0]`** (resolved; A contributes no light to B's source) is
  what a resolved pair means, and gives the recovery results above.
- **falling through to the bodies' own `flux_G`** gives
  `f_A = flux_G(A)/flux_G(B) = 50` for a secondary 4 mag fainter, so B's
  "photocentre" is dominated by A. The test measures the two spellings
  differing by far more than 1 in `ln_like` at every separation.

The tier-3 escape (the bodies' own fluxes) is the documented, non-positional,
recommended path everywhere else in `G23HObs`, and here it is the wrong answer.
The two ways to override it are both awkward:

- a **system-level `fluxratio` vector**, which is shared by name across every
  `G23HObs` in the system with the same companion count — fine for a symmetric
  two-source pair, a foot-gun otherwise;
- an **observation-level `variables=` block**, which *replaces* `G23HObs`'s
  default block outright, so you must hand-write `σ_AL`, `σ_att`, `σ_calib`,
  `transit_priorities`, `transits`, `transits_dr2` yourself.

The clean fix is a small `G23HObs` addition — a constant `fluxratio=` /
`fluxratio_hip=` keyword appended to the *defaults* rather than replacing them,
so `fluxratio=[0.0]` can be written on the observation. That is a (small)
feature change to a release-gating likelihood, so it is also yours to call.

## Gap 3 — the overhead is not in the channels

`channels=` restricts the *datum*. It does not restrict the cost, the model
dimension, or the data you must have.

**Model dimension.** `transit_priorities` marginalizes the epoch selection over
the whole forecast pool and is declared regardless of surviving channels. On the
62-transit fixture, a two-channel `G23HObs` has model dimension **62** with every
body variable fixed — one sampled dimension per pooled transit. A real DR3 source
has ~800–1000 forecast transits, so a "DR3 position and proper motion only"
observation is a ~1000-dimensional nuisance block. `freeze_epochs=true` removes it
and is the lever that actually matters.

**Likelihood cost.** Measured on the fixture, full log density including the
trajectory solve:

| configuration | channels | epochs | dim | log density |
|---|---|---|---|---|
| companion source, DR3 pos/pm only | 2 | 66 | 8 | 0.026 ms |
| companion source, + DR2/DR32 | 6 | 66 | 8 | 0.026 ms |
| companion source, all available | 7 | 66 | 11 | 0.021 ms |
| full `G23HObs` on the primary | 12 | 213 | 17 | 0.035 ms |
| full `G23HObs` + Hipparcos IAD | 12 | 213 | 17 | 0.040 ms |
| `HGCAObs`-equivalent | 6 | 213 | 8 | 0.051 ms |

(`freeze_epochs=true`; 500 evaluations × 5 batches, minimum taken.) Everything is
20–50 μs and the spread between configurations is within run-to-run noise on a
shared machine — the DR3-only configuration is at best ~1.3–1.6× the full one, not
the order of magnitude "overhead" suggests. Dropping the DR2/DR32 channels costs
nothing because the DR2 five-parameter refit runs unconditionally inside
`if any_active`; only the *comparison* is dropped.

**Data.** For a DR3-position/PM-only observation, a greedy minimisation of the
catalog row finds **37 of 67 fields still required** — including every Hipparcos
and DR2/DR32 proper motion, error and correlation that the selected channels never
touch (`pmra_hip`, `pmra_hg`, `pmra_dr2`, `pmra_dr32`, all six `epoch_*`, …). Plus
`_g23h_merge_dr2_sidecar` is called whenever `gaia_id` is given, so the 300 MB
`G23H_DR2Transits` DataDep is needed even with no DR2 channel unless you supply
`astrometric_matched_observations_dr2` inline. The 14 GB `G23H_Catalog` DataDep is
needed unless you hand-build the row.

**So the honest statement of the overhead is: the catalog row and the epoch
marginalization, not the channel set.** A `GaiaDR3CompanionObs` type would avoid
both — but by re-deriving the scan-law refit rather than reusing it.

## Convenience constructor — sketched, not shipped

A `HGCAObs`-style helper would be:

```julia
GaiaDR3SourceObs(; source, about, gaia_id, kwargs...) =
    G23HObs(; host=source, companions=(about,), gaia_id,
              channels=(:ra_dr3, :dec_dr3), ueva_mode=:none,
              include_iad=false, include_rv=false, frame_shift=false, kwargs...)
```

It is **not trivially thin**, and shipping it would be actively harmful: it
cannot set `fluxratio = [0.0]`, so it papers over Gap 1 while leaving Gap 2's
wrong default in place — a one-line call that quietly models a resolved
secondary as its own bright primary. It becomes worth shipping the moment Gap 2
is fixed, and not before.

## Recommendation

**Do not port `GaiaDR3CompanionObs`.** Reasons, in order:

1. **v9 removed its reason to exist.** `host=`/`companions=`/`ref=` plus a
   shared frame already expresses "this catalog source is that body". The
   feature works today; what is missing is ergonomics, not machinery.
2. **DR3 has no public epoch astrometry**, so any DR3 companion likelihood must
   reproduce the five-parameter reduction from a scan-law forecast. `G23HObs`
   already does that, validated. A second implementation is a second thing to
   keep correct.
3. **DR4 supersedes the whole question.** When epoch astrometry is published,
   the companion case is `GaiaDR4AstromObs(scans; target=B, ref=Barycentre)` —
   which already builds and evaluates today against the `dr4-like-scans.csv`
   fixture (`refspecs == ((:B,), ())`, model dimension **4**: its own four
   nuisance parameters and nothing else). No catalog row, no 14 GB DataDep, no
   epoch marginalization, no five-parameter refit. Building a DR3-only type now
   buys a few months of a strictly worse answer.
4. **The two real defects are in `G23HObs`, not in a missing type.** Gaps 1 and
   2 would still be there with a new type beside them, and would still bite
   anyone modelling a wide pair with two `G23HObs`.

Suggested sequencing: fix Gap 1 (silent), then Gap 2 (`fluxratio=` keyword),
then ship the convenience constructor, then close #104 as superseded — closing
is your call, not the bot's.
