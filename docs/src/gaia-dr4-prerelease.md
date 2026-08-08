# Fitting Gaia DR4 Pre-Release Epoch Astrometry

In June 2026 ESA pre-released Gaia DR4 epoch astrometry for **12 sources**. Unlike the
[Gaia DR4 Simulation](@ref data-simulation-dr4) tutorial — which generates *synthetic*
epoch astrometry from a scan law — this tutorial fits **real Gaia epoch astrometry**.

We will fit **Gaia-4b**, a ~11.8 M<sub>jup</sub> companion on a ~571 day orbit
(Stefansson et al. 2025), from its along-scan measurements alone.

!!! note "Data provenance and credit"
    The data used here is derived from
    [`gaia-dr4-prerelease-epoch-astrometry_2026-06-26.zip`](https://www.cosmos.esa.int/web/gaia/dr4-prerelease).

    * Detailed information: <https://www.cosmos.esa.int/web/gaia/dr4-prerelease>
    * Data license: <https://www.cosmos.esa.int/web/gaia-users/license>
    * Acknowledgement and citation instructions:
      <https://gea.esac.esa.int/archive/documentation/GDR3/Miscellaneous/sec_credit_and_citation_instructions/>

    If you use this data in a publication, please follow ESA/DPAC's acknowledgement
    and citation instructions above.

## Setup

```julia
using Octofitter
using Distributions
using CairoMakie
using CSV, DataFrames
using Statistics
using Pigeons
```

!!! note "Two things to know before you start"
    [`GaiaDR4AstromObs`](@ref) takes epoch astrometry, references and variables only —
    it does **not** take `gaia_id=` and does not query the archive, so where you need the
    DR3 solution you query it yourself (shown below). And every mass is in **solar
    masses**.

## Loading the data

The pre-release ships as a VOTABLE, which needs a little massaging before it can be
used (see [Where this file came from](@ref dr4-prerelease-votable) at the bottom of this
page if you're curious). To keep this tutorial focused on the modelling, we ship a
ready-to-use CSV of the Gaia-4 measurements with the documentation:
[`gaia4_epoch_astrometry.csv`](https://github.com/sefffal/Octofitter.jl/blob/main/docs/src/gaia4_epoch_astrometry.csv).

```julia
const GAIA4_SOURCE_ID = 1457486023639239296
const REF_EPOCH_MJD   = 57936.375   # Gaia DR4 reference epoch, J2017.5

df = CSV.read(joinpath(@__DIR__, "gaia4_epoch_astrometry.csv"), DataFrame, comment="#")
df = mapcols(collect, df)   # see the note below
println("$(nrow(df)) CCD-level observations, $(length(unique(df.transit_id))) transits")
```

This gives one row per **CCD observation**: 1077 rows spread over 109 transits.
The columns are:

| Column | Units | Meaning |
|---|---|---|
| `transit_id` | — | Groups the CCD observations belonging to one field-of-view transit |
| `epoch` | MJD | Observation time (converted from TCB nanoseconds — see below) |
| `scan_pos_angle` | **radians** | Position angle of the scan, ψ |
| `centroid_pos_al` | mas | Along-scan centroid offset from the reference point |
| `centroid_pos_error_al` | mas | Formal along-scan centroid uncertainty |
| `ipd_error_al` | mas | Image parameter determination uncertainty |
| `parallax_factor_al` | — | Along-scan parallax factor (one value per transit) |
| `used_by_agis_al` | bool | Whether this CCD observation was used in the astrometric solution |
| `outlier_flag` | 0/1 | `0` where `used_by_agis_al` is true; `GaiaDR4AstromObs` skips rows with `outlier_flag > 0` |

!!! warning "Scan angles must be in radians"
    `GaiaDR4AstromObs` takes `sincos(scan_pos_angle)` directly, so this column **must be in
    radians**. The raw VOTABLE ships degrees; the conversion has already been applied in
    the CSV above. If you prepare your own table, remember `deg2rad`.

    A degree-valued column does not error — it produces a wrong, quietly plausible fit —
    so this is on the list to make safe by construction rather than by warning once the
    DR4 column conventions are final.

!!! note "`mapcols(collect, df)`"
    When Julia is started with multiple threads, `CSV.read` can return
    `SentinelArrays.ChainedVector` columns. These currently break `GaiaDR4AstromObs`'s
    per-epoch indexing (`AssertionError: wrong ChainedVectorIndex`). Materializing each
    column with `mapcols(collect, df)` makes the tutorial independent of thread count.

## From CCD observations to transit-level data

Gaia records roughly 9 usable CCD observations per field-of-view transit (SM, then
AF1–AF9). These are taken seconds apart, so for orbit fitting they carry essentially the
same astrometric information but are *not* statistically independent — their errors share
attitude and calibration systematics.

We therefore collapse each transit down to a single measurement. The reduction we
currently recommend is to drop the CCD observations not used by AGIS, then take the
per-transit median row:

```julia
function median_row(sdf; key=:centroid_pos_error_al)
    s = sort(sdf, key)
    n = nrow(s)
    lo, hi = (n + 1) ÷ 2, n ÷ 2 + 1     # odd: lo == hi (one row); even: the two middle rows
    out = DataFrame()
    for col in names(s)
        col == "transit_id" && continue  # grouping key; `combine` re-adds it
        v = s[!, col]
        if eltype(v) <: Real && !(eltype(v) <: Bool)
            out[!, col] = [(v[lo] + v[hi]) / 2]   # midpoint, exactly like median
        else
            out[!, col] = [v[lo]]                 # non-numeric (e.g. flags): lower-middle
        end
    end
    return out
end

gdf = groupby(subset(df, :used_by_agis_al), :transit_id)
transit_level_data = combine(median_row, gdf)

println("$(nrow(transit_level_data)) transits over ",
        round((maximum(transit_level_data.epoch) - minimum(transit_level_data.epoch))/365.25, digits=2),
        " yr")
```

For Gaia-4 this leaves **93 transits** spanning MJD 57038.6 – 58843.2 (4.94 yr, a bit
over three orbits at the published 571 day period).

!!! note "What this binning does and does not model"
    The median row keeps the *formal* single-CCD uncertainty. That does not credit the
    averaging of ~9 CCD observations, and it does not model the correlated
    attitude/calibration component shared within a transit — two omissions that push the
    resulting error bar in opposite directions. We recommend this reduction as the
    starting point today because it is simple and transparent about what it assumes.

    A refinement worth exploring is to take the median with bootstrap uncertainties,
    floored at the IPD uncertainty divided by `sqrt(N_CCD)`; the `ipd_error_al` column is
    carried in the CSV above so you can try it without regenerating the data. As always
    with pre-release data, check that your conclusions are not sensitive to the choice of
    reduction before publishing them.

## Building the likelihood

The astrometric nuisance parameters live in the likelihood's own `@variables` block. The
frame zero-point offsets (`ra_offset_mas`, `dec_offset_mas`) absorb the difference between
the DR3 catalogue position and the DR4 frame, and `pmra`/`pmdec` are the barycentre's
proper motion.

```julia
gaia_obs = GaiaDR4AstromObs(
    transit_level_data,
    target = Photocentre,   # what the catalogue source is: the blended, flux-weighted point
    ref    = Barycentre,    # ...measured against the system barycentre
    name="GaiaDR4",
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)   # mas
        ra_offset_mas  ~ Normal(0, 1000)               # frame zero-point (absorbs DR3<->DR4 offset)
        dec_offset_mas ~ Normal(0, 1000)
        pmra  ~ Uniform(-1000, 1000)                   # mas/yr
        pmdec ~ Uniform(-1000, 1000)
        ref_epoch = $REF_EPOCH_MJD
    end
)

# The observation carries data, references and variables only, so query the DR3
# solution explicitly when you need it:
sol = Octofitter._query_gaia_dr3(gaia_id=GAIA4_SOURCE_ID)
println("DR3 solution: plx=$(sol.parallax)  pmra=$(sol.pmra)  pmdec=$(sol.pmdec)")
```

The parallax the likelihood uses comes from the system's own `plx`; there is nothing to
forward into the observation's block.

!!! warning "`plx` is the barycentre's parallax, not the photocentre's"
    The frame's `plx` places the system's *barycentre*, while what the epoch astrometry
    actually constrains is the parallax of the **photocentre**. For a compact orbit the
    two are the same to well within the errors. For a wide one they are not, and a model
    that pins `plx` tightly to a catalog value — which was itself fitted as a single-star
    photocentre solution — can produce a posterior with a long curved ridge that samplers
    struggle to explore. If your orbit is wide compared to the mission baseline, give
    `plx` room and check the `plx`–`a`–`mass` corner for that ridge.

!!! note "`_query_gaia_dr3` is internal for now"
    The leading underscore is deliberate: the archive query layer will be replaced once
    the DR4 TAP endpoints are published, at which point the intent is that you pass a
    `gaia_id` and Octofitter fetches everything it needs. Until then, treat this call as
    a convenience whose signature may change.

## The model

```julia
orbit_ref_epoch = mean(gaia_obs.table.epoch)

A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(0.644, 0.02), lower=0.1)   # host mass [Msol] (Stefansson et al. 2025)
        flux = 1.0                           # the host sets the photocentre's flux scale
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(0.3mjup, 100mjup)  # companion mass [Msol]; `mjup` is a constant
        flux = 0.0                           # dark companion => the photocentre is the host
        a ~ Uniform(0, 10)                   # AU  (Gaia-4b is at ~1.17 AU)
        e ~ Uniform(0, 0.99)
        ω ~ UniformCircular()
        i ~ Sine()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()                # position angle at `epoch`
        # `θ` + `epoch` is a phase parametrization the orbit constructor accepts
        # directly, and an orbit's total mass comes from the hierarchy.
        epoch = $orbit_ref_epoch
    end
)

sys = System(
    name="Gaia4",
    bodies=[A, b],
    observations=[gaia_obs],
    variables=@variables begin
        plx ~ truncated(Normal(sol.parallax, sol.parallax_error), lower=1)
    end
)

model = Octofitter.LogDensityModel(sys; verbosity=4)
```

There is no system-level `M` any more: each body carries its own `mass`, and an orbit's
total mass comes from the hierarchy. A `Photocentre` target needs at least one body to
declare a flux, which is what `flux = 1.0` on `A` is for; `flux = 0.0` on `b` makes the
photocentre exactly the host, which is the right model for a dark companion.

## Initializing and sampling

```julia
init_chain = initialize!(model, (;
    plx = sol.parallax,
    bodies = (
        A = (; mass = 0.644),
        b = (; mass = 11.8mjup, a = 1.17, e = 0.1, i = 1.0),
    ),
    observations = (GaiaDR4 = (
        astrometric_jitter = 0.1,
        ra_offset_mas = 0.0,
        dec_offset_mas = 0.0,
        pmra = sol.pmra,
        pmdec = sol.pmdec,
    ),),
))
octoplot(model, init_chain)
```

Starting values nest under `bodies=` now, not `planets=`, and the host star is a body like
any other, so its mass is initialized in the same place.

Astrometry-only orbit fits are strongly multi-modal, so parallel tempering is the safer
choice here:

```julia
chain, pt = octofit_pigeons(
    model,
    n_chains=16,
    n_rounds=8,
    n_chains_variational=0,
    variational=nothing,
)
```

!!! tip "An astrometry-only fit wants parallel tempering"
    This fit runs under [`octofit_pigeons`](@ref) with a non-variational reference
    (`n_chains_variational=0`, `variational=nothing`), which is what an astrometry-only
    posterior wants. Dropping to HMC here is a real downgrade: seed it well
    ([`initialize!`](@ref) / [`startingpoints!`](@ref)) and check carefully for missed
    modes if you do.

## Results

```julia
q(v) = round.(quantile(v, (0.16, 0.5, 0.84)), sigdigits=5)

a    = vec(chain["b_a"])
e    = vec(chain["b_e"])
inc  = vec(chain["b_i"])
mb   = vec(chain["b_mass"]) ./ mjup           # Msol -> Mjup, for comparison with the paper
Mtot = vec(chain["A_mass"]) .+ vec(chain["b_mass"])
Pday = sqrt.(a.^3 ./ Mtot) .* 365.25

println("period  [day] : ", q(Pday), "   (Stefansson et al. 2025: 571.3 ± 1.4)")
println("a       [AU]  : ", q(a))
println("e             : ", q(e))
println("i       [deg] : ", round.(rad2deg.(quantile(inc, (0.16, 0.5, 0.84))), digits=1))
println("mass_b  [Mjup]: ", q(mb),   "   (Stefansson et al. 2025: 11.8 ± 0.7)")
```

The chain columns follow the bodies: `b_mass` and `A_mass`, both in solar masses.

And the usual plots:

```julia
octoplot(model, chain)
```

```julia
# A single posterior draw, plotted against the Gaia along-scan data.
# As with RV, this only works for individual draws: the Gaia points are "detrended"
# using the parameters of that particular draw.
#
# `GaiaDR4AstromObs` declares an `:along_scan` plot channel, so the generic panel covers
# the abscissa-vs-time view — restrict the PosteriorSeries to the draw you want:
idx = rand(1:size(chain, 1))
octoplot(Octofitter.PosteriorSeries(model, chain; ii=[idx]))

# `gaiastarplot` is the sky-plane version of the same draw: the reflex track with each
# transit's residual re-projected along its own scan angle.
Octofitter.gaiastarplot(model, chain, idx)
```

```julia
using PairPlots
octocorner(model, chain, small=true)
```

## [Where this file came from](@id dr4-prerelease-votable)

You do not need this section to follow the tutorial — it documents how
`gaia4_epoch_astrometry.csv` was derived from the official VOTABLE, so the reduction is
reproducible. We expect the eventual DR4 release to be reachable through a friendlier
query interface, at which point this step should disappear.

The pre-release VOTABLE stores **one row per transit**, with several columns holding
*arrays* of 10 values — the SM sample followed by AF1–AF9. `obs_time_tcb`,
`scan_pos_angle`, `centroid_pos_al`, `centroid_pos_error_al`, `ipd_error_al` and
`used_by_agis_al` are all array-valued, while `parallax_factor_al`, `ra0` and `dec0` are
per-transit scalars.

To get the flat, one-row-per-CCD-observation table Octofitter expects:

1. Select the rows for your `source_id`.
2. Flatten the array-valued columns, broadcasting the per-transit scalars across the
   10 CCD entries. Drop any CCD entry where a required quantity is masked.
3. Convert the time: `obs_time_tcb` is in **nanoseconds since JD 2455197.5 (TCB)**, so
   `jd = 2455197.5 + obs_time_tcb / 1e9 / 86400` and `mjd = jd - 2400000.5`.
   (Note this differs from the Gaia BH3 paper's published table, where the equivalent
   column is already in days.)
4. Convert `scan_pos_angle` from **degrees to radians**.
5. Set `outlier_flag = 0` where `used_by_agis_al` is true, else `1`.

For Gaia-4 (`source_id` 1457486023639239296) this yields 1077 CCD observations over 109
transits, of which 824 are flagged as used by AGIS. The reference point for
`centroid_pos_al` is `ra0 = 209.506326888 deg`, `dec0 = 31.695499700 deg`; these values,
along with the full conversion, are recorded in the comment header of the CSV itself.

!!! note "The pre-release includes Gaia BH3"
    The 12 pre-released sources are `2237987199365376`, `2309425390592896`,
    `10973744521070720`, `20694084440761600`, `60730287810150016`, `435469040545191680`,
    `1457486023639239296`, `1663617687609809280`, `3926186255616949504`,
    `3937211745905473024`, `4181040337841125632` and `4318465066420528000` — the last of
    which is **Gaia BH3** (α = 294.82786°, δ = +14.93098°, G = 11.23), with 77 transits.

    The [Gaia BH 3](@ref) example nonetheless still uses `astrom.dat`, the table published
    alongside the BH3 paper, which is already in the flat one-row-per-CCD-observation form
    rather than the array-valued pre-release form. The two agree where they overlap: for a
    shared `transit_id`, the paper table's scan angles reproduce the pre-release
    `scan_pos_angle` (declared in the VOTABLE as `unit="deg"`) to five significant figures.
    Re-basing that example onto the pre-release file has not been done.
