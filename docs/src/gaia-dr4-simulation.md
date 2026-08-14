# [Simulating and Fitting Gaia DR4 Data](@id data-simulation-dr4)

This tutorial demonstrates how to use Octofitter to simulate Gaia DR4 data for a particular target. 

**Before reading this, take a look at the general [Generating and Fitting Simulated Data](@ref data-simulation) tutorial.**


----

To be as realistic as possible, we will select a particular star as the basis of our simulation. This will allow us to give reasonable estimates of
* the measurement epochs and scan angles (retrieved from GOST)
* the measurement noise, from that source's *measured* Gaia performance

> *"But wait, I just want to simulate a hypothetical planet, why do I have to pick a real star?"*

Gaia's sensitivity is far from uniform! If you're not sure what star to use, just pick a the Gaia ID of your favourite star that has a representative magnitude and colour.

## Setup

```@example 1
using Octofitter, Distributions
using CairoMakie, PairPlots, Pigeons
using Statistics
```

## Prepare Star Info and Noise

A realistic simulation needs three things: **when** Gaia looked at the star and at what
scan angle, the **parallax factor** of each of those looks, and **how precisely** each one
was measured. The first two come from GOST, Gaia's scan-law forecast service; the third
comes from the G23H catalog (Thompson et al. 2026), which carries a per-source, measured
along-scan noise calibration — the same one [`G23HObs`](@ref) uses.

```@example 1
gaia_id = 5064625130502952704
dr3 = gaia_dr3_solution(; gaia_id)   # α, δ, ϖ, G … from the DR3 archive (cached)
(; dr3.ra, dr3.dec, dr3.parallax, dr3.phot_g_mean_mag)
```

### The noise model

[`g23h_scan_uncertainty`](@ref) reads this source's calibrated along-scan uncertainties
out of the catalog:

```@example 1
# The docs build must not touch the 14 GB G23H catalog, so it substitutes the # hide
# handful of columns of this source's row that the uncertainty model reads. # hide
# Drop the `catalog=` keyword to read the real catalog (downloaded on first use). # hide
catalog = (; # hide
    gaia_source_id = gaia_id, # hide
    phot_g_mean_mag_dr3 = 6.941057, # hide
    sig_AL = 0.05813228279803717, # hide
    sig_att_radec = 0.07765105456802447, # hide
    sig_cal = 0.14901214069992602, # hide
    astrometric_n_good_obs_al_dr3 = 882, # hide
    astrometric_matched_transits_dr3 = 100, # hide
) # hide
σ = g23h_scan_uncertainty(;
    gaia_id,
    catalog, # hide
)
```

Gaia measures a source about nine times per field-of-view transit — the sky-mapper sample,
then AF1 through AF9 — and it is the **transit**, not the individual CCD observation, that
DR4 publishes an abscissa for. The three calibrated terms split along exactly that line:

* `σ_AL` and `σ_att` are per-CCD-observation and independent, so they *do* average down
  within a transit. Their quadrature sum `σ_formal` = 0.097 mas is the per-CCD formal
  error, and `σ_transit_formal` = `σ_formal/√n_ccd` = 0.0327 mas is the formal error of one
  transit-level abscissa (`n_ccd` = 8.82 for this source, from its own DR3 counts).
* `σ_calib` = 0.149 mas is a calibration error shared by all the CCD observations of a
  transit, so it does **not** average down — and it is not included in Gaia's formal
  uncertainties. For most sources it dominates the *actual* per-transit scatter,
  `σ_transit_true` = 0.153 mas.

!!! note "Why not just pick a number"
    These are calibrated per source against real Gaia performance and vary by a factor of a
    few from star to star. As a check on the scale: `σ_formal` reproduces the median
    *published per-CCD* `centroid_pos_error_al` of the three
    [Gaia DR4 pre-release](@ref dr4-prerelease-others) sources to 4–13%, and
    `σ_transit_formal` lands 20–35% below the per-transit error bars that the recommended
    bootstrap-median reduction assigns to those same real transits — about what the
    inefficiency of a median of nine CCD samples, relative to an optimal combination of
    them, predicts.

    If your target has no G23H calibration — the catalog does not cover everything — pass
    the three σ yourself and say where they came from.

### The scan geometry, and the table

[`gaia_dr4_transit_template`](@ref) queries GOST for the transits Gaia is forecast to make
of this position over the DR4 baseline and returns them in the format
[`GaiaDR4AstromObs`](@ref) reads:

```@example 1
transits = gaia_dr4_transit_template(;
    ra = dr3.ra,
    dec = dr3.dec,
    # We *simulate* the true scatter; see the note below on what a real DR4
    # table would instead quote.
    σ_al = σ.σ_transit_true,
    baseline = :dr4,
)
```

One row per transit, with `epoch` in MJD, `scan_pos_angle` ψ in **degrees** (the unit the
Gaia archive publishes it in, and the unit `GaiaDR4AstromObs` ingests), GOST's own
`parallax_factor_al`, and `centroid_pos_al` left at zero — Octofitter's data simulation
fills those in below. You can of course put your own measurements there instead.

A forecast is optimistic about *how many* transits you get: it lists every scheduled one,
while real DR4 loses some to dead time and more to AGIS's outlier rejection, which is
harshest for bright stars. Gaia-4 goes 122 forecast → 109 in the pre-release → 93 used by
AGIS. If that matters for your question, drop rows.

!!! note "Transit level, not CCD level"
    An earlier version of this tutorial expanded each visibility window into ~9 synthetic
    CCD observations and shrank the error bars to match. Don't do that: the CCD
    observations within one transit are taken seconds apart and share their attitude and
    calibration errors, so treating them as independent measurements overstates the
    astrometric information by up to √9. DR4 publishes one abscissa per transit, so
    simulating at transit level is both simpler and closer to the data you will actually
    fit. The [pre-release tutorial](@ref dr4-prerelease-reduction) shows the matching
    reduction for real CCD-level tables.

!!! note "The parallax factors are GOST's, not the Earth's"
    How good is a forecast? For Gaia-4, whose real DR4 transits we have, GOST's forecast
    covers the same MJD 56890–58843 span, matches 108 of the 109 real transits, and
    reproduces their published `scan_pos_angle` to 0.004° rms and their published
    `parallax_factor_al` to 5 × 10⁻⁵ rms.

    That last number is why `gaia_dr4_transit_template` carries GOST's own parallax factors
    through instead of recomputing them from an Earth ephemeris, as an earlier version of
    this tutorial did: Gaia observes from L2, about 0.01 AU from the geocentre, so an
    Earth-centred factor is off by 0.005 rms — a hundred times worse, and 0.07 mas at
    ϖ = 13 mas, about twice the per-transit precision.

!!! note "What a real DR4 table would quote"
    We put the *true* per-transit scatter in `centroid_pos_error_al`, so the data are
    scattered by exactly the uncertainty the likelihood is told about — the simplest
    self-consistent simulation. A real DR4 table instead quotes the formal error
    (`σ_transit_formal`), leaving the calibration term for the observation's
    `astrometric_jitter` to absorb. To simulate *that*, generate the data as below with
    `σ_al = σ.σ_transit_true`, then rebuild the observation from the simulated table with
    `centroid_pos_error_al` overwritten by `σ.σ_transit_formal`, and keep
    `astrometric_jitter` free — it should come back at ≈ `σ.σ_calib`.

```@example 1
# The Gaia DR4 reference epoch, J2017.5.
ref_epoch_mjd = 57936.375

gaiaIADobs = GaiaDR4AstromObs(transits;
    target = Photocentre,
    ref = Barycentre,
    name = "GaiaDR4",
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10) # mas
        ra_offset_mas  ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000) # mas/yr
        pmdec ~  Uniform(-1000, 1000) # mas/yr
        ref_epoch = $ref_epoch_mjd
    end
)

nothing # hide
```

!!! note "What the observation carries"
    `GaiaDR4AstromObs` takes data, references and variables and nothing else — no
    `gaia_id=`, and no archive query. It declares `target=Photocentre` (the system's
    flux-weighted point, the default) measured against `ref=Barycentre`, and reads its
    parallax from the system's own `plx`.

## Define the Model

Now, we define a model that incorporates this data. Bodies are model nodes in their own
right — the host star included — and every mass is in **solar masses**:

```@example 1
orbit_ref_epoch = mean(gaiaIADobs.table.epoch)

A = Body(
    name="A",
    variables=@variables begin
        mass = 1.0     # Msol
        flux = 1.0     # sets the flux scale the photocentre is weighted by
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(0.01mjup, 1000mjup)   # Msol
        flux = 0.0                              # dark companion
        a ~ LogUniform(0.01, 100)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        epoch = $orbit_ref_epoch
    end
)

sys = System(
    name="target_1",
    bodies=[A, b],
    observations=[gaiaIADobs],
    variables=@variables begin
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```

A `Photocentre` target needs at least one body to declare a flux, so `A` carries
`flux = 1.0`; every other body's `flux` is then a contrast ratio against it. With
`b.flux = 0.0` the photocentre is exactly the host, which is the right model for a dark
companion.

## Simulating a Planet

Now we have defined our template model! Normally, we would just fit this...
Here, we will hook into Octofitter's simulation capabilities to generate a new model for us with a simulated planet.



## Generate Synthetic Data

We have two choices for generating simulated data:
1. Draw values from the priors
2. Specifying values manually

We will look at each.

### 1. Draw values from the priors

We can draw a value from the priors like so:
```@example 1
params_to_simulate = Octofitter.drawfrompriors(model.system)
```


### 2. Specifying values manually
We can also specify all values for the simulation manually. This process is a bit more involved. 


!!! warning
    Note that the output below is just an example, you must generate your own template from your model and modify it as needed. The exact structure is not garuanteed to be stable between versions of Octofitter.

Note the shape: system variables at the top level, then **`bodies`** — the host star
included — then `observations`.

```@example 1
params_to_simulate = (
    plx = 16.138978209522527,
    bodies = (
        A = (
            mass = 1.0,
            flux = 1.0,
        ),
        b = (
            mass = 10.0 * Octofitter.mjup,
            a = 2.0,
            e = 0.01,
            ω = 0.01,
            i = 0.01,
            Ω = 0.6,
            θ = 3.5226970272017826,
            epoch = orbit_ref_epoch,
            flux = 0.0,
        ),
    ),
    observations = (
        GaiaDR4 = (
            astrometric_jitter = 0.015590772157762368,
            ra_offset_mas = 0.05413791355838311,
            dec_offset_mas = 0.0889816167388366,
            pmra = 5.301074192615374,
            pmdec = -24.188882826919325,
            ref_epoch = 57936.375,
        ),
    ),
)
nothing # hide
```


## Generate synthetic system with simulated data

We will call `Octofitter.generate_from_params` to generate a model with new synthetic observations.
If you set `add_noise = true`, the generated data points will have scatter according to the `centroid_pos_error_al` specified above. 


```@example 1
sim_system = Octofitter.generate_from_params(model.system, params_to_simulate; add_noise=true)
sim_model = Octofitter.LogDensityModel(sim_system)
```


## Fit simulated data

Follow the usual flow documented elsewhere.

Find starting point for MCMC (via variational approximation)
```@example 1
init_chain = initialize!(sim_model)
octoplot(sim_model, init_chain)
```


Fit:
```@example 1
chain, pt = octofit_pigeons(sim_model, n_rounds=9)
```

!!! note "Why parallel tempering here"
    A single epoch-astrometry time series constrains the orbit only through the
    along-scan abscissa, and that leaves genuinely separated solutions — most
    conspicuously the ±180° ambiguity in Ω and the degeneracy between period and
    the fraction of an orbit the scan window covers. [`octofit_pigeons`](@ref) runs
    a ladder of tempered chains down to the prior, where those solutions merge, so
    replicas can move between them and each stays represented in the posterior.

    [`octofit`](@ref) will sample this model too and is considerably cheaper. If you
    use it, seed it with [`initialize!`](@ref) as above and run a few chains from
    different starting points, so a mode one chain missed shows up as disagreement
    rather than as false confidence.


```@example 1
octoplot(sim_model, chain)
```

`GaiaDR4AstromObs` declares one plot channel (the along-scan abscissa), so `octoplot` draws
the data, the modelled abscissae, a residual strip and a residual histogram automatically.
Restrict the [`PosteriorSeries`](@ref) to get the "one particular draw" version:

```@example 1
octoplot(Octofitter.PosteriorSeries(sim_model, chain; ii=[1]))
```

[`gaiastarplot`](@ref) shows that single draw in the sky plane instead — the reflex track
with each transit's along-scan residual re-projected along its own scan angle, which is
where the one-dimensional nature of a Gaia measurement becomes visible:

```@example 1
Octofitter.gaiastarplot(sim_model, chain, 1)
```

Finally, compare the recovered orbit against the truth: rebuild the
`PlanetOrbits.System` for each draw with `construct_system` and plot the tracks directly,
which is what the sky panel does internally:

```@example 1
ts = range(minimum(transits.epoch), minimum(transits.epoch) + 4000, length=300)

fig = Figure(size=(500,500))
ax = Axis(fig[1,1], xlabel="Δα⋆ [mas]", ylabel="Δδ [mas]",
          aspect=DataAspect(), xreversed=true)

for i in rand(1:size(chain,1), 50)
    posys = construct_system(sim_model, chain, i)
    traj = orbitsolve(posys, ts)
    lines!(ax, [raoff(traj[k], :b, :A) for k in eachindex(ts)],
               [decoff(traj[k], :b, :A) for k in eachindex(ts)],
           color=(:black, 0.1))
end

# the true orbit, in red
true_sys = construct_system(sim_model, params_to_simulate)
true_traj = orbitsolve(true_sys, ts)
lines!(ax, [raoff(true_traj[k], :b, :A) for k in eachindex(ts)],
           [decoff(true_traj[k], :b, :A) for k in eachindex(ts)],
       color=:red, linewidth=3, label="truth")
scatter!(ax, [0], [0], color=:orange, markersize=12)
axislegend(ax)
fig
```

```@example 1
octocorner(sim_model, chain, small=true)
```
