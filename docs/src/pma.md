# [Fit Proper Motion Anomaly](@id fit-pma)

Octofitter.jl supports fitting orbit models to astrometric motion in the form of GAIA-Hipparcos proper motion anomaly (HGCA; [https://arxiv.org/abs/2105.11662](https://arxiv.org/abs/2105.11662)).
These data points are calculated by finding the difference between a long term proper motion of a star between the Hipparcos and GAIA catalogs, and their proper motion calculated within the windows of each catalog. This gives four data points that can constrain the dynamical mass & orbits of planetary companions (assuming we subtract out the net trend).

If your star of interest is in the HGCA, all you need is it's GAIA DR3 ID number. You can find this number by searching for your target on [SIMBAD](http://simbad.cds.unistra.fr).

For this tutorial, we will examine the star and companion [HD 91312 A & B](https://arxiv.org/abs/2109.12124) discovered by SCExAO. We will use their published astrometry and proper motion anomaly extracted from the HGCA.

The first step is to find the GAIA source ID for your object. For HD 91312, SIMBAD tells us the GAIA DR3 ID is `756291174721509376`.

!!! tip "`freeze_epochs=true` for a quick approximation"
    The Gaia epoch *selections* are sampled parameters by default. We know what epochs and scan angles Gaia could possibly observe the target, and we know how many observations were ultimately included by the DR3 AGIS solution, but we don't know which subset of observations were actually used. When `freeze_epochs=false` (the default), we add a bunch of nuisance variables that let us marginalize over this uncertainty. Passing `freeze_epochs=true` fixes that selection and makes fitting
    dramatically faster, at the cost of an approximation---good for quick docs examples. 

## Fitting Astrometric Motion Only

Initial setup:
```@example 1
using Octofitter, Distributions, Random
```

We begin by finding orbits that are consistent with the astrometric motion. Later, we will
add in relative astrometry to the fit from direct imaging to further constrain the planet's
orbit and mass — see [Astrometry, PMA, and RV](@ref astrom-pma-rv).


### The bodies

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.61, 0.1), lower=0.1)   # Msol
        flux_G  = 1.0     # the host sets the contrast scale in each band
        flux_Hp = 1.0
    end
)

planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(0.5mjup, 1000mjup)   # Msol
        a ~ LogUniform(0.1, 100.0)             # au
        e ~ Uniform(0, 0.9)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()                             # the Sine() distribution is defined by Octofitter
        Ω ~ Uniform(0, 2pi)
        tp ~ Uniform(50000, 60000)             # time of periastron [MJD]
        flux_G  = 0.0     # dark to Gaia…
        flux_Hp = 0.0     # …and to Hipparcos
    end
)
nothing # hide
```

`flux_G` and `flux_Hp` are the *contrast ratios against the host* (hence `flux_G = 1.0` on
`A`), and they are what the code uses to model the Gaia photocentre and to modulate the
Hipparcos abscissa. A luminous, unresolved companion gets a real prior here — `flux_G ~
Uniform(0, 1)`. If you ommit them, the companion is assumed to be dark.

### Retrieving the HGCA

```@example 1
using Arrow, DataFrames, CSV # hide
# The docs build must touch neither the 14 GB G23H catalog nor the network, so it # hide
# substitutes a one-row catalog subset and a cached GOST scan forecast. Drop the # hide
# `catalog=`/`forecast_table=` keywords below to fetch both for real. # hide
catalog = DataFrame(Arrow.Table(joinpath(@__DIR__, "..", "..", "test", "G23H-test-subset.feather"))) # hide
gost = CSV.read(joinpath(@__DIR__, "GOST-158.30707896392835-40.42555422701387-dr3.csv"), Table, normalizenames=true) # hide
forecast = Table( # hide
    epoch = Octofitter.jd2mjd.(gost.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_), # hide
    scanAngle_rad = gost.scanAngle_rad_, # hide
    parallaxFactorAlongScan = gost.parallaxFactorAlongScan, # hide
) # hide
hgca_obs = HGCAObs(
    gaia_id = 756291174721509376,
    host = A,
    companions = (planet_b,),
    ref = Barycentre,
    freeze_epochs = true,       # fast but approximate; drop for a production fit
    catalog = catalog,          # hide
    forecast_table = forecast,  # hide
)
```


### System Model & Specifying Proper Motion Anomaly


We specify priors on `plx` as usual, using the `gaia_plx` helper to read the parallax and
uncertainty from the Gaia DR3 catalog by source ID.

We also add parameters for the star's long term proper motion. This is usually close to the
long term trend between the Hipparcos and GAIA measurements. Use wide priors on `pmra` and `pmdec` because the data are going to constrain them.


```@example 1
sys = System(
    name="HD91312_pma",
    bodies=[A, planet_b],
    observations=[hgca_obs],
    variables=@variables begin
        plx ~ gaia_plx(gaia_id=756291174721509376)

        # Priors on the center of mass proper motion
        pmra ~ Uniform(-137 - 100, -137 + 100)
        pmdec ~ Uniform(2 - 100, 2 + 100)

        # The rest of the absolute frame
        ra = $hgca_obs.catalog.ra
        dec = $hgca_obs.catalog.dec
        rv = 0.0   # barycentre's RV in [m/s] used perspective acceleration correction etc.
        ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
    end
)

model_pma = Octofitter.LogDensityModel(sys)
```

### Sampling from the posterior (PMA only)

Because proper motion anomaly data is quite sparse, it can often produce multi-modal
posteriors. If your orbit already has several relative astrometry or RV data points, this is
less of an issue.

```@example 1
chain_pma, pt = octofit_pigeons(model_pma, n_rounds=10)
display(chain_pma)
```

### Analysis

The first step is to look at the table output above generated by MCMCChains.jl.
The `rhat` column gives a convergence measure. Each parameter should have an `rhat` very close to 1.000.
If not, you may need to run the model for more iterations or tweak the parameterization of the model to improve sampling.
The `ess` column gives an estimate of the effective sample size.
The `mean` and `std` columns give the mean and standard deviation of each parameter.

The second table summarizes the 2.5, 25, 50, 75, and 97.5 percentiles of each parameter in the model.

### Pair Plot

If we wish to examine the covariance between parameters in more detail, we can construct a
pair-plot (aka. corner plot).

```@example 1
# Create a corner plot / pair plot.
# We can access any property from the chain specified in Variables
using CairoMakie: Makie
using PairPlots
octocorner(model_pma, chain_pma, small=true)
```

### Proper Motion Panels

`G23HObs` declares `pmra` and `pmdec` plot channels, so [`octoplot`](@ref) draws the five
catalog proper motions — Hipparcos, Hipparcos–Gaia, DR2, DR3−DR2 and DR3 — against the
model's reflex proper-motion curve, with a residual strip below each:

```@example 1
octoplot(model_pma, chain_pma)
```

The horizontal bar on each point is its *averaging window*, not an epoch uncertainty: a catalog
proper motion is a value computed over a full mission span.

### Posterior Mass vs. Semi-Major Axis

Given that this posterior is quite unconstrained, it is useful to make a simplified plot
marginalizing over all orbital parameters besides separation. [`dotplot`](@ref) is that
summary — mass against separation (or period, with `mode=:period`), coloured by
eccentricity, with marginal histograms:

```@example 1
Octofitter.dotplot(model_pma, chain_pma)
```


## See Also

- [Astrometry, PMA, and RV](@ref astrom-pma-rv) — the same target with relative astrometry
  and radial velocities added, plus a model comparison across data subsets.
- [Joint Gaia-Hipparcos (G23H)](@ref fit-g23h) — the full channel set this page restricts.
- [Hipparcos IAD](@ref fit-hipparcos) — the Hipparcos abscissae on their own.
