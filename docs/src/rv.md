# [Fit RV and Proper Motion Anomaly](@id fit-rv-pma)

In this example, we will fit an orbit model to a combination of radial velocity and Hipparcos-GAIA proper motion anomaly for the star $\epsilon$ Eridani. We will use some of the radial velocity data collated in [Mawet et al 2019](https://iopscience.iop.org/article/10.3847/1538-3881/aaef8a).

!!! note
    The public RV archive loaders and the marginalized RV likelihood are supplied by the extension package OctofitterRadialVelocity.
    While v9 is an unregistered prerelease it must be installed from the `v2` branch, in the same `Pkg.add` as Octofitter itself —
    see [Installation](@ref install).

Datasets from radial velocity instruments are modelled together with separate jitters and instrumental offsets.


```@example 1
using Octofitter, OctofitterRadialVelocity, Distributions, PlanetOrbits, CairoMakie
using Pigeons

gaia_id = 5164707970261890560
```

We build the bodies first. The host is an
ordinary body, with a mass and-0-because Gaia and Hipparcos see a *photocentre*---
a flux in each of their bands. Setting the host's flux to 1 makes
every other body's flux a contrast ratio. Set dark companions to 0.

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(0.82, 0.02), lower=0.5, upper=1.5) # M⊙ (Baines & Armstrong 2011)
        flux_G  = 1.0   # Gaia G
        flux_Hp = 1.0   # Hipparcos Hp
    end
)

b = Body(
    name="b",
    about=A,
    # No relative astrometry is included since the planet has not yet been
    # directly detected.
    variables=@variables begin
        # For speed of example, we are fitting a circular orbit only.
        e = 0
        ω = 0.0
        # Masses are solar masses; `mjup` is a plain multiplicative constant.
        mass ~ Uniform(0, 3mjup)      # M⊙
        a ~ Uniform(3, 10)            # AU
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        # Phase: the mean anomaly at a reference epoch.
        M0 ~ Uniform(0, 2pi)
        epoch = 58849.0
        flux_G  = 0.0   # dark to Gaia and Hipparcos
        flux_Hp = 0.0
    end
)
nothing # hide
```

We load in data from one RV instrument.
We use `MarginalizedRVObs` instead of
[`RadialVelocityObs`](@ref) to analytically marginalize out
the radial velocity zero point of each instrument, saving one parameter.
```@example 1
hires_data = OctofitterRadialVelocity.HIRES_rvs("HD22049")
rvlike_hires = MarginalizedRVObs(
    hires_data;
    target=A, ref=Barycentre,   # the star's reflex motion against the barycentre
    name="HIRES",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)
nothing # hide
```

!!! warning "No offset or jitter is added for you"
    Octofitter never invents a prior. With `MarginalizedRVObs` the offset is
    integrated out analytically (and declaring one is an error), but the `jitter`
    line above is required — without it the model fits with no white-noise term
    at all.

We load the G23H data for this target. `target=` and `blends=` declare which
bodies this source is modelled from:
```@example 1
using Arrow, DataFrames, CSV # hide
# The docs build must touch neither the 14 GB G23H DataDep nor the network, so it # hide
# substitutes a catalog subset and a cached GOST scan forecast. Drop the # hide
# `catalog=`/`forecast_table=` keywords to fetch both for real. # hide
catalog = DataFrame(Arrow.Table(joinpath(@__DIR__, "..", "..", "test", "G23H-test-subset.feather"))) # hide
gost = CSV.read(joinpath(@__DIR__, "GOST-53.22829341517546--9.458168216292322-dr3.csv"), Table, normalizenames=true) # hide
forecast = Table( # hide
    epoch = Octofitter.jd2mjd.(gost.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_), # hide
    scanAngle_rad = gost.scanAngle_rad_, # hide
    parallaxFactorAlongScan = gost.parallaxFactorAlongScan, # hide
) # hide
pma = G23HObs(; gaia_id, target=A, blends=(b,), ref=Barycentre, freeze_epochs=true,
    catalog = catalog, # hide
    forecast_table = forecast, # hide
)
nothing # hide
```

`freeze_epochs=true` fixes which of Gaia's forecast scans were actually used,
instead of marginalizing over that selection. It is the fast-but-approximate
setting (the same one the [G23H tutorial](@ref fit-g23h) recommends for
exploration). Without it, the model gains one free parameter per forecast scan —
31 for this target. Set `freeze_epochs=false` (the default) for production fits.

The catalog row it fetched is available as `pma.catalog`, which is convenient for
centring the frame priors below:
```@example 1
cat = pma.catalog
(cat.parallax, cat.pmra_dr3, cat.pmdec_dr3)
```

Now the system. Absolute astrometry needs a **full absolute frame**: `plx`, `ra`,
`dec`, `pmra`, `pmdec`, `rv` and `ref_epoch` must all be declared together in the
system block (a partial frame is rejected at model-build time). `pmra`/`pmdec`
are the reference point's proper motion, which is what `G23HObs` compares its
catalog proper motions against.

```@example 1
sys = System(
    name="ϵEri",
    bodies=[A, b],
    observations=[pma, rvlike_hires],
    variables=@variables begin
        plx ~ gaia_plx(; gaia_id)
        ra  = $(cat.ra)
        dec = $(cat.dec)
        pmra  ~ Normal(cat.pmra_dr3, 10)
        pmdec ~ Normal(cat.pmdec_dr3, 10)
        rv = 0.0             # systemic radial velocity [km/s]
        ref_epoch = 57388.5  # Gaia DR3 reference epoch, J2016.0 (MJD)
    end
)
# Build model
model = Octofitter.LogDensityModel(sys)
```

!!! warning "Use wide priors for `pmra`/`pmdec` with absolute astrometry"
    The catalog proper motions already contain the companion's signal, so pinning
    the frame's proper motion tightly to a catalog value double-counts it. Give
    `pmra` and `pmdec` room to move, as above — and note that they are *frame*
    variables: they must be declared together with the rest of the absolute
    frame, and they cannot be declared alone.

!!! tip "Interpolating fitted values into a model"
    `ra = $(cat.ra)` uses `$` interpolation, which splices a value computed
    outside the model into a **derived** (`=`) variable. Prior (`~`) lines are
    evaluated in your own scope already, so they take plain expressions:
    `pmra ~ Normal(cat.pmra_dr3, 10)` — no `$` (and `$` on a `~` line is a
    syntax error).

Find good starting points and visualize the starting position + data:
```@example 1
init_chain = initialize!(model)
octoplot(model, init_chain)
```


```@example 1
using Pigeons
results, pt = octofit_pigeons(model, n_rounds=10, n_chains=10, n_chains_variational=0, explorer=SliceSampler());
nothing # hide
```

We can now plot the results with a multi-panel plot. [`octoplot`](@ref) derives
its panels from the model, so the sky orbit, the RV time series and its
phase-folded panel all appear without being asked for:
```@example 1
octoplot(model, results)
```

[`rvplot`](@ref) is the complementary view. Where `octoplot` draws many samples
from the posterior — one RV curve per draw, each carrying its own instrument
offsets, with the measurements left exactly as reported — `rvplot` renders a
**single** draw, which is what lets it put every instrument back on one axis:
```@example 1
rvplot(model, results)
```

We can see what the orbit looks like for the maximum a-posteriori sample (note, we would need to run an optimizer to get the true MAP value; this is just the MCMC sample with highest posterior density). Slicing the chain is how you plot a single draw:
```@example 1
i_max = argmax(vec(results[:logpost]))
res = octoplot(model, results[i_max:i_max, :, :])
Label(res.figure[0, 1], "Maximum a-posteriori orbit sample")
Makie.resize_to_layout!(res.figure)
res.figure
```


And a corner plot:
```@example 1
using PairPlots
octocorner(model, results, small=true)
```
