# [Astrometry, PMA, and RV](@id astrom-pma-rv)

Octofitter.jl supports fitting orbit models to astrometric motion in the form of GAIA-Hipparcos proper motion anomaly (HGCA; [https://arxiv.org/abs/2105.11662](https://arxiv.org/abs/2105.11662)).
These data points are calculated by finding the difference between a long term proper motion of a star between the Hipparcos and GAIA catalogs, and their proper motion calculated within the windows of each catalog. This gives four data points that can constrain the dynamical mass & orbits of planetary companions (assuming we subtract out the net trend).

If your star of interest is in the HGCA, all you need is it's GAIA DR3 ID number. You can find this number by searching for your target on [SIMBAD](http://simbad.cds.unistra.fr).

For this tutorial, we will examine the star and companion [HD 91312 A & B](https://arxiv.org/abs/2109.12124) discovered by SCExAO. We will use their published astrometry and proper motion anomaly extracted from the HGCA.

We will also perform a model comparison: we will fit the same model to four different subsets of data to see how each dataset are impacting the final constraints. This is an important consistency check, especially with proper motion / absolute astrometry data which be susceptible to systematic errors.

The first step is to find the GAIA source ID for your object. For HD 91312, SIMBAD tells us the GAIA DR3 ID is `756291174721509376`.

```@contents
Pages = ["astrom-pma-rv.md"]
Depth = 5
```

!!! note 
    Work through the the [Proper Motion Anomaly](@ref fit-pma) tutorial before attempting a joint fit.

## Model: PMA Only

Initial setup:
```@example 1
using Octofitter, Distributions, Random
using CairoMakie
using PairPlots
using Pigeons
```

We begin by finding orbits that are consistent with the astrometric motion. Later, we will add in relative astrometry to the fit from direct imaging to further constrain the planet's orbit and mass.

Compared to previous tutorials, we now add a prior on the mass of the companion, called
`mass`. **All masses are in solar masses** — `mjup` is a plain multiplicative constant, so
a Jupiter-mass prior reads `LogUniform(0.5mjup, 1000mjup)`. There is also nothing to do to
"place a prior on the host mass rather than the system total mass": every body carries its
own `mass`, and an orbit's total mass comes from the hierarchy.

### The bodies

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.61, 0.1), lower=0.1)   # Msol
        flux_G  = 1.0    # the host sets the contrast scale
        flux_Hp = 1.0
    end
)

planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(0.5mjup, 1000mjup)   # Msol
        a ~ LogUniform(0.1, 20)
        e ~ Uniform(0, 0.999)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()                             # the Sine() distribution is defined by Octofitter
        Ω ~ Uniform(0, 2pi)
        # `θ` (position angle) + `epoch` is a phase parametrization the orbit
        # constructor understands directly.
        θ ~ Uniform(0, 2pi)
        epoch = 57423.0                        # epoch of the GAIA measurement
        flux_G  = 0.0                          # dark to Gaia…
        flux_Hp = 0.0                          # …and to Hipparcos
    end
)
nothing # hide
```

The two `flux_*` lines are contrast ratios against the host — set the host to `1.0` and
every other body's flux is a ratio — and they are how the model places the Gaia photocentre
and modulates the Hipparcos abscissa. Give a luminous, unresolved companion a real prior
here (`flux_G ~ Uniform(0, 1)`); omitting all four lines is legal and equivalent to `0.0`.

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
    target = A,
    blends = (planet_b,),
    ref = Barycentre,
    freeze_epochs = true,       # fast but approximate; drop for a production fit
    catalog = catalog,          # hide
    forecast_table = forecast,  # hide
)
hgca_obs.table.kind
```

### System Model & Specifying Proper Motion Anomaly

Now that we have our bodies, we create a system model to contain them. Observations are no
longer attached to a planet — they are a flat list on the system, and each one names what it
is a measurement *of* (`target`/`host`) and what it is measured *against* (`ref`).

We specify priors on `plx` as usual, but here we use the `gaia_plx` helper function to read
the parallax and uncertainty directly from the Gaia DR3 catalog using its source ID.

We also add parameters for the star's long term proper motion. This is usually close to the long term trend between the Hipparcos and GAIA measurements.

!!! warning "Use wide priors for pmra/pmdec with HGCA data"
    When fitting HGCA data, use **wide, uninformative priors** for `pmra` and `pmdec`, such as `Normal(0, 1000)` (0 ± 1000 mas/yr). **Do not** use Gaia DR3 proper motion values as informative priors—this would double-count the information since the HGCA already incorporates Gaia astrometry and will constrain the system's proper motion through the likelihood. The pmra/pmdec parameters represent the center-of-mass proper motion, which the HGCA measurements help determine.


```@example 1
ra_deg  = 158.30707896392835
dec_deg = 40.42555422701387

sys = System(
    name="HD91312_pma",
    bodies=[A, planet_b],
    observations=[hgca_obs],
    variables=@variables begin
        plx ~ gaia_plx(gaia_id=756291174721509376)

        # Priors on the center of mass proper motion
        pmra ~ Normal(-137, 10)
        pmdec ~ Normal(2, 10)

        # The rest of the absolute frame
        ra = $ra_deg
        dec = $dec_deg
        rv = 0.0                        # m/s
        ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
    end
)

model_pma = Octofitter.LogDensityModel(sys)
```

### Sampling from the posterior (PMA only)

Because proper motion anomaly data is quite sparse, it can often produce multi-modal posteriors. If your orbit already has several relative astrometry or RV data points, this is less of an issue.

```@example 1
chain_pma, pt = octofit_pigeons(model_pma, n_rounds=8)
display(chain_pma)
```

### Pair Plot
If we wish to examine the covariance between parameters in more detail, we can construct a pair-plot (aka. corner plot).

```@example 1
# Create a corner plot / pair plot.
# We can access any property from the chain specified in Variables
octocorner(model_pma, chain_pma, small=true)
```

## Model: PMA & Relative Astrometry

The first orbit fit to only Hipparcos/GAIA data was very unconstrained. We will now add six epochs of
relative astrometry (measured from direct images) gathered from the [discovery paper](https://arxiv.org/abs/2109.12124).

```@example 1
astrom_dat = Table(;
    epoch = [mjd("2016-12-15"), mjd("2017-03-12"), mjd("2017-03-13"), mjd("2018-02-08"), mjd("2018-11-28"), mjd("2018-12-15")],
    ra    = [133., 126., 127., 083., 058., 056.],
    dec   = [-174., -176., -172., -133., -122., -104.],
    σ_ra  = [07.0, 04.0, 04.0, 10.0, 10.0, 08.0],
    σ_dec = [07.0, 04.0, 04.0, 10.0, 20.0, 08.0],
    cor   = [0.2, 0.3, 0.1, 0.4, 0.3, 0.2]
)

astrom_obs = RelAstromObs(
    astrom_dat,
    target = planet_b,   # the body the offsets are of
    ref    = A,          # …and the body they are measured against
    name = "SCExAO",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)
scatter(astrom_obs.table.ra, astrom_obs.table.dec)
```

We use the same bodies as before, but now condition the model on the astrometry by adding it
to the system's `observations` list. We wrap it in
[`ObsPriorONeil2019`](@ref) — the observable-based prior of O'Neil et al. (2019), which
reweights the orbit prior by the Jacobian of the observables.

!!! warning "Attach only the wrapper"
    List **`ObsPriorONeil2019(astrom_obs)`** in `observations=`, not the wrapper *and*
    `astrom_obs`. Listing both counts the astrometry twice.

    The wrapper needs to know which orbit the prior applies to. It defaults to the wrapped
    observation's `target`, which is correct here and for relative RVs. A stellar-reflex
    `RadialVelocityObs(…; target=A, ref=Barycentre)` has no orbit of its own, so wrapping
    one requires `ObsPriorONeil2019(rvs; orbit=planet_b)` explicitly (`orbit=(b, c)` sums
    over several orbits).

```@example 1
sys_astrom = System(
    name="HD91312_pma_astrom",
    bodies=[A, planet_b],
    observations=[hgca_obs, ObsPriorONeil2019(astrom_obs)],
    variables=@variables begin
        plx ~ gaia_plx(gaia_id=756291174721509376)

        pmra ~ Normal(-137, 10)
        pmdec ~ Normal(2, 10)

        ra = $ra_deg
        dec = $dec_deg
        rv = 0.0
        ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
    end
)

model_pma_astrom = Octofitter.LogDensityModel(sys_astrom, verbosity=4)

chain_pma_astrom, pt = octofit_pigeons(model_pma_astrom, n_rounds=7, explorer=SliceSampler())
nothing # hide
```

```@example 1
octoplot(model_pma_astrom, chain_pma_astrom)
```

`octoplot` draws one time-series panel per plottable data channel: `RelAstromObs`
contributes separation and position angle (or Δα⋆/Δδ), and `G23HObs` contributes `pmra`
and `pmdec` — the five catalog proper motions against the model's reflex curve, with each
point's mission averaging window drawn as a horizontal bar.

Wrapping an observation changes nothing about how it is drawn: `ObsPriorONeil2019`
reweights the prior, and the astrometry inside it still contributes its points to the sky
panel and its own time-series panels, calibrated by the fitted `platescale`/`northangle`
exactly as an unwrapped dataset is. The panels are labelled with the wrapper's name
(`obspri_SCExAO` rather than `SCExAO`), which is also the name its calibration parameters
take in the chain.

## Model: PMA & Relative Astrometry & RVs

We now add in three additional epochs of stellar RVs.


```@example 1
rv_dat_abs = Table(;
    epoch = [mjd("2008-05-01"), mjd("2010-02-15"), mjd("2016-03-01")],
    rv    = [1300, 700, -2700],
    σ_rv  = [150, 150, 150]
)

rvlike = RadialVelocityObs(
    rv_dat_abs,
    target = A,            # the star's own reflex motion…
    ref = Barycentre,      # …against the system barycentre
    name = "SOPHIE",
    # Per dataset, and the default: with a full absolute frame like this one,
    # the secular (perspective) acceleration term is non-zero and modelled.
    # Pass `:data_corrected` for a series whose pipeline already removed it.
    secular_acceleration = :model,
    variables = @variables begin
        jitter ~ truncated(Normal(10, 5), lower=0)  # m/s
        offset ~ Normal(0, 1000)                    # m/s
    end
)
nothing # hide
```

`RadialVelocityObs` lives in core Octofitter and covers both cases: `target=A,
ref=Barycentre` is the star's reflex motion against the barycentre, and `target=b, ref=A`
is a companion's velocity relative to its host. No `using OctofitterRadialVelocity` is
needed for either — that package is for Celerite, `MarginalizedRVObs`, and the archive
loaders.

```@example 1
sys_rv_astrom = System(
    name="HD91312_pma_rv_astrom",
    bodies=[A, planet_b],
    observations=[hgca_obs, rvlike, ObsPriorONeil2019(astrom_obs)],
    variables=@variables begin
        plx ~ gaia_plx(gaia_id=756291174721509376)

        # Priors on the centre of mass proper motion
        pmra ~ Normal(-137, 10)
        pmdec ~ Normal(2, 10)

        ra = $ra_deg
        dec = $dec_deg
        rv = 0.0           # m/s
        ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
    end
)

model_pma_rv_astrom = Octofitter.LogDensityModel(sys_rv_astrom, verbosity=4)
chain_pma_rv_astrom, pt = octofit_pigeons(model_pma_rv_astrom, n_rounds=9, explorer=SliceSampler())
display(chain_pma_rv_astrom)
```

The mass vs. semi-major axis posterior is now much more constrained:
```@example 1
pairplot(
    (; a=chain_pma_rv_astrom["b_a"][:], mass=chain_pma_rv_astrom["b_mass"][:] ./ mjup) =>
        (
            PairPlots.Scatter(color=:red, markersize=5),
            PairPlots.MarginHist(),
            PairPlots.MarginQuantileText()
        ),
    labels=Dict(:mass=>"mass [Mⱼᵤₚ]", :a=>"sma. [au]"),
)
```

It is now useful to display the orbits projected onto the plane of the sky using `octoplot`.
This produces a sky panel plus one time-series panel per data channel — here, separation and
position angle from the astrometry, and radial velocity (with a phase-folded panel) from the
RVs.
```@example 1
octoplot(model_pma_rv_astrom, chain_pma_rv_astrom)
```

## Model:  Relative Astrometry & RVs (no PMA)

There is a final model we should consider: one using the RV and astrometry data, but not the
proper motion anomaly. With no absolute astrometry there is no absolute frame to declare, so
the system block shrinks to `plx` alone:
```@example 1
sys_final = System(
    name="HD91312_rv_astrom",
    bodies=[A, planet_b],
    observations=[rvlike, ObsPriorONeil2019(astrom_obs)],
    variables=@variables begin
        plx ~ gaia_plx(gaia_id=756291174721509376)
    end
)

model_rv_astrom = Octofitter.LogDensityModel(sys_final, verbosity=4)

chain_rv_astrom, pt = octofit_pigeons(model_rv_astrom, n_rounds=11)
nothing # hide
```

```@example 1
octoplot(model_rv_astrom, chain_rv_astrom)
```

## Model Comparison
Let's now display the constraints provided by each data set in a single corner plot
```@example 1
octocorner(
    model_pma,
    chain_pma,
    chain_pma_astrom,
    chain_rv_astrom,
    chain_pma_rv_astrom,
    small=true,
    axis=(;
        b_a = (;lims=(low=0, high=25))
    ),
    viz=(
        PairPlots.MarginDensity(),
        PairPlots.Scatter()
    )
)
```

We see that the constraints provided by the PMA, the astrometry, and the radial velocity data all individually overlap, and agree with the joint model constraint. 
This is means that none of the datasets are in tension with each other, which might suggest an issue with the data or with the modelling assumptions (e.g. single planet). 
