# [Detection Limits](@id detection-limits)

!!! warning
    This tutorial is a work in progress.

This guide shows how to calculate detection limits, in mass, or in photometry, as a function of orbital parameters for different combinations of data.

There are a few use cases for this:

* Mass limit vs semi-major axis given one or more images and/or contrast curves
* Mass limit vs semi-major axis given an RV non-detection
* Mass limit vs semi-major axis given proper motion anomaly from the Hipparcos-Gaia astrometry
* Any combination of the above

The recipe is the same in every case: build a model in which the companion's *brightness* is derived from its **mass** through an evolutionary model, fit the data, and read off which masses survive.

!!! note "About the target and the data on this page"
    The star is Gaia DR3 `6166183842771027328` (HIP 65808), a ~32 pc solar-type
    dwarf. It is chosen because its proper motion is a **non-detection**: the HGCA
    χ² is 17.8 and every Hipparcos/Hipparcos–Gaia/DR3 channel sits within ~3σ of
    the long-term trend, so the PMA posterior below is a genuine *upper limit* —
    the broad envelope the other data sets are compared against. Do not retarget
    this page onto a star with a real astrometric companion (HD 91312, used by
    [Fit Proper Motion Anomaly](@ref fit-pma), has a 16σ anomaly): the PMA
    posterior then pins itself to the top of the mass prior and the page stops
    demonstrating anything.

    So that the page builds offline it reads the Hipparcos–Gaia catalog subset and
    the cached scan-law forecast that ship with Octofitter's tests, together with
    the example L′ contrast image from [Fitting Images](@ref fit-images) — which is
    a stand-in, not a real observation of this star. For a real target you simply
    drop the `catalog=` and `forecast_table=` keywords and Octofitter fetches both
    itself.

```@example 1
using Octofitter
using OctofitterImages
using Distributions
using AstroImages
using CairoMakie
using PairPlots
using Statistics
using Pigeons
using Arrow, DataFrames, CSV # hide
```


## Photometry Model

We will need to decide on an atmosphere model to map image intensities into mass. Here we use the Sonora Bobcat cooling and atmosphere models which will be auto-downloaded by Octofitter:
```@example 1
const cooling_tracks = Octofitter.sonora_cooling_interpolator()
const sonora_temp_mass_L = Octofitter.sonora_photometry_interpolator(:Keck_L′)
```

!!! warning "Masses are in solar masses"
    There is a single mass unit throughout, M⊙, and these interpolators take M⊙ by
    default. `mjup` is a plain multiplicative constant, so a Jupiter-mass threshold is
    written `10mjup`, **not** `10` — a bare `10` means ten solar masses, which lands off
    the end of the grid and quietly returns `NaN`. Pass `mass_unit=:Mjup` to the
    constructors if you would rather keep a script in Jupiter masses.

## Proper Motion Anomaly Data

We start by defining and sampling from a model that only includes the Hipparcos-Gaia proper motion anomaly.

`HGCAObs` is a helper over the joint Gaia-Hipparcos likelihood [`G23HObs`](@ref), restricted to the six HGCA channels. It models the actual Gaia and Hipparcos scan epochs; `freeze_epochs=true` below fixes the epoch selection, which is much faster and is what makes a detection-limit sweep practical.

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass = system.M_pri
        flux_G  = 1.0     # the host defines the contrast scale in each band
        flux_Hp = 1.0
    end
)

B = Body(
    name="B",
    about=A,
    variables=@variables begin
        mass = system.M_sec       # M⊙
        a ~ LogUniform(1, 65)
        e ~ Uniform(0,0.9)
        ω ~ Uniform(0,2pi)
        i ~ Sine()                # The Sine() distribution is defined by Octofitter
        Ω ~ Uniform(0,pi)
        θ ~ Uniform(0,2pi)
        epoch = 57423.0           # epoch of the Gaia measurement
        flux_G  = 0.0             # dark to Gaia...
        flux_Hp = 0.0             # ...and to Hipparcos
    end
)

gaia_id = 6166183842771027328   # HIP 65808 -- see the note above
catalog = DataFrame(Arrow.Table(joinpath(@__DIR__, "..", "..", "test", "G23H-test-subset.feather"))) # hide
gost = CSV.read(joinpath(@__DIR__, "GOST-202.33684537818345--35.571314001279916-dr3.csv"), Table, normalizenames=true) # hide
forecast = Table( # hide
    epoch = Octofitter.jd2mjd.(gost.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_), # hide
    scanAngle_rad = gost.scanAngle_rad_, # hide
    parallaxFactorAlongScan = gost.parallaxFactorAlongScan, # hide
) # hide
nothing # hide

pma = HGCAObs(;
    gaia_id = gaia_id,
    target = A,
    blends = (B,),
    ref = Barycentre,
    freeze_epochs = true,       # faster but approximate; set false for real use
    catalog = catalog,          # hide
    forecast_table = forecast,  # hide
)
nothing # hide
```

The system block supplies the frame. `G23HObs` (and therefore `HGCAObs`) needs the reference point's proper motion, which means a **full** absolute frame: `plx, ra, dec, pmra, pmdec, rv, ref_epoch`. Declaring only some of them is a build-time error: Octofitter will not guess at half a frame.

```@example 1
cat = pma.catalog

HD_pma = System(
    name="limits_pma",
    bodies=[A, B],
    observations=[pma],
    variables=@variables begin
        M_pri ~ truncated(Normal(0.95, 0.05), lower=0.1)    # M⊙
        M_sec ~ LogUniform(0.2mjup, 65mjup)                # M⊙ -- note the units!

        plx ~ truncated(Normal(cat.parallax, cat.parallax_error), lower=0.1)
        ra  = $(cat.ra)
        dec = $(cat.dec)
        # Wide priors on the centre-of-mass proper motion: the data constrain these.
        pmra  ~ Uniform(cat.pmra_dr3 - 100, cat.pmra_dr3 + 100)
        pmdec ~ Uniform(cat.pmdec_dr3 - 100, cat.pmdec_dr3 + 100)
        rv = $(isnan(cat.radial_velocity) ? 0.0 : cat.radial_velocity * 1e3)   # m/s
        ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
    end
)
model_pma = Octofitter.LogDensityModel(HD_pma)
```

!!! warning "`pmra` and `pmdec` are frame variables"
    Use wide priors on `pmra`/`pmdec`. the data will constrain them.

Sample:
```@example 1
init_pma = initialize!(model_pma)
chain_pma, pt = octofit_pigeons(model_pma, n_chains=16, n_chains_variational=16, n_rounds=12)
nothing # hide
```

Plot the marginal mass vs. semi-major axis posterior with contours using PairPlots.jl. Note that `B_mass` is in solar masses, so we convert for the plot:
```@example 1
pairplot(
    PairPlots.Series(
        (;
            sma=log.(chain_pma[:B_a][:],),
            mass=log.(chain_pma[:B_mass][:] ./ mjup),
        ),
        label="PMA",
        color=Makie.wong_colors()[1],
    )=>(
        PairPlots.Scatter(markersize=3,alpha=0.35),
        PairPlots.Contour(sigmas=[1,3]),
        PairPlots.MarginStepHist(),
    ),
    labels=Dict(
        :sma=>"log Semi-major axis [au]",
        :mass=>"log Mass [Mⱼᵤₚ]"
    )
)
```


## Image Data

Now the same star, with an L′-band contrast image instead. This is where the flux/band unification pays off: the companion's brightness is a **body** variable, `flux_L`, derived from its mass — so the image constrains the mass directly.

```@example 1
download(
    "https://github.com/sefffal/Octofitter.jl/raw/main/docs/image-examples-1.fits",
    "image-examples-1.fits"
)
image = AstroImages.load("image-examples-1.fits").*2e-7 # units of contrast
img_dat_table = Table([
     (image=AstroImages.recenter(image), platescale=4.0, epoch=57423.6),
])
nothing # hide
```

```@example 1
B_img = Body(
    name="B",
    about=A,
    variables=@variables begin
        mass = system.M_sec

        # Calculate companion temperature from the cooling track and its mass
        tempK = $cooling_tracks(system.age, mass)
        # Calculate absolute magnitude
        abs_mag_L = $sonora_temp_mass_L(tempK, mass)
        # Deal with out-of-grid values by clamping to grid max and min.
        # NB: the threshold is `10mjup`, not `10` -- masses are in M⊙.
        abs_mag_L′ = if isfinite(abs_mag_L)
            abs_mag_L
        elseif mass > 10mjup
            8.2 # jump to absurdly bright
        else
            16.7 # jump to absurdly dim
        end
        # Calculate relative magnitude
        rel_mag_L = abs_mag_L′ - system.rel_mag + 5log10(1000/system.plx)
        # Convert to contrast -- the same units as the image. This
        # 10^(-0.4 Δmag) step is what makes `flux_L` a *linear* flux.
        flux_L = 10.0^(rel_mag_L/-2.5)

        a ~ LogUniform(1, 65)
        e ~ Uniform(0,0.9)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,pi)
        θ ~ Uniform(0,2pi)
        epoch = 57423.6
    end
)

image_data = ImageObs(
    img_dat_table,
    targets = (B_img,),
    ref     = A,
    band    = :L,
    name    = "imgdat-sim",
    variables=@variables begin
        # Optional parameters for marginalizing over instrument systematics:
        # Platescale uncertainty multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        platescale = 1.0
        # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
        northangle = 0.0
    end
)
nothing # hide
```



```@example 1
A_img = Body(
    name="A",
    variables=@variables begin
        mass = system.M_pri
    end
)

HD_img = System(
    name="limits_img",
    bodies=[A_img, B_img],
    observations=[image_data],
    variables=@variables begin
        # age ~ truncated(Normal(40, 15),lower=0, upper=200)
        age = 10                                           # Myr
        M_pri ~ truncated(Normal(0.95, 0.05), lower=0.1)    # M⊙
        # Mass of the secondary.
        # Make sure to pick only a mass range that is covered by your models.
        M_sec ~ LogUniform(0.55mjup, 65mjup)               # M⊙
        plx ~ truncated(Normal(cat.parallax, cat.parallax_error), lower=0.1)
        rel_mag = 5.65
    end
)
model_img = Octofitter.LogDensityModel(HD_img)
```

```@example 1
init_img = initialize!(model_img)
chain_img, pt = octofit_pigeons(model_img, n_chains=5, n_chains_variational=5, n_rounds=7)
nothing # hide
```

Plot mass vs. semi-major axis posterior:
```@example 1
vis_layers = (
    PairPlots.Contour(sigmas=[1,3]),
    PairPlots.MarginStepHist(),
)
pairplot(
    PairPlots.Series(
        (;
            sma=log.(chain_pma[:B_a][:],),
            mass=log.(chain_pma[:B_mass][:] ./ mjup),
        ),
        label="PMA",
        color=Makie.wong_colors()[1],
    )=>vis_layers,
    PairPlots.Series(
        (;
            sma=log.(chain_img[:B_a][:],),
            mass=log.(chain_img[:B_mass][:] ./ mjup),
        ),
        label="IMG",
        color=Makie.wong_colors()[2],
    )=>vis_layers,
    labels=Dict(
        :sma=>"log Semi-major axis [au]",
        :mass=>"log Mass [Mⱼᵤₚ]"
    )
)
```


## Image and PMA data

Combining the two is just a matter of listing both observations on the system. Observations are a flat list and each names its own references, so there is nothing else to wire up.

```@example 1
pma_joint = HGCAObs(;
    gaia_id = gaia_id,
    target = A,
    blends = (B_img,),
    ref = Barycentre,
    freeze_epochs = true,       # faster but approximate; set false for real use
    catalog = catalog,          # hide
    forecast_table = forecast,  # hide
)

HD_both = System(
    name="limits_both",
    bodies=[A, B_img],
    observations=[image_data, pma_joint],
    variables=@variables begin
        age = 10
        M_pri ~ truncated(Normal(0.95, 0.05), lower=0.1)
        M_sec ~ LogUniform(0.55mjup, 65mjup)
        rel_mag = 5.65

        plx ~ truncated(Normal(cat.parallax, cat.parallax_error), lower=0.1)
        ra  = $(cat.ra)
        dec = $(cat.dec)
        pmra  ~ Uniform(cat.pmra_dr3 - 100, cat.pmra_dr3 + 100)
        pmdec ~ Uniform(cat.pmdec_dr3 - 100, cat.pmdec_dr3 + 100)
        rv = $(isnan(cat.radial_velocity) ? 0.0 : cat.radial_velocity * 1e3)
        ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
    end
)
model_both = Octofitter.LogDensityModel(HD_both)

init_both = initialize!(model_both)
chain_both, pt = octofit_pigeons(model_both, n_chains=5, n_chains_variational=5, n_rounds=10)
nothing # hide
```

Compare all three posteriors to see limits:
```@example 1
vis_layers = (
    PairPlots.Contour(sigmas=[1,3]),
    PairPlots.MarginStepHist(),
)
pairplot(
    PairPlots.Series(
        (;
            sma=log.(chain_pma[:B_a][:],),
            mass=log.(chain_pma[:B_mass][:] ./ mjup),
        ),
        label="PMA",
        color=Makie.wong_colors()[1],
    )=>vis_layers,
    PairPlots.Series(
        (;
            sma=log.(chain_img[:B_a][:],),
            mass=log.(chain_img[:B_mass][:] ./ mjup),
        ),
        label="IMG",
        color=Makie.wong_colors()[2],
    )=>vis_layers,
        PairPlots.Series(
        (;
            sma=log.(chain_both[:B_a][:],),
            mass=log.(chain_both[:B_mass][:] ./ mjup),
        ),
        label="IMG + PMA",
        color=Makie.wong_colors()[3],
    )=>vis_layers,
    labels=Dict(
        :sma=>"log Semi-major axis [au]",
        :mass=>"log Mass [Mⱼᵤₚ]"
    )
)
```

For a systematic version of this analysis — injecting companions on a grid of mass and separation and measuring what fraction are recovered — see [Detection Completeness Mapping](@ref completeness).
