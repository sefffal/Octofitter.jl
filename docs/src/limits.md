## Detection Limits

!!! warning
    This tutorial is a work in progress.

This guide shows how to calculate detection limits, in mass, or in photometry, as a function of orbital parameters for different combinations of data.

There are a few use cases for this:

* Mass limit vs semi-major axis given one or more images and/or contrast curves
* Mass limit vs semi-major axis given an RV non-detection
* Mass limit vs semi-major axis given proper motion anomaly from the GAIA-Hipparcos Catalog of Accelerations
* Any combination of the above

We will once more use some sample data from the system [HD 91312 A & B](https://arxiv.org/abs/2109.12124) discovered by SCExAO. 


```@example 1
using Octofitter
using OctofitterImages
using OctofitterRadialVelocity
using Distributions
using Pigeons
using CairoMakie
using PairPlots
```

## Photometry Model

We will need to decide on an atmosphere model to map image intensities into mass. Here we use the Sonora Bobcat cooling and atmosphere models which will be auto-downloaded by Octofitter:
```@example 1
const cooling_tracks = Octofitter.sonora_cooling_interpolator()
const sonora_temp_mass_L = Octofitter.sonora_photometry_interpolator(:Keck_L′)
```

## Proper Motion Anomaly Data
We start by defining and sampling from a model that only includes proper motion anomaly data from the HGCA:
```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(1, 65)
    e ~ Uniform(0,0.9)
    ω ~ Uniform(0,2pi)
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)# ~ UniformCircular()
    mass = system.M_sec
    θ ~ Uniform(0,2pi)
    tp = θ_at_epoch_to_tperi(system,B,57423.0) # epoch of GAIA measurement
end
@system HD91312_pma begin
    M_pri ~ truncated(Normal(0.95, 0.05), lower=0) # Msol
    M_sec ~ LogUniform(0.2, 65) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol

    plx ~ gaia_plx(gaia_id=6166183842771027328)
            
    # Priors on the center of mass proper motion
    pmra ~ Normal(0, 1000)
    pmdec ~ Normal(0,  1000)
end HGCALikelihood(gaia_id=6166183842771027328) B
model_pma = Octofitter.LogDensityModel(HD91312_pma)
```

Sample:
```@example 1
using Pigeons
chain_pma, pt = octofit_pigeons(model_pma, n_chains=16, n_chains_variational=16, n_rounds=12);
nothing # hide
```

Plot the marginal mass vs. semi-major axis posterior with contours using PairPlots.jl:
```@example 1
pairplot(
    PairPlots.Series(
        (;
            sma=log.(chain_pma[:B_a][:],),
            mass=log.(chain_pma[:B_mass][:]),
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


```@example 1
using AstroImages
download(
    "https://github.com/sefffal/Octofitter.jl/raw/main/docs/image-examples-1.fits",
    "image-examples-1.fits"
)

# Or multi-extension FITS (this example)
image = AstroImages.load("image-examples-1.fits").*2e-7 # units of contrast

image_data = ImageLikelihood(
    (band=:L, image=AstroImages.recenter(image), platescale=4.0, epoch=57423.6),
)
```



```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(1, 65)
    e ~ Uniform(0,0.9)
    ω ~ Uniform(0,2pi)
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)
    mass = system.M_sec

    # Calculate planet temperature from cooling track and planet mass variable
    tempK = cooling_tracks(system.age, B.mass)
    # Calculate absolute magnitude
    abs_mag_L = sonora_temp_mass_L(B.tempK, B.mass)
    # Deal with out-of-grid values by clamping to grid max and min
    abs_mal_L′ = if isfinite(B.abs_mag_L)
        B.abs_mag_L
    elseif B.mass > 10 
        8.2 # jump to absurdly bright
    else
        16.7 # jump to absurdly dim
    end
    # Calculate relative magnitude
    rel_mag_L = B.abs_mal_L′ - system.rel_mag + 5log10(1000/system.plx)
    # Convert to contrast (same units as image)
    L = 10.0^(B.rel_mag_L/-2.5)

    θ ~ Uniform(0,2pi)
    tp = θ_at_epoch_to_tperi(system,B,57423.6)
end image_data 

@system HD91312_img begin
    # age ~ truncated(Normal(40, 15),lower=0, upper=200)
    age = 10
    M_pri ~ truncated(Normal(0.95, 0.05), lower=0) # Msol
    # Mass of secondary
    # Make sure to pick only a mass range that is covered by your models
    M_sec ~ LogUniform(0.55, 65) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol
    plx ~ gaia_plx(gaia_id=6166183842771027328)
    # Priors on the center of mass proper motion
    # pmra ~ Normal(0, 1000)
    # pmdec ~ Normal(0,  1000)
    rel_mag = 5.65
end  B
model_img = Octofitter.LogDensityModel(HD91312_img)
```

```@example 1
using Pigeons
chain_img, pt = octofit_pigeons(model_img, n_chains=5, n_chains_variational=5, n_rounds=7)
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
            mass=log.(chain_pma[:B_mass][:]),
        ),
        label="PMA",
        color=Makie.wong_colors()[1],
    )=>vis_layers,
    PairPlots.Series(
        (;
            sma=log.(chain_img[:B_a][:],),
            mass=log.(chain_img[:B_mass][:]),
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


# Image and PMA data

```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(1, 65)
    e ~ Uniform(0,0.9)
    ω ~ Uniform(0,2pi)
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)
    mass = system.M_sec

    # Calculate planet temperature from cooling track and planet mass variable
    tempK = cooling_tracks(system.age, B.mass)
    # Calculate absolute magnitude
    abs_mag_L = sonora_temp_mass_L(B.tempK, B.mass)
    # Deal with out-of-grid values by clamping to grid max and min
    abs_mal_L′ = if isfinite(B.abs_mag_L)
        B.abs_mag_L
    elseif B.mass > 10 
        8.2 # jump to absurdly bright
    else
        16.7 # jump to absurdly dim
    end
    # Calculate relative magnitude
    rel_mag_L = B.abs_mal_L′ - system.rel_mag + 5log10(1000/system.plx)
    # Convert to contrast (same units as image)
    L = 10.0^(B.rel_mag_L/-2.5)

    # L ~ Uniform(0,1)

    θ ~ Uniform(0,2pi)
    tp = θ_at_epoch_to_tperi(system,B,57423.6)
end image_data 

@system HD91312_both begin
    # age ~ truncated(Normal(40, 15),lower=0, upper=200)
    age = 10
    M_pri ~ truncated(Normal(0.95, 0.05), lower=0) # Msol
    # Mass of secondary
    # Make sure to pick only a mass range that is covered by your models
    M_sec ~ LogUniform(0.55, 65) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol
    plx ~ gaia_plx(gaia_id=6166183842771027328)
    # Priors on the center of mass proper motion
    pmra ~ Normal(0, 1000)
    pmdec ~ Normal(0,  1000)
    rel_mag = 5.65
end HGCALikelihood(gaia_id=6166183842771027328) B
model_both = Octofitter.LogDensityModel(HD91312_both)

using Pigeons
chain_both, pt = octofit_pigeons(model_both,n_chains=5,n_chains_variational=5,n_rounds=10)
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
            mass=log.(chain_pma[:B_mass][:]),
        ),
        label="PMA",
        color=Makie.wong_colors()[1],
    )=>vis_layers,
    PairPlots.Series(
        (;
            sma=log.(chain_img[:B_a][:],),
            mass=log.(chain_img[:B_mass][:]),
        ),
        label="IMG",
        color=Makie.wong_colors()[2],
    )=>vis_layers,
        PairPlots.Series(
        (;
            sma=log.(chain_both[:B_a][:],),
            mass=log.(chain_both[:B_mass][:]),
        ),
        label="IMG + PMA",
        color=Makie.wong_colors()[3],
    )=>vis_layers,
    labels=Dict(
        :sma=>"log₂ Semi-major axis [au]",
        :mass=>"log₂ Mass [Mⱼᵤₚ]"
    )
)
```
