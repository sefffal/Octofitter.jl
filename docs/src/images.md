# [Fitting Images](@id fit-images)

One of the key features of Octofitter.jl is the ability to search for planets directly from images of the system. Sampling from images is much more computationally demanding than sampling from astrometry, but it allows for a few very powerful results:

1. You can search for a planet that is not well detected in a single image. By this, we mean you can feed in images of a system with no clear detections, and see if a planet is hiding in the noise based off of its Kepelerian motion.

2. Not detecting a planet in a given image can be almost as useful as a detection for constraining its orbit.  If you have a clear detection in one epoch, but no detection in another, Octofitter can use the image from the second epoch to rule out large swathes of possible orbits.

Sampling from images can be freely combined with any known astrometry points, as well as astrometric acceleration. See advanced models for more details.

!!! note
    Image modelling is supported in Octofitter via the extension package OctofitterImages.
    While v9 is an unregistered prerelease it must be installed from the `v2` branch, in the same `Pkg.add` as Octofitter itself —
    see [Installation](@ref install).

## Preparing images
The first step will be to load your images. For this, we will use our AstroImages.jl package.

Start by loading your images:
```@example 1
using Octofitter
using OctofitterImages
using Distributions
using AstroImages
using CairoMakie
using Pigeons

# Load individual iamges
# image1 = load("image1.fits")
# image2 = load("image2.fits")

# Or slices from a cube:
# cube = load("cube1.fits")
# image1 = cube[:,:,1] 

# Download sample images from GitHub
download(
    "https://github.com/sefffal/Octofitter.jl/raw/main/docs/image-examples-1.fits",
    "image-examples-1.fits"
)

# Or multi-extension FITS (this example)
images = AstroImages.load("image-examples-1.fits",:)
```

You can preview the image using `imview` from AstroImages:
```@example 1
# imshow2(image1, cmap=:magma) # for a single image
hcat(imview.(images, clims=(-1.0, 4.0))...)
```

Your images should either be convolved with a gaussian of diameter one λ/D, or be matched filtered. This is so that the values of the pixels in the image represent the photometry at that location. 

If you want to perform the convolution in Julia, see ImageFiltering.jl.

## Build the model

First, we create a table of our image observations:

```@example 1
image_dat = Table(;
    epoch = 56000 .+ [1238.6, 1584.7, 3220.0, 7495.9, 7610.4],
    image = [
        AstroImages.recenter(images[1]),
        AstroImages.recenter(images[2]), 
        AstroImages.recenter(images[3]),
        AstroImages.recenter(images[4]),
        AstroImages.recenter(images[5])
    ],
    platescale = [10.0, 10.0, 10.0, 10.0, 10.0]
)
nothing # hide
```

Provide one entry for each image you want to sample from. `platescale` should be the pixel scale of your images, in milliarseconds / pixel. `epoch` should be the Modified Julian Day (MJD) that your image was taken; you can use the `mjd("2021-09-09")` function to calculate this for you. (The sample images above ship with arbitrary epoch labels spanning about 17 years, so we offset them onto plausible MJDs. Only the *spacing* of the epochs affects the fit.)

Areas of the image where there is no data should be filled with `NaN` and will not contribute to the likelihood of your model.

Each image must be re-centered so that index `[0,0]` falls on the position of `ref` — the body the images are measured against, usually the host star. That is what `AstroImages.recenter` is doing above.

Now the bodies. The host star `A` is a body like any other, and the planet `b` orbits it (`about=A`):

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(2.0, 0.1), lower=0.1) # M⊙
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        # Brightness of this planet in the H band, in image units.
        flux_H ~ Normal(3.8, 0.5)

        a ~ truncated(Normal(13, 4), lower=0.1, upper=100)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()       # position angle at `epoch`
        epoch = 57238.6             # reference epoch for θ; near the first image
    end
)
nothing # hide
```

See [Fit Relative Astrometry](@ref fit-astrometry) for a description of the different orbital parameters, and conventions used.

Finally the observation itself, and the system:

```@example 1
image_obs = ImageObs(
    image_dat,
    targets = (b,),        # every source these images are modelled to contain
    ref     = A,           # the point index [0,0] of each image sits on
    band    = :H,          # read each target's `flux_H` variable
    name    = "SPHERE",
    variables=@variables begin
        # The following are optional parameters for marginalizing over instrument systematics:
        # Platescale uncertainty multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        platescale = 1.0
        # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
        northangle = 0.0
    end
)

sys = System(
    name="HD82134",
    bodies=[A, b],
    observations=[image_obs],
    variables=@variables begin
        plx ~ truncated(Normal(45., 0.02), lower=0.1)
    end
)
nothing # hide
```

### Planet Flux Variables

Each body you want to model should have a `flux_<band>` variable, e.g. `flux_H`.

```julia
c = Body(name="c", about=A, variables=@variables begin
    flux_H ~ Normal(1.2, 0.5)
    a ~ truncated(Normal(30, 8), lower=0.1, upper=100)
    # ... the rest of c's orbit ...
end)

image_obs = ImageObs(image_dat; targets=(b, c), ref=A, band=:H, name="SPHERE")
```

!!! warning "Two companions can land on top of each other"
    Nothing in the likelihood stops two `targets` from occupying the same pixels in the
    same epoch: the model would then explain one real source twice and leave the other
    unconstrained, and the posterior can genuinely go there. Multi-companion image
    fitting is usable but not yet guarded. Add an [`OrbitOrderPrior`](@ref), a
    [`NonCrossingPrior`](@ref), or explicit separation priors that keep the sources
    apart, and check the posterior for the degenerate mode before believing it.

A few consequences worth knowing:

* `band=:H` selects which `flux_<band>` variable each target is read from. You may omit `band=` only when the bodies declare exactly one band (a bare `flux = ...`); with several bands defined and no `band=`, you get an error listing them.
* The units of `flux_H` are the units of your image pixels — contrast, magnitudes, Jy, or arbitrary, as long as they are consistent. Setting the host's `flux_H = 1.0` makes every other body's flux a contrast ratio, but for images the host is usually not in `targets` at all (it is behind the coronagraph), so it needs no flux variable.
* `targets` is a *structural* statement about which sources the forward model contains. It is deliberately not inferred from "every body with a flux in this band": leaving `c` out is a different model from including `c` with a flux prior pushed near zero.
* `ref` is what the image is centred against. `ref=Barycentre` and `ref=b` (one companion measured against another) are both legal, and useful for hierarchical stellar systems.

By default, the contrast of the images is calculated automatically, but you can supply your own contrast curve as well by also passing a `contrast` column: `contrast=OctofitterImages.contrast_interp(AstroImages.recenter(my_image))`. That callable maps separation **in pixels** to the 1σ flux uncertainty there. A 2-D `contrastmap` column is honoured too.

You can freely mix and match images from different instruments as long as you specify the correct platescale — one `ImageObs` per instrument. 
You can also provide images from multiple bands and they will be sampled independently. If you wish to tie them together, see [Connecting Mass with Photometry](@ref mass-photometry).

You can also do some very clever things like searching for planets that are co-planar and/or have a specific resonance between their periods. To do this, put the period of the system or base period in the system variables and derive the body variables from those values of the system. 

## Sampling
Sampling from images is much more challenging than relative astrometry or proper motion anomaly, so the fitting process tends to take longer.

This is because the posterior is much "bumpier" with images.
One way this manifests is very high tree depths. You might see a sampling report that says `max_tree_depth_frac = 0.9` or even `1.0`.
To encourage the sampler to take larger steps and explore the images,
it's recommended to lower the target acceptance ratio to around 0.5±0.2 and also increase the number of adapataion steps.

```@example 1
model = Octofitter.LogDensityModel(sys)

chain, pt = octofit_pigeons(model, n_rounds=10)
display(chain)
```

!!! note
    `octofit_pigeons` scales very well across multiple cores. Start julia with `julia --threads=auto` to make sure you have multiple threads available for sampling.

!!! tip "Why parallel tempering here"
    Image posteriors are strongly multi-modal: a randomly drawn orbit almost always puts the planet on empty sky, where the likelihood is flat, and distinct orbits can explain the same set of point sources. Parallel tempering explores that far better than HMC does.

    If you do sample an image model with [`octofit`](@ref) instead, treat a single chain as exploring *one* mode: run [`initialize!`](@ref) first, then several chains from different starting points (`octofit(model, MCMCThreads(), ...)`), and check that they agree before believing a single-peaked posterior.


## Diagnostics
The first thing you should do with your results is check a few diagnostics to make sure the sampler converged as intended.

The acceptance rate should be somewhat lower than when fitting just astrometry, e.g. around the 0.6 target.

You can make a trace plot:
```@example 1
lines(
    chain["b_a"][:],
    axis=(;
        xlabel="iteration",
        ylabel="semi-major axis (aU)"
    )
)
```

And an auto-correlation plot:
```@example 1
using StatsBase
lines(
    autocor(chain["b_e"][:], 1:500),
    axis=(;
        xlabel="lag",
        ylabel="autocorrelation",
    )
)
```
For this model, there is somewhat higher correlation between samples. Some thinning to remove this correlation is recommended.


## Analysis

We can now view the orbit fit:
```@example 1
octoplot(model, chain)
```

`octoplot` returns an [`OctoPlotResult`](@ref), which carries the `figure` itself plus a named tuple of the `axes` it built. That is what lets us reach into the sky panel and draw underneath the orbits:

```@example 1
res = octoplot(model, chain)
fig = res.figure
ax = res.axes.sky.sky   # the sky-plane axis

# We have to do some annoying work to get the image orientated correctly,
# since we want the RA axis increasing to the left.
image_idx = 2
platescale = image_dat.platescale[image_idx]
img = AstroImages.recenter(AstroImage(collect(image_dat.image[image_idx])[end:-1:begin,:]))
imgax1 = dims(img,1) .* platescale
imgax2 = dims(img,2) .* platescale
h = heatmap!(ax, imgax1, imgax2, collect(img), colormap=:greys)
Makie.translate!(h, 0,0,-1) # Send heatmap to back of the plot

# Add colorbar for image
Colorbar(fig[1,2], h, label="image flux")

Makie.resize_to_layout!(fig)
fig
```

Another useful view would be the orbits over a stack of the maximum pixel values of all images.
```@example 1
res = octoplot(model, chain)
fig = res.figure
ax = res.axes.sky.sky

# We have to do some annoying work to get the image orientated correctly
# since we want the RA axis increasing to the left.
platescale = image_dat.platescale[image_idx]
imgs = maximum(stack(image_dat.image),dims=3)[:,:]
img = AstroImages.recenter(AstroImage(imgs[end:-1:begin,:]))
imgax1 = dims(img,1) .* platescale
imgax2 = dims(img,2) .* platescale
h = heatmap!(ax, imgax1, imgax2, collect(img), colormap=:greys)
Makie.translate!(h, 0,0,-1) # Send heatmap to back of the plot
Makie.resize_to_layout!(fig)
fig
```


## Pair Plot
We can show the relationships between variables on a pair plot (aka corner plot):

```@example 1
using CairoMakie, PairPlots
octocorner(model, chain, small=true)
```
Note that this time, we also show the recovered photometry in the corner plot.


## Assessing Detections
To assess a detection, we can treat all the orbital variables as nuisance parameters. 
We start by plotting the marginal distribution of the flux parameter, `flux_H`:

Body variables are named `<body>_<variable>` in the chain, so the planet's H-band flux is `b_flux_H`.

```@example 1
hist(chain["b_flux_H"][:], axis=(xlabel="flux", ylabel="counts"))
```


We can calculate an analog of the traditional signal to noise ratio (SNR) using that same histogram:
```@example 1
flux = chain["b_flux_H"]
snr = mean(flux)/std(flux)
```
It might be better to consider a related measure, like the median flux over the interquartile distance. This will depend on your application.
