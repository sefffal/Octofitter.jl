# [Fitting Images](@id fit-images)

One of the key features of Octofitter.jl is the ability to search for planets directly from images of the system. Sampling from images is much more computationally demanding than sampling from astrometry, but it allows for a few very powerful results:

1. You can search for a planet that is not well detected in a single image
By this, we mean you can feed in images of a system with no clear detections, and see if a planet is hiding in the noise based off of its Kepelerian motion.

2. Not detecting a planet in a given image can be almost as useful as a detection for constraining its orbit. 
If you have a clear detection in one epoch, but no detection in another, Octofitter can use the image from the second epoch to rule out large swathes of possible orbits.

Sampling from images can be freely combined with any known astrometry points, as well as astrometric acceleration. See advanced models for more details.

!!! note
    Image modelling is supported in Octofitter via the extension package OctofitterImages. To install it, run 
    `pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterImages`

## Preparing images
The first step will be to load your images. For this, we will use our AstroImages.jl package.

Start by loading your images:
```@example 1
using Octofitter
using OctofitterImages
using Distributions
using Pigeons
using AstroImages
using CairoMakie

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

First, we create a table of our image data that will be attached to the `Planet`:

```@example 1
image_data = ImageLikelihood(
    (band=:H, image=AstroImages.recenter(images[1]), platescale=10.0, epoch=1238.6),
    (band=:H, image=AstroImages.recenter(images[2]), platescale=10.0, epoch=1584.7),
    (band=:H, image=AstroImages.recenter(images[3]), platescale=10.0, epoch=3220.0),
    (band=:H, image=AstroImages.recenter(images[4]), platescale=10.0, epoch=7495.9),
    (band=:H, image=AstroImages.recenter(images[5]), platescale=10.0, epoch=7610.4),
)
```
Provide one entry for each image you want to sample from. Ensure that each image has been re-centered so that index `[0,0]` is the position of the star. Areas of the image where there is no data should be filled with `NaN` and will not contribute to the likelihood of your model. `platescale` should be the pixel scale of your images, in milliarseconds / pixel. `epoch` should be the Modified Julian Day (MJD) that your image was taken. You can use the `mjd("2021-09-09")` function to calculate this for you.
`band` should be a symbol that matches the name you supplied when you created the `Planet`.

By default, the contrast of the images is calculated automatically, but you can supply your own contrast curve as well by also passing `contrast=OctofitterImages.contrast_interp(AstroImages.recenter(my_image))`.

You can freely mix and match images from different instruments as long as you specify the correct platescale. 
You can also provide images from multiple bands and they will be sampled independently. If you wish to tie them together, see [Connecting Mass with Photometry](@ref mass-photometry).


Now specify the planet:
```@example 1
@planet b Visual{KepOrbit} begin
    H ~ Normal(3.8, 0.5)
    a ~ truncated(Normal(13, 4), lower=0, upper=100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,1238.6)
end image_data
```
Note how we also provided a prior on the photometry called `H`. We can put any name we want here, as long as it's used consistently throughout the model specification.

See [Fit PlanetRelAstromLikelihood](@ref fit-astrometry) for a description of the different orbital parameters, and conventions used.


Finally, create the system and pass in the planet.
```@example 1
@system HD82134 begin
    M ~ truncated(Normal(2.0, 0.1),lower=0)
    plx ~ truncated(Normal(45., 0.02),lower=0)
end b
```

If you want to search for two or more planets in the same images, just create multiple `Planet`s and pass the same images to each. You'll need to adjust the priors in some way to prevent overlap.

You can also do some very clever things like searching for planets that are co-planar and/or have a specific resonance between their periods. To do this, put the planet of the system or base period as variables of the system and derive the planet variables from those values of the system. 

## Sampling
Sampling from images is much more challenging than relative astrometry or proper motion anomaly, so the fitting process tends to take longer.

This is because the posterior is much "bumpier" with images.
One way this manifests is very high tree depths. You might see a sampling report that says `max_tree_depth_frac = 0.9` or even `1.0`.
To encourage the sampler to take larger steps and explore the images,
it's recommended to lower the target acceptance ratio to around 0.5±0.2 and also increase the number of adapataion steps.

```@example 1
model = Octofitter.LogDensityModel(HD82134)

chain, pt = octofit_pigeons(model, n_rounds=10)
display(chain)
```

!!! note
    `octofit_pigeons` scales very well across multiple cores. Start julia with `julia --threads=auto` to make sure you have multiple threads available for sampling.


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
fig = octoplot(model, chain)
```

With a bit of work, we can plot one of the images under the orbit.
```@example 1
fig = octoplot(model, chain)
ax = fig.content[1] # grap first axis in the figure

# We have to do some annoying work to get the image orientated correctly
# since we want the RA axis increasing to the left.
image_idx = 2
platescale = image_data.table.platescale[image_idx]
img = AstroImages.recenter(AstroImage(collect(image_data.table.image[image_idx])[end:-1:begin,:]))
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
fig = octoplot(model, chain)
ax = fig.content[1] # grap first axis in the figure

# We have to do some annoying work to get the image orientated correctly
# since we want the RA axis increasing to the left.
platescale = image_data.table.platescale[image_idx]
imgs = maximum(stack(image_data.table.image),dims=3)[:,:]
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
We start by plotting the marginal distribution of the flux parameter, `H`:


```@example 1
hist(chain["b_H"][:], axis=(xlabel="H", ylabel="counts"))
```


We can calculate an analog of the traditional signal to noise ratio (SNR) using that same histogram:
```@example 1
flux = chain["b_H"]
snr = mean(flux)/std(flux)
```
It might be better to consider a related measure, like the median flux over the interquartile distance. This will depend on your application.
