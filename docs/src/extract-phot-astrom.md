# Extracting Traditional Photometry and Astrometry

Though not its primary purpose, you can use Octofitter to extract traditional astrometry and photometry from one or more images. This uses the functionality in the [Fit Orbits to Images tutorial](@ref fit-images), but with a much simpler model. 

Instead of fitting an entire orbit, we will simply fit an X / Y position and brightness.


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
    "https://zenodo.org/records/6823071/files/HR8799.2021.fits?download=1",
    "HR8799-2021.fits"
)

# Or multi-extension FITS (this example)
image = AstroImages.load("HR8799-2021.fits")
```

You can preview the image using `imview` from AstroImages:
```@example 1
imview(image)
```

Note that to accurately extract astrometry and photometry, the input image should have already been convolved with the star or planet point spread function. If this isn't available, a convolution by a Gaussian or Airy disk might be an acceptable approximation.

## Build the model

First, we create a table of our image data that will be attached to the `Planet`:

```@example 1
imglike = ImageLikelihood(
    Table(
        image=[AstroImages.recenter(image)],
        platescale=[9.971],
        epoch=[mjd("2021")]
    ),
    variables=@variables begin
        # Planet flux in image units -- could be contrast, mags, Jy, or arb. as long as it's consistent with the units of the data you provide
        flux ~ Uniform(0, 1)
        # The following are optional parameters for marginalizing over instrument systematics:
        # Platescale uncertainty multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        platescale = 1.0
        # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
        northangle = 0.0
    end
)
```
Note that you can also supply a contrast curve or map directly. If not provided, a simple contrast curve will be calculated directly from the data.

Next create the simplest possible model of 2D position, plus a contrast variable matching the band name used in the `ImageLikelihood` above:
```@example 1
planet_b = Planet(
    name="b",
    basis=Visual{Octofitter.FixedPosition},
    likelihoods=[imglike],
    variables=@variables begin
        sep ~ Uniform(0, 2000)
        pa ~ Uniform(0,2pi)
        # Contrast ratio
    end
)

sys = System(
    name="sys",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        plx = 24.4620
    end
)

model = Octofitter.LogDensityModel(sys, verbosity=4)
```

## Sample from the model (locally)

If you already know where the planet is and you only want to extract astrometry from that known location, you can specify a starting point and use hamiltonian monte carlo as follows. This will be very very fast.
```@example 1
initialize!(model, (;
    planets=(;
        b=(;
            sep=1704,
            pa=deg2rad(70.63),
        )
    )
))
chain = octofit(model, iterations=10000)
```

## Sample from the model (globally)

You could also try sampling across the entire image, without necessarily specifying a starting position.
Note that if there are multiple candidates, taking the naive mean and standard deviation will average across all planets.
```@example 1
using Pigeons
initialize!(model)
chain, pt = octofit_pigeons(model, n_rounds=11)
```

## Access results
```@example 1
samples_sep = chain[:b_sep]
samples_pa = chain[:b_pa]
println("The median separation is ", median(samples_sep))

flux = chain[:b_images_flux]
println("The flux is ", mean(flux), " ± ", std(flux))
println("The \"SNR\" is ", mean(flux)/std(flux))
```

## Visualize
```@example 1
using CairoMakie, PairPlots
octocorner(model,chain)
```
