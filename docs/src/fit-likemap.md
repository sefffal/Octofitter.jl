# [Fitting Likelihood Maps](@id fit-likemap)

There are circumstances where you might have a 2D map of planet likelihood vs. position in the plane of the sky ($\Delta$ R.A. and Dec.). These could originate from:
* cleaned interferometer data 
* some kind of spectroscopic cube model fitting to detect planets
* some other imaginative modelling process I haven't thought of!

You can feed such 2D likelihood maps in to Octofitter. Simply pass in a list of maps and a platescale mapping the pixel size to a separation in milliarcseconds. 
You can of course also mix these likelihood maps with relative astrometry, radial velocity, proper motion, images, etc.

If your likelihood map is not centered on the star, you can specify offset dimensions as shown below.

!!! note
    Image modelling is supported in Octofitter via the extension package OctofitterImages. To install it, run 
    `pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterImages`

!!! note
    For simple models of interferometer data, OctofitterInterferometry.jl can already handle fitting point sources directly to visibilities.

```@example 1
using Octofitter
using Distributions
using OctofitterImages
using Pigeons
using AstroImages
```

Typically one would load your likelihood maps from eg. FITS files like so:
```julia
# If centred at the star:
image1 = AstroImages.recenter(AstroImages.load("image-example-1.fits",1))

# If not centered at the star:
image1 = AstroImages.load("image-example-1.fits")
image1_offset = AstroImage(
    image1,
    # Specify coordinates here:
    (
        # X coordinates should go from negative to positive.
        # The image should have +RA at the left.
        X(-4.85878653527304:1.0:95.14121346472696),
        Y(-69.0877222942365:1.0:30.9122777057635)
    )
    # Below, there is a platescale option. `platescale` multiplies
    # these values by a scaling factor. It can be 1 if the coordinates
    # above are already in milliarcseconds.
)
imview(image1_offset)
```
If you're using a FITS file, make sure to store your data in 64-bit floating point format.



For this demonstration, however, we will construct two synthetic likelihood maps using a template orbit. We will create three peaks in two epochs.
```@example 1

orbit_template = orbit(
    a = 1.0,
    e = 0.1,
    i = 0.0,
    ω = 0.0,
    Ω = 0.5,
    plx = 50.0,
    M = 1
)
epochs = [
    mjd("2024-01-30"),
    mjd("2024-02-29"),
]

## Create simulated data with three likelihood peaks at both epochs
x1,y1 = raoff(orbit_template,epochs[1]), decoff(orbit_template,epochs[1])
# The three peaks in our likelihood map
d1 = MvNormal([x1, y1], [
    5.0 0.2
    0.2 8.0
])
d2 = MvNormal([x1+8.0,y1+4], [
    5.0 0.6
    0.6 5.0
])
d3 = MvNormal([x1+9.0, y1-10.0], [
    6.0 0.6
    0.6 6.0
])
d = MixtureModel([d1, d2, d3], [0.5, 0.3, 0.2])

# Calculate a log-likelihood map over a +-50 mas patch around (x1, x2)
lm1 = broadcast(x1 .+ (-50:50), y1 .+ (-50:1:50)') do x, y
    logpdf(d, [x, y])
end

# Place in an AstroImage with appropriate offset coordinates
image1_offset = AstroImage(
    lm1,
    # Specify coordinates here:
    (
        # X coordinates should go from negative to positive.
        # The image should have +RA at the left.
        X(x1 .+ (-50:50)),
        Y(y1 .+ (-50:1:50))
    )
    # Below, there is a platescale option. `platescale` multiplies
    # these values by a scaling factor. It can be 1 if the coordinates
    # above are already in milliarcseconds.
)
imview(10 .^ image1_offset)
```

That was the first epoch. We now generate data for the second epoch:
```@example 1

x2,y2 = raoff(orbit_template,epochs[2]), decoff(orbit_template,epochs[2])
# The three peaks in our likelihood map
d1 = MvNormal([x2, y2], [
    5.0 0.2
    0.2 8.0
])
d2 = MvNormal([x2+10.0,y2], [
    5.0 0.6
    0.6 5.0
])
d3 = MvNormal([x2-4.0, y2-10.0], [
    6.0 0.6
    0.6 6.0
])
d = MixtureModel([d1, d2, d3], [0.5, 0.3, 0.2])



lm2 = broadcast(x2 .+ (-50:50), y2 .+ (-50:1:50)') do x, y
    logpdf(d, [x, y])
end
# Place in an AstroImage with appropriate offset coordinates
image2_offset = AstroImage(
    lm2,
    # Specify coordinates here:
    (
        # X coordinates should go from negative to positive.
        # The image should have +RA at the left.
        X(x2 .+ (-50:50)),
        Y(y2 .+ (-50:1:50))
    )
    # Below, there is a platescale option. `platescale` multiplies
    # these values by a scaling factor. It can be 1 if the coordinates
    # above are already in milliarcseconds.
)
imview(10 .^ image2_offset)
```


Okay, we have our synthetic data. We now set up a `LogLikelihoodMap` object to contain our matrices of log likelihood values:
```@example 1
loglikemap = LogLikelihoodMap(
    (;
        epoch=epochs[1],
        map=image1_offset,
        platescale=1.0 # milliarcseconds/pixel of the map
    ),
    (;
        epoch=epochs[2],
        map=image2_offset,
        platescale=1.0 # milliarcseconds/pixel of the map
    )
);
```

!!! note
    The likelihood maps will be interpolated using a simple bi-linear interpolation. 

We now create a one-planet model and run the fit using `octofit_pigeons`. This parallel-tempered sampler is slower than the regular `octofit`, but is recommended over the default Hamiltonian Monte Carlo sampler due to the multi-modal nature of the data. 
```@example 1

@planet b Visual{KepOrbit} begin
    a ~ Uniform(0, 10)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,60339.0)  # reference epoch for θ. Choose an MJD date near your data.
end loglikemap
@system Tutoria begin
    M ~ truncated(Normal(1.0, 0.1), lower=0.1)
    plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
end b 
model = Octofitter.LogDensityModel(Tutoria)
chain, pt = octofit_pigeons(model, n_rounds=10) # increase n_rounds until log(Z₁/Z₀) converges.
display(chain)
```

!!! note
    `octofit_pigeons` scales very well across multiple cores. Start julia with `julia --threads=auto` to make sure you have multiple threads available for sampling.

Display the results:
```@example 1
octoplot(model, chain)
```

Corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model,chain,small=true,)
```

And finally let's look at the posterior predictive distributions at both epochs:
```@example 1
els = Octofitter.construct_elements(chain,:b, :)
x = raoff.(els, loglikemap.table.epoch[1])
y = decoff.(els, loglikemap.table.epoch[1])
pairplot(
    (;x,y),
    axis=(
        x = (;
            lims=(low=100,high=-100)
        ),
        y = (;
            lims=(low=-100,high=100)
        )
    )
)
```

```@example 1
x = raoff.(els, loglikemap.table.epoch[2])
y = decoff.(els, loglikemap.table.epoch[2])
pairplot(
    (;x,y),
    axis=(
        x = (;
            lims=(low=100,high=-100)
        ),
        y = (;
            lims=(low=-100,high=100)
        )
    )
)
```


## Resume sampling for additional rounds

If you would like to add additional rounds of sampling, you may do the following:
```@example 1
pt = increment_n_rounds!(pt, 2)
chain, pt = octofit_pigeons(pt)
```

Updated corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model,chain,small=false,)
```
