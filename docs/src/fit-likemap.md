# [Fitting Likelihood Maps](@id fit-likemap)

There are circumstances where you might have a 2D map of planet likelihood vs. position in the plane of the sky ($\Delta$ R.A. and Dec.). These could originate from:
* cleaned interferometer data 
* some kind of spectroscopic cube model fitting to detect planets
* some other imaginative modelling process I haven't thought of!

You can feed such 2D likelihood maps in to Octofitter. Simply pass in a list of maps and a platescale mapping the pixel size to a separation in milliarcseconds. 
You can of course also mix these likelihood maps with relative astrometry, radial velocity, proper motion, images, etc.

If your likelihood map is not centered on the star, you can specify offset dimensions as shown below.

!!! note
    Image modelling is supported in Octofitter via the extension package OctofitterImages.
    While v9 is an unregistered prerelease it must be installed from the `v2` branch, in the same `Pkg.add` as Octofitter itself —
    see [Installation](@ref install).

!!! note
    For simple models of interferometer data, OctofitterInterferometry.jl can already handle fitting point sources directly to visibilities.

```@example 1
using Octofitter
using Distributions
using OctofitterImages
using Pigeons
using AstroImages
using PlanetOrbits
using CairoMakie
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

The template is a PlanetOrbits.jl system built directly, rather than an Octofitter model — we only want somewhere to evaluate a position, not a fit. Note that `Body`, `Orbit` and `System` have to be qualified with `PlanetOrbits.` here, because Octofitter exports model-building types of the same names.

```@example 1

star_template = PlanetOrbits.Body(mass=1.0, name=:A)
planet_template = PlanetOrbits.Body(mass=0.0, name=:b)
orbit_template = PlanetOrbits.System(
    (star_template, planet_template),
    (PlanetOrbits.Orbit(planet_template; about=star_template,
        a = 1.0,
        e = 0.1,
        i = 0.0,
        ω = 0.0,
        Ω = 0.5,
        tp = 0.0,
    ),);
    plx = 50.0
)
epochs = [
    mjd("2024-01-30"),
    mjd("2024-02-29"),
]

## Create simulated data with three likelihood peaks at both epochs
sol1 = orbitsolve(orbit_template, epochs[1])
x1, y1 = raoff(sol1, :b, :A), decoff(sol1, :b, :A)
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

sol2 = orbitsolve(orbit_template, epochs[2])
x2, y2 = raoff(sol2, :b, :A), decoff(sol2, :b, :A)
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


Okay, we have our synthetic data. We now set up a `LogLikelihoodMapObs` object to contain our matrices of log likelihood values:

First, create a table of our likelihood map observations:
```@example 1
likemap_dat = Table(;
    epoch = epochs,
    map = [image1_offset, image2_offset],
    platescale = [1.0, 1.0] # milliarcseconds/pixel of the map
)

loglikemap = LogLikelihoodMapObs(
    likemap_dat,
    target = :b,       # the companion the map describes
    ref    = :A,       # the point index [0,0] of the map sits on
    name   = "GRAVITY",
    variables=@variables begin
        platescale = 1.0               # Platescale multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        northangle = 0.0               # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
    end
);
nothing # hide
```

`target` and `ref` take the full reference grammar: a `Body` model node, a `Symbol` naming one, or `Barycentre`. `ref` is the point each map is centred on (index `[0,0]`), which need not be the star.

Only a single `target` body is accepted per observation object.

!!! note
    The likelihood maps will be interpolated using a simple bi-linear interpolation. 

We now create a one-planet model and run the fit using `octofit_pigeons`. This parallel-tempered sampler is slower than the regular `octofit`, but is recommended over the default Hamiltonian Monte Carlo sampler due to the multi-modal nature of the data.

```@example 1

A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.1), lower=0.1)   # M⊙
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        a ~ Uniform(0, 10)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 60339.0    # reference epoch for θ. Choose an MJD date near your data.
    end
)

sys = System(
    name="Tutoria",
    bodies=[A, b],
    observations=[loglikemap],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
init_chain = initialize!(model)
chain, pt = octofit_pigeons(model, n_rounds=10) # increase n_rounds until log(Z₁/Z₀) converges.
display(chain)
```

Display the results:
```@example 1
using CairoMakie
octoplot(model, chain)
```

Corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model,chain,small=true,)
```

And finally let's look at the posterior predictive distributions at both epochs. `construct_system(model, chain)` rebuilds one PlanetOrbits system per posterior draw:
```@example 1
posteriors = construct_system(model, chain)

x = [raoff(orbitsolve(p, loglikemap.table.epoch[1]), :b, :A) for p in posteriors]
y = [decoff(orbitsolve(p, loglikemap.table.epoch[1]), :b, :A) for p in posteriors]
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
x = [raoff(orbitsolve(p, loglikemap.table.epoch[2]), :b, :A) for p in posteriors]
y = [decoff(orbitsolve(p, loglikemap.table.epoch[2]), :b, :A) for p in posteriors]
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

This picks up where the previous run left off — no need to start over — which is the
usual way to keep adding rounds until `log(Z₁/Z₀)` stops drifting.

Updated corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model,chain,small=false,)
```
