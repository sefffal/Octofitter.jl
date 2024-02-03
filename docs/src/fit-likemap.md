# [Fitting Images](@id fit-images)

```@example 1
using Octofitter
using Distributions
using OctofitterImages
using AstroImages
```

Typically one would load the likelihood maps from eg. FITS files like so:
```julia
images = AstroImages.recenter(load("image-example-1.fits",1))
```

However for this demonstration we will construct a synthetic likelihood map using 
a template orbit. We will create three peaks in two epochs.
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


# For calculations we save the full map in log likelihood
llmap_epoch1 = broadcast(-500:1:500, (-500:1:500)') do x, y
    logpdf(d, [x, y])
end

# For display we show just a subset in linear units:
lm1 = broadcast(x1 .+ (-50:50), y1 .+ (-50:1:50)') do x, y
    pdf(d, [x, y])
end
imview(lm1,clims=extrema)
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


# For calculations we save the full map in log likelihood
llmap_epoch2 = broadcast(-500:1:500, (-500:1:500)') do x, y
    logpdf(d, [x, y])
end

# For display we show just a subset in linear units:
lm2 = broadcast(x2 .+ (-50:50), y2 .+ (-50:1:50)') do x, y
    pdf(d, [x, y])
end
imview(lm2,clims=extrema)
```


Okay, we have our synthetic data. We now set up a `LogLikelihoodMap` object to contain our matrices of log likelihood values:
```@example 1
loglikemap = LogLikelihoodMap(
    (;
        epoch=epochs[1],
        map=AstroImages.recenter(AstroImage(llmap_epoch1)),
        platescale=1.0 # milliarcseconds/pixel of the map
    ),
    (;
        epoch=epochs[2],
        map=AstroImages.recenter(AstroImage(llmap_epoch2)),
        platescale=1.0 # milliarcseconds/pixel of the map
    )
);
```

We now create a one-planet model and run the fit using `octofit_pigeons`.
```@example 1

@planet b Visual{KepOrbit} begin
    a ~ Uniform(0, 10)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P*365.25 + 50420 # reference epoch for τ. Choose an MJD date near your data.
end loglikemap
@system Tutoria begin
    M ~ truncated(Normal(1.0, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b 
model = Octofitter.LogDensityModel(Tutoria)
chain, pt = octofit_pigeons(model, n_rounds=13)
```

Display the results:
```@example 1
using Plots : Plots
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
    (;x,y)
)
```

```@example 1
x = raoff.(els, loglikemap.table.epoch[2])
y = decoff.(els, loglikemap.table.epoch[2])
pairplot(
    (;x,y)
)
```
