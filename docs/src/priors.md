# [Priors](@id priors)

All parameters to your model must have a prior defined.
You may provide any continuous, univariate distribution from the Distributions.jl.
A few useful distributions include:

* `Normal`
* `Uniform`
* `LogNormal`
* `LogUniform`
* `TrucatedNormal`
* `VonMises`

This pacakge also defines the `Sine()` distribution for e.g. inclination priors and `UniformCircular()` for periodic variables.
Internally, `UniformCircular()` creates two standard normal variables and finds the angle between them using `arctan`. This allows the sampler to smoothly cycle past the ends of the domain. You can specify a different circular domain than (0,2pi) by passing the size of the domain e.g. `τ = UniformCircular(1.0)`.

The VonMise distribution is notable but not commonly used. It is the analog of a normal distribution defined on a circular domain (-π, +π). If you have a Gaussian prior on an angular parameter, a Von Mises distribution is probably more appropriate.

Behind the scenes, Octofitter remaps your parameters to unconstrained domains using the Bijectors.jl (and corrects the priors accordingly). This is essential for good sampling efficiency with HMC based samplers.

This means that e.g. if you define the eccentricity prior as `e=Uniform(0,0.5)`, the sampler will actually generate values across the whole real line and transform them back into the `[0,0.5]` range before evaluating the orbit.
**It is therefore essential that your priors do not include invalid domains.**

For example, setting `a=Normal(3,2)` will result in poor sampling efficiency as sometimes negative values for semi-major axis will be drawn (especially if you're using the parallel tempered sampler).

Instead, for parameters like semi-major axis, eccentricity, parallax, and masses, you should truncate any distributions that have negative tails.
This can easily be accomplished with `TrauncatedNormal` or `Trunacted(dist, low, high)` for any arbitrary distribution.


## Kernel Density Estimate Priors

Octofitter has support for sampling from smoothed kernel density estimate priors. These are non-parametric distributions fit to a 1D dataset consisting of random draws. This is one way to include the output of a different model as the prior to a new model. That said, it's usually best to try and incorporate the model directly into the code. There are a few examples on GitHub of this, including atmosphere model grids, cooling tracks, etc.

### Using a KDE
First, we will generate some data. In the real world, you would load this data eg. from a CSV file.
```@example 1
using Octofitter, Distributions

# create a smoothed KDE estimate of the samples from a 10+-1 gaussian
kde = Octofitter.KDEDist(randn(1000).+10)
```

Note that in Octofitter the KDE will have its support truncated to the minimum and maximum values that occur in your dataset, ie. it doesn't allow for infinite long tails.

Now add it to your model as a prior:
```@example 1
@planet b Visual{KepOrbit} begin
    a ~ kde # Sample from the KDE here
    e ~ Uniform(0.0, 0.99)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50000)
end 
@system Tutoria begin 
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b
model = Octofitter.LogDensityModel(Tutoria)
chain = octofit(model)
```

We now examine the posterior and verify that it matches our KDE prior:
```@example 1
dat = chain[:b_a][:]
@show mean(dat) std(dat)
```

## Observable Based Priors

Octofitter implements observable-based priors from O'Neil 2019 for relative astrometry. You can fit a model to astrometry using observable-based priors using the following recipe:


```@example 1
using Octofitter, Distributions

astrom_like = PlanetRelAstromLikelihood(
    (;epoch=mjd("2020-12-20"), ra=400.0, σ_ra=5.0, dec=400.0, σ_dec=5.0)
)

@planet b Visual{KepOrbit} begin
    # For using with ObsPriors:
	P ~ Uniform(0.001, 1000)
    a = cbrt(system.M * b.P^2)

    e ~ Uniform(0.0, 1.0)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    mass ~ LogUniform(0.01, 100)

    τ ~ UniformCircular(1.0)
    tp =  b.τ*b.P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
end astrom_like ObsPriorAstromONeil2019(astrom_like);

@system System1 begin
    plx ~ Normal(21.219, 0.060)
	M ~ truncated(Normal(1.1, 0.2),lower=0)
end b

model = Octofitter.LogDensityModel(System1, verbosity=4)

```

```@example 1
octofit(
    model, 0.95;
    adaptation = 1000,
    iterations = 1000,
    verbosity = 4,
    max_depth = 10,
)
```