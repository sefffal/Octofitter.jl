# Fit Relative RV Data

Octofitter includes support for fitting relative radial velocity data. Currently this is only tested with a single companion. Please open an issue if you would like to fit multiple companions simultaneously.

The convention we adopt is that positive relative radial velocity is the velocity of the companion (exoplanets) minus the velocity of the host (star).

To fit relative RV data, start by creating a likelihood object:
```@example 1
using Octofitter
using OctofitterRadialVelocity
using CairoMakie
using Distributions

rel_rv_like = PlanetRelativeRVLikelihood(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11),
    )
```
The columns in the table should be . See the standard radial velocity tutorial for examples on how this data can be loaded from a CSV file.

The relative RV likelihood does not incorporate an instrument-specific RV offset. A jitter parameter called `jitter` can still be specified in the `@planet` block, as can parameters for a gaussian process model of stellar noise. Currently only a single instrument jitter parameter is supported. If you need to model relative radial velocities from multiple instruments with different jitters, please open an issue on GitHub.

Next, create a planet and system model, attaching the relative rv likelihood to the planet. Make sure to add a `jitter` parameter (optionally jitter=0) to the planet.

```@example 1
@planet b Visual{KepOrbit} begin
    a ~ truncated(Normal(10, 4), lower=0, upper=100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50420)

    jitter ~ LogUniform(1, 1000) # m/s
end rel_rv_like

@system ExampleSystem begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b

model = Octofitter.LogDensityModel(ExampleSystem)
```


```@example 1
using Random
rng = Random.Xoshiro(1234)

chain = octofit(rng, model)
```


```@example 1
octoplot(model, chain)
```