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
    (epoch=5000, rv=-24022.74287804528, σ_rv=1500.0),
    (epoch=5100, rv=-18571.333891168735, σ_rv=1500.0),
    (epoch=5200, rv= 14221.562142944855, σ_rv=1500.0),
    (epoch=5300, rv= 26076.885281031347, σ_rv=1500.0),
    (epoch=5400, rv=  -459.2622916989299, σ_rv=1500.0),
    (epoch=5500, rv=-26319.264894263993, σ_rv=1500.0),
    (epoch=5600, rv=-13430.95547916007, σ_rv=1500.0),
    (epoch=5700, rv= 19230.962951723584, σ_rv=1500.0),
    (epoch=5800, rv= 23580.261108170227, σ_rv=1500.0),
    (epoch=5900, rv= -6786.277919597756, σ_rv=1500.0),
    (epoch=6000, rv=-27161.777481651112, σ_rv=1500.0),
    (epoch=6100, rv= -7548.583094927461, σ_rv=1500.0),
    (epoch=6200, rv= 23177.948014103342, σ_rv=1500.0),
    (epoch=6300, rv= 19780.94394128632, σ_rv=1500.0),
    (epoch=6400, rv=-12738.38520102873, σ_rv=1500.0),
    (epoch=6500, rv=-26503.73597982596, σ_rv=1500.0),
    (epoch=6600, rv= -1249.188767321913, σ_rv=1500.0),
    (epoch=6700, rv= 25844.465894406647, σ_rv=1500.0),
    (epoch=6800, rv= 14888.827293969505, σ_rv=1500.0),
    (epoch=6900, rv=-17986.75959839915, σ_rv=1500.0),
    (epoch=7000, rv=-24381.49393255423, σ_rv=1500.0),
    (epoch=7100, rv=  5119.21707156116, σ_rv=1500.0),
    (epoch=7200, rv= 27083.2046462065, σ_rv=1500.0),
    (epoch=7300, rv=  9174.176455190982, σ_rv=1500.0),
    (epoch=7400, rv=-22241.45434114139, σ_rv=1500.0),
    )
```
The columns in the table should be . See the standard radial velocity tutorial for examples on how this data can be loaded from a CSV file.

The relative RV likelihood does not incorporate an instrument-specific RV offset. A jitter parameter called `jitter` can still be specified in the `@planet` block, as can parameters for a gaussian process model of stellar noise. Currently only a single instrument jitter parameter is supported. If you need to model relative radial velocities from multiple instruments with different jitters, please open an issue on GitHub.

Next, create a planet and system model, attaching the relative rv likelihood to the planet. Make sure to add a `jitter` parameter (optionally jitter=0) to the planet.

```@example 1
@planet b KepOrbit begin
    a ~ Uniform(0,10)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P*365.25 + 6000 # reference epoch for τ. Choose an MJD date near your data.

    jitter ~ LogUniform(0.1, 1000) # m/s
end rel_rv_like

@system ExampleSystem begin
    M ~ truncated(Normal(1.2, 0.1), lower=0.1)
end b

model = Octofitter.LogDensityModel(ExampleSystem)
```


```@example 1
using Random
rng = Random.Xoshiro(123)
chain = octofit(rng, model)
```


```@example 1
octoplot(model, chain, ts=range(4500, 8000, length=200), show_physical_orbit=true)
```