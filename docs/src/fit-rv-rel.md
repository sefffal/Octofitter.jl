# Fit Relative RV Data

Octofitter includes support for fitting relative radial velocity data. Currently this is only tested with a single companion. Please open an issue if you would like to fit multiple companions simultaneously.

The convention we adopt is that positive relative radial velocity is the velocity of the companion (exoplanets) minus the velocity of the host (star).

To fit relative RV data, start by creating a likelihood object:
```@example 1
using Octofitter
using OctofitterRadialVelocity
using CairoMakie
using Distributions

rv_dat_1 = Table(
    epoch=55000:100:57400,
    rv = [
         -24022.74
        -18571.33
        14221.56
        26076.89
        -459.26
        -26319.26
        -13430.96
        19230.96
        23580.26
        -6786.28
        -27161.78
        -7548.58
        23177.95
        19780.94
        -12738.39
        -26503.74
        -1249.19
        25844.47
        14888.83
        -17986.76
        -24381.49
        5119.22
        27083.2
        9174.18
        -22241.45
    ],
    # Hint! Type as \sigma + <TAB>
    σ_rv= fill(15000.0, 25),
)


rel_rv_like = PlanetRelativeRVLikelihood(
    rv_dat_1, 
    instrument_name="simulated data",
    variables = @variables begin
        jitter ~ LogUniform(0.1, 1000) # m/s
    end
)
```
See the standard radial velocity tutorial for examples on how this data can be loaded from a CSV file.

The relative RV likelihood does not incorporate an instrument-specific RV offset. A jitter parameter called can still be specified in the `@planet` block, as can parameters for a gaussian process model of stellar noise. Currently only a single instrument jitter parameter is supported. If you need to model relative radial velocities from multiple instruments with different jitters, please open an issue on GitHub.

Next, create a planet and system model, attaching the relative rv likelihood to the planet. Make sure to add a `jitter` parameter (optionally jitter=0) to the planet.

```@example 1
planet_1 = Planet(
    name="b",
    basis=RadialVelocityOrbit,
    likelihoods=[rel_rv_like],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1) # total mass in solar masses
        a ~ Uniform(0,10)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        τ ~ UniformCircular(1.0)
        P = √(this.a^3/this.M)
        tp =  this.τ*this.P*365.25 + 60000 # reference epoch for τ. Choose an MJD date near your data.

    end
)
sys = System(
    name = "Example-System",
    companions=[planet_1],
    likelihoods=[],
    variables=@variables begin
    end
)

model = Octofitter.LogDensityModel(sys)
```


### Initialize the model and verify starting point

```@example 1
init_chain = initialize!(model)

octoplot(model, init_chain)
```


```@example 1
using Random
rng = Random.Xoshiro(123)
chain = octofit(rng, model)
```


```@example 1
octoplot(model, chain, show_physical_orbit=true, mark_epochs_mjd=[mjd("2015-07-15")])
```