# [Fit Radial Velocity](@id fit-rv)

You can use Octofitter as a basic tool for fitting radial velocity data, by itself, or in combination with other kinds of data.
Multiple instruments (up to five) are supported.

!!! tip
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. To install it, run 
    `pkg> add OctofitterRadialVelocity`

Load the packages we'll need:
```@example 1
using Octofitter, OctofitterRadialVelocity, Distributions, PlanetOrbits, Plots
```

We can specify a table of radial velocity data manually by creating a [`RadialVelocityLikelihood`](@ref). An example of this is in [the other RV tutorial](@ref fit-rv-pma).

We can also directly load in data from the HARPS RVBank dataset:
```@example 1
rvs = OctofitterRadialVelocity.HARPS_rvs("GJ436")
```

Now, create a planet. Since we're only fitting radial velocity data, we
fix some of these parameters
```@example 1
@planet b Visual{KepOrbit} begin
    τ ~ UniformCircular(1.0)
    mass ~ truncated(Normal(21.3*0.00314558, 8*0.00314558),lower=0)
    a ~ truncated(Normal(0.028, 0.02), lower=0)
    i = pi/2
    e ~ Uniform(0, 0.5)
    ω ~ UniformCircular()
    Ω = 0
end

# Nor create the system
@system HD82134 begin
    M ~ truncated(Normal(0.425,0.009),lower=0)
    plx ~ 0.425

    # RV zero point of the system for instrument 1
    rv0_1 ~ Normal(0,10)
    # Jitter term for instrument 1
    jitter_1 ~ truncated(Normal(0,10),lower=0)
end rvs b
```

Build model:
```@example 1
model = Octofitter.LogDensityModel(HD82134; autodiff=:ForwardDiff, verbosity=4) # defaults are ForwardDiff, and verbosity=0
```


Sample from chains:
```@example 1
results = octofit(
    model, 0.8;
    adaptation =  1000,
    iterations =  1000,
    verbosity = 4,
    max_depth = 14
)
```

```@example 1
using Plots
octoplot(model, results)
```
![](HD82134-plot-grid.png)