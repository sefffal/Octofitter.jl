# [Fit Radial Velocity](@id fit-rv)

You can use Octofitter as a basic tool for fitting radial velocity data, by itself, or in combination with other kinds of data.
Multiple instruments (up to five) are supported.

!!! tip
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. You'll need
    to add both packaged to continue.

Load the packages we'll need:
```julia
using Octofitter, OctofitterRadialVelocity, Distributions, PlanetOrbits, Plots
```

We can specify a table of radial velocity data manually by creating a [`RadialVelocity`](@ref). An example of this is in [the other RV tutorial](@ref (@id fit-rv-pma).

We can also directly load in data from the HARPS RVBank dataset:
```julia
rvs = OctofitterRadialVelocity.HARPS_rvs("GJ436")
```

Now, create a planet. Since we're only fitting radial velocity data, we
fix some of these parameters
```julia
@named b = Planet{VisualOrbit}(
    Variables(
        τ = UniformCircular(1.0),
        mass = truncated(Normal(21.3*0.00314558, 8*0.00314558),lower=0),
        a=truncated(Normal(	0.028, 0.02), lower=0),
        i=pi/2,
        e = Uniform(0, 0.5),
        ω=UniformCircular(),
        Ω=0,
    ),
)

# Nor create the system
@named GJ436 = System(
    Variables(
        M = truncated(Normal(0.425,0.009),lower=0),
        plx = 0.425,

        # RV zero point of the system for instrument 1
        rv0_1 = Normal(0,10),
        # Jitter term for instrument 1
        jitter_1 = truncated(Normal(0,10),lower=0),
    ),
    rvs,
    b
)
```

Build model:
```julia
model = Octofitter.LogDensityModel(GJ436; autodiff=:ForwardDiff, verbosity=4) # defaults are ForwardDiff, and verbosity=0
```


Sample from chains:
```julia
results = Octofitter.advancedhmc(
    model, 0.5;
    adaptation =  600,
    iterations =  600,
    verbosity = 4,
    tree_depth = 8
)


timeplot(results, :b, "b_mass", :rv)
```