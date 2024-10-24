# [Fit Radial Velocity and Astrometry](@id fit-rv)

You can use Octofitter to jointly fit relative astrometry data and radial velocity data. 
Below is an example. For more information on these functions, see previous guides.


Import required packages
```@example 1
using Octofitter
using OctofitterRadialVelocity
using CairoMakie
using PairPlots
using Distributions
using PlanetOrbits
```


We now use PlanetOrbits.jl to create sample data. We start with a template orbit and record it's positon and velocity at a few epochs.
```@example 1
orb_template = orbit(
    a = 1.0,
    e = 0.7,
    i= pi/4,
    Ω = 0.1,
    ω = 1π/4,
    M = 1.0,
    plx=100.0,
    m =0,
    tp =58829
)
Makie.lines(orb_template)
```


Sample position and store as relative astrometry measurements:
```@example 1
epochs = [58849,58852,58858,58890]
astrom = PlanetRelAstromLikelihood(Table(
    epoch=epochs,
    ra=raoff.(orb_template, epochs),
    dec=decoff.(orb_template, epochs),
    σ_ra=fill(1.0, size(epochs)),
    σ_dec=fill(1.0, size(epochs)),
))
```

And plot our simulated astrometry measurments:
```@example 1
fig = Makie.lines(orb_template,)
Makie.scatter!(astrom.table.ra, astrom.table.dec)
fig
```


Generate a simulated RV curve from the same orbit:
```@example 1
using Random
Random.seed!(1)

epochs = 58849 .+ (20:20:660)
planet_sim_mass = 0.001 # solar masses here
rvlike = StarAbsoluteRVLikelihood(Table(
    epoch=epochs,
    rv=radvel.(orb_template, epochs, planet_sim_mass) .+ randn.() .+ 100randn(),
    σ_rv=fill(5.0, size(epochs)),
))
Makie.scatter(rvlike.table.epoch[:], rvlike.table.rv[:])
```


Now specify model and fit it:
```@example 1

@planet b Visual{KepOrbit} begin
    e ~ Uniform(0,0.999999)
    a ~ truncated(Normal(1, 1),lower=0)
    mass ~ truncated(Normal(1, 1), lower=0)
    i ~ Sine()
    Ω ~ UniformCircular()
    ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,58849.0)  # reference epoch for θ. Choose an MJD date near your data.
end astrom



@system test begin
    M ~ truncated(Normal(1, 0.04),lower=0) # (Baines & Armstrong 2011).
    plx = 100.0
    jitter ~ truncated(Normal(0,10),lower=0)
    rv0 ~ Normal(0,10)
end rvlike b

model = Octofitter.LogDensityModel(test; autodiff=:ForwardDiff,verbosity=4)

using Random
rng = Xoshiro(0) # seed the random number generator for reproducible results

results = octofit(rng, model)
```


Display results as a corner plot:
```@example 1
octocorner(model, results, small=true)
```

Plot RV curve, phase folded curve, and binned residuals:
```@example 1
OctofitterRadialVelocity.rvpostplot(model, results,)
```

Display RV, PMA, astrometry, relative separation, position angle, and 3D projected views:
```@example 1
octoplot(model, results)
```

