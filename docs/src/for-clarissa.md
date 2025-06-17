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


epochs = 58849 .+ (20:20:660)
planet_sim_mass = 0.001 # solar masses here


orb_template = orbit(
    a = 1.0,
    e = 0.7,
    # i= pi/4, # You can remove I think
    # Ω = 0.1, # You can remove I think
    ω = 1π/4, # radians
    M = 1.0, # Total mass, not stellar mass FYI
    plx=100.0,
    tp =58829 # Epoch of periastron passage. 
)
# Makie.lines(orb_template)


rvlike = StarAbsoluteRVLikelihood(
    Table(
        epoch=epochs,
        rv=radvel.(orb_template, epochs, planet_sim_mass),
        σ_rv=fill(5.0, size(epochs)),
    ),
    name=["simulated"]
)


fg,ax,pl = Makie.scatter(rvlike.table.epoch[:], rvlike.table.rv[:])
Makie.errorbars!(rvlike.table.epoch[:], rvlike.table.rv[:], rvlike.table.σ_rv[:])
fg
```


Now specify model and fit it:
```@example 1
first_epoch_for_tp_ref = first(epochs)
@planet b RadialVelocityOrbit begin
    e ~ Uniform(0,0.999999)
    a ~ truncated(Normal(1, 1),lower=0)
    mass ~ truncated(Normal(1, 1), lower=0)

    # Remove these, we don't need 'em
    # i ~ Sine()
    # Ω ~ UniformCircular()
    # ω ~ UniformCircular()
    # θ ~ UniformCircular()

    ω ~ Uniform(0,2pi)

    # τ ~ UniformCircular(1.0)
    τ ~ Uniform(0.0, 1.0)
    tp =  b.τ*√(b.a^3/system.M)*365.25 + $first_epoch_for_tp_ref # reference epoch for τ. Choose to be near data
end 

@system SimualtedSystem begin
    M ~ truncated(Normal(1, 0.04),lower=0) # (Baines & Armstrong 2011).
    plx = 100.0
    jitter ~ truncated(Normal(0,10),lower=0)
    rv0 ~ Normal(0, 100)
end rvlike b

model = Octofitter.LogDensityModel(SimualtedSystem)

using Random
rng = Xoshiro(0) # seed the random number generator for reproducible results

results = octofit(rng, Octofitter.prior_o model, iterations=5000)
```


Display results as a corner plot:
```@example 1
octocorner(model, results, small=true)
```

Plot RV curve, phase folded curve, and binned residuals:
```@example 1
OctofitterRadialVelocity.rvpostplot(model, results,1)
```

Display RV, PMA, astrometry, relative separation, position angle, and 3D projected views:
```@example 1
octoplot(model, results)
```

