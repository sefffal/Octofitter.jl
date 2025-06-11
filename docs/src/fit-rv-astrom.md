# [Fit Radial Velocity and Astrometry](@id fit-rv-astrom)

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
    tp =58829-40
)
Makie.lines(orb_template,axis=(;autolimitaspect=1))
```


Sample position and store as relative astrometry measurements:
```@example 1
epochs = [58849,58852,58858,58890]
astrom_dat = Table(
    epoch=epochs,
    ra=raoff.(orb_template, epochs),
    dec=decoff.(orb_template, epochs),
    σ_ra=fill(1.0, size(epochs)),
    σ_dec=fill(1.0, size(epochs)),
    cor=fill(0.0, size(epochs))
)

astrom = PlanetRelAstromLikelihood(
    astrom_dat,
    name = "simulated",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)
```

And plot our simulated astrometry measurments:
```@example 1
fig = Makie.lines(orb_template,axis=(;autolimitaspect=1))
Makie.scatter!(astrom.table.ra, astrom.table.dec)
fig
```


Generate a simulated RV curve from the same orbit:
```@example 1
using Random
Random.seed!(1)

epochs = 58849 .+ range(0,step=1.5, length=20)
planet_sim_mass = 0.001 # solar masses here


rvlike = MarginalizedStarAbsoluteRVLikelihood(
    Table(
        epoch=epochs,
        rv=radvel.(orb_template, epochs, planet_sim_mass) .+ 150,
        σ_rv=fill(5.0, size(epochs)),
    ),
    name="inst1",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)

epochs = 58949 .+ range(0,step=1.5, length=20)

rvlike2 = MarginalizedStarAbsoluteRVLikelihood(
    Table(
        epoch=epochs,
        rv=radvel.(orb_template, epochs, planet_sim_mass) .- 150,
        σ_rv=fill(5.0, size(epochs)),
    ),
    name="inst2",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)

fap = Makie.scatter(rvlike.table.epoch[:], rvlike.table.rv[:])
Makie.scatter!(rvlike2.table.epoch[:], rvlike2.table.rv[:])
fap
```


Now specify model and fit it:
```@example 1

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom],
    variables=@variables begin
        e ~ Uniform(0,0.999999)
        a ~ truncated(Normal(1, 1),lower=0.1)
        mass ~ truncated(Normal(1, 1), lower=0.)
        i ~ Sine()
        M = super.M
        Ω_x ~ Normal()
        Ω_y ~ Normal()
        Ω = atan(Ω_y, Ω_x)
        ω_x ~ Normal()
        ω_y ~ Normal()
        ω = atan(ω_y, ω_x)
        θ_x ~ Normal()
        θ_y ~ Normal()
        θ = atan(θ_y, θ_x)
        tp = θ_at_epoch_to_tperi(θ, 58849.0; M, e, a, i, ω, Ω)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

sys = System(
    name="test",
    companions=[planet_b],
    likelihoods=[rvlike, rvlike2],
    variables=@variables begin
        M ~ truncated(Normal(1, 0.04),lower=0.1) # (Baines & Armstrong 2011).
        plx = 100.0
    end
)

model = Octofitter.LogDensityModel(sys)

using Random
rng = Xoshiro(0) # seed the random number generator for reproducible results

results = octofit(rng, model, max_depth=9, adaptation=300, iterations=400)
```

Display results as a corner plot:
```@example 1
octocorner(model,results, small=true)
```

Plot RV curve, phase folded curve, and binned residuals:
```@example 1
Octofitter.rvpostplot(model, results)
```

Display RV, PMA, astrometry, relative separation, position angle, and 3D projected views:
```@example 1
octoplot(model, results)
```

