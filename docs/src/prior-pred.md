# Prior Predictive Checks

The prior predictive distributin of a Bayesian model what you get by sampling parameters directly from the priors and calculating where the model would place the data.
For example, if sampling from relative astrometry, the prior predictive model is the distribution of (simulated) astrometry points corresponding to orbits drawn from the prior. For radial velocity data, these would be simulated RV points based on an RV curve drawn from the priors.

To generate a prior predictive distribution, one first needs to create a model. We will use the model and sample data from the [Fit Astrometry](@ref fit-astrometry) tutorial:


```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions

astrom_like = PlanetRelAstromLikelihood(Table(;
    epoch= [50000,50120,50240,50360,50480,50600,50720,50840,],
    ra = [-505.764,-502.57,-498.209,-492.678,-485.977,-478.11,-469.08,-458.896,],
    dec = [-66.9298,-37.4722,-7.92755,21.6356, 51.1472, 80.5359, 109.729, 138.651,],
    σ_ra = fill(10.0, 8),
    σ_dec = fill(10.0, 8),
    cor = fill(0.0, 8)
))
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom_like],
    variables=@variables begin
        M = super.M
        a ~ truncated(Normal(10, 4), lower=0, upper=100)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω_x ~ Normal()
        ω_y ~ Normal()
        ω = atan(ω_y, ω_x)
        Ω_x ~ Normal()
        Ω_y ~ Normal()
        Ω = atan(Ω_y, Ω_x)
        θ_x ~ Normal()
        θ_y ~ Normal()
        θ = atan(θ_y, θ_x)
        tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

Tutoria = System(
    name="Tutoria",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)
nothing #hide
```

We can now draw one sample from the prior:
```@example 1
prior_draw_system = generate_from_params(Tutoria)
prior_draw_astrometry = prior_draw_system.planets.b.observations[4]
```

And plot the generated astrometry:
```@example 1
Makie.scatter(prior_draw_astrometry.table.ra, prior_draw_astrometry.table.dec,color=:black, axis=(;autolimitaspect=1,xreversed=true))

```

We can repeat this many times to get a feel for our chosen priors in the domain of our data:
```@example 1
using Random
Random.seed!(1)


fig = Figure()
ax = Axis(
    fig[1,1], xlabel="ra offset [mas]", ylabel="dec offset [mas]",
    xreversed=true,
    aspect=1
)
for i in 1:50
    prior_draw_system = generate_from_params(Tutoria)
    prior_draw_astrometry = prior_draw_system.planets.b.observations[4]
    Makie.scatter!(
        ax,
        prior_draw_astrometry.table.ra,
        prior_draw_astrometry.table.dec,
        color=Makie.cgrad(:turbo)[i/50],
    )
end

fig
```

The heavy black crosses are our actual data, while the colored ones are simulations drawn from our priors. Notice that our real data lies at a greater separation than most draws from the prior? That might mean the priors could be tweaked.
