# Posterior Predictive Checks

A posterior predictive check compares our true data with simulated data drawn from the posterior. This allows us to evaluate if the model is able to reproduce our observations appropriately. Samples drawn from the posterior predictive distribution should match the locations of the original data.

To demonstrate, we will fit a model to relative astrometry data:

```@example 1
using Octofitter
using Distributions

astrom_dat = Table(
    epoch= [50000,50120,50240,50360,50480,50600,50720,50840,],
    ra = [-505.764,-502.57,-498.209,-492.678,-485.977,-478.11,-469.08,-458.896,],
    dec = [-66.9298,-37.4722,-7.92755,21.6356, 51.1472, 80.5359, 109.729, 138.651,],
    σ_ra = fill(10.0, 8),
    σ_dec = fill(10.0, 8),
    cor = fill(0.0, 8)
)
astrom_like = PlanetRelAstromLikelihood(astrom_dat, name="simulated astrom")

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom_like],
    variables=@variables begin
        a ~ truncated(Normal(10, 4), lower=0.1, upper=100)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ,50420; super.M, a, e, i, ω, Ω)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

sys = System(
    name="Tutoria",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)

using Random
Random.seed!(0)
chain = octofit(model)
```

We now have our posterior as approximated by the MCMC chain. Convert these posterior samples into orbit objects:
```@example 1
# Instead of creating orbit objects for all rows in the chain, just pick
# every twentieth row.
ii = 1:20:1000
orbits = Octofitter.construct_elements(model, chain, :b, ii)
```

Calculate and plot the location the planet would be at each observation epoch:
```@example 1
using CairoMakie 

epochs = astrom_like.table.epoch' # transpose

x = raoff.(orbits, epochs)[:]
y = decoff.(orbits, epochs)[:]

fig = Figure()
ax = Axis(
    fig[1,1], xlabel="ra offset [mas]", ylabel="dec offset [mas]",
    xreversed=true,
    aspect=1
)
for orbit in orbits
    Makie.lines!(ax, orbit, color=:lightgrey)
end

Makie.scatter!(
    ax,
    x, y,
    markersize=3,
)
fig

Makie.scatter!(ax, astrom_like.table.ra, astrom_like.table.dec,color=:black, label="observed")
Makie.scatter!(ax, [0],[0], marker='⋆', color=:black, markersize=20,label="")
Makie.xlims!(400,-700)
Makie.ylims!(-200,200)
fig

```

Looks like a great match to the data! Notice how the uncertainty around the middle point is lower than the ends. That's because the orbit's posterior location at that epoch is also constrained by the surrounding data points. We can know the location of the planet in hindsight better than we could measure it!

You can follow this same procedure for any kind of data modelled with Octofitter.