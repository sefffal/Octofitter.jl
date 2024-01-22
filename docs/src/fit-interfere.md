# Fitting Interferometric Observables

In this tutorial, we fit a planet & orbit model to a sequence of interferometric observations.
Closure phases and squared visibilities are supported.

We load the observations in [OI-FITS format](https://github.com/emmt/OIFITS.jl) and model them as a point source orbiting a star.


!!! note
    Interferometer modelling is supported in Octofitter via the extension package OctofitterInterferometry. To install it, run 
    `pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterInterferometry`


```@example 1
using Octofitter
using OctofitterInterferometry
using Distributions
using CairoMakie
using PairPlots
```

Download simulated JWST AMI observations from our examples folder on GitHub:
```@example 1
download("https://github.com/sefffal/Octofitter.jl/raw/main/examples/AMI_data/Sim_data_2023_1_.oifits", "Sim_data_2023_1_.oifits")
download("https://github.com/sefffal/Octofitter.jl/raw/main/examples/AMI_data/Sim_data_2023_2_.oifits", "Sim_data_2023_2_.oifits")
download("https://github.com/sefffal/Octofitter.jl/raw/main/examples/AMI_data/Sim_data_2024_1_.oifits", "Sim_data_2024_1_.oifits")
```

Create the likelihood object:
```@example 1
vis_like = InterferometryLikelihood(
    (; filename="Sim_data_2023_1_.oifits", epoch=mjd("2023-06-01"), band=:F480M, use_vis2=false),
    (; filename="Sim_data_2023_2_.oifits", epoch=mjd("2023-08-15"), band=:F480M, use_vis2=false),
    (; filename="Sim_data_2024_1_.oifits", epoch=mjd("2024-06-01"), band=:F480M, use_vis2=false),
)
```

```@example 1
@planet b Visual{KepOrbit} begin
    a ~ Uniform(0, 100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()

    # Our prior on the planet's photometry
    F480M ~ Normal(0., 10)

    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
end

@system Tutoria begin
    M ~ truncated(Normal(1.0, 1), lower=0)
    plx ~ truncated(Normal(100., 50), lower=0)
end vis_like b
```

Create the model object and run `octofit`:
```@example 1
model = Octofitter.LogDensityModel(Tutoria)
results = octofit(model)
```

Plot the resulting orbit:
```@example 1
using Plots: Plots
plotchains(results,:b, kind=:astrometry)
```


Plot position at each epoch:
```@example 1
using PlanetOrbits
els = Octofitter.construct_elements(results,:b,:);
fig = Figure()
ax = Axis(
    fig[1,1],
    aspect=1,
    xreversed=true,
    xlabel="ΔR.A. (mas)",
    ylabel="ΔDec. (mas)",
)
for epoch in vis_like.table.epoch
    scatter!(
        ax,
        raoff.(els, epoch)[:],
        decoff.(els, epoch)[:],
        label=string(mjd2date(epoch)),
        markersize=3,
    )
end
Legend(fig[1,2], ax, "date")
fig
```


Examine the recovered photometry posterior:
```@example 1
hist(results[:b_F480M][:], axis=(;xlabel="F480M"))
```

And create a corner plot examining contrast vs. separation at the epoch of the first observation:
```@example 1
pairplot(
    (;
        sep=projectedseparation.(els, vis_like.table.epoch[1]),
        F480M=results["b_F480M"][:]
    ),
    labels=Dict(:sep=>"sep [mas]")
)
```