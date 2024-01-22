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
    a ~ truncated(Normal(2,0.1), lower=0)
    e ~ truncated(Normal(0, 0.05),lower=0, upper=1.0)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()

    # Our prior on the planet's photometry
    # 0 +- 10% of stars brightness (assuming this is unit of data files)
    F480M ~ truncated(Normal(0, 0.1),lower=0)

    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
end

@system Tutoria begin
    M ~ truncated(Normal(1.5, 0.01), lower=0)
    plx ~ truncated(Normal(100., 0.1), lower=0)
end vis_like b
```

Create the model object and run `octofit`:
```@example 1
model = Octofitter.LogDensityModel(Tutoria)

using Random
rng = Xoshiro(0)
results = octofit(rng, model)
```


Examine the recovered photometry posterior:
```@example 1
hist(results[:b_F480M][:], axis=(;xlabel="F480M"))
```

Determine the significance of the detection:
```@example 1
using Statistics
phot = results[:b_F480M][:]
snr = mean(phot)/std(phot)
```

Plot the resulting orbit:
```@example 1
using Plots: Plots
plotchains(results,:b, kind=:astrometry, color=:b_F480M)
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
        markersize=1.5,
    )
end
Legend(fig[1,2], ax, "date")
fig
```


We can use PairPlots.jl to create a contour plot of positions at all three epochs:
```@example 1
els = Octofitter.construct_elements(results,:b,:);
fig = pairplot(
    [
        (;
                ra=raoff.(els, epoch)[:],
                dec=decoff.(els, epoch)[:],
        )=>(PairPlots.Contourf(),)
        for epoch in vis_like.table.epoch
    ]...,
    bodyaxis=(;width=400,height=400),
    axis=(;
        ra=(;reversed=true, lims=(;low=250,high=-250,)),
        dec=(;lims=(;low=-250,high=250,)),
    ),
    labels=Dict(:ra=>"ra offset [mas]", :dec=>"dec offset [mas]"),
)
scatter!(fig.content[1], [0],[0],marker='⭐', markersize=30, color=:black)
fig
```

Finally we can examine the joint photometry and orbit posterior as a corner plot:
```@example 1
octocorner(model, results)
```