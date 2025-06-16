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
data = Table([
    (; filename="Sim_data_2023_1_.oifits", epoch=mjd("2023-06-01"), use_vis2=false),
    (; filename="Sim_data_2023_2_.oifits", epoch=mjd("2023-08-15"), use_vis2=false),
    (; filename="Sim_data_2024_1_.oifits", epoch=mjd("2024-06-01"), use_vis2=false),
])
vis_like = InterferometryLikelihood(
    data,
    name="NIRISS-AMI",
    variables=@variables begin
        # For single planet:
        flux ~ truncated(Normal(0, 0.1), lower=0)  # Planet flux/contrast (array with one element)
        
        # For multiple planets (array - one per planet):
        # flux ~ Product([truncated(Normal(0, 0.1), lower=0), truncated(Normal(0, 0.1), lower=0)])
        
        # Optional calibration parameters:
        platescale = 1.0               # Platescale multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        northangle = 0.0               # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
        σ_cp_jitter = 0.0  # closure phase jitter [could use ~ LogUniform(0.1, 100))
    end
)
```

!!! note
    If you want to include multiple bands, group these into different `InterferometryLikelihood` objects
    with different instrument names (i.e. include the band in the name for the sake of bookkeeping)


Plot the closure phases:
```@example 1
fig = Makie.Figure()
ax = Axis(
    fig[1,1],
    xlabel="index",
    ylabel="closure phase",
)
Makie.stem!(
    vis_like.table.cps_data[1][:],
    label="epoch 1",
)
Makie.stem!(
    vis_like.table.cps_data[2][:],
    label="epoch 2"
)
Makie.stem!(
    vis_like.table.cps_data[3][:],
    label="epoch 3"
)
Makie.Legend(fig[1,2], ax)
fig
```

```@example 1
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[],
    variables=@variables begin
        M = super.M
        a ~ truncated(Normal(2,0.1), lower=0.1)
        e ~ truncated(Normal(0, 0.05),lower=0, upper=0.90)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()

        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 60171; M, e, a, i, ω, Ω)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

sys = System(
    name="Tutoria",
    companions=[planet_b],
    likelihoods=[vis_like],
    variables=@variables begin
        M ~ truncated(Normal(1.5, 0.01), lower=0.1)
        plx ~ truncated(Normal(100., 0.1), lower=0.1)
    end
)
```

Create the model object and run `octofit_pigeons`:
```@example 1
model = Octofitter.LogDensityModel(sys)

using Pigeons
results,pt = octofit_pigeons(model, n_rounds=10);
nothing # hide
```

Note that we use Pigeons paralell tempered sampling (`octofit_pigeons`) instead of HMC (`octofit`) because interferometry data is almost always multi-modal (or more precisely non-convex, there is often still a single mode that dominates).


Examine the recovered photometry posterior:
```@example 1
hist(results[:NIRISS_AMI_flux][:], axis=(;xlabel="flux"))
```

Determine the significance of the detection:
```@example 1
using Statistics
phot = results[:NIRISS_AMI_flux][:]
snr = mean(phot)/std(phot)
```

Plot the resulting orbit:
```@example 1
octoplot(model, results)
```


Plot only the position at each epoch:
```@example 1
using PlanetOrbits
els = Octofitter.construct_elements(model, results,:b,:);
fig = Makie.Figure()
ax = Makie.Axis(
    fig[1,1],
    autolimitaspect = 1,
    xreversed=true,
    xlabel="ΔR.A. (mas)",
    ylabel="ΔDec. (mas)",
)
for epoch in vis_like.table.epoch
    Makie.scatter!(
        ax,
        raoff.(els, epoch)[:],
        decoff.(els, epoch)[:],
        label=string(mjd2date(epoch)),
        markersize=1.5,
    )
end
Makie.Legend(fig[1,2], ax, "date")
fig
```


Finally we can examine the joint photometry and orbit posterior as a corner plot:
```@example 1
using PairPlots
using CairoMakie: Makie
octocorner(model, results)
```