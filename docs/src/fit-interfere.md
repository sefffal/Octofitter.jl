# Fitting Interferometric Observables

In this tutorial, we fit a planet & orbit model to a sequence of interferometric observations.
Closure phases and squared visibilities are supported.

We load the observations in OI-FITS format and model them as a collection of point sources.


!!! note
    Interferometer modelling is supported in Octofitter via the extension package OctofitterInterferometry.
    It is unregistered — install it with a `PackageSpec` pointing at the Octofitter
    repository; see [Installation](@ref install).


```@example 1
using Octofitter
using OctofitterInterferometry
using Distributions
using CairoMakie
using PairPlots
using Statistics
using Pigeons
```

Download simulated JWST AMI observations from our examples folder on GitHub:
```@example 1
download("https://github.com/sefffal/Octofitter.jl/raw/main/examples/AMI_data/Sim_data_2023_1_.oifits", "Sim_data_2023_1_.oifits")
download("https://github.com/sefffal/Octofitter.jl/raw/main/examples/AMI_data/Sim_data_2023_2_.oifits", "Sim_data_2023_2_.oifits")
download("https://github.com/sefffal/Octofitter.jl/raw/main/examples/AMI_data/Sim_data_2024_1_.oifits", "Sim_data_2024_1_.oifits")
```

## The forward model

An interferometric observation in Octofitter is a set of point sources whose complex visibility is

```math
V(u,v) = \frac{\sum_j f_j \, e^{-2\pi i (u\,\Delta\alpha^*_j + v\,\Delta\delta_j)}}{\sum_j f_j}
```

`targets` names the bodies in that sum, `ref` is the phase centre the offsets are measured from, and each `f_j` is that body's `flux_<band>` variable — **the host included**. There is no privileged primary: a source may orbit any body, so a moon or the wide component of a hierarchical system is expressible.

## Build the model

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.5, 0.01), lower=0.1)   # M⊙
        flux_K = 1.0        # the host is an ordinary source now
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0
        flux_K ~ truncated(Normal(0, 0.1), lower=0)   # contrast ratio against the host

        a ~ truncated(Normal(2,0.1), lower=0.1)
        e ~ truncated(Normal(0, 0.05),lower=0, upper=0.90)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 60171.0   # reference epoch for θ. Choose an MJD date near your data.
    end
)
nothing # hide
```

Now the observation itself:

```@example 1
data = Table([
    (; filename="Sim_data_2023_1_.oifits", epoch=mjd("2023-06-01"), use_vis2=false),
    (; filename="Sim_data_2023_2_.oifits", epoch=mjd("2023-08-15"), use_vis2=false),
    (; filename="Sim_data_2024_1_.oifits", epoch=mjd("2024-06-01"), use_vis2=false),
])
vis_obs = InterferometryObs(
    data,
    targets = (A, b),      # every source in the visibility sum, host included
    ref     = A,           # phase centre
    band    = :K,          # which `flux_<band>` to read
    name    = "NIRISS-AMI",
    variables=@variables begin
        # Optional calibration parameters:
        platescale = 1.0               # Platescale divisor [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        northangle = 0.0               # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
        σ_cp_jitter = 0.0  # closure phase jitter, degrees [could use ~ LogUniform(0.1, 100))
    end
)

sys = System(
    name="Tutoria",
    bodies=[A, b],
    observations=[vis_obs],
    variables=@variables begin
        plx ~ truncated(Normal(100., 0.1), lower=0.1)
    end
)
nothing # hide
```

A second companion is a third entry in `targets` plus its own `flux_K` on its own body — not a second element of a `Product` distribution.

!!! note
    If you want to include multiple bands, group these into different `InterferometryObs` objects
    with different instrument names (i.e. include the band in the name for the sake of bookkeeping),
    and give the bodies one `flux_<band>` variable per band.

!!! note "Choosing `ref`"
    Closure phases, kernel phases and squared visibilities are all invariant to the
    phase centre, so `ref` is a free choice — but only modulo 360°: baseline phases
    are folded into (−180°, 180°] and the triangle sum is not. Keep `ref` near the flux
    centroid; `Barycentre` (the default) is fine for a faint companion, and `A` is the
    conventional choice for a bright one.

Plot the closure phases:
```@example 1
fig = Makie.Figure()
ax = Axis(
    fig[1,1],
    xlabel="index",
    ylabel="closure phase",
)
Makie.stem!(
    vis_obs.table.cps_data[1][:],
    label="epoch 1",
)
Makie.stem!(
    vis_obs.table.cps_data[2][:],
    label="epoch 2"
)
Makie.stem!(
    vis_obs.table.cps_data[3][:],
    label="epoch 3"
)
Makie.Legend(fig[1,2], ax)
fig
```

Create the model object and sample:
```@example 1
model = Octofitter.LogDensityModel(sys)

init_chain = initialize!(model)
results, pt = octofit_pigeons(model, n_rounds=10)
nothing # hide
```

!!! tip "Why parallel tempering here"
    Interferometry data is almost always multi-modal (or more precisely non-convex —
    there is often still a single dominant mode), which is why this page samples with
    [`octofit_pigeons`](@ref) rather than HMC.

    If you do use [`octofit`](@ref) instead, treat a single chain as an exploration of
    one mode: run several chains from different starting points and compare them before
    believing a single-peaked posterior.

Examine the recovered photometry posterior. The contrast is a body variable, so its chain key is `b_flux_K`:
```@example 1
hist(results[:b_flux_K][:], axis=(;xlabel="flux (K band, relative to host)"))
```

Determine the significance of the detection:
```@example 1
phot = results[:b_flux_K][:]
snr = mean(phot)/std(phot)
```

Plot the resulting orbit:
```@example 1
octoplot(model, results)
```


Plot only the position at each epoch. `construct_system(model, chain)` rebuilds one PlanetOrbits system per posterior draw:
```@example 1
using PlanetOrbits
posteriors = construct_system(model, results)
fig = Makie.Figure()
ax = Makie.Axis(
    fig[1,1],
    autolimitaspect = 1,
    xreversed=true,
    xlabel="ΔR.A. (mas)",
    ylabel="ΔDec. (mas)",
)
for epoch in vis_obs.table.epoch
    sols = [orbitsolve(p, epoch) for p in posteriors]
    Makie.scatter!(
        ax,
        [raoff(s, :b, :A) for s in sols],
        [decoff(s, :b, :A) for s in sols],
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
