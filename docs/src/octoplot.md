# Orbit Visualization with `octoplot`

`octoplot` is a versatile visualization function that creates publication-quality figures showing orbit fits and data. It assembles a single figure containing:

* Projected orbits in the plane of the sky (or in AU, if the model has no parallax)
* One time-series panel per data channel, each with a residual strip and a marginal residual histogram
* Phase-folded panels for radial velocity channels
* Any bespoke panels an observation type provides

Here's a basic example showing how to create a plot from your MCMC chain:
```julia
using Octofitter
using CairoMakie      # a Makie backend must be loaded
# After running your MCMC fit...
res = octoplot(model, chain)
```

`octoplot` returns an [`OctoPlotResult`](@ref), which displays as its figure and also
carries the pieces you need to customize it:

```julia
res.figure     # the Makie Figure
res.axes       # a nested NamedTuple of the axes, e.g. res.axes.sky.sky, res.axes.rv.main
res.series     # the underlying PosteriorSeries — the solved orbits behind every panel
```

## What gets drawn

`octoplot` derives its panels from the model, not from a list of flags. Every observation
declares its *plot channels* (`plotchannels(obs)`), each of which is one measured
quantity with an associated [`ObservableQuery`](@ref) — for example
`radvel(:A, Barycentre)` for a stellar reflex RV instrument, or `raoff(:b, :A)` and
`decoff(:b, :A)` for relative astrometry. Channels measuring the same quantity share a
panel, so two RV instruments overlay on one axis, with each instrument's own offset and
jitter applied.

The sky panel is drawn whenever the model has angular observables (i.e. the system block
defines `plx`) and at least one orbit. It shows one phase-coloured track per **hierarchy
row** — each row plotted as its exterior side relative to its interior side, which is
exactly the relationship that row parametrizes. For a star with planets that is each
planet about the star; for a hierarchical system it generalizes with no special cases (a
moon about its planet, an inner pair's barycentre about the outer body).

The two panel families can be turned off explicitly:

```julia
octoplot(model, chain;
    show_sky = false,     # skip the sky-plane panel
    show_phase = false,   # skip phase-folded RV panels
)
```

!!! note "Not every observation type has plot channels yet"
    `RelAstromObs`, `RadialVelocityObs`, `MarginalizedRVObs` and `GaiaDR4AstromObs`
    declare plot channels and appear in `octoplot`. `ImageObs`, `InterferometryObs`,
    `PhotometryObs`, `G23HObs` and `HipparcosIADObs` do not yet: `octoplot` still draws
    the orbits, but those observations' data are not overlaid. This is a gap in the
    plotting layer rather than a deliberate omission.

## Controlling the draws

By default, `octoplot` draws 250 orbits selected without replacement from your posterior,
with a fixed seed so that reruns and separate panels agree:
```julia
# Plot fewer orbits for faster rendering
octoplot(model, chain, N=50)

# A different (but still reproducible) subset
octoplot(model, chain, N=250, seed=42)
```

`ndraws` limits how many of the selected draws each *panel* actually renders, which is
useful when you want a dense set of orbits computed but a light plot:
```julia
octoplot(model, chain, N=250, ndraws=25)
```

To choose the draws yourself — for example to plot only the maximum-posterior sample —
build the [`PosteriorSeries`](@ref) explicitly and pass it in:
```julia
idx_MAP = argmax(chain[:logpost][:])
series = Octofitter.PosteriorSeries(model, chain; ii=[idx_MAP])
octoplot(series)
```

## Other keywords

| Keyword | Meaning |
|---|---|
| `N=250` | Number of posterior draws to solve |
| `seed=0` | Seed for the draw selection |
| `ndraws=nothing` | Cap on the draws each panel renders |
| `show_sky=nothing` | Force the sky panel on or off (default: on when the model is angular) |
| `show_phase=nothing` | Force phase-folded panels on or off (default: on under a Keplerian propagator) |
| `whiten=nothing` | Divide residuals by their uncertainty (default: on when several draws are shown) |
| `figscale=1.0` | Scale the whole figure |
| `fname=nothing` | If set, save the figure to this path |

Nothing is written to disk unless you pass `fname`.

## Building your own figure

The panels `octoplot` assembles are public, so you can lay out your own figure. They all
take a Makie grid position and a `PosteriorSeries`:

```julia
using CairoMakie
series = Octofitter.PosteriorSeries(model, chain; N=100)

fig = Figure(size=(700, 900))
skypanel!(fig[1,1], series)
timeseriespanel!(fig[2,1], series, entries)   # entries: (observation, channel) pairs
phasefoldpanel!(fig[3,1], series, entries; row=1)
fig
```

To get the model curves as plain numbers rather than a plot, evaluate an
[`ObservableQuery`](@ref) over the series:

```julia
q = ObservableQuery(radvel, :A, Barycentre)
curves = modelcurves(series, q)   # one vector per draw, over series.ts
best   = mapcurve(series, q)      # the maximum-posterior draw
```

## Predicting positions at future epochs

Solve the system for a draw and query it at whatever epochs you like:

```julia
posys = construct_system(model, chain, 1)
sols  = orbitsolve(posys, [mjd("2028-01-01"), mjd("2029-01-01")])

raoff(sols[1], :b, :A)                 # ΔRA  [mas]
decoff(sols[1], :b, :A)                # ΔDec [mas]
projectedseparation(sols[1], :b, :A)   # [mas]
posangle(sols[1], :b, :A)              # [rad]
```

Loop over draws to build a predicted distribution, and plot it however you like.

# Customizing Appearance

## Post-Creation Customization
Since `octoplot` returns a result carrying named axes, you can customize any panel after
creation using ordinary Makie calls:
```julia
res = octoplot(model, chain)

# Named access, rather than guessing an index
ax = res.axes.sky.sky
xlims!(ax, -100, 100)
ylims!(ax, -100, 100)
ax.title = "HD 12345 Orbital Fit"

res.figure
```

The keys of `res.axes` name the panels: `sky` for the sky panel, and one entry per
channel group (e.g. `raoff`, `decoff`, `rv`, `rv_phase_b`), each of which is itself a NamedTuple of
`main`, `resid` and `hist` axes. `println(keys(res.axes))` lists them for your model.

!!! tip
    Take care when modifying time-based panels as they share synchronized x-axes.
    Modifying the time limits of one panel will affect all others.


# Saving Figures

When using CairoMakie (recommended for publication-quality outputs), you can save figures in several formats:

```julia
res = octoplot(model, chain)

# Or pass fname= to octoplot and skip this step
save("orbit_plot.png", res.figure)
save("orbit_plot.pdf", res.figure)   # great for publications
save("orbit_plot.svg", res.figure)   # good for further editing

# Increase PNG resolution
save("orbit_plot.png", res.figure, px_per_unit=5)
```

!!! note
    Many plot elements are internally rasterized for performance, so extremely high
    `px_per_unit` values may not improve quality.

If using GLMakie instead of CairoMakie, you get interactive figures that can be zoomed and panned before saving. However, only PNG output is supported. GLMakie is great for exploration, while CairoMakie is preferred for final publication figures.

# Corner Plots

[`octocorner`](@ref) draws the posterior pair plot, using PairPlots.jl:

```julia
using CairoMakie, PairPlots
octocorner(model, chain, small=true)
```

`small=true` keeps only each body's `a`, `e`, `i` and `mass`; `includecols` and
`excludecols` add or remove columns by name. `UniformCircular` helper pairs and fixed
values are dropped automatically. Passing several chains overlays them, which is a
convenient way to compare a prior draw with a posterior, or two models with each other.
