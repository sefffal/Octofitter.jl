# RV Posterior Plots

[`rvpostplot`](@ref) is the radial-velocity summary figure for a **single** posterior
draw: an RV time series carrying every instrument at once, with a residual strip and a
marginal histogram, then one phase-folded panel per planet.

```julia
res = rvpostplot(model, chain)          # the highest-posterior-density draw
res = rvpostplot(model, chain, 42)      # a draw you pick
```

The single draw is what makes the shared panel possible. A calibrated RV series is only
defined *per draw* — each instrument's offset, its jitter, and the other planets' signals
that are subtracted before folding all move from sample to sample — so several draws
cannot share one panel without misrepresenting at least one of them.

For the many-draws view, ask [`octoplot`](@ref) for the RV channels directly:

```julia
octoplot(model, chain; show_sky=false, channels=PlanetOrbits.radvel)
```

It groups panels by *observable*, not by instrument, so several instruments measuring
`radvel(A, Barycentre)` still share one panel, coloured per instrument. The difference is
what the panel is consistent with: the curves span every draw, while the data on them can
only be calibrated once, so they are calibrated with the **MAP** draw's offsets and trend.
For a fit where the offsets are well determined that is a distinction without much
difference; where they are not, `rvpostplot`'s single draw is the honest picture.

`rvpostplot` is built out of the same panel functions as the rest of the figure set, so
the error-bar conventions and the returned axes match everywhere else.

While [`octoplot`](@ref) gives a broad overview of *all* the data in a fit, this
page covers the parts of that figure that concern radial velocity: which panels
you get, what the error bars mean, and how to build a custom RV figure out of
the same pieces.

The panel anatomy — a time series with a residual strip, then one phase-folded panel per
planet:

![](assets/rv-postplot-1.png)

## Basic usage

```julia
using Octofitter, CairoMakie

# After running your fit...
res = octoplot(model, chain)
res.figure     # the Makie Figure
res.axes       # named axes, for annotating with ordinary Makie calls
res.series     # the underlying PosteriorSeries
```

`octoplot` returns an [`OctoPlotResult`](@ref), which displays as its figure. The
`axes` field is a nested NamedTuple keyed by panel: for a model with one
`RadialVelocityObs` and one planet `b` you get `res.axes.rv.main`,
`res.axes.rv.resid`, and `res.axes.rv_phase_b.main`. Pass `fname="rv.png"` to
write the figure out.

## Which panels appear, and why

Panels are derived from the model, not requested by flag. Each observation type
declares its plot channels (`Octofitter.plotchannels`) and how to compute
residuals (`Octofitter.residuals`), and `octoplot` assembles the stack:

* **Sky panel** — drawn when the system has angular observables (i.e. it defines
  `plx`). A pure RV model has none, so an RV-only fit gets no sky panel. Override
  with `show_sky=true/false`.
* **One time-series panel per channel group.** Channels are grouped by the
  *observable query* behind them, so every instrument measuring
  `radvel(star, barycentre)` shares one panel, with a different colour per
  instrument, while a relative-RV observation (`radvel(b, A)`) gets its own
  panel. This is what makes reflex and relative RV coexist in one figure.
* **One phase-folded panel per (RV channel group × foldable hierarchy row)**,
  under the default Keplerian propagator. In a two-planet fit you get a panel
  folded on `b` and a panel folded on `c`; each subtracts the other rows' signals
  from the data before folding. Override with `show_phase=true/false`.
* Any bespoke panels an observation declares through
  `Octofitter.defaultpanels`.

All epoch axes are linked, so panning one time-series panel pans them all.

## Reading the panels

### Time series panel
* RV measurements from each instrument, in different colours.
* One model curve per posterior draw, including any Gaussian process stellar
  activity model.
* Error bars: the coloured bar is the raw measurement uncertainty; the grey
  extension is measurement plus that instrument's jitter.
* The data are *calibrated* before plotting — the fitted `offset` and the
  `trend_function` contribution are subtracted from both the data and the model,
  so the points lie on the curve that `radvel` predicts.

### Residual strip
Directly under each time-series panel, with a marginal histogram of the
residuals beside it. When more than one draw is plotted, residuals are whitened
(divided by their effective uncertainty) by default, so a heteroscedastic
dataset still reads as a single scatter; pass `whiten=false` for residuals in
m/s.

### Phase-folded panels
For each planet: the data folded on that planet's orbital period with the other
planets' signals removed, noise-weighted binned means, and the isolated
single-planet model curve per draw, plus its own residual strip.

## Options

```julia
octoplot(model, chain;
    N = 250,            # posterior draws to take from the chain (default 250)
    seed = 0,           # RNG seed for choosing those draws
    ndraws = 50,        # how many of them to actually draw as curves
    show_sky = nothing, # nothing = automatic; true/false to force
    show_phase = nothing,
    whiten = nothing,   # nothing = automatic; false for un-whitened residuals
    figscale = 1.0,     # scale the whole figure
    fname = nothing,    # write the figure to this path
)
```

### Choosing which draw

`rvpostplot(model, chain, i)` renders sample `i`; with no index it picks the
highest-posterior-density sample. Slicing the chain works too, and works for every
plotting function rather than just this one:

```julia
rvpostplot(model, chain, 42)            # an arbitrary draw
rvpostplot(model, chain[42:42, :, :])   # the same thing, by slicing
```

### Animations

[`rvpostplot_animated`](@ref) sweeps through the posterior: it rebuilds the figure for `N`
single-draw slices and records them, so the phase folds move with each draw's own period.

```julia
Octofitter.rvpostplot_animated(model, chain; N=50, framerate=4, fname="rv-posterior.mp4")
```

Any extension Makie can write from a `VideoStream` works (`.mp4`, `.gif`, `.mkv`). For a
different figure per frame, loop over single-draw slices yourself:

```julia
draws = rand(1:size(chain, 1), 50)
for (k, i) in enumerate(draws)
    octoplot(model, chain[i:i, :, :]; fname="frame-$(lpad(k, 3, '0')).png")
end
```

## Building a custom RV figure

The panels are public functions, so you can assemble your own layout from a
[`PosteriorSeries`](@ref):

```julia
using Octofitter, CairoMakie

series = Octofitter.PosteriorSeries(model, chain; N=200)

# `entries` is the (observation, channel) list for the panel you want.
entries = [(rvlike, ch) for ch in Octofitter.plotchannels(rvlike)]

fig = Figure(size=(700, 700))
timeseriespanel!(fig[1, 1], series, entries)         # data vs model + residuals
phasefoldpanel!(fig[2, 1], series, entries; row=1)   # folded on hierarchy row 1
fig
```

Listing several observations' channels in one `entries` vector is what puts
several instruments in one panel. See [`timeseriespanel!`](@ref),
[`phasefoldpanel!`](@ref) and [`skypanel!`](@ref).

`octoplot`'s `channels=` keyword does the same restriction declaratively — it is exactly
what `rvpostplot` is — and takes an observable function, a channel or observable name, or
a collection of either:

```julia
octoplot(model, chain; channels=radvel)              # RV panels only
octoplot(model, chain; channels=(:sep, :pa))         # relative astrometry only
```

For the underlying data rather than a figure, `Octofitter.residuals(obs, ctx)`
returns, per channel, the epochs, the calibrated data, the model, the residual
and the effective uncertainty; [`obscontext`](@ref) builds the `ctx` from a
`PosteriorSeries`. Both `RadialVelocityObs` and `MarginalizedRVObs` implement it.
