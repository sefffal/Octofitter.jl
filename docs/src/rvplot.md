# Radial Velocity Figures

Two figures cover radial velocity, and they answer different questions.

[`rvplot`](@ref) is the **single-draw summary**: one time series carrying every
instrument at once, with a residual strip and a marginal histogram, then one
phase-folded panel per planet with noise-weighted binned means. This is the
publication figure.

[`octoplot`](@ref) is the **posterior view**: the spread of orbits the fit
allows, with each instrument on its own panel.

```julia
res = rvplot(model, chain)          # the highest-posterior-density draw
res = rvplot(model, chain, 42)      # a draw you pick
```

![](assets/rv-postplot-1.png)

!!! note "Renamed in v9"
    This figure was called `rvpostplot` before v9. It never plotted a
    posterior — it plots one draw — so it is now `rvplot`, and
    `rvpostplot_animated` is `rvplot_animated`. The old names still work and
    forward to the new ones with a deprecation warning.

## Why one draw, and why that lets instruments share an axis

A calibrated RV series is only defined **per draw**. Each instrument's zero
point, its jitter, the long-term trend and the other planets' signals that are
subtracted before folding all move from sample to sample. Fix a draw and every
one of those is a number: the data can be shifted onto a single model curve and
several spectrographs can be compared on one axis, honestly, because everything
in the frame belongs to the same sample.

Show many draws and that stops being true — there is no one calibration to put
the data in. So `octoplot` does the opposite: it leaves each instrument's
measurements exactly as reported, gives each instrument its own panel, and lets
**each draw's model curve carry that draw's own offset and trend**. The curves
move; the data stay put. That is the only arrangement in which a posterior
cloud and a dataset are both drawn correctly.

```julia
octoplot(model, chain; show_sky=false, channels=radvel)
```

This is a policy, not a special case for RV. Each observation type declares
whether its calibration is tight enough to share a panel through
[`sharepanel`](@ref); relative astrometry says yes (a `platescale` pinned to a
fraction of a percent moves nothing visibly), radial velocity says no.

## Reading the panels

### Time series
* Every instrument, one colour and one marker shape each — the shapes matter in
  print and for readers who cannot rely on colour.
* One model curve per draw. In `rvplot` there is exactly one.
* Error bars: the black bar is the raw measurement uncertainty; the grey
  extension adds that instrument's fitted jitter (and its Gaussian-process
  predictive variance, if there is one).
* Calendar dates label the top of the figure; MJD numbers label the bottom.

### Residual strip
Directly under the time series, with a marginal histogram beside it.

With **one draw** the residuals are in m/s, and each is the measurement minus
everything the model accounts for — orbit, offset, trend, and the
Gaussian-process activity prediction where one was fitted.

With **several draws** they are whitened: per point, the median and 16–84 %
interval of `(data − model)/σ_eff` across the draws, with the histogram drawn
against a unit normal. This is not a display preference. A raw residual is not
comparable between draws — each draw has its own jitter, so the same point sits
a different number of σ from the model in each — which is why
`octoplot(...; whiten=false)` is an *error* unless you also ask for `ndraws=1`.
`σ_eff` is measurement uncertainty, jitter and GP variance combined, so the
z-scores are standard normal for a fit that is working.

### Phase-folded panels
One per planet: the data folded on that planet's period with the other planets'
signals removed exactly (through the reference grammar, not plot-side
subtraction), the isolated single-planet model curve, and noise-weighted binned
means in red. Phase zero is the signal's upward zero crossing, as before.

### Correlated noise
Where a `gaussian_process` was fitted, its prediction is drawn as a band around
the model curve, subtracted from the residuals, and folded into `σ_eff`. See
[Fit Gaussian Process](@ref fit-rv-gp).

## Options

```julia
rvplot(model, chain;
    show_phase = true,          # phase-folded panels
    show_hist = true,           # marginal residual histogram
    show_legend = true,         # instrument + symbol legend below the figure
    show_perspective = false,   # the frame's own RV drift; see below
    gpband = true,              # the correlated-noise band
    figscale = 1.0,
    fname = nothing,            # write the figure here
)
```

!!! warning "`show_perspective` is a diagnostic, not part of the fit"
    v9's `radvel` is a strictly relative observable, and `RadialVelocityObs`
    does not currently add the system barycentre's own secular radial-velocity
    drift (`PlanetOrbits.frame_rv`) to its forward model. `show_perspective=true`
    overlays that drift as a dashed orange line for a model with an absolute
    frame, so you can see how large it is — but it is not in the fitted curve,
    which is why it is off by default. v8 included the term in the model.

### Choosing a draw

```julia
rvplot(model, chain, 42)            # an arbitrary draw
rvplot(model, chain[42:42, :, :])   # the same thing, by slicing
```

Slicing works for every plotting function, not just this one.

### Animations

[`rvplot_animated`](@ref) sweeps through the posterior: it rebuilds the figure
for `N` single-draw slices and records them, so the phase folds move with each
draw's own period.

```julia
Octofitter.rvplot_animated(model, chain; N=50, framerate=4, fname="rv-posterior.mp4")
```

Any extension Makie can write from a `VideoStream` works (`.mp4`, `.gif`,
`.mkv`). For a different figure per frame, loop over single-draw slices
yourself:

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
several instruments in one panel — and `calibrate=:data` is what makes that
legitimate, which is why you should only do it for a single-draw series. See
[`timeseriespanel!`](@ref), [`phasefoldpanel!`](@ref) and [`skypanel!`](@ref).

For the underlying data rather than a figure, `Octofitter.residuals(obs, ctx)`
returns, per channel, the epochs, the calibrated data, the model, the residual
and the effective uncertainty; [`obscontext`](@ref) builds the `ctx` from a
`PosteriorSeries`. Both `RadialVelocityObs` and `MarginalizedRVObs` implement
it.
