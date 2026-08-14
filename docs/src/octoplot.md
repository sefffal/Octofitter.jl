# [Orbit plots with `octoplot`](@id octoplot-page)

`octoplot(model, chain)` is the one-call summary of a fit: the orbits on the sky,
your data with the model drawn through it, and residuals underneath. This page
starts from that figure and then works through the changes people most often want
to make to it — different panels, a different time range, fewer orbits, one panel
on its own for a paper.

You do not need to know anything about Makie to follow it. Everything up to
[Advanced topics](@ref octoplot-advanced) is a short recipe you can paste next to
your own `octoplot` call.

## The example fit

Every example below uses the fit built here — one planet, relative astrometry from
one instrument, and radial velocities from two spectrographs. **If you already have
your own `model` and `chain`, skip this block**; everything on the page applies
unchanged.

```@example 1
using Octofitter, Distributions, CairoMakie, PairPlots

A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.2, 0.1), lower=0.1)   # [M⊙]
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(1mjup, 100mjup)   # [M⊙]
        a ~ Uniform(0.5, 20)                # [AU]
        e ~ Uniform(0.0, 0.6)
        i ~ Sine()                          # [rad]
        ω ~ Uniform(0, 2pi)                 # [rad]
        Ω ~ Uniform(0, 2pi)                 # [rad]
        θ ~ Uniform(0, 2pi)                 # [rad]
        epoch = 57000.0                     # [MJD]
    end
)

astrom = RelAstromObs(
    Table(
        epoch = [55927.0, 56444.0, 57023.0, 57754.0, 58484.0, 59366.0],
        ra    = [-102.40, -36.35, 137.31, -101.39, 22.67, 23.57],
        dec   = [29.74, 165.11, 12.63, 9.38, 152.21, -99.09],
        σ_ra  = fill(5.0, 6),
        σ_dec = fill(5.0, 6),
        cor   = fill(0.0, 6),
    );
    target=b, ref=A, name="GPI"
)

harps = RadialVelocityObs(
    Table(
        epoch = [55927.0, 56191.5, 56456.1, 56720.6, 56985.2, 57249.7, 57514.2,
                 57778.8, 58043.3, 58307.8, 58572.4, 58836.9, 59101.5, 59366.0],
        rv    = [411.6, 595.4, 427.9, 108.1, -304.2, -707.8, -523.4,
                 379.2, 600.0, 419.4, 134.6, -289.9, -693.1, -562.2],
        σ_rv  = fill(6.0, 14),
    );
    target=A, ref=Barycentre, name="HARPS",
    variables=@variables begin
        offset ~ Uniform(-500, 500)
        jitter ~ LogUniform(0.1, 30)
    end
)

hires = RadialVelocityObs(
    Table(
        epoch = [56809.0, 57076.3, 57343.7, 57611.0, 57878.3,
                 58145.7, 58413.0, 58680.3, 58947.7, 59215.0],
        rv    = [84.4, -339.0, -653.4, -50.1, 660.3, 684.6, 426.0, 74.1, -355.1, -641.2],
        σ_rv  = fill(9.0, 10),
    );
    target=A, ref=Barycentre, name="HIRES",
    variables=@variables begin
        offset ~ Uniform(-500, 500)
        jitter ~ LogUniform(0.1, 30)
    end
)

sys = System(
    name="HD12345",
    bodies=[A, b],
    observations=[astrom, harps, hires],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
init_chain = initialize!(model, (; bodies=(; b=(; a=3.2, e=0.25),),))
chain = octofit(model, iterations=1000)
nothing # hide
```

## The default figure

```@example 1
using CairoMakie      # a Makie backend must be loaded before plotting

octoplot(model, chain)
```

Reading it from the top:

* **The sky panel.** One orbit track per planet, drawn once for each posterior
  sample, coloured by orbital phase. The star sits at the origin, and right
  ascension increases to the **left**, as it does on the sky. Your relative
  astrometry is overlaid, labelled by instrument.
* **One panel per measured quantity, per instrument.** Here that is Δα\* and Δδ
  from the GPI astrometry, then radial velocity from HARPS and from HIRES on
  separate panels. All of these share one horizontal axis: calendar dates along
  the top of the figure, MJD numbers along the bottom.
* **A residual strip under each panel**, with a marginal histogram at its right.

A few words are used throughout the rest of the page, and they all mean something
small and concrete:

| Word | What it means |
|---|---|
| **panel** | One of the stacked plots — the sky panel, the HARPS panel, and so on. |
| **draw** | One sample from your chain. `octoplot` overplots many, which is why each panel shows a bundle of curves rather than a single line. |
| **channel** | One measured quantity from one observation: `sep`, `pa`, `raoff`, `decoff`, `rv`, … Panels are built one per channel. |
| **row** | One orbit in your model — one `Body(..., about=...)`. This model has a single row, planet `b` about star `A`. Rows only matter for phase folds and for building your own figure. |

## What you get back

`octoplot` returns an [`OctoPlotResult`](@ref). It displays as the figure, and it
also carries the pieces you need in order to change anything:

```@example 1
res = octoplot(model, chain)

res.figure     # the Makie Figure — pass this to `save`
res.axes       # the panels, by name
res.series     # the solved orbits behind every panel
nothing # hide
```

`res.axes` names every panel. Print the keys to see what your own model produced:

```@example 1
keys(res.axes)
```

`sky` is the sky panel; the rest are named after the channel they show. Where
several instruments measure the same quantity, the instrument name is appended —
`rv_HARPS`, `rv_HIRES` — because `rv` alone would no longer identify one panel.
Each entry is itself a small named group:

```@example 1
keys(res.axes.rv_HARPS)      # main plot, residual strip, marginal histogram
```

So `res.axes.rv_HARPS.main` is an ordinary Makie `Axis` that you can set a title
or limits on. The sky panel's axis is `res.axes.sky.sky`.

## Draw fewer orbits

By default `octoplot` solves 250 draws from your posterior, chosen without
replacement with a fixed seed so that reruns and separate panels agree.

```@example 1
octoplot(model, chain, N=50)     # faster, lighter figure
nothing # hide
```

```@example 1
octoplot(model, chain, N=250, seed=42)   # a different, still reproducible subset
nothing # hide
```

`N` sets how many draws are *computed*. `ndraws` caps how many each panel
actually *renders*, which is what you want if a light-looking plot should still be
backed by a well-sampled posterior:

```@example 1
octoplot(model, chain, N=250, ndraws=25)
nothing # hide
```

`ndraws=1` is a special case worth knowing: with a single draw the residuals are
shown in data units rather than whitened, phase-folded panels appear, and
correlated-noise bands are drawn. See [Show a single orbit](@ref).

## Choose which panels appear

Two keywords do almost all of this.

`show_sky=false` drops the sky panel:

```@example 1
octoplot(model, chain; show_sky=false)
nothing # hide
```

`channels=` restricts the figure to some of the measurements. It takes a channel
name as a `Symbol` (`:sep`, `:pa`, `:raoff`, `:decoff`, `:rv`), an observable
function (`radvel`, `pmra`, `projectedseparation`, …), or a collection of either:

```@example 1
octoplot(model, chain; channels=radvel, show_sky=false)
```

```@example 1
octoplot(model, chain; channels=(:raoff, :decoff))   # astrometry only, sky kept
nothing # hide
```

Relative astrometry is special in a useful way: whether your table holds `ra`/`dec`
or `sep`/`pa`, **both** parametrizations are available, because one is a
parameter-free rotation of the other. The figure above was fitted to `ra`/`dec`,
and asking for separation works anyway:

```@example 1
octoplot(model, chain; channels=:sep, show_sky=false)
```

Nothing is written to disk by any of these calls unless you pass `fname=`.

## Plot a quantity you didn't measure

Naming an **observable** that your model has no data for draws it as a
*prediction*: the model curves alone, with nothing to compare them to. This fit
contains no proper-motion data, but the orbit implies a reflex proper motion:

```@example 1
octoplot(model, chain; channels=pmra, show_sky=false)
```

This is how you ask a fit what signal it implies before spending telescope time on
it — `channels=radvel` on an astrometry-only fit is the common case. Note that this
only happens for an *observable* (`radvel`, `pmra`); naming a *channel*
(`channels=:rv`) means "the RV data", and asking for data a fit does not have is
not a request to invent it.

## Zoom the time axis

Every time panel shares one x axis, so set the limits on any of them and the whole
stack follows. Epochs are MJD numbers; [`mjd`](@ref) converts a date string:

```@example 1
res = octoplot(model, chain; N=50)
xlims!(res.axes.raoff.main, mjd("2012-01-01"), mjd("2018-01-01"))
res.figure
```

The model curves exist over the epoch grid `octoplot` chose, which you can inspect:

```@example 1
mjd2date.(extrema(res.series.ts))
```

That grid covers your data, padded slightly and widened to at least one orbital
period, so zooming out past it will show empty axes rather than more curve. There
is no keyword to extend it — to evaluate the orbit at epochs of your own choosing,
see [Positions at future epochs](@ref octoplot-future-epochs).

!!! tip
    `ylims!` is *not* linked, so you can set a vertical range on one panel without
    touching the others.

## Change titles, labels and limits

`res.axes` gives you the panels by name, so anything Makie can do to an `Axis`, you
can do here:

```@example 1
res = octoplot(model, chain; N=50)

res.axes.sky.sky.title = "HD 12345 b"
res.axes.rv_HARPS.main.ylabel = "RV [m/s]"
ylims!(res.axes.rv_HIRES.main, -900, 900)

res.figure
```

`legend=false` removes the instrument key from the sky and shared data panels, and
`figscale=` scales the whole figure at once.

## Save the figure

Either pass `fname=` and let `octoplot` write it, or save `res.figure` yourself:

```julia
octoplot(model, chain; fname="orbit.png")

res = octoplot(model, chain)
save("orbit.png", res.figure)
save("orbit.pdf", res.figure)   # vector output, good for a journal
save("orbit.svg", res.figure)   # vector output, good for further editing

save("orbit.png", res.figure, px_per_unit=4)   # higher-resolution raster
```

With CairoMakie (the usual choice for publication figures) all three formats work.
GLMakie gives you an interactive window you can zoom and pan before saving, but
only writes PNG.

!!! note
    Some plot elements are rasterized internally for performance, so very large
    `px_per_unit` values stop improving quality. If you need a truly resolution-free
    figure, save PDF or SVG.

## One panel on its own

For a paper you often want one panel by itself, at your own size and aspect ratio.
Build a [`PosteriorSeries`](@ref) — the solved orbits, computed once — and hand it
to the panel function directly.

`entries` says which observation and which channel the panel should draw.
[`plotchannels`](@ref) lists the channels an observation exposes:

```@example 1
series = Octofitter.PosteriorSeries(model, chain; N=100)

entries = [(astrom, ch) for ch in Octofitter.plotchannels(astrom) if ch.name === :raoff]

fig = Figure(size=(600, 320))
timeseriespanel!(fig[1, 1], series, entries)
fig
```

The panel brings its own residual strip and marginal histogram, as it does inside
`octoplot`; `show_hist=false` drops the histogram. The sky panel needs no entries at
all:

```@example 1
fig = Figure(size=(500, 500))
skypanel!(fig[1, 1], series)
fig
```

## Lay out your own figure

The panels are ordinary Makie recipes taking a grid position, so you can arrange
them however you like — side by side, for instance:

```@example 1
fig = Figure(size=(1100, 480))
skypanel!(fig[1, 1], series; colorbar=false)
timeseriespanel!(fig[1, 2], series, entries)
colsize!(fig.layout, 1, Aspect(1, 1.0))
fig
```

The three panel builders are [`skypanel!`](@ref), [`timeseriespanel!`](@ref) and
[`phasefoldpanel!`](@ref). They all take a grid position first and the same
`PosteriorSeries` second, so every panel in your figure shows the same draws.

## Phase-folded panels

A phase fold collapses the epoch axis onto one orbital period. `octoplot` leaves
these **off by default** — see [Why phase folds are off by default](@ref) — and
`show_phase=true` asks for them:

```@example 1
octoplot(model, chain; channels=radvel, show_sky=false, show_phase=true, N=100)
```

Each fold gets its own panel, named after the channel and the row it was folded on
(`rv_HARPS_phase_b`). `show_phase=true` folds *every* channel that can be folded,
not only radial velocity, so a well-sampled astrometric orbit can be folded too:

```julia
octoplot(model, chain; channels=posangle, show_sky=false, show_phase=true)
```

For the conventional single-draw radial-velocity figure with phase folds and binned
means, use [`rvplot`](@ref) instead — see [Radial Velocity Figures](@ref).

## Show a single orbit

Passing `ndraws=1` renders one draw. Because there is then a single set of
parameters, residuals are shown in data units, phase folds appear, and any
correlated-noise band is drawn:

```@example 1
octoplot(model, chain; ndraws=1)
nothing # hide
```

To choose *which* draw, either slice the chain or build the series yourself:

```@example 1
i_map = argmax(chain[:logpost][:])          # the maximum-posterior sample

octoplot(model, chain[i_map:i_map, :, :])   # by slicing

series_map = Octofitter.PosteriorSeries(model, chain; ii=[i_map])
octoplot(series_map)                        # by building the series
nothing # hide
```

Slicing works for every plotting function, not just this one.

## Colours and marker styles

`octoplot` picks its colours itself and has no colour keywords: each orbit row gets
an accent colour that its sky track and its panels share, and each instrument on a
shared panel gets its own colour *and* its own marker shape (so datasets stay
separable in greyscale print). `legend=false` drops the instrument key.

To choose the colours, build the panel yourself — [`timeseriespanel!`](@ref) and
[`phasefoldpanel!`](@ref) accept `curvecolor=` and a `datastyle=` NamedTuple:

```@example 1
fig = Figure(size=(600, 320))
timeseriespanel!(fig[1, 1], series, entries;
    curvecolor = :firebrick,
    datastyle  = (; markersize=10, σ_color=:instrument))
fig
```

`datastyle` accepts `markersize`, `resid_markersize`, `strokewidth`, `σ_linewidth`,
`σ_color`, `σeff_linewidth` and `σeff_color`; an unrecognised key is an error rather
than a silent no-op. A Makie theme (`set_theme!`, `with_theme`) still governs fonts,
sizes and axis decorations, but the curve and marker colours are set explicitly by
the panels and are not read from the theme.

## All keywords

| Keyword | Meaning |
|---|---|
| `N=250` | Number of posterior draws to solve |
| `seed=0` | Seed for the draw selection |
| `ndraws=nothing` | Cap on the draws each panel renders |
| `show_sky=nothing` | Force the sky panel on or off (default: on when the model has angular observables and at least one orbit) |
| `show_phase=nothing` | `nothing` folds only when a single draw is shown; `true` folds every foldable channel; `false` folds nothing |
| `channels=nothing` | Restrict the figure to some channels or observables |
| `whiten=nothing` | Divide residuals by their uncertainty (default: on when several draws are shown, and required then) |
| `boxwidth=nothing` | Width of the residual boxes in x units (default: from the epoch spacing) |
| `gpband=nothing` | Correlated-noise bands (default: on only for a single draw) |
| `legend=true` | Instrument keys on the sky and shared data panels |
| `figscale=1.0` | Scale the whole figure |
| `figure=nothing` | Draw into an existing (emptied) `Figure` rather than a new one |
| `fname=nothing` | If set, save the figure to this path |

`N` and `seed` belong to the `octoplot(model, chain; ...)` form; if you build a
[`PosteriorSeries`](@ref) yourself, they are its keywords instead.

## Related figures

* [`rvplot`](@ref) — the single-draw radial-velocity figure, with every instrument
  on one axis and phase folds by default. See [Radial Velocity Figures](@ref).
* [`octocorner`](@ref) — the posterior corner plot, below.
* [`dotplot`](@ref) — mass against separation or period for every body, coloured by
  eccentricity.
* [`gaiastarplot`](@ref), [`gaiatimeplot`](@ref), [`hipparcosplot`](@ref),
  [`skytrackplot`](@ref) — absolute-astrometry figures in their own geometry.

## Corner plots

[`octocorner`](@ref) draws the posterior pair plot, using PairPlots.jl:

```@example 1
using CairoMakie, PairPlots
octocorner(model, chain, small=true)
```

`small=true` keeps only each body's `a`, `e`, `i` and `mass`; `includecols` and
`excludecols` add or remove columns by name. `UniformCircular` helper pairs and
fixed values are dropped automatically. Passing several chains overlays them, which
is a convenient way to compare a prior draw with a posterior, or two models with
each other.

---

# [Advanced topics](@id octoplot-advanced)

Everything below explains *why* the default figure looks the way it does, and how to
reach the machinery underneath it. None of it is needed to make a plot.

## How the panels are chosen

`octoplot` derives its panels from the model, not from a list of flags. Every
observation declares its plot channels ([`plotchannels`](@ref)), each of which is
one measured quantity with an associated [`ObservableQuery`](@ref) — for example
`radvel(:A, Barycentre)` for a stellar reflex RV instrument, or `raoff(:b, :A)` and
`decoff(:b, :A)` for relative astrometry.

The sky panel is drawn whenever the model has angular observables (i.e. the system
block defines `plx`) and at least one orbit. It shows one phase-coloured track per
**hierarchy row** — each row plotted as its exterior side relative to its interior
side, which is exactly the relationship that row parametrizes. For a star with
planets that is each planet about the star; for a hierarchical system it generalizes
with no special cases (a moon about its planet, an inner pair's barycentre about the
outer body).

The queries used for a *predicted* channel are chosen by [`default_queries`](@ref):
one per hierarchy row for separations and angles, and the host body against the
system barycentre for reflex quantities like `radvel` and `pmra`.

!!! note "Not every observation type has plot channels yet"
    `RelAstromObs`, `RadialVelocityObs`, `MarginalizedRVObs` and `GaiaDR4AstromObs`
    declare plot channels and appear in `octoplot`. `ImageObs`, `InterferometryObs`,
    `PhotometryObs`, `G23HObs` and `HipparcosIADObs` do not yet: `octoplot` still
    draws the orbits, but those observations' data are not overlaid. This is a gap in
    the plotting layer rather than a deliberate omission.

## One panel per instrument — except where sharing is honest

Channels group by the quantity they measure **and** by which observation reported
it, unless the observation type opts into sharing via [`sharepanel`](@ref). The
question that setting answers is whether one draw's calibration can stand in for all
of them:

* **Relative astrometry shares.** `platescale` and `northangle` are calibration
  constants pinned to a fraction of a percent, so every instrument's points land in
  the same place under any draw. All of them merge onto one separation panel, one
  position-angle panel, one Δα\* panel and one Δδ panel, calibrated by the
  maximum-posterior draw — the same behaviour Octofitter has always had.
* **Radial velocity does not.** An instrument zero point is a free parameter of order
  the data range and moves visibly between draws. Each spectrograph gets its own
  panel, its measurements are drawn **exactly as reported**, and each draw's model
  curve carries that draw's own offset and trend. The data stay put; the curves move.
  That is the only way a posterior cloud and a dataset can share a panel without
  misrepresenting one of them.

[`rvplot`](@ref) is the deliberate exception on the other side: it puts every RV
instrument back on one axis precisely because it shows a single draw.

## Every measurement on every panel it belongs on

`(sep, pa)` and `(Δα*, Δδ)` are the same measurement in two bases, related by a
rotation with no free parameter in between. A campaign that switched conventions
partway through is still one dataset, so **all** of its astrometry appears on all
four panels — the non-native pair projected deterministically, with first-order
error propagation. The likelihood still scores only the pair the table actually
carries; the projected channels are marked `derived` and exist for the figure.

A derived channel only opens a panel of its own if you ask for it by name
(`channels=:sep`, as above); otherwise it rides along with a panel some other
instrument declared natively.

## Why residuals are whitened, and why they are boxplots

When more than one draw is shown the residuals are **whitened**, and this is not
optional: a raw residual is not comparable between draws — the jitter, the offsets
and the other bodies' subtracted signals all move — so `whiten=false` throws unless
you also pass `ndraws=1`. `σ_eff` includes fitted jitter and, where an observation
has a Gaussian-process noise model, that model's predictive variance; its mean is
subtracted from the residual too, so the strip shows what the fit is actually left
with.

The same argument applies once more inside the strip: those nuisance terms move from
draw to draw, so a measurement's z-score is not one number but a distribution, and a
single mark per point cannot say how much of the scatter is the fit's own
uncertainty. Each point is therefore drawn as a **boxplot** of `(data − model)/σ_eff`
over the draws — median, quartiles, 1.5 IQR whiskers — with the marginal histogram
pooling every draw against a unit normal. Boxes are sized automatically from the
spacing of the epochs and are narrow by design; `boxwidth=` (in x units — days on a
time panel, cycles on a phase panel) sets them by hand when a dataset defeats that.
With `ndraws=1` there is one z-score per point and it is drawn as a marker instead.

Correlated-noise bands ([`noisemodel`](@ref)) are likewise off unless a single draw
is shown: 250 draws of "orbit" and 250 of "orbit + activity" on one axis is twice
the ink for an envelope that is a different draw's activity model at every epoch.
`gpband=true` overrides that.

## Why phase folds are off by default

Folding collapses the epoch axis through one ephemeris. For a single draw that is
its own period and its own phase zero, and the fold is exact. For many draws the
data can still only be folded once — on the maximum-posterior period — so a chain
whose period is broad, or multi-modal, has its measurements placed at phases most of
the drawn orbits disagree with, and the panel reads as scatter about a curve instead
of as the ambiguity it actually is. Plenty of posteriors are tight enough for it to
look fine, which is the trap: whether the picture is honest is a property of the fit,
not of the figure.

A channel is foldable when one hierarchy row's contribution to it can be isolated
exactly — always true when a single row moves the quantity, and true for observables
linear in the separation otherwise. Each draw's *curve* is folded on its own period
even in the many-draw case, so the spread of the curves is the fold uncertainty
drawn honestly; it is the single folding of the data underneath them that the
default is cautious about.

## Model curves as numbers

To get the model curves as plain numbers rather than a plot, evaluate an
[`ObservableQuery`](@ref) over the series:

```@example 1
q = ObservableQuery(radvel, :A, Barycentre)

curves = modelcurves(series, q)   # one vector per draw, over series.ts
best   = mapcurve(series, q)      # the maximum-posterior draw

length(curves), length(best)
```

For the data rather than the curves, `Octofitter.residuals(obs, ctx)` returns, per
channel, the epochs, the calibrated data, the model, the residual and the effective
uncertainty; [`obscontext`](@ref) builds the `ctx` from a `PosteriorSeries`.

## [Positions at future epochs](@id octoplot-future-epochs)

Solve the system for a draw and query it at whatever epochs you like — this is the
route to any epoch outside the figure's own grid:

```@example 1
posys = construct_system(model, chain, 1)
sols  = orbitsolve(posys, [mjd("2028-01-01"), mjd("2029-01-01")])

raoff(sols[1], :b, :A)                 # ΔRA  [mas]
decoff(sols[1], :b, :A)                # ΔDec [mas]
projectedseparation(sols[1], :b, :A)   # [mas]
posangle(sols[1], :b, :A)              # [rad]
```

Loop over draws to build a predicted distribution, and plot it however you like.
