# [Cross-Validation](@id cross-validation)

Cross-validation asks how well a fitted model predicts data it has not seen. In Octofitter that question is answered by *subsetting observations*: every observation type that carries a data table knows how to produce a copy of itself restricted to a chosen set of rows, and the machinery below builds derived models out of those copies.

There are two granularities:

* **Pointwise** — one held-out data row at a time. This is what leave-one-out cross-validation and PSIS-LOO need, and [`Octofitter.pointwise_like`](@ref) computes it from a single chain with no refitting.
* **Whole likelihoods, or groups of epochs** — drop an entire instrument, keep only one, or feed the data in cumulatively. These produce new `System`s that you refit.

We will use a model with two relative-astrometry instruments so that both granularities have something to say:

```@example 1
using Octofitter
using Distributions
using CairoMakie
using Random

astrom_dat_gpi = Table(;
    epoch = [50000., 50120, 50240, 50360],
    ra    = [-505.764, -502.57, -498.209, -492.678],
    dec   = [-66.9298, -37.4722, -7.92755, 21.6356],
    σ_ra  = fill(10.0, 4),
    σ_dec = fill(10.0, 4),
    cor   = fill(0.0, 4),
)
astrom_dat_sphere = Table(;
    epoch = [50480., 50600, 50720, 50840],
    ra    = [-485.977, -478.11, -469.08, -458.896],
    dec   = [51.1472, 80.5359, 109.729, 138.651],
    σ_ra  = fill(10.0, 4),
    σ_dec = fill(10.0, 4),
    cor   = fill(0.0, 4),
)

A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.2, 0.1), lower=0.1)  # M⊙
    end
)
b = Body(
    name="b",
    about=A,
    variables=@variables begin
        a ~ truncated(Normal(10, 4), lower=0.1, upper=100)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 50420.0
    end
)

gpi    = RelAstromObs(astrom_dat_gpi;    target=b, ref=A, name="GPI")
sphere = RelAstromObs(astrom_dat_sphere; target=b, ref=A, name="SPHERE")

sys = System(
    name="Tutoria",
    bodies=[A, b],
    observations=[gpi, sphere],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
Random.seed!(0)
chain = octofit(model)
display(chain)
nothing # hide
```

## Calculating Pointwise Likelihoods

After you have defined a model and sampled from its posterior (e.g. via `octofit`), you can see how each datapoint is influencing the posterior:

```@example 1
@time likelihood_mat, epochs = Octofitter.pointwise_like(model, chain)
size(likelihood_mat)
```

`likelihood_mat` is an `N_sample × N_data` matrix, and `epochs` labels each column with the epoch of the data row it came from.

The columns are ordered exactly as the data are defined in the model: observation by observation in the order you listed them in `System(...; observations=[...])`, and within each observation, row by row in table order.

```@example 1
epochs
```


!!! note "Prior-shaped terms are excluded"
    A `~` line written inside a `@variables` block, an `LL += ...` line, and the
    `UnitLengthPrior` that sits behind every `UniformCircular` are all *prior* terms:
    they reshape the prior rather than adding data, so they do not participate in this
    machinery and get no column.

    `ObsPriorONeil2019` is *not* a prior-shaped term in this sense — it wraps a real
    likelihood — so it still contributes one column per row of the observation it wraps.

The consequence is that the columns sum to the model's log-likelihood **minus** those prior terms — which is the quantity used by PSIS-LOO:

```@example 1
θ = Octofitter.mcmcchain2result(model, chain, 1)
lnlike = Octofitter.make_ln_like(model.system)
(sum(likelihood_mat[1, :]), lnlike(model.system, θ))
```

The difference is exactly the three `UnitLengthPrior` terms this model's three `UniformCircular` variables contribute.

Observations with no epochs get one column, labelled `NaN`.
[`PhotometryObs`](@ref) carries data but has no epoch column; it contributes one
column per photometry row, and those columns are labelled `NaN` in `epochs`.

## Pareto-Smoothed Importance Sampling

Once you have `likelihood_mat` you can use the Julia package [ParetoSmooth.jl](https://turinglang.org/ParetoSmooth.jl/dev/) to efficiently calculate a leave-one-out cross-validation score. This technique takes a single posterior chain and, using the pointwise likelihoods, estimates what the posterior would have been had each datapoint in turn been held out.

In broad terms, one might say that this test verifies that no individual datapoints are overly skewing the results.

`psis_loo` comes from ParetoSmooth.jl, not from Octofitter; Octofitter supplies the matrix it consumes. Note the transpose — ParetoSmooth wants `N_data × N_sample`:

```@example 1
using ParetoSmooth
result = psis_loo(
    collect(likelihood_mat'),
    chain_index=ones(Int, size(chain, 1))
)
display(result)
nothing # hide
```

`chain_index` tells ParetoSmooth which chain each row of the matrix came from. `pointwise_like` flattens all chains into the sample dimension in order, so for a single chain `ones(Int, size(chain,1))` is right; for several chains use `repeat(1:size(chain,3), inner=size(chain,1))`.

`result.pointwise` is a `KeyedArray` of per-row diagnostics, indexed by name:

```@example 1
result.pointwise(:pareto_k)
```

The available statistics are `:cv_elpd` (the leave-one-out expected log predictive
density of that row), `:naive_lpd` (the same quantity without holding the row out),
`:p_eff` (the difference between them — how much of that row the fit has effectively
absorbed), `:mcse`, and `:pareto_k`.

The diagnostic to look at is the Pareto shape parameter ``\hat{k}`` for each point. Values above about 0.7 mean the importance-sampling approximation is unreliable for that point — usually because it is highly influential — and that point deserves a genuine refit with the row held out (see [Holding out whole observations](@ref cv-holdout) below).

```@example 1
using CairoMakie

nrows = size(likelihood_mat, 2)
fig = Figure(size=(650, 620))

ax = Axis(
    fig[1,1],
    xlabel="data row",
    ylabel="Pareto k̂",
    xticks=1:nrows
)
scatter!(ax, collect(result.pointwise(:pareto_k)))
hlines!(ax, [0.7], color=:red, linestyle=:dash)
ylims!(ax, -0.1, 0.9)   # keep the threshold visible even when nothing is near it


ax = Axis(
    fig[2,1],
    xlabel="data row",
    ylabel="MCSE",
    xticks=1:nrows
)
scatter!(ax, collect(result.pointwise(:mcse)))


ax = Axis(
    fig[3,1],
    xlabel="data row",
    ylabel="p_eff",
    xticks=1:nrows
)
scatter!(ax, collect(result.pointwise(:p_eff)))


fig
```

Every ``\hat{k}`` here is comfortably below the line, which is what a healthy fit looks
like: no single astrometric epoch is holding the posterior up on its own, so the
LOO score above can be trusted as it stands.

### [Finding an influential point](@id cv-influential)

The diagnostic is only interesting when something *is* wrong, so here is the same model
with one SPHERE measurement displaced by 60 mas — six times its quoted uncertainty —
standing in for a mis-registered frame or a background object mistaken for the companion:

```@example 1
astrom_dat_sphere_bad = Table(;
    epoch = astrom_dat_sphere.epoch,
    ra    = astrom_dat_sphere.ra,
    dec   = [51.1472, 80.5359 + 60.0, 109.729, 138.651],   # row 2 nudged
    σ_ra  = fill(10.0, 4),
    σ_dec = fill(10.0, 4),
    cor   = fill(0.0, 4),
)

sys_bad = System(
    name="Tutoria_outlier",
    bodies=[A, b],
    observations=[
        gpi,
        RelAstromObs(astrom_dat_sphere_bad; target=b, ref=A, name="SPHERE"),
    ],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model_bad = Octofitter.LogDensityModel(sys_bad)
Random.seed!(0)
chain_bad = octofit(model_bad, verbosity=0)
mat_bad, epochs_bad = Octofitter.pointwise_like(model_bad, chain_bad)
result_bad = psis_loo(collect(mat_bad'), chain_index=ones(Int, size(chain_bad, 1)))
display(result_bad)
nothing # hide
```

Two things changed. `p_eff` — the number of parameters the data are effectively
buying — has jumped, because one row is now doing work no smooth orbit can absorb; and
ParetoSmooth has warned that some ``\hat{k}`` crossed 0.7, which it does on its own
whenever the importance-sampling approximation has failed somewhere.

To find *which* rows those are, threshold the vector — and, because the columns are in
model order, map them back to the observation and epoch they came from with the same walk
over the observations used above:

```@example 1
k̂ = collect(result_bad.pointwise(:pareto_k))
flagged = findall(>(0.7), k̂)

data_obs_bad = filter(!Octofitter._isprior, sys_bad.observations)
rowinst = reduce(vcat, [fill(Octofitter.likelihoodname(o), max(Octofitter._nrows(o), 1))
                        for o in data_obs_bad])

[(row=i, instrument=rowinst[i], epoch=epochs_bad[i], k̂=round(k̂[i], digits=2))
 for i in flagged]
```

Now plot it. The top panel is the diagnostic — one ``\hat{k}`` per data row, with the 0.7
threshold drawn and the rows above it picked out and labelled by epoch. The bottom panel
is a [`skypanel!`](@ref) of the same fit, zoomed to the data, so you can see *where* the
flagged epochs sit relative to the posterior orbits:

```@example 1
fig = Figure(size=(700, 900))

ax = Axis(fig[1,1], xlabel="data row", ylabel="Pareto k̂",
          xticks=1:length(k̂), title="PSIS-LOO influence diagnostic")
hlines!(ax, [0.7], color=:red, linestyle=:dash)
hlines!(ax, [0.5], color=:gray, linestyle=:dot)
fine = setdiff(eachindex(k̂), flagged)
scatter!(ax, fine, k̂[fine], color=:black, markersize=12)
scatter!(ax, flagged, k̂[flagged], color=:red, marker=:diamond, markersize=17)
for i in flagged
    right = i > length(k̂) / 2      # keep the label inside the axis
    text!(ax, i, k̂[i];
        text = "epoch $(round(Int, epochs_bad[i]))\n($(rowinst[i]))",
        align = (right ? :right : :left, :center),
        offset = (right ? -12 : 12, 0),
        color = :red, fontsize = 12)
end
ylims!(ax, -0.1, 1.15)
rowsize!(fig.layout, 1, Relative(0.3))

# The same draws Octofitter's own figures use, restricted to the data region.
series = PosteriorSeries(model_bad, chain_bad; N=50)
skyax = skypanel!(fig[2,1], series; colorbar=false).sky
ra_all  = vcat(astrom_dat_gpi.ra,  astrom_dat_sphere_bad.ra)
dec_all = vcat(astrom_dat_gpi.dec, astrom_dat_sphere_bad.dec)
scatter!(skyax, ra_all[flagged], dec_all[flagged],
    marker=:circle, markersize=26, strokewidth=3, strokecolor=:red,
    color=(:white, 0.0))
# Equal spans, so the sky panel's DataAspect gives a square view.
ra_lo, ra_hi = extrema(ra_all)
dec_lo, dec_hi = extrema(dec_all)
half = max(ra_hi - ra_lo, dec_hi - dec_lo) / 2 + 60
cx, cy = (ra_lo + ra_hi) / 2, (dec_lo + dec_hi) / 2
limits!(skyax, cx - half, cx + half, cy - half, cy + half)

fig
```

The corrupted epoch is the one the top panel picks out, and the ring around it in the
bottom panel shows why: it sits off the track that every posterior draw wants to follow,
so the posterior is being pulled by that single row and holding it out would move the
answer.

It is not usually alone. The last epoch of the campaign crosses the line here too, and
that is worth understanding, because nothing was done to it: the ends of a baseline are
the rows an orbit fit leans on hardest, which is already visible in the `p_eff` panel of
the healthy fit above — rows 1 and 8 are several times any interior row. Tilting the
orbit to accommodate a displaced point therefore lands on whichever rows were carrying
the fit anyway. Read ``\hat{k}`` as pointing at a *region* of the dataset that has become
load-bearing, not as a list of bad measurements.

What to do about a flagged row is a judgement call, not a rule. PSIS-LOO is telling you
that its approximation failed there, not that the measurement is wrong. Refit with the
row genuinely held out (below) and compare; if the posterior moves materially, the
measurement deserves a look before the orbit does.

## [Holding out whole observations](@id cv-holdout)

The functions below return new `System` objects, which you then wrap in a `LogDensityModel` and refit. Every one of them keeps the model's prior-shaped terms intact — dropping the `UnitLengthPrior` behind a `UniformCircular` would change the *prior*, not the data set.

Leave one instrument out at a time:

```@example 1
kfold_systems = Octofitter.generate_kfold_systems(sys)
map(s -> Octofitter.likelihoodname.(s.observations), kfold_systems)
```

Keep only one instrument at a time — useful for checking that two instruments agree before combining them:

```@example 1
per_like_systems = Octofitter.generate_systems_per_like(sys)
map(s -> Octofitter.likelihoodname.(s.observations), per_like_systems)
```

Or select by an arbitrary predicate on the observation:

```@example 1
filtered = Octofitter.generate_system_filtered_like(
    o -> Octofitter.likelihoodname(o) == "GPI", sys)
Octofitter.likelihoodname.(filtered.observations)
```

Refit any of them the usual way:

```@example 1
gpi_only_model = Octofitter.LogDensityModel(filtered)
gpi_only_chain = octofit(gpi_only_model, verbosity=0)
nothing # hide
```

### Scoring the held-out dataset

Refitting without a dataset is only half the exercise; the other half is asking how well
the *reduced* posterior predicts the data you removed. Evaluate the full model's
per-row likelihoods at the reduced model's draws, and keep the columns belonging to the
held-out observation:

```@example 1
# Which columns belong to each observation, in the order described above:
# observation by observation, prior-shaped terms skipped, and one column per
# table row (or a single column for an observation that carries no table).
data_obs = filter(!Octofitter._isprior, sys.observations)
widths = [max(Octofitter._nrows(o), 1) for o in data_obs]
bounds = cumsum(widths)
starts = [1; bounds[1:end-1] .+ 1]
heldout = [i for (o, s, e) in zip(data_obs, starts, bounds)
             if Octofitter.likelihoodname(o) != "GPI" for i in s:e]

# Per-row likelihoods of *all* the data, under the GPI-only posterior
mat_all, _ = Octofitter.pointwise_like(model, gpi_only_chain)
mat_heldout = mat_all[:, heldout]

# Log pointwise predictive density of the held-out rows: log mean exp over draws,
# summed over rows.
using StatsBase: mean
lppd = sum(log(mean(exp, col)) for col in eachcol(mat_heldout))
```

`lppd` is on the same scale as a log Bayes factor per held-out dataset: larger is better,
and comparing it across folds says which instruments the model is and is not able to
predict from the others. It is not as directly interpretable as PSIS-LOO — there is no
importance-sampling diagnostic to warn you when it is unreliable — but it needs no
approximation either, because each fold is a genuine refit. Note that an instrument
carrying its own `offset`/`jitter` variables is being predicted *with those nuisances free*,
so a poor score there can mean a bad calibration rather than a bad orbit.

## Epoch-level folds

Data rows are numbered globally — observation by observation, then table row by table row — so row 5 means "the fifth data row in the model" regardless of which observation it lives in.

One system per data row, each containing only that row:

```@example 1
per_epoch_systems, per_epoch_epochs = Octofitter.generate_system_per_epoch(sys)
(length(per_epoch_systems), per_epoch_epochs)
```

One system per data row, each containing rows 1 through *i* — for watching a posterior tighten as data accumulate:

```@example 1
cumulative_systems, cumulative_epochs = Octofitter.generate_cumulative_system_per_epoch(sys)
length.(cumulative_epochs)
```

And the general form, which takes the row groups you want:

```@example 1
grouped, grouped_epochs = Octofitter.generate_systems_with_epoch_groups(
    sys,
    [[1, 2], [3, 4, 5], [6, 7, 8]],
    g -> "_group_$g",
)
grouped_epochs
```

## Observations that cannot be subset

Not every likelihood decomposes into independent per-row terms. The canonical case is `MarginalizedRVObs` from OctofitterRadialVelocity: it integrates the instrument's zero point out analytically, which couples every point in that instrument, so "the likelihood of row 7" does not exist as a quantity. Asking for it produces an error naming the offending observation and its type, rather than a bare failure from inside the machinery:

```julia
Octofitter.pointwise_like(model_with_marginalized_rv, chain)
# ERROR: ArgumentError: Cannot subset the data of observation "HIRES" (a MarginalizedRVObs), …
```

Use a plain `RadialVelocityObs` with an explicit `offset` variable if you need to cross-validate radial velocities.

Holding out rows from a `RadialVelocityObs` that carries a Gaussian process works with **both** GP backends — **Celerite** and **AbstractGPs** — and does all the behind the scenes work to do cross-validation correctly in the presence of a Gaussian process: the process is conditioned on the retained rows only, and each held-out point is then scored against the predictive mean and variance there, with that point's own measurement error and jitter added.

## See also

* [Posterior Predictive Checks](@ref) — the qualitative counterpart: does the fitted model reproduce the data at all?
* [Bayesian evidence](@ref bayesian-evidence) — comparing whole models rather than individual points.
* [`prior_only_model`](@ref) — a copy of the model with its data likelihoods replaced by no-ops.
