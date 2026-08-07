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
nothing # hide
```

## Calculating Pointwise Likelihoods

After you have defined a model and sampled from its posterior (e.g. via `octofit`), you can see how each datapoint is influencing the posterior:

```@example 1
likelihood_mat, epochs = Octofitter.pointwise_like(model, chain)
size(likelihood_mat)
```

`likelihood_mat` is an `N_sample × N_data` matrix, and `epochs` labels each column with the epoch of the data row it came from.

The columns are ordered exactly as the data are defined in the model: observation by observation in the order you listed them in `System(...; observations=[...])`, and within each observation, row by row in table order.

```@example 1
epochs
```

Two properties are worth knowing about, because they are what makes the matrix usable as a PSIS-LOO input:

!!! note "Prior-shaped terms are excluded"
    A `~` line written inside a `@variables` block, an `LL += ...` line, and the
    `UnitLengthPrior` that sits behind every `UniformCircular` are all *prior* terms:
    they reshape the prior rather than adding data, so they do not participate in this
    machinery and get no column.

    `ObsPriorONeil2019` is *not* a prior-shaped term in this sense — it wraps a real
    likelihood — so it still contributes one column per row of the observation it wraps.

The consequence is that the columns sum to the model's log-likelihood **minus** those prior terms — which is the quantity PSIS-LOO wants:

```@example 1
θ = Octofitter.mcmcchain2result(model, chain, 1)
lnlike = Octofitter.make_ln_like(model.system)
(sum(likelihood_mat[1, :]), lnlike(model.system, θ))
```

The difference is exactly the three `UnitLengthPrior` terms this model's three `UniformCircular` variables contribute.

!!! note "Observations with no epochs get one column, labelled `NaN`"
    [`PhotometryObs`](@ref) carries data but has no epoch column; it contributes one
    column per photometry row, and those columns are labelled `NaN` in `epochs`.

## Pareto-Smoothed Importance Sampling

Once you have `likelihood_mat` you can use the Julia package [ParetoSmooth.jl](https://turinglang.org/ParetoSmooth.jl/dev/) to efficiently calculate a leave-one-out cross-validation score. This technique takes a single posterior chain and, using the pointwise likelihoods, estimates what the posterior would have been had each datapoint in turn been held out.

In broad terms, one might say that this test verifies that no individual datapoints are overly skewing the results.

`psis_loo` comes from ParetoSmooth.jl, not from Octofitter; Octofitter supplies the matrix it consumes. Note the transpose — ParetoSmooth wants `N_data × N_sample`:

```julia
using ParetoSmooth
result = psis_loo(
    collect(likelihood_mat'),
    chain_index=ones(Int,size(chain,1))
)
```
```
Results of PSIS-LOO-CV with 1000 Monte Carlo samples and 8 data points.
┌───────────┬────────┬──────────┬───────┬─────────┐
│           │  total │ se_total │  mean │ se_mean │
├───────────┼────────┼──────────┼───────┼─────────┤
│   cv_elpd │ -53.92 │     0.44 │ -6.74 │    0.06 │
│ naive_lpd │ -53.33 │     0.26 │ -6.67 │    0.03 │
│     p_eff │   0.60 │     0.18 │  0.07 │    0.02 │
└───────────┴────────┴──────────┴───────┴─────────┘
```

`chain_index` tells ParetoSmooth which chain each row of the matrix came from. `pointwise_like` flattens all chains into the sample dimension in order, so for a single chain `ones(Int, size(chain,1))` is right; for several chains use `repeat(1:size(chain,3), inner=size(chain,1))`.

The diagnostic to look at is the Pareto shape parameter `k` for each point. Values above about 0.7 mean the importance-sampling approximation is unreliable for that point — usually because it is highly influential — and that point deserves a genuine refit with the row held out (see [Holding out whole observations](@ref cv-holdout) below).

Plot like so:
```julia
using CairoMakie

fig = Figure()

ax = Axis(
    fig[1,1],
    xlabel="data #",
    ylabel="Pareto K"
)
scatter!(ax, result.pointwise(:pareto_k))
hlines!(ax, [0.7], color=:red, linestyle=:dash)


ax = Axis(
    fig[2,1],
    xlabel="data #",
    ylabel="MCSE"
)
scatter!(ax, result.pointwise(:mcse))


ax = Axis(
    fig[3,1],
    xlabel="data #",
    ylabel="P_EFF"
)
scatter!(ax, result.pointwise(:p_eff))


fig
```

!!! note
    ParetoSmooth.jl is not a dependency of Octofitter, and is not installed as part of
    the documentation build, which is why the two blocks above are not executed here.
    Install it with `] add ParetoSmooth` to run them.

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

In addition, holding out rows from a `RadialVelocityObs` with an `AbstractGPs` based GP is not supported. It works correctly with the **Celerite** backend, and does all the behind the scenes work to do cross-validation correctly in the presence of a Gaussian process.

## See also

* [Posterior Predictive Checks](@ref) — the qualitative counterpart: does the fitted model reproduce the data at all?
* [Bayesian evidence](@ref bayesian-evidence) — comparing whole models rather than individual points.
* [`prior_only_model`](@ref) — a copy of the model with its data likelihoods replaced by no-ops.
