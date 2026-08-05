# Adding a Custom Observation Type

It's fairly straightforward to add support for a new kind of observation to Octofitter.jl.
You can also follow the same workflow if you want to handle an existing kind of observation in a new way—say, tweaking a calculation, or using Gaussian processes to better model noise in radial velocity data.

All the existing observation types are listed in [`src/likelihoods`](https://github.com/sefffal/Octofitter.jl/tree/main/src/likelihoods) and can be used as examples. Two are worth reading first:

* [`PhotometryObs`](@ref) (`src/likelihoods/photometry.jl`) — the smallest complete type. It touches no epoch and no orbit, so it is nothing but the interface.
* [`RelAstromObs`](@ref) (`src/likelihoods/relative-astrometry.jl`) — the canonical one. Everything below is a stripped-down version of it.

## First: is it an observation at all?

This is the question to settle before writing any code, and the answer follows one rule:

!!! note "Observation, or prior term?"
    **If it carries data you might want to hold out or simulate, it is an observation.
    If it only reshapes the prior, it is a `~` line in a `@variables` block.**

A `~` line whose left-hand side is an expression rather than a fresh variable name is already a fully-supported likelihood contribution. It needs no type, no file, and no method:

```@example 1
using Octofitter, Distributions

A = Body(name="A", variables=@variables begin
    mass ~ truncated(Normal(1.2, 0.1), lower=0.1)
end)

b = Body(name="b", about=A, variables=@variables begin
    a ~ LogUniform(1.0, 100.0)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ Uniform(0, 2pi)
    Ω ~ Uniform(0, 2pi)
    θ ~ Uniform(0, 2pi)
    epoch = 50420.0
    # An ad-hoc constraint: a published dynamical-stability argument puts the
    # apoapsis near 12 au and rules out anything beyond 20 au. This is a `~`
    # line, not an observation.
    a * (1 + e) ~ truncated(Normal(12.0, 3.0), upper=20.0)
end)
nothing # hide
```

Under the hood that becomes an `Octofitter.UserLikelihood` with `_isprior = true`: it reshapes the prior rather than adding data, so it is reported under the chain's `logprior` rather than its `loglike`, has no table to subset, and has no `generate_from_params`. It is not held out by cross-validation, and it cannot be used to simulate data.

That is the trade. Write a real observation type when:

* the numbers are **measurements**, so you want [Cross-Validation](@ref cross-validation) to be able to hold rows out; or
* you want [`generate_from_params`](@ref) to replicate them, for [Posterior Predictive Checks](@ref), [Simulation Based Calibration](@ref sbc), or completeness mapping; or
* the comparison needs the *solved trajectory* — a position, a velocity, an epoch — rather than the parameters alone.

`PhotometryObs` is the borderline case that was decided this way on purpose: `flux_H ~ Normal(15.0, 3.0)` in a body's block expresses exactly the same arithmetic, but only `PhotometryObs(phot_data; target=b, band=:H, name="NIRC2")` can be cross-validated or simulated.

## The interface

An observation is **data + a `target` reference + a `ref` reference + its own variables**. There is one context type, [`ObsContext`](@ref).

| method | required? | what it does |
|---|---|---|
| `Octofitter.ln_like(obs, ctx::ObsContext)` | **yes** | the log-likelihood |
| `Octofitter.refspecs(obs)` | if it touches an orbit | tuple of reference specs, validated at model-build time |
| `Octofitter.epochs(obs)` | if `obs.table.epoch` is not the whole story | epochs [MJD] to solve at, in table order |
| `Octofitter.likelihoodname(obs)` | no (defaults to `obs.name`) | labels this observation's variables in the chain |
| `Octofitter.likeobj_from_epoch_subset(obs, inds)` | for cross-validation | a copy restricted to rows `inds` |
| `Octofitter.simulate(obs, ctx)` | no | the forward model, as data |
| `Octofitter.generate_from_params(obs, ctx; add_noise)` | for simulation/SBC | a copy carrying simulated data |
| `Octofitter._isprior(obs)` | no (defaults `false`) | `true` for terms that carry no data |
| `Octofitter._refdesc(obs)` | no | one-line description for `show` |
| `Octofitter.plotchannels(obs)` / `Octofitter.residuals(obs, ctx)` | no | participate in `octoplot` |

## Creating an observation type

The first step is to create a new data type to hold the observations. We will write `SepObs`: the projected separation (in mas) between a target and a reference, at a set of epochs.

```@example 1
using Octofitter: AbstractObs, ObsContext, Priors, Derived
using Octofitter.Tables: columnnames

"""
    SepObs(data; target, ref, name, variables=@variables begin end)

Projected separation [mas] of `target` from `ref`.

`data` needs `:epoch` [MJD], `:sep` [mas] and `:σ_sep` [mas] columns.

# Variables
  - `jitter` [mas] — added in quadrature to `σ_sep`.
"""
struct SepObs{TTable<:Table,TT,TR} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    target::TT
    ref::TR
    name::String
end

function SepObs(observations;
                target,
                ref,
                name,
                variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    table = Table(observations)
    Octofitter.equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    expected = (:epoch, :sep, :σ_sep)
    issubset(expected, columnnames(table)) ||
        error("Expected columns $expected")

    # Store *reference specs*, not the model nodes the caller passed. A spec is
    # a singleton carrying its content in type parameters, so resolving it
    # against a sample's system constant-folds and allocates nothing.
    t, r = Octofitter.refspec(target), Octofitter.refspec(ref)
    return SepObs{typeof(table),typeof(t),typeof(r)}(
        table, priors, derived, t, r, String(name))
end
nothing # hide
```

The observation type carries:

- **`table`** — the observational data, as a TypedTables.Table
- **`name`** — used for variable naming in MCMC chains, and unique within a system
- **`priors`** and **`derived`** — observation-specific variables from the `@variables` block
- **`target`** and **`ref`** — the two references this observation measures between

Try to follow the advice in the Julia Manual's performance tips section to ensure you've created a fully "concrete" type. This won't affect correctness, but will be important for performance down the road.

!!! note "References do all the work"
    `raoff(sol, target, ref)` is the whole positional model, and it is defined for every
    reference the grammar can express: a body, another companion, a barycentre, a
    photocentre. A likelihood never reconstructs a companion's position by hand.

## Declaring the references and the epochs

```@example 1
# Validated against the system's body list at model-build time, and printed by `show`.
Octofitter.refspecs(obs::SepObs) = (obs.target, obs.ref)

# The default `epochs` already reads `obs.table.epoch`, so this method is
# redundant here — it is written out to show what the contract is. Return
# `Float64[]` if your observation needs no trajectory at all.
Octofitter.epochs(obs::SepObs) = collect(Float64, obs.table.epoch)
nothing # hide
```

## Create the likelihood function

Now extend `Octofitter.ln_like` for your type. Inside it you receive an [`ObsContext`](@ref):

- **`ctx.θ_system`** — system-level parameters (`ctx.θ_system.plx`) and every body's namespace (`ctx.θ_system.bodies.b.a`)
- **`ctx.θ_obs`** — this observation's own variables, from its `@variables` block
- **`ctx.system`** — the `PlanetOrbits.System` built from this sample's parameters
- **`ctx.traj`** — the trajectory solved *once* for the whole model, over the sorted, deduplicated union of every observation's epochs
- **`ctx.buf`** — a per-sample scratch arena (Bumper.jl) for temporaries
- **`solutionat(ctx, i)`** — the per-epoch solution for **row `i` of this observation's table**

```@example 1
function Octofitter.ln_like(obs::SepObs, ctx::ObsContext)
    # Never assert `Float64`: under ForwardDiff these values are `Dual`s.
    T = Octofitter._system_number_type(ctx.θ_system)

    θ_obs = ctx.θ_obs
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)

    # Resolve the references ONCE, outside the epoch loop. `Octofitter.ref`
    # returns a `BodyRef` or a `WeightedPoint`; `resolverefs(ctx, specs)` does
    # a whole tuple at a time.
    target = Octofitter.ref(ctx, obs.target)
    reference = Octofitter.ref(ctx, obs.ref)

    ll = zero(T)
    for i in eachindex(obs.table.epoch)
        # Use the pre-solved trajectory. Index it by *table row*, not by epoch:
        # `solutionat` maps through this observation's own row → column map.
        sol = solutionat(ctx, i)

        Δα = raoff(sol, target, reference)
        Δδ = decoff(sol, target, reference)
        sep_model = hypot(Δα, Δδ)

        resid = sep_model - obs.table.sep[i]
        σ² = obs.table.σ_sep[i]^2 + jitter^2
        ll += -(1 / 2) * resid^2 / σ² - log(sqrt(2π * σ²))
    end
    return ll
end
nothing # hide
```

Points worth copying:

* **Resolve references outside the loop.** `Octofitter.ref(ctx, spec)` constant-folds; calling it per epoch does not.
* **Use `Octofitter._system_number_type(ctx.θ_system)`, never `Float64`.** Anything typed to `Float64` breaks ForwardDiff, and Octofitter's default sampler is gradient-based.
* **Use `solutionat(ctx, i)`, not `ctx.traj[i]`.** The trajectory is over the deduplicated, sorted epoch *union* of the whole model — several observations share it — so table row `i` is generally not trajectory column `i`.
* **For temporaries, use `ctx.buf`** rather than allocating: `@no_escape ctx.buf begin; v = @alloc(T, n); …; end` (with `using Bumper`). One arena, sized to the model at build time, covers the whole evaluation.

If your likelihood requires new parameters, define them in the `variables` block of the constructor using the `@variables begin; end` syntax. These will be accessible via `ctx.θ_obs`. If any parameter has a restricted domain, make sure the prior is truncated with `Distributions.truncated()` — Octofitter remaps the variable using Bijectors.jl to keep the sampler inside the valid region.

## Using your custom observation

```@example 1
sep_dat = Table(;
    epoch = [50000., 50120, 50240, 50360, 50480, 50600],
    sep   = [509.9, 503.9, 498.3, 493.4, 488.7, 484.9],
    σ_sep = fill(5.0, 6),
)

sep_obs = SepObs(
    sep_dat;
    target = b,
    ref = A,
    name = "MY_INSTRUMENT",
    variables = @variables begin
        jitter ~ LogUniform(0.01, 10.0)   # mas
    end
)
```

Note that references take the full grammar — `ref=Barycentre`, `ref=c`, `target=Photocentre(:H, (b, c))` are all legal, and are checked against the model's body list when the `System` is built.

Observations go in the `System`'s `observations=` list; they are never attached to a body:

```@example 1
sys = System(
    name = "Tutoria",
    bodies = [A, b],
    observations = [sep_obs],
    variables = @variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
```

The observation's own variables appear in the chain as `<observation name>_<variable>`, with the name normalized to a valid identifier — so `jitter` above is `MY_INSTRUMENT_jitter`:

```@example 1
using Random
chain = octofit(model, verbosity=0, iterations=1000, adaptation=1000)
chain[:MY_INSTRUMENT_jitter][1:3]
```

!!! note "How the columns are named"
    An observation's variables are named `<observation>_<variable>`, with no body prefix —
    an observation belongs to the system, not to a body.
    [`Octofitter.checkchain`](@ref) reports the mismatch if you load a chain written
    against a differently-named model.

## Bonus: the generative model

The above is sufficient to start sampling from the posterior. Ideally, you will also add a function that does the reverse: generate observations from a set of parameters. This is what makes your type usable for posterior predictive checks, simulation-based calibration, and injection–recovery.

It is conventional to split it in two — a `simulate` method returning the forward model, and `generate_from_params` wrapping it in the noise model — so that plotting code and the likelihood can share the forward model:

```@example 1
function Octofitter.simulate(obs::SepObs, ctx::ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    target = Octofitter.ref(ctx, obs.target)
    reference = Octofitter.ref(ctx, obs.ref)
    sep_model = map(eachindex(obs.table.epoch)) do i
        sol = solutionat(ctx, i)
        hypot(raoff(sol, target, reference), decoff(sol, target, reference))
    end
    return (; sep_model, epochs=obs.table.epoch)
end

function Octofitter.generate_from_params(obs::SepObs, ctx::ObsContext; add_noise)
    sim = Octofitter.simulate(obs, ctx)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : 0.0
    sep = sim.sep_model
    if add_noise
        sep = sep .+ randn.() .* hypot.(obs.table.σ_sep, jitter)
    end
    # Return a NEW observation of the same type, with everything but the
    # measured values carried over.
    return SepObs(
        Table(; epoch=obs.table.epoch, sep, σ_sep=obs.table.σ_sep);
        target = obs.target,
        ref = obs.ref,
        obs.name,
        variables = (obs.priors, obs.derived),
    )
end
nothing # hide
```

Now the whole simulation machinery works on your type:

```@example 1
Random.seed!(0)
sim_sys = generate_from_params(sys; add_noise=true)
sim_sys.observations[1].table.sep
```

## Bonus: cross-validation

Implement `likeobj_from_epoch_subset` and your type can be held out row by row, which is what [`Octofitter.pointwise_like`](@ref) and PSIS-LOO need:

```@example 1
function Octofitter.likeobj_from_epoch_subset(obs::SepObs, inds)
    return SepObs(
        obs.table[inds, :, 1];
        target = obs.target,
        ref = obs.ref,
        obs.name,
        variables = (obs.priors, obs.derived),
    )
end

likelihood_mat, epochs_out = Octofitter.pointwise_like(model, chain, verbosity=0)
size(likelihood_mat)
```

!!! note
    `likeobj_from_epoch_subset(obs, inds)` **keeps** the rows in `inds`; it does not drop
    them.

    Make sure every keyword survives the round trip. If you drop `name` or `variables`
    here, every derived model gets a differently-named observation and its chain will not
    line up with the original.

If your likelihood genuinely does not decompose into independent per-row terms — as with `MarginalizedRVObs`, whose analytic marginalization couples every point in the instrument — do **not** implement this method. The default raises an informative error, and Octofitter's cross-validation machinery turns that into a message naming your observation and the analysis that wanted it.

## Bonus: appearing in `octoplot`

Implement `Octofitter.plotchannels(obs)` and `Octofitter.residuals(obs, ctx)` to get time-series and residual panels for free. See `src/plotting-api.jl` for the `PlotChannel` type, and `RelAstromObs`/`RadialVelocityObs` for worked examples. Without them, `octoplot` still draws the model's orbit panels — it simply does not overlay your data.

## Where to look next

* `src/model/obs.jl` — the interface, and the three "membership tiers" for observations of blended sources.
* `src/model/refs.jl` — the reference grammar: `BodyRefSpec`, `Barycentre`, `Photocentre`, `refspec`, `resolveref`.
* `src/likelihoods/sky-offset.jl` — [`sky_offset`](@ref) and [`sky_calibration`](@ref), the shared front-end for any likelihood that wants a sky-plane offset with a plate scale and a north angle. Use it rather than rolling your own rotate-and-scale block; the sign convention is stated there in full.
