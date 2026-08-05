# ---------------------------------------------------
# Simulation-based calibration   (agent G)
#
# Draws from the prior, simulates data with `generate_from_params`, refits,
# and checks rank uniformity.
#
# This file also holds the system-level `generate_from_params` driver that SBC,
# completeness mapping and the data-simulation tutorials all call. The per-type
# `generate_from_params(obs, ctx; add_noise)` methods live with their
# observation types; this is the thing that solves the trajectory once and
# hands each of them its context.
# ---------------------------------------------------

"""
    generate_from_params(obs::AbstractObs, ctx::ObsContext; add_noise)

Fallback: a prior-shaped term has no data to regenerate, so it passes through
unchanged. An observation that *does* carry data and has not implemented
forward simulation is an error rather than a silent pass-through — otherwise
a "simulated" system would quietly retain the real measurements, and every
downstream check (SBC rank statistics, completeness) would be measuring the
wrong thing.
"""
function generate_from_params(@nospecialize(obs::AbstractObs), ctx::ObsContext; add_noise)
    _isprior(obs) && return obs
    error("""
    Observation "$(likelihoodname(obs))" (a $(typeof(obs).name.name)) does not implement
    `generate_from_params`, so data cannot be simulated from this model.

    Implement `generate_from_params(obs::$(typeof(obs).name.name), ctx::ObsContext; add_noise)`
    for it, or drop it from the system before simulating.
    """)
end

"""
    generate_from_params(system, θ=drawfrompriors(system); add_noise=false)

A copy of `system` whose observations hold data simulated at the parameters
`θ` — the model's own generative direction, used for prior-predictive checks,
posterior-predictive checks, simulation-based calibration and
injection-recovery.

`θ` is a nested parameter NamedTuple, as produced by
[`drawfrompriors`](@ref), `model.arr2nt(θ_flat)`, or
[`mcmcchain2result`](@ref).

With `add_noise=true` each observation adds a draw from its own noise model
(including any jitter `θ` specifies); with `add_noise=false` the data are the
noiseless model prediction, which is what you want to check that a fit
recovers a known truth exactly.

The system is solved once over the union of every observation's epochs and
each observation is handed the same `ObsContext` its likelihood would see, so
simulated data and the likelihood that will be fit to them cannot drift apart.
"""
function generate_from_params(sys::System, θ=drawfrompriors(sys); add_noise::Bool=false)
    unique_ep, maps = epoch_plan(sys)
    posys = _build_po_system(sys, θ)
    traj = isempty(unique_ep) ? nothing :
           orbitsolve(posys, unique_ep; _solvekw(sys)...)
    obsns = hasproperty(θ, :observations) ? θ.observations : (;)

    newobs = map(sys.observations) do obs
        key = normalizename(likelihoodname(obs))
        θ_obs = hasproperty(obsns, key) ? getproperty(obsns, key) : (;)
        return generate_from_params(obs, ObsContext(θ, θ_obs, posys, traj, maps[obs]);
            add_noise)
    end

    return _rebuild_system(sys; observations=newobs)
end
export generate_from_params

"""
    _with_seed(f, seed)

Run `f()` with the task's default RNG seeded to `seed`, then put the caller's
stream back exactly as it was.

The per-type `generate_from_params(obs, ctx; add_noise)` methods draw their
noise straight from the default RNG — the interface has no `rng` argument to
thread one through. That leaves a "reproducible from its seed" trial only half
reproducible: the parameter draw is seeded, the noise is not. Seeding around
the call closes that gap without the trial quietly consuming (or resetting)
the caller's stream, which a script that seeded before calling would notice.

Remove this once `generate_from_params` takes an `rng`.
"""
function _with_seed(f, seed)
    rng = Random.default_rng()
    saved = copy(rng)
    try
        Random.seed!(rng, seed)
        return f()
    finally
        copy!(rng, saved)
    end
end

# ---------------------------------------------------
# One SBC trial
# ---------------------------------------------------

"""
    calibrationhmc(system; θ=sample_priors(rng, system), target_accept=0.85, add_noise=true, kwargs...)

One simulation-based-calibration trial (Talts et al. 2018): draw parameters
from the prior, simulate a data set at them, refit, and report where the true
value falls in the resulting posterior.

Returns `(priorsampledict, rdict, chains)` — the drawn parameter values, the
rank of each within the posterior as a percentage, and the chain itself.
Repeated over many trials, a well-specified model and an unbiased sampler give
a flat histogram of ranks for every parameter.

# Keywords
  - `θ` — the flat prior draw to calibrate against; defaults to a fresh one.
  - `add_noise` — whether the simulated data get a draw from the noise model.
  - `target_accept`, `verbosity` and anything else are forwarded to
    [`octofit`](@ref).

!!! warning "Changed from v1"
    `add_noise` now defaults to `true`. v1 simulated *noiseless* data here
    (it took `generate_from_params`'s `add_noise=false` default), which is not
    a draw from the likelihood: ranks computed against it are not the SBC
    statistic and the histograms it produces are not the diagnostic they look
    like. Pass `add_noise=false` to reproduce the old behaviour.
"""
function calibrationhmc(
    system::System;
    rng=Random.default_rng(),
    verbosity=0,
    target_accept=0.85,
    θ=sample_priors(rng, system),
    add_noise::Bool=true,
    kwargs...
)
    θ_flat = θ
    arr2nt = make_arr2nt(system)
    θ_nt = arr2nt(θ_flat)

    verbosity > 0 && println("= SBC " * "="^70)
    verbosity > 0 && println("Drew parameters: ", stringify_nested_named_tuple(θ_nt))

    newsystem = generate_from_params(system, θ_nt; add_noise)

    model = LogDensityModel(newsystem; verbosity=max(0, verbosity - 1))

    chains = octofit(rng, model, target_accept; verbosity, kwargs...)

    # The log posterior and log likelihood *at the truth*, which Modrák et al.
    # (2022) use as a summary statistic alongside the per-parameter ranks: it
    # catches miscalibration that no single marginal shows.
    loglike = make_ln_like(newsystem, θ_nt)(newsystem, θ_nt)
    # `sampled=false` asks for the plain prior density rather than the
    # transformed one — no change-of-variables Jacobian, because θ_flat is in
    # the natural domain. (v1 called `make_ln_prior`, which v2 folded into
    # `make_ln_prior_transformed`'s second argument.)
    logprior = make_ln_prior_transformed(newsystem)(θ_flat, false)
    logpost = logprior + loglike

    θ_array = result2mcmcchain([(; loglike, logpost, θ_nt...)])

    # Rank of the true value within the posterior, for every parameter the
    # chain and the truth have in common.
    chainkeys = string.(intersect(keys(θ_array), keys(chains)))
    priorsampledict = OrderedDict()
    rdict = OrderedDict()
    for key in chainkeys
        posteriorsamples = vec(chains[key])
        priorsample = only(θ_array[key])
        r = count(<(priorsample), posteriorsamples)
        priorsampledict[key] = priorsample
        rdict[key] = r / length(posteriorsamples) * 100   # on a 0-100 scale
    end

    return priorsampledict, rdict, chains
end

"""
    sbctrial(system, chainparams, saveas)

Run one [`calibrationhmc`](@ref) trial and write it to disk as four files:

  - `<saveas>_sampler_parameters.toml` — the sampler settings used
  - `<saveas>_parameters.toml` — the parameter values drawn from the prior
  - `<saveas>_rank_stats.toml` — the rank of each true value in the posterior
  - `<saveas>_chains.fits` — the full chain

`chainparams` is a NamedTuple or Dict of keyword arguments; `θ` (a flat prior
draw) and any sampler settings go in it. Run this script many times with
different seeds, then histogram the rank statistics — see the
"Simulation Based Calibration" tutorial.
"""
function sbctrial(system::System, chainparams, saveas::AbstractString)
    chainparams = Dict{Symbol,Any}(pairs(chainparams))

    verbosity = chainparams[:verbosity] = get(chainparams, :verbosity, 2)
    target_accept = get(chainparams, :target_accept, 0.85)
    delete!(chainparams, :target_accept)

    priorsampledict, rdict, chains =
        calibrationhmc(system; target_accept, chainparams...)

    # TOML cannot represent a distribution or a function, and the sampler
    # settings are the one thing here a user might have put either in.
    open("$(saveas)_sampler_parameters.toml", "w") do f
        TOML.print(f, Dict(string(k) => string(v) for (k, v) in chainparams))
    end
    open("$(saveas)_parameters.toml", "w") do f
        TOML.print(f, Dict(string(k) => v for (k, v) in priorsampledict))
    end
    open("$(saveas)_rank_stats.toml", "w") do f
        TOML.print(f, Dict(string(k) => v for (k, v) in rdict))
    end

    savechain("$(saveas)_chains.fits", chains)

    return nothing
end

# v1 also carried `getplotdata` and `calibrationplots` here. Neither is ported:
# both were dead code that could not have run (an undefined `saveas` inside
# `getplotdata`, `Plots.jl` calls in a package that does not depend on Plots,
# and DataFrames indexing on a `TypedTables.Table`). The SBC tutorial's own
# analysis script — load the per-trial TOML files, filter on R̂, histogram the
# ranks with Makie — is the supported path.
