# ---------------------------------------------------
# Cross-validation and PSIS-LOO
#
# Requires `likeobj_from_epoch_subset` on every observation type it is asked
# to hold out; types that structurally cannot provide pointwise likelihoods
# (`MarginalizedRVObs`, whose analytic offset marginalization couples every
# point in the instrument) must fail loudly rather than silently.
#
# Observations are a flat list on the system, so every walk here is one loop
# and no function has to ask which body a likelihood belongs to.
# ---------------------------------------------------

# ---------------------------------------------------
# Rebuilding a system with a different observation list
# ---------------------------------------------------

# The `@variables` block of an observation, as `System`/`BlankLikelihood`
# want it. Prior-shaped terms own no variables and have no such fields.
_obs_variables(@nospecialize obs::AbstractObs) =
    hasproperty(obs, :priors) ? (obs.priors, obs.derived) : (Priors(), Derived())

# The prior-shaped terms `@variables` emitted in the *system* block. Node-owned
# terms travel with their node, so they are not listed here.
_system_priorterms(sys::System) = Tuple(t for (owner, t) in sys.priorterms if owner === :system)

# A copy of a node with its `@variables`-emitted prior terms dropped. Only
# `prior_only_model(; exclude_all=true)` wants this — see its docstring for
# why dropping a `UnitLengthPrior` is a deliberate, narrow choice.
_strip_priorterms(n::Body) = Body(; name=n.name, about=n.about,
    variables=(n.priors, n.derived))
_strip_priorterms(n::Orbit) = Orbit(; name=n.name, exterior=n.exterior, about=n.about,
    variables=(n.priors, n.derived))

"""
    _rebuild_system(sys; observations, nodes=sys.nodes, sysextra=…, name=sys.name)

`sys` with a different observation list, keeping its variables, its solve
configuration and its hierarchy. Everything in this file that "drops a
likelihood" goes through here, so the reconstruction stays in one place and
the full `System` validation (unique names, references that exist, bands that
are declared) runs on every derived model rather than only on the original.
"""
function _rebuild_system(sys::System;
                         observations,
                         nodes=sys.nodes,
                         sysextra=_system_priorterms(sys),
                         name=sys.name)
    return System(;
        name,
        bodies=nodes,
        observations=Tuple(observations),
        variables=(sys.priors, sys.derived, sysextra...),
        method=sys.method,
        observing_geometry=sys.observing_geometry,
        barycentric_lighttime=sys.barycentric_lighttime,
    )
end

# ---------------------------------------------------
# Prior-only models
# ---------------------------------------------------

"""
    prior_only_model(system; exclude_all=false)

A copy of `system` whose data likelihoods have been stripped out, so that
sampling it samples the prior. Used for prior-predictive checks, for
tempering, and — with `exclude_all=true` — for the reference `log_Z0` that
turns a Pigeons log evidence *ratio* into a log evidence.

Each real observation is replaced by a `BlankLikelihood` carrying the
same name and the same `@variables` block, so the parameter vector keeps its
shape and a chain from the prior-only model lines up column for column with
one from the full model.

`exclude_all=false` (the default) keeps prior-shaped terms — `lhs ~ dist` and
`LL +=` lines, and the `UnitLengthPrior` behind each `UniformCircular`. They
reshape the prior rather than adding data, so a "prior only" model is the one
that still has them.

`exclude_all=true` drops those as well, leaving nothing but the declared
prior distributions. That is what you want when the point is to normalize:

```julia
prior_model = Octofitter.LogDensityModel(Octofitter.prior_only_model(model.system, exclude_all=true))
_, pt_prior = octofit_pigeons(prior_model, n_rounds=10)
log_Z0 = stepping_stone(pt_prior)
```

!!! note
    With `exclude_all=true` a `UniformCircular` variable loses the prior that
    keeps its (x, y) pair off the origin, where the angle is undefined. That
    is intentional — it is exactly the term whose normalization you are
    measuring — but it means such a model is for evidence bookkeeping, not
    for inference.
"""
function prior_only_model(system::System; exclude_all::Bool=false)
    newobs = map(system.observations) do obs
        # A prior-shaped observation carries no data, so there is nothing to
        # strip out of it; keep it whole unless the caller asked for a bare
        # model.
        if !exclude_all && _isprior(obs)
            return obs
        end
        return BlankLikelihood(_obs_variables(obs), likelihoodname(obs))
    end
    nodes = exclude_all ? map(_strip_priorterms, system.nodes) : system.nodes
    sysextra = exclude_all ? () : _system_priorterms(system)
    return _rebuild_system(system; observations=newobs, nodes, sysextra)
end
export prior_only_model

# ---------------------------------------------------
# Subsetting, and failing usefully when it is impossible
# ---------------------------------------------------

"""
    _count_likeobj(system)

Number of real likelihood objects — prior-shaped terms do not count, because
holding one out is not a cross-validation fold.
"""
_count_likeobj(sys::System)::Int = count(o -> !_isprior(o), sys.observations)

"""
    _subset(obs, inds; what)

`likeobj_from_epoch_subset(obs, inds)`, with the failure turned into a
message that names the observation and says what wanted it.

Some observation types cannot supply pointwise likelihoods at all. The
canonical case is `MarginalizedRVObs`: its analytic marginalization over the
instrument's zero point couples every point in that instrument, so "the
likelihood of row 7" does not exist as a quantity, and no amount of
refactoring makes it exist. The generic fallback in `model/obs.jl` throws for
types that simply have not implemented subsetting yet. Both arrive here, and
both should tell the user which observation in *their* model is the problem
rather than surfacing a bare type name from inside a `map`.
"""
function _subset(@nospecialize(obs::AbstractObs), inds; what::AbstractString="this analysis")
    try
        return likeobj_from_epoch_subset(obs, inds)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("""
        Cannot subset the data of observation "$(likelihoodname(obs))" (a $(typeof(obs).name.name)), \
        which $what requires.

        The underlying error was:
        $(sprint(showerror, err))

        Either drop this observation from the model before running $what, or — if it is a \
        type whose likelihood genuinely does not decompose into independent per-epoch terms, \
        such as a marginalized-offset radial velocity likelihood — use the non-marginalized \
        equivalent, which does.
        """))
    end
end

# Global numbering of the model's data rows: one entry per row, in
# (observation order, table order). This single plan is what makes the
# epoch-group functions and `pointwise_like` agree on what "epoch 7" means.
function _row_plan(sys::System)
    plan = Tuple{Int,Int}[]
    for (k, o) in enumerate(sys.observations)
        for r in 1:_nrows(o)
            push!(plan, (k, r))
        end
    end
    return plan
end

# Rows of data an observation can be split on. Prior-shaped terms and
# table-less observations have none, and are carried whole into every subset.
_nrows(@nospecialize obs::AbstractObs)::Int =
    (_isprior(obs) || !hasproperty(obs, :table)) ? 0 : Tables.rowcount(obs.table)

# Representative epoch of one data row, or NaN when the type has no epoch
# column (photometry, say).
function _row_epoch(@nospecialize(obs::AbstractObs), r::Int)
    hasproperty(obs, :table) || return NaN
    hasproperty(obs.table, :epoch) || return NaN
    return Float64(obs.table.epoch[r])
end

# ---------------------------------------------------
# Whole-likelihood folds
# ---------------------------------------------------

"""
    generate_kfold_systems(system)

One copy of `system` per real likelihood object, each with that likelihood
dropped. Prior-shaped terms are kept in every copy.
"""
function generate_kfold_systems(sys::System)
    n = _count_likeobj(sys)
    return map(1:n) do i_obs
        j = 0
        keep = AbstractObs[]
        for o in sys.observations
            if _isprior(o)
                push!(keep, o)
                continue
            end
            j += 1
            j == i_obs || push!(keep, _subset(o, (:); what="k-fold cross-validation"))
        end
        _rebuild_system(sys; observations=keep, name=Symbol(sys.name, "_kfold_", i_obs))
    end
end

"""
    generate_systems_per_like(system)

One copy of `system` per real likelihood object, each keeping only that
likelihood (plus every prior-shaped term).
"""
function generate_systems_per_like(sys::System)
    n = _count_likeobj(sys)
    return map(1:n) do i_obs
        j = 0
        keep = AbstractObs[]
        for o in sys.observations
            if _isprior(o)
                push!(keep, o)
                continue
            end
            j += 1
            j == i_obs && push!(keep, _subset(o, (:); what="per-likelihood model generation"))
        end
        _rebuild_system(sys; observations=keep, name=Symbol(sys.name, "_like_", i_obs))
    end
end

"""
    generate_system_filtered_like(predicate, system)

A copy of `system` keeping only the observations for which
`predicate(obs)` is true. Prior-shaped terms are always kept: they carry no
data, so filtering them changes the prior rather than the data set. (v8
applied the predicate to them too, which silently dropped the
`UnitLengthPrior` behind every `UniformCircular` unless the user's predicate
happened to allow for it.)
"""
function generate_system_filtered_like(predicate, sys::System)
    keep = AbstractObs[]
    kept_j = Int[]
    j = 0
    for o in sys.observations
        if _isprior(o)
            push!(keep, o)
            continue
        end
        j += 1
        if predicate(o)
            push!(keep, _subset(o, (:); what="filtered model generation"))
            push!(kept_j, j)
        end
    end
    return _rebuild_system(sys; observations=keep,
        name=Symbol(sys.name, "_filt_obs_", join(kept_j, '_')))
end

# ---------------------------------------------------
# Epoch-level folds
# ---------------------------------------------------

"""
    generate_systems_with_epoch_groups(system, epoch_groups, name_suffix_func)

One copy of `system` per entry of `epoch_groups`, each containing exactly the
data rows that entry names. Rows are numbered globally, in observation order
then table order, so `3` means "the third data row in the model" regardless
of which observation it lives in. Observations left with no rows are dropped
from that copy; observations with no rows to split (prior-shaped terms,
tables-less types) are included in all of them.

Returns `(systems, epochs)`, where `epochs[i]` lists the epochs actually
included in `systems[i]`.

!!! warning "Changed from v8"
    v8 passed the *complement* of the wanted rows to
    `likeobj_from_epoch_subset`, on the strength of a comment claiming that
    function drops the indices it is given. It does not — every
    implementation of it returns `table[inds]`, i.e. it *keeps* them. So a
    group asking for row `i` got every row except `i`. That inverted
    `generate_system_per_epoch`, and through it `pointwise_like`: each
    "pointwise" column held the likelihood of all the data bar one point.
"""
function generate_systems_with_epoch_groups(sys::System, epoch_groups, name_suffix_func)
    plan = _row_plan(sys)
    isempty(plan) && return System[], Vector{Float64}[]

    epochs_out = Vector{Float64}[]
    systems = map(enumerate(epoch_groups)) do (g, wanted)
        wanted_set = Set{Int}(wanted)
        rows_for = Dict{Int,Vector{Int}}()
        for (gi, (k, r)) in enumerate(plan)
            gi in wanted_set && push!(get!(rows_for, k, Int[]), r)
        end

        keep = AbstractObs[]
        ep = Float64[]
        for (k, o) in enumerate(sys.observations)
            if _nrows(o) == 0
                push!(keep, o)
                continue
            end
            rows = get(rows_for, k, Int[])
            isempty(rows) && continue
            push!(keep, _subset(o, rows; what="epoch-group model generation"))
            append!(ep, (_row_epoch(o, r) for r in rows))
        end
        push!(epochs_out, ep)
        return _rebuild_system(sys; observations=keep,
            name=Symbol(sys.name, name_suffix_func(g)))
    end

    return systems, epochs_out
end

"""
    generate_system_per_epoch(system) -> (systems, epochs)

One copy of `system` per data row, each containing only that row.
"""
function generate_system_per_epoch(sys::System)
    total = length(_row_plan(sys))
    total == 0 && return System[], Float64[]
    systems, epoch_vectors = generate_systems_with_epoch_groups(
        sys, [[i] for i in 1:total], i -> "_epoch_$i")
    return systems, reduce(vcat, epoch_vectors; init=Float64[])
end

"""
    generate_cumulative_system_per_epoch(system) -> (systems, epochs)

One copy of `system` per data row, each containing rows `1` through `i` —
for watching a posterior tighten as data accumulate.
"""
function generate_cumulative_system_per_epoch(sys::System)
    total = length(_row_plan(sys))
    total == 0 && return System[], Vector{Float64}[]
    return generate_systems_with_epoch_groups(
        sys, [collect(1:i) for i in 1:total], i -> "_cumulative_epoch_$i")
end

# ---------------------------------------------------
# Pointwise likelihoods
# ---------------------------------------------------

# One column of the pointwise-likelihood matrix: an observation restricted to
# a single data row (or a whole table-less observation), plus the namespace
# its variables live in and the epoch to label the column with.
struct PointwiseTerm{TO}
    obs::TO
    key::Symbol
    epoch::Float64
end

function _pointwise_terms(sys::System)
    terms = PointwiseTerm[]
    for o in sys.observations
        _isprior(o) && continue
        key = normalizename(likelihoodname(o))
        n = _nrows(o)
        if n == 0
            # No table to split: the whole term is one column.
            push!(terms, PointwiseTerm(o, key, NaN))
        else
            for r in 1:n
                push!(terms, PointwiseTerm(
                    _subset(o, [r]; what="pointwise likelihood evaluation"),
                    key, _row_epoch(o, r)))
            end
        end
    end
    return terms
end

"""
    pointwise_like(model, chain; verbosity=1) -> (likelihood_mat, epochs)

The log-likelihood of each individual data point under each posterior sample:
an `N_sample × N_data` matrix, plus the epoch labelling each column.

Columns are in model order — observation by observation, then table row by
table row. Prior-shaped terms (`lhs ~ dist` lines, `LL +=` lines, the
`UnitLengthPrior` behind a `UniformCircular`) are **not** columns: they carry
no data, so including them would both invent data points and add the same
constant to every column. Observations that carry data but have no epoch
table contribute one column each, labelled `NaN`.

Feed the result to [ParetoSmooth.jl](https://turinglang.org/ParetoSmooth.jl)
for leave-one-out cross-validation:

```julia
likelihood_mat, epochs = Octofitter.pointwise_like(model, chain)

using ParetoSmooth
result = psis_loo(collect(likelihood_mat'), chain_index=ones(Int, size(chain, 1)))
```

# Implementation
v8 built one whole `System` per data point and compiled a fresh
`RuntimeGeneratedFunction` likelihood for each — `N_data` model builds before
the first number came out. v9 needs one: the per-row observations are
evaluated directly against a single trajectory solved over the union of
their epochs, exactly as `make_ln_like` does it, so the cost is one solve per
sample rather than one model compile per data point.

The consequence worth knowing about is that the columns sum to the model's
total log-likelihood *minus* its prior-shaped terms. That is the quantity
PSIS-LOO wants; v8's columns did not have that property.
"""
function pointwise_like(model, chain; verbosity::Int=1)
    sys = model.system

    verbosity >= 1 && @info "Resolving chain"
    sample_nts = mcmcchain2result(model, chain)

    verbosity >= 1 && @info "Planning pointwise terms"
    terms = _pointwise_terms(sys)
    isempty(terms) && return zeros(length(sample_nts), 0), Float64[]

    # The epoch union these per-row observations need, and each one's map from
    # its (single) table row into that union.
    all_ep = Float64[]
    for t in terms
        append!(all_ep, epochs(t.obs))
    end
    unique_ep = sort!(unique!(all_ep))
    index_of = Dict(e => k for (k, e) in enumerate(unique_ep))
    maps = [Int[index_of[e] for e in epochs(t.obs)] for t in terms]

    kw = _solvekw(sys)
    LL_out = zeros(length(sample_nts), length(terms))

    verbosity >= 1 && @info "Calculating pointwise likelihoods" samples=length(sample_nts) points=length(terms)
    Threads.@threads for i_sample in eachindex(sample_nts)
        θ = sample_nts[i_sample]
        posys = _build_po_system(sys, θ)
        traj = isempty(unique_ep) ? nothing : orbitsolve(posys, unique_ep; kw...)
        obsns = hasproperty(θ, :observations) ? θ.observations : (;)
        for (j, t) in enumerate(terms)
            θ_obs = hasproperty(obsns, t.key) ? getproperty(obsns, t.key) : (;)
            LL_out[i_sample, j] = ln_like(t.obs, ObsContext(θ, θ_obs, posys, traj, maps[j]))
        end
    end

    return LL_out, Float64[t.epoch for t in terms]
end
