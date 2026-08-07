# ---------------------------------------------------
# Code generation
#
# Everything the sampler touches is unrolled at model-build time into a
# RuntimeGeneratedFunction: the parameter vector → nested NamedTuple map, the
# prior density, the prior sampler, and the log-likelihood. Unrolling is what
# lets user `@variables` expressions inline into one tight numeric function
# instead of looping over a heterogeneous container with runtime dispatch.
#
# Parameter order — used identically by every generated function here:
#     system priors → node priors (node order) → observation priors (obs order)
#
# Evaluation order (§8.4): bodies look up and outward; the system block looks
# down and inward, after the fact; siblings never see each other.
#     1. system block, non-deferred lines        → `system.*`
#     2. every node's block, in declaration order (sees `system.*`)
#     3. system block, deferred lines            (sees `system.*` and `b.*`)
#     4. observation blocks                      (sees `system.*`)
#     5. build the PlanetOrbits.System, solve once, evaluate likelihoods
#
# Note what is *not* here: no `M` injection into user blocks. Both of its uses
# in v1 — `a = cbrt(M*P^2)` and the `θ`-to-`tp` conversion — moved into
# `PlanetOrbits.Orbit`'s constructor, which already holds the row's
# gravitating mass because `about=` carries the Body values.
# ---------------------------------------------------

# --- one variable namespace ---------------------------------------------------

# Emits the prior-unpacking assignments and the chain of derived-variable
# steps for a single namespace. `i` is the running parameter index; `bindings`
# are extra `name = value` statements made visible to every derived
# expression (e.g. `system = …` inside a body block).
function _emit_namespace(tag::Symbol, priors::Priors, derived::Derived,
                         i::Base.RefValue{Int}, bindings::Vector;
                         derived_keys=collect(keys(derived.variables)),
                         preexisting::Vector{Symbol}=Symbol[])
    prior_exprs = Expr[]
    for key in keys(priors.priors)
        d = priors.priors[key]
        if length(d) > 1
            parts = Expr[]
            for _ in 1:length(d)
                i[] += 1
                push!(parts, :(arr[$(i[])]))
            end
            push!(prior_exprs, :($key = ($(parts...),)))
        else
            i[] += 1
            push!(prior_exprs, :($key = arr[$(i[])]))
        end
    end

    # `$x` interpolations from the `@variables` block, re-bound as locals so
    # they are compile-time constants rather than global lookups.
    #
    # The value goes in as a `QuoteNode`, not spliced directly: splicing a
    # `Symbol` into an expression produces an *identifier*, so a captured
    # `$key` where `key = :mass_b` used to compile to `k = mass_b` and fail
    # with an `UndefVarError` naming a variable the user never wrote.
    # QuoteNode is correct for every value type, including numbers and arrays.
    captured = Expr[Expr(:(=), derived.captured_names[j],
                         QuoteNode(derived.captured_vals[j]))
                    for j in eachindex(derived.captured_names)]
    prior_keys = collect(keys(priors.priors))

    steps = Expr[]
    for (j, key) in enumerate(derived_keys)
        expr = derived.variables[key]
        # Everything already in this namespace is in scope, in order.
        avail = unique(vcat(preexisting, prior_keys, derived_keys[1:j-1]))
        prev = Symbol(tag, j - 1)
        push!(steps, :(
            $(Symbol(tag, j)) = let _prev = $prev
                (; _prev..., $key = let
                    $(captured...)
                    $(bindings...)
                    $((:($k = _prev.$k) for k in avail)...)
                    $expr
                end)
            end
        ))
    end
    return prior_exprs, steps
end

# `tag0 = (; priors...); steps...; tagN` as a single expression.
function _namespace_block(tag::Symbol, prior_exprs, steps)
    return :(begin
        $(Symbol(tag, 0)) = (; $(prior_exprs...))
        $(steps...)
        $(Symbol(tag, length(steps)))
    end)
end

_nstag(name::Symbol) = Symbol("_ns_", name, "_")

# --- arr2nt -------------------------------------------------------------------

"""
    make_arr2nt(system)

Build the function mapping a flat parameter vector (drawn from priors,
sampled from the posterior, or carrying Duals) to the nested NamedTuple the
likelihoods read:

    (; <system variables>...,
       bodies = (; A = (; mass = …), b = (; mass = …, a = …, e = …)),
       observations = (; GaiaDR4 = (; jitter = …)))

All derived variables are evaluated in the process, in the order documented
at the top of this file.
"""
make_arr2nt(sys::System) = @RuntimeGeneratedFunction(arr2nt_expr(sys))

# Split out so the generated code can be inspected:
#   Base.remove_linenums!(Octofitter.arr2nt_expr(sys)) |> println
function arr2nt_expr(sys::System)
    i = Ref(0)

    # 1. system block, non-deferred lines
    nondef = [k for k in keys(sys.derived.variables) if !(k in sys.deferred)]
    sys_priors, sys_steps = _emit_namespace(:_sys, sys.priors, sys.derived, i, Any[];
        derived_keys=nondef)
    sysvar = Symbol(:_sys, length(sys_steps))

    # 2. every node's block
    node_blocks = Expr[]
    for n in sys.nodes
        tag = _nstag(n.name)
        p, s = _emit_namespace(tag, n.priors, n.derived, i, Any[:(system = $sysvar)])
        push!(node_blocks, :($(n.name) = $(_namespace_block(tag, p, s))))
    end

    # 3. system block, deferred lines — node namespaces bound by name.
    # No priors are consumed here (they were all taken in step 1), so the
    # index counter is a throwaway.
    node_bindings = Any[:($(n.name) = _bodies.$(n.name)) for n in sys.nodes]
    _, def_steps = _emit_namespace(:_sysd, Priors(), sys.derived, Ref(0), node_bindings;
        derived_keys=sys.deferred,
        preexisting=Symbol[collect(keys(sys.priors.priors))..., nondef...])
    deferred_block = isempty(def_steps) ? Expr[] : Expr[:(_sysd0 = $sysvar), def_steps...]
    sysfull = isempty(def_steps) ? sysvar : Symbol(:_sysd, length(def_steps))

    # 4. observation blocks
    obs_blocks = Expr[]
    for o in sys.observations
        hasproperty(o, :priors) || continue
        tag = _nstag(Symbol(normalizename(likelihoodname(o))))
        p, s = _emit_namespace(tag, o.priors, o.derived, i, Any[:(system = $sysfull)])
        (isempty(p) && isempty(s)) && continue
        push!(obs_blocks, :($(normalizename(likelihoodname(o))) =
            $(_namespace_block(tag, p, s))))
    end

    n = i[]
    return :(function (arr)
        @boundscheck if length(arr) != $n
            error("Expected exactly $($n) elements in the parameter vector, got ", length(arr))
        end
        _sys0 = (; $(sys_priors...))
        $(sys_steps...)
        _bodies = (; $(node_blocks...))
        $(deferred_block...)
        return (; $sysfull..., bodies=_bodies, observations=(; $(obs_blocks...)))
    end)
end

# --- priors -------------------------------------------------------------------

"""
    _list_priors(system)

Flat vector of prior distributions, in the same order as [`make_arr2nt`](@ref).
"""
function _list_priors(sys::System)
    out = Any[]
    append!(out, values(sys.priors.priors))
    for n in sys.nodes
        append!(out, values(n.priors.priors))
    end
    for o in sys.observations
        hasproperty(o, :priors) || continue
        append!(out, values(o.priors.priors))
    end
    return map(identity, out)
end

# Every namespace's priors, in parameter order, as (name, distribution) pairs.
function _prior_pairs(sys::System)
    out = Pair{Symbol,Any}[]
    for (k, d) in sys.priors.priors
        push!(out, k => d)
    end
    for n in sys.nodes, (k, d) in n.priors.priors
        push!(out, k => d)
    end
    for o in sys.observations
        hasproperty(o, :priors) || continue
        for (k, d) in o.priors.priors
            push!(out, k => d)
        end
    end
    return out
end

"""
    make_ln_prior_transformed(system)

Log prior density with the Bijectors change-of-variables correction applied,
as an unrolled function of the flat parameter vector in the natural domain.
"""
function make_ln_prior_transformed(sys::System)
    i = 0
    evals = Expr[]
    for (_, d) in _prior_pairs(sys)
        if length(d) > 1
            samples = Expr[]
            for _ in 1:length(d)
                i += 1
                push!(samples, :(arr[$i]))
            end
            sample_ex = :(SVector($(samples...),))
        else
            i += 1
            sample_ex = :(arr[$i])
        end
        push!(evals, :(
            p = $logpdf_with_trans($d, $sample_ex, sampled);
            # Out-of-bounds values here come from finite precision in the
            # unconstrained space, not from a bad model; clamp rather than
            # propagate a NaN into the sampler.
            if !isfinite(p) && $(eltype(d) <: AbstractFloat)
                if sign(p) > 1
                    return prevfloat(typemax($(eltype(d))))
                else
                    return nextfloat(typemin($(eltype(d))))
                end
            end;
            lp += p
        ))
    end
    isempty(evals) && error("Model includes no free variables")
    return @RuntimeGeneratedFunction(:(function (arr, sampled)
        @boundscheck if length(arr) != $i
            error("Expected exactly $($i) elements in the parameter vector (got $(length(arr)))")
        end
        lp = zero(first(arr))
        @inbounds begin
            $(evals...)
        end
        return lp
    end))
end

"""
    make_prior_sampler(system)

Draw one IID sample from the model's priors, as a flat tuple in parameter
order.
"""
function make_prior_sampler(sys::System)
    exprs = Expr[]
    for (_, d) in _prior_pairs(sys)
        push!(exprs, :(sample = $rand(rng, $d)))
        for j in 1:length(d)
            push!(exprs, :(prior_samples = (prior_samples..., sample[$j])))
        end
    end
    return @RuntimeGeneratedFunction(:(function (rng)
        prior_samples = ()
        @inbounds begin
            $(exprs...)
        end
        return prior_samples
    end))
end

"""
    make_Bijector_invlinkvec(priors_vec)

Unrolled `Bijectors.invlink` over the prior list — the transformed →
natural map on the sampler's hot path.
"""
function make_Bijector_invlinkvec(priors_vec)
    i = 0
    transforms = Expr[]
    for d in priors_vec
        if length(d) > 1
            access = Expr[]
            for _ in 1:length(d)
                i += 1
                push!(access, :(arr[$i]))
            end
            push!(transforms, :($(Bijectors.invlink)($d, ($(access...),))...))
        else
            i += 1
            push!(transforms, :($(Bijectors.invlink)($d, arr[$i])))
        end
    end
    return @RuntimeGeneratedFunction(:(function (arr)
        @boundscheck if length(arr) != $i
            error("Expected exactly $($i) elements in the parameter vector (got $(length(arr)))")
        end
        @inbounds theta_out = tuple($(transforms...))
        return theta_out
    end))
end

"""
    drawfrompriors(system)
    drawfrompriors(system; overrides, rng=Random.default_rng())

One draw from the model's priors, already expanded into the nested NamedTuple
structure — system variables at the top level, then `bodies`, then
`observations`.

`overrides` pins chosen **free** (`~`) variables to values you supply, in the
same nested shape, and is the right way to build a parameter set for
[`generate_from_params`](@ref):

    θ = drawfrompriors(sys; overrides=(;
        plx = 24.5,
        bodies = (; b = (; mass = 85mjup, a = 45.0, e = 0.15)),
    ))

The values are written into the *flat* parameter vector before it is expanded,
so every derived (`=`) variable is recomputed from them. That is what
`merge`ing into a template cannot do: `merge` replaces an entry in the
already-expanded structure, leaving any derived variable that was computed from
it stale — and `generate_from_params` reads the derived elements, not the
sampled ones they came from. Overriding a derived variable is an error naming
the model's free variables, rather than a silent no-op.
"""
drawfrompriors(sys::System; overrides=(;), rng::Random.AbstractRNG=Random.default_rng()) =
    make_arr2nt(sys)(_apply_overrides!(collect(make_prior_sampler(sys)(rng)), sys, overrides))
export drawfrompriors

# --- the likelihood -----------------------------------------------------------

# v1 threaded the per-planet Kepler solve from here; v2 batches it inside
# PlanetOrbits (SIMD across epochs, and across the Dual chunk), so Octofitter
# owns no inner threading at all. The flag survives only as the guard that
# stops the initializer's outer loops from nesting threads inside a threaded
# solve — it is now always false, and the outer loops are always free to
# thread.
const _kepsolve_use_threads = Ref(false)

"""
    _obsctx(sys, θ, θ_obs, posys, traj, map)

`ObsContext` for one observation of a **built** `System`, carrying that
system's solve configuration.

Every path that evaluates a likelihood outside the generated function —
cross-validation, `generate_from_params`, the plotting API — goes through this
rather than calling `ObsContext` directly, so that a setting like
`einstein_rv=:off` cannot be honoured by the sampler and quietly ignored by
everything that reads the result back.
"""
@inline _obsctx(sys::System, θ, θ_obs, posys, traj, map) =
    ObsContext(θ, θ_obs, posys, traj, map, Bumper.default_buffer(), sys.einstein_rv)

"""
    epoch_plan(system)

The sorted, deduplicated union of every observation's epochs, plus a per
observation map from table row to column of the solved `Trajectory`.

Deduplication is what makes the single-solve design pay off: two likelihoods
sharing an epoch — an RV point and an astrometric point on the same night —
cost one solve, not two, and the plan is the "union of states needed" the
reference declarations give for free.
"""
function epoch_plan(sys::System)
    all_ep = Float64[]
    for o in sys.observations
        append!(all_ep, epochs(o))
    end
    unique_ep = sort!(unique(all_ep))
    index_of = Dict(e => k for (k, e) in enumerate(unique_ep))
    maps = Dict{Any,Vector{Int}}()
    for o in sys.observations
        maps[o] = Int[index_of[e] for e in epochs(o)]
    end
    return unique_ep, maps
end

"""
    BumpAlloc(buf)

Column allocator handed to `PlanetOrbits.Trajectory`: `alloc(T, dims...)`
carves the column out of a Bumper buffer.

This is a named callable rather than the closure it looks like it should be.
A `@generated` function body may not contain a closure, and
`RuntimeGeneratedFunction`s are generated functions: written inline, the
allocator stayed uninferred, `Trajectory`'s column type parameters came out
abstract, and every `traj[k]` then boxed a `TrajectorySolution` — 220 epochs
× ~64 bytes, i.e. 14 kB and ~250 allocations per likelihood evaluation, on a
path whose whole point is to allocate nothing.
"""
struct BumpAlloc{B}
    buf::B
end
@inline (a::BumpAlloc)(::Type{S}, dims::Vararg{Integer,N}) where {S,N} =
    Bumper.alloc!(a.buf, S, dims...)

# ---------------------------------------------------
# Per-sample scratch arena
#
# Bumper's default buffer is a `SlabBuffer{1 MB}`, and its checkpoint restore
# frees **every slab past the first** when a `@no_escape` block exits. Slab #1
# is therefore the only storage that survives between evaluations: a model
# whose per-sample scratch exceeds it hands the whole excess back to the
# allocator every sample and pays to fault it in again on the next one.
#
# That is invisible on small models and brutal on large ones. Measured on a
# 2-body RV model at ForwardDiff `Dual{10}` (∇ℓπ, µs), against the same
# workload run on a buffer whose single slab covers the working set:
#
#   epochs |  500   1000   1500   2000   3000   4000   8000
#   default|  393    797   2118   3915   5583   7557  15083
#   sized  |  394    793   1206   1608   2421   3239   6723
#
# Linear up to ~1000 epochs, then 1.8-2.4x and staying there. The sized column
# tracks a plain heap-allocated trajectory to within 4%, so this recovers the
# whole gap rather than trading it somewhere else.
#
# So: size a slab to the model at build time and keep it in task-local
# storage, exactly where Bumper keeps its own default buffer (which is what
# makes the default thread-safe — one arena per task, not one per process).
# Models that fit in the default 1 MB keep using it, so nothing changes for
# them and no second arena is held.
# ---------------------------------------------------

# Task-local, keyed by the slab size, which is a build-time constant baked into
# the generated likelihood. `Val{N}` rather than a tuple or an interpolated
# symbol because task-local storage is an `IdDict`: the key has to be `===`
# across calls, and types are interned while freshly-built tuples are not.
@inline function _scratch_buffer(::Val{SlabSize}) where {SlabSize}
    return get!(() -> Bumper.SlabBuffer{SlabSize}(),
        task_local_storage(), Val{SlabSize})::Bumper.SlabBuffer{SlabSize}
end
@inline _scratch_buffer(::Val{nothing}) = Bumper.default_buffer()

const _DEFAULT_SLAB = 1 << 20
# A model wanting more than this per sample gets extra slabs and the cost above.
# The cap exists because the arena is held for the life of the task: one 256 MB
# buffer per sampling chain is already more than any real model needs, and
# beyond it correctness (extra slabs) is preferable to the memory.
const _MAX_SLAB = 1 << 28

"""
    _slab_size(sys, unique_ep, θ_example=nothing) -> Int or nothing

Slab size for this model's scratch arena, or `nothing` to keep Bumper's
default buffer.

Sized for the widest single ForwardDiff chunk the model can ask for — one
`Dual` partial per free parameter — since that is what `LogDensityModel`
defaults to and it is the largest element type on the ordinary gradient path.
Nested `Dual`s (Hessians) exceed it and fall back to extra slabs: slower, but
correct, which is the right way round for a case nobody's inner loop is.

Pass `θ_example` (a nested NamedTuple, as `make_ln_like` already receives) and
`build` (the system constructor it already compiled) to avoid redoing either
just to price the trajectory.
"""
function _slab_size(sys::System, unique_ep, θ_example=nothing, build=nothing)
    D = length(_list_priors(sys))
    T = ForwardDiff.Dual{Nothing,Float64,D}
    posys = try
        _example_posys(sys, θ_example, build)
    catch err
        err isa InterruptException && rethrow()
        return nothing   # can't price it; the default buffer is always correct
    end
    isnothing(posys) && return nothing
    bytes = PlanetOrbits.trajectory_storage(T, posys, unique_ep)
    # Plus the likelihoods' own `@alloc` scratch, which shares this arena:
    # no current one takes more than three columns of its own table.
    bytes += 3 * sizeof(T) * sum(o -> length(epochs(o)), sys.observations; init=0)
    # Headroom for the allocator's alignment padding, then round up: a slab
    # that is a hair too small costs the whole penalty above.
    want = nextpow(2, bytes + bytes ÷ 4 + 65536)
    want <= _DEFAULT_SLAB && return nothing
    return min(want, _MAX_SLAB)
end

# One PlanetOrbits.System, purely to price its trajectory. Returns `nothing`
# if it cannot be built — sizing is an optimization, so it must never be the
# thing that fails model construction.
#
# The RNG is local and fixed-seeded, not `Random.default_rng()`: model
# construction must not consume the global stream, or a script that seeds
# before building its model gets different samples depending on whether the
# arena happened to need sizing.
function _example_posys(sys::System, θ_example=nothing, build=nothing)
    θ = isnothing(θ_example) ?
        make_arr2nt(sys)(collect(make_prior_sampler(sys)(Random.Xoshiro(0)))) :
        θ_example
    f = isnothing(build) ? @RuntimeGeneratedFunction(:(function (θ)
        T = $_system_number_type(θ)
        $(_build_system_expr(sys)...)
    end)) : build
    return f(θ)
end

# Expression building the PlanetOrbits.Body / Orbit / System for one sample.
function _build_system_expr(sys::System)
    stmts = Expr[]
    bodyvar = Dict{Symbol,Symbol}()
    for n in bodynodes(sys)
        v = Symbol("_pob_", n.name)
        bodyvar[n.name] = v
        names = varnames(n)
        massex = :mass in names ? :(θ.bodies.$(n.name).mass) : :(zero(T))
        fluxvars = _flux_vars(n)
        fluxex = isempty(fluxvars) ? :((;)) :
                 Expr(:tuple, Expr(:parameters,
                     (Expr(:kw, _flux_band(fv), :(θ.bodies.$(n.name).$fv)) for fv in fluxvars)...))
        push!(stmts, :($v = PlanetOrbits.Body(mass=$massex, flux=$fluxex,
            name=$(QuoteNode(n.name)))))
    end

    _spec(members) = length(members) == 1 ? bodyvar[members[1]] :
                     Expr(:tuple, (bodyvar[m] for m in members)...)

    orbvars = Symbol[]
    for (k, (owner, ext, int)) in enumerate(sys.rows)
        node = _node(sys, owner)
        v = Symbol("_poo_", k)
        push!(orbvars, v)
        kws = [Expr(:kw, e, :(θ.bodies.$(owner).$e)) for e in _element_vars(node)]
        push!(stmts, :($v = PlanetOrbits.Orbit($(_spec(ext));
            about=$(_spec(int)), $(kws...))))
    end

    framekws = [Expr(:kw, f, :(θ.$f)) for f in sys.framevars]
    push!(stmts, :(return PlanetOrbits.System(
        $(Expr(:tuple, (bodyvar[nm] for nm in sys.bodynames)...)),
        $(Expr(:tuple, orbvars...));
        $(framekws...))))
    return stmts
end

"""
    make_ln_like(system, θ_example; include_priors=true)

Build the model's log-likelihood: one `PlanetOrbits.System` per sample, one
`orbitsolve!` over the whole epoch union into Bumper-allocated storage, then
every observation evaluated against index views into that one trajectory.

A likelihood asks for `raoff(sol, b, A)` and gets the exact relative offset for
whatever hierarchy the model describes, under whichever propagator is
configured — there is nothing for it to reconstruct by hand.

`include_priors=false` builds the same function over the **data** terms only,
skipping every `_isprior` observation (`UserLikelihood` from a `~` line in a
`@variables` block, `UnitLengthPrior`, the dynamical stability priors) and the
`@variables`-emitted prior terms. That is what belongs in a chain's `loglike`
column; the terms it skips reshape the prior and are reported as `logprior`.
The full version — priors included — is what the sampler's log posterior needs,
and is unchanged.
"""
function make_ln_like(sys::System, θ_example=nothing; include_priors::Bool=true)
    unique_ep, maps = epoch_plan(sys)

    build = @RuntimeGeneratedFunction(:(function (θ)
        T = $_system_number_type(θ)
        $(_build_system_expr(sys)...)
    end))

    # Observation calls, in order: real observations first, then the
    # prior-shaped terms `@variables` emitted (which read their owner's
    # namespace rather than one of their own).
    ein = sys.einstein_rv
    calls = Expr[]
    j = 0
    for (k, o) in enumerate(sys.observations)
        !include_priors && _isprior(o) && continue
        nm = normalizename(likelihoodname(o))
        θ_obs = :(hasproperty(θ.observations, $(QuoteNode(nm))) ? θ.observations.$nm : (;))
        push!(calls, :($(Symbol(:ll, j + 1)) = $(Symbol(:ll, j)) + ln_like(
            system.observations[$k],
            ObsContext(θ, $θ_obs, posys, traj, $(maps[o]), buf, $ein))))
        j += 1
    end
    if include_priors
        for (k, (owner, _)) in enumerate(sys.priorterms)
            θ_obs = owner === :system ? :(θ) : :(θ.bodies.$owner)
            push!(calls, :($(Symbol(:ll, j + 1)) = $(Symbol(:ll, j)) + ln_like(
                system.priorterms[$k][2],
                ObsContext(θ, $θ_obs, posys, traj, $(Int[]), buf, $ein))))
            j += 1
        end
    end

    method = sys.method
    og = sys.observing_geometry
    blt = sys.barycentric_lighttime
    slab = _slab_size(sys, unique_ep, θ_example, build)

    evaluate = @RuntimeGeneratedFunction(:(function (system, θ, posys)
        T = $_system_number_type(θ)
        # Sized to this model at build time; see `_slab_size`. The likelihoods
        # get it through the context so their own `@alloc` scratch shares one
        # arena rather than falling back to the default buffer's 1 MB slab.
        buf = $_scratch_buffer($(Val(slab)))
        return @no_escape buf begin
            # One trajectory for the whole model. Its columns come from the
            # bump allocator through PlanetOrbits' caller-allocated
            # constructor, so the column *set* stays PlanetOrbits' business —
            # it has already grown once (the observing-geometry pass) and
            # would otherwise have to be tracked in two places.
            traj = PlanetOrbits.Trajectory(BumpAlloc(buf), T, posys, $unique_ep)
            PlanetOrbits.orbitsolve!(traj, posys; method=$method,
                observing_geometry=$og, barycentric_lighttime=$blt)
            ll0 = zero(T)
            $(calls...)
            $(Symbol(:ll, j))
        end
    end))

    return GeneratedLnLike(build, evaluate)
end

"""
    GeneratedLnLike

The compiled log-likelihood: `build` maps parameters to a
`PlanetOrbits.System`, `evaluate` solves it once and sums every
observation's contribution.

Two layers, deliberately. Orbit construction can fail on a proposal the
priors admit but the elements do not (`e == 1` exactly, a non-positive `a`
out of a derived expression), and Bumper's `@no_escape` is not
exception-safe — a throw inside it would leak the bump buffer. Catching
around construction only, *outside* the block, is also what keeps `posys`
concretely typed: the catch branch returns rather than falling through.

Both halves are kept as fields rather than closed over so a model can be
inspected: `Base.remove_linenums!(m.evaluate.body)` for the generated code.
"""
struct GeneratedLnLike{B,E}
    build::B
    evaluate::E
end

function (f::GeneratedLnLike)(system::System, θ)
    T = _system_number_type(θ)
    posys = try
        f.build(θ)
    catch err
        err isa InterruptException && rethrow()
        @warn "Failed to construct the orbital system from these parameters" exception = (err, catch_backtrace()) maxlog = 1
        return convert(T, -Inf)
    end
    return f.evaluate(system, θ, posys)::T
end

# Helper for determining the number type to use in a likelihood function.
# Walks the nested NamedTuple and promotes every numeric leaf.
_system_number_type(T::NamedTuple) = _system_number_type(typeof(T))
@generated function _system_number_type(::Type{NamedTuple{Keys,Vals}}) where {Keys,Vals}
    T = Bool
    for V in fieldtypes(Vals)
        if V <: Number
            T = promote_type(T, V)
        elseif V <: NamedTuple
            # Only recurse into concretely-typed nested tuples: values
            # reconstructed outside the sampler hot path can carry an
            # abstractly-typed field, whose contents cannot be inspected here.
            isconcretetype(V) && (T = promote_type(T, _system_number_type(V)))
        elseif V <: Tuple
            T = promote_type(T, eltype(V))
        end
    end
    return float(T)
end
