# ---------------------------------------------------
# Dynamical configuration priors
#
# `OrbitOrderPrior`, `NonCrossingPrior`, `LimitClosestApproachAUPrior` and
# `HillStabilityPrior`. All carry no data — they only reshape the prior — so
# they are `_isprior = true` terms that read `ctx.system` and nothing else.
# All are written over explicit body refs and fixed-width rows, and allocate
# nothing per evaluation.
# ---------------------------------------------------

# ---------------------------------------------------
# Shared row plumbing
#
# The dynamics live in the system's hierarchy rows, and a body names a row
# only indirectly — through the row that places it. These three helpers are
# that indirection, plus a non-allocating ordered walk over the rows.
# ---------------------------------------------------

"""
    _placing_row(posys, j) -> Int

Index of the hierarchy row that places body `j`, or `0` if none does (the
hierarchy root).

"The row that places a body" is the most *local* relationship in which it is
an exterior member — the one binding the fewest bodies. In a Jacobi or
astrocentric chain each body is the exterior of exactly one row and the choice
is not a choice; in a 2+2 quadruple `Ba` is *also* a member of the wide row's
set exterior, and it is the tight row that describes `Ba`'s own motion. Ties
(two rows binding equally many bodies, both with `j` on the exterior side)
take the first, which is the innermost by row order.

Runs over masks rather than index tuples, so it neither allocates nor branches
on body count; `NR` and `NB` are type parameters, so it is a handful of
boolean loads.
"""
@inline function _placing_row(posys::PlanetOrbits.System{NB,NR}, j::Int) where {NB,NR}
    best = 0
    bestn = NB + 1
    @inbounds for k in 1:NR
        s = posys.specs[k]
        s.ext[j] || continue
        n = sum(s.int) + sum(s.ext)
        if n < bestn
            bestn = n
            best = k
        end
    end
    return best
end

@inline _allrows(::PlanetOrbits.System{NB,NR}) where {NB,NR} = ntuple(identity, Val(NR))

# The solved system's own scalar type. Everything these priors read comes off
# the rows and the masses, so the bookkeeping is done in *that* type rather
# than in `_system_number_type(ctx.θ_system)`: the two are usually equal, but a
# model whose only sampled parameters live in an observation's namespace has
# Duals in θ and a plain `Float64` system, and mixing the two would infer a
# `Union` for the accumulator.
@inline _rowtype(::PlanetOrbits.System{NB,NR,Ts}) where {NB,NR,Ts} = Ts

"""
    _selected_rows(ctx, specs::Tuple) -> NTuple{N,Int}

The hierarchy rows a dynamical prior applies to: the rows placing the named
bodies, or — given no names — every row in the system, which is what v8 always
did (it looped over `orbits`).

Restricting matters in v9 in a way it could not in v8: a system's rows now
include things like the wide orbit of a 2+2 quadruple, and "does the inner
pair's orbit cross the wide binary's" is not a question the apsidal test
answers meaningfully.
"""
@inline _selected_rows(ctx::ObsContext, ::Tuple{}) = _allrows(ctx.system)
@inline function _selected_rows(ctx::ObsContext, specs::Tuple)
    refs = resolverefs(ctx, specs)
    return map(refs, specs) do r, s
        k = _placing_row(ctx.system, r.idx)
        k == 0 && _err_unplaced(s)
        return k
    end
end

@noinline _err_unplaced(spec) = error(
    "dynamical prior: body $(_refstr(spec)) is not placed by any orbit — it is " *
    "the hierarchy root, so it has no orbit of its own. Name the companions " *
    "whose orbits the prior constrains instead.")

"""
    _next_by_sma(posys, ks, a_prev, k_prev, has_prev) -> (k, a)

Next row of `ks` in ascending semi-major axis after `(a_prev, k_prev)`, or
`(0, …)` when there is none. Ordering is lexicographic in `(a, k)`, so equal
semi-major axes are broken by row index and the walk is a total order (and
therefore terminates) even for a system with degenerate rows.

This is the non-allocating replacement for v8's
`sort(collect(orbits), by=semimajoraxis)`, which built a `Vector` on **every**
likelihood evaluation and carried a `# TODO: would be nice to make this
non-allocating` for it. Rows are fixed-width and few, so the O(n²) walk is
cheaper than any sort that touches the heap, and it stays type-stable under
ForwardDiff (a permutation vector would not be — the `by` values are Duals).
"""
@inline function _next_by_sma(posys::PlanetOrbits.System, ks::Tuple,
                              a_prev, k_prev::Int, has_prev::Bool)
    k_best = 0
    a_best = a_prev
    @inbounds for k in ks
        a = posys.rows[k].a
        if has_prev && !(a > a_prev || (a == a_prev && k > k_prev))
            continue
        end
        if k_best == 0 || a < a_best || (a == a_best && k < k_best)
            k_best = k
            a_best = a
        end
    end
    return k_best, a_best
end

# v1's `periapsis`/`apoapsis`, which PlanetOrbits v2 does not carry: they are
# two lines of algebra over the row, and their hyperbolic branch (an "apoapsis"
# for an unbound orbit) is a v1 convention rather than a physical quantity, so
# it belongs with the priors that consume it rather than in the orbit library.
@inline _row_periapsis(row) = row.e < 1 ? row.a * (1 - row.e) : -row.a * (row.e - 1)
@inline _row_apoapsis(row) = row.e < 1 ? row.a * (1 + row.e) : -row.a * (1 + row.e)

# Total mass of a masked subset of bodies. (PlanetOrbits builds row masses the
# same way; going through its helper keeps the two definitions from drifting.)
@inline _submass(posys::PlanetOrbits.System, mask) =
    PlanetOrbits._maskmass(posys.masses, mask)

# ---------------------------------------------------
# Ordering
# ---------------------------------------------------

"""
    OrbitOrderPrior(b, c, d, …; name="OrbitOrderPrior")

Assert that the named bodies' orbits are ordered in semi-major axis:
`a_b < a_c < a_d < …`. Anything else gets `-Inf`.

    System(…, observations=[astrom_b, astrom_c, OrbitOrderPrior(b, c)], …)

This breaks the labelling degeneracy of a multi-companion fit. Without it a
sampler is free to swap which companion is called `b` and which `c`, and the
posterior is a superposition of every permutation — technically correct, and
useless to summarize.

Each argument is a `Body` model node or a `Symbol` naming one; the orbit meant
is the row that *places* that body (see `_placing_row`), so a body may be
named whatever hierarchy convention the model uses. v8's `PlanetOrderPrior`
took `Planet`s and looked their orbits up positionally with `findfirst` over
`keys(θ_system.planets)`.

Carries no data, so it cannot be held out or simulated: it is a prior term
(`_isprior = true`), not an observation.
"""
struct OrbitOrderPrior{TB<:Tuple} <: AbstractObs
    bodies::TB
    name::String
end

function OrbitOrderPrior(bodies...; name="OrbitOrderPrior")
    specs = map(refspec, Tuple(bodies))
    all(s -> s isa BodyRefSpec, specs) || error(
        "OrbitOrderPrior orders the orbits of individual bodies, so each argument " *
        "must be a `Body` model node or a `Symbol` naming one.")
    length(specs) >= 2 || error(
        "OrbitOrderPrior needs at least two bodies to order; got $(length(specs)).")
    return OrbitOrderPrior{typeof(specs)}(specs, String(name))
end

export OrbitOrderPrior

# v1 name, kept so existing scripts keep parsing. The v1 constructor took
# `Planet` objects, which no longer exist; the replacement takes `Body` nodes.
"""
    PlanetOrderPrior

Deprecated alias for [`OrbitOrderPrior`](@ref).
"""
const PlanetOrderPrior = OrbitOrderPrior
export PlanetOrderPrior

_isprior(::OrbitOrderPrior) = true
likeobj_from_epoch_subset(p::OrbitOrderPrior, _) = p
TypedTables.Table(::OrbitOrderPrior) = nothing
refspecs(p::OrbitOrderPrior) = p.bodies
epochs(::OrbitOrderPrior) = Float64[]
_refdesc(p::OrbitOrderPrior) = join(map(_refstr, p.bodies), " < ")

function ln_like(prior::OrbitOrderPrior, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    ks = _selected_rows(ctx, prior.bodies)
    posys = ctx.system
    @inbounds for i in 1:(length(ks)-1)
        if posys.rows[ks[i]].a >= posys.rows[ks[i+1]].a
            return convert(T, -Inf)
        end
    end
    return zero(T)
end

# ---------------------------------------------------
# Apsidal separation
# ---------------------------------------------------

"""
    LimitClosestApproachAUPrior(soft_closest_approach_au; bodies=(), name=…)
    LimitClosestApproachAUPrior(hard_closest_approach_au, soft_closest_approach_au; …)
    NonCrossingPrior(; bodies=(), name=…)

Require successive orbits to stay apart: sorting the selected orbits by
semi-major axis, the gap between one's apoapsis and the next one's periapsis
must exceed `hard_closest_approach_au` [AU], and is softly penalised (a `1/x²`
repulsive term) below `soft_closest_approach_au`.

`NonCrossingPrior()` is the both-zero case — orbits may not cross — and is the
usual thing to reach for.

`bodies=` restricts the prior to the rows placing those bodies; with no list,
every hierarchy row in the system is included, which is what v8 did. Name the
bodies in any system where the rows are not all planetary orbits of one star —
comparing an inner pair's apsides with a wide binary's is not a meaningful
test.

Carries no data: a prior term (`_isprior = true`), not an observation.

!!! note "This is an apsidal test, not a stability test"
    Non-crossing is necessary for a coplanar system to be long-term stable and
    neither necessary nor sufficient in general (mutually inclined orbits may
    cross in projection and never approach; resonant pairs cross and survive).
    See [`HillStabilityPrior`](@ref) for a mass-aware criterion.
"""
struct LimitClosestApproachAUPrior{TB<:Tuple} <: AbstractObs
    hard_closest_approach_au::Float64
    soft_closest_approach_au::Float64
    bodies::TB
    name::String
end

function LimitClosestApproachAUPrior(hard_closest_approach_au::Real,
                                     soft_closest_approach_au::Real;
                                     bodies=(),
                                     name="LimitClosestApproachAUPrior")
    specs = _bodyspecs(bodies, "LimitClosestApproachAUPrior")
    return LimitClosestApproachAUPrior{typeof(specs)}(
        Float64(hard_closest_approach_au), Float64(soft_closest_approach_au),
        specs, String(name))
end
LimitClosestApproachAUPrior(soft_closest_approach_au::Real; kwargs...) =
    LimitClosestApproachAUPrior(0.0, soft_closest_approach_au; kwargs...)

"""
    NonCrossingPrior(; bodies=(), name="NonCrossingPrior")

Forbid successive orbits from crossing: the zero-threshold case of
[`LimitClosestApproachAUPrior`](@ref), whose docstring documents both.
"""
NonCrossingPrior(; bodies=(), name="NonCrossingPrior") =
    LimitClosestApproachAUPrior(0.0, 0.0; bodies, name)

export NonCrossingPrior, LimitClosestApproachAUPrior

_isprior(::LimitClosestApproachAUPrior) = true
likeobj_from_epoch_subset(p::LimitClosestApproachAUPrior, _) = p
TypedTables.Table(::LimitClosestApproachAUPrior) = nothing
refspecs(p::LimitClosestApproachAUPrior) = p.bodies
epochs(::LimitClosestApproachAUPrior) = Float64[]
_refdesc(p::LimitClosestApproachAUPrior) = _rowsetdesc(p.bodies)

function ln_like(prior::LimitClosestApproachAUPrior, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    posys = ctx.system
    ks = _selected_rows(ctx, prior.bodies)
    length(ks) < 2 && return zero(T)

    ll = zero(T)
    k_in, a_in = _next_by_sma(posys, ks, zero(_rowtype(posys)), 0, false)
    while true
        k_out, a_out = _next_by_sma(posys, ks, a_in, k_in, true)
        k_out == 0 && break
        @inbounds begin
            # k_out has sma >= k_in by construction.
            sep_farthest_inner = _row_apoapsis(posys.rows[k_in])
            sep_closest_outer = _row_periapsis(posys.rows[k_out])
        end
        closest_approach = sep_closest_outer - sep_farthest_inner
        if closest_approach <= prior.hard_closest_approach_au
            return convert(T, -Inf)
        end
        # Repulsive potential, only inside the soft radius (and only if the
        # soft radius is the wider of the two).
        if closest_approach < prior.soft_closest_approach_au
            ll -= 1 / (closest_approach - prior.soft_closest_approach_au)^2
        end
        k_in, a_in = k_out, a_out
    end
    return ll
end

# ---------------------------------------------------
# Hill stability
# ---------------------------------------------------

"""
    HillStabilityPrior(; bodies=(), name="HillStabilityPrior")

Require successive orbits to be separated by more than `2√3` mutual Hill
radii — the Gladman (1993) criterion for a two-planet system on nearly
circular, nearly coplanar orbits. Anything closer gets `-Inf`.

Sorting the selected orbits by semi-major axis, each adjacent pair
`(inner, outer)` must satisfy

    a_out − a_in  >  2√3 · R_H,     R_H = a_out · ((m_in + m_out) / 3M★)^(1/3)

where `m_in`/`m_out` are the total masses of the two rows' exterior bodies and
`M★` is the mass interior to the outer row, excluding the inner row's bodies.

`bodies=` restricts the prior as for [`NonCrossingPrior`](@ref). Carries no
data: a prior term (`_isprior = true`), not an observation.

!!! warning "Not bit-identical to v8"
    v8's `HillStabilityPrior` had a copy-paste bug — it assigned
    `θ_planet_b = θ_system.planets[idx_a]` immediately after assigning it from
    `idx_b`, so *both* masses in `R_H` and in `M★` were the inner planet's, and
    the outer planet's mass never entered the criterion at all. This port
    implements the intended formula, so it will reject configurations v8
    accepted (and vice versa) whenever the two masses differ.

    The definition of `M★` also had to change: v8 read `θ_planet.M`, the
    per-planet total mass v9 has no equivalent of. Reading the mass interior to
    the outer row instead gives `M★ = M_A` for both the Jacobi and the
    astrocentric spelling of the same physical system, which v8's expression
    did not.
"""
struct HillStabilityPrior{TB<:Tuple} <: AbstractObs
    bodies::TB
    name::String
end

function HillStabilityPrior(; bodies=(), name="HillStabilityPrior")
    specs = _bodyspecs(bodies, "HillStabilityPrior")
    return HillStabilityPrior{typeof(specs)}(specs, String(name))
end

export HillStabilityPrior

_isprior(::HillStabilityPrior) = true
likeobj_from_epoch_subset(p::HillStabilityPrior, _) = p
TypedTables.Table(::HillStabilityPrior) = nothing
refspecs(p::HillStabilityPrior) = p.bodies
epochs(::HillStabilityPrior) = Float64[]
_refdesc(p::HillStabilityPrior) = _rowsetdesc(p.bodies)

function ln_like(prior::HillStabilityPrior, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    posys = ctx.system
    ks = _selected_rows(ctx, prior.bodies)
    length(ks) < 2 && return zero(T)

    k_in, a_in = _next_by_sma(posys, ks, zero(_rowtype(posys)), 0, false)
    while true
        k_out, a_out = _next_by_sma(posys, ks, a_in, k_in, true)
        k_out == 0 && break
        s_in = @inbounds posys.specs[k_in]
        s_out = @inbounds posys.specs[k_out]
        m_in = _submass(posys, s_in.ext)
        m_out = _submass(posys, s_out.ext)
        # Everything the outer row orbits, minus whatever the inner row places.
        # For a Jacobi chain the outer row's interior already contains the inner
        # planet, and for an astrocentric set it does not; masking it out makes
        # the two spellings of the same system give the same M★.
        M_star = max(zero(_rowtype(posys)), _submass(posys, s_out.int .& .!s_in.ext))
        R_H = a_out * cbrt((m_in + m_out) / (3M_star))
        if a_out - a_in <= 2 * sqrt(3) * R_H
            return convert(T, -Inf)
        end
        k_in, a_in = k_out, a_out
    end
    return zero(T)
end

# ---------------------------------------------------

# `bodies=` normalization, shared by the two priors that take an optional list.
function _bodyspecs(bodies, who::AbstractString)
    specs = map(refspec, Tuple(bodies))
    all(s -> s isa BodyRefSpec, specs) || error(
        "$who: `bodies=` names the bodies whose orbits the prior applies to, so " *
        "each entry must be a `Body` model node or a `Symbol` naming one.")
    return specs
end

_rowsetdesc(::Tuple{}) = "every orbit"
_rowsetdesc(specs::Tuple) = "orbits of " * join(map(_refstr, specs), ", ")
