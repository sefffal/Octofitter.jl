# ---------------------------------------------------
# Observations
#
# An observation is *data + references + its own variables*. It never lives
# under a companion: with the references explicit there is nothing left for
# per-body attachment to express, so the old
# PlanetObservationContext/SystemObservationContext split collapses into one
# context type.
# ---------------------------------------------------

"""
    AbstractObs

Supertype for observations. Concrete types bundle a data table with the
nuisance parameters they add to the model, and declare the references they
measure via [`refspecs`](@ref).

The interface is:

  - `ln_like(obs, ctx::ObsContext)` — required.
  - `likelihoodname(obs)` — defaults to `obs.name`.
  - `refspecs(obs)` — the reference specs this observation resolves; `()` for
    observations that touch no orbit (priors, for instance).
  - `epochs(obs)` — the epochs it needs solved; `Float64[]` if none.
  - `likeobj_from_epoch_subset(obs, inds)` — optional, for cross-validation.

# What an observation observes

An observation names a *query* against the solved system — a target
reference and a reference reference. Only *which question is asked* is
static in the type; the answer is recomputed from that sample's masses and
fluxes every draw.

For a blended catalog source that means deciding **when** its membership is
decided. There are three tiers, and they all land on the same primitive, a
`PlanetOrbits.WeightedPoint`:

 1. **Build time.** `target=Photocentre(:G, (Aa, Ab))` — the subset spec.
    For sources whose membership is structurally certain. It constant-folds,
    appears in `refspecs(obs)`, is validated against the body list at
    model-build time, and shows up in `show`.

 2. **Sample time.** The observation computes effective weights itself —
    `w_j ∝ member_j · f_j`, with `member_j` from sampled or derived
    variables (a resolved-flag, a source-assignment latent) — and builds the
    `WeightedPoint` once per evaluation from
    `PlanetOrbits.fluxes(ctx.system, band)`. Membership latents live at
    *system* level; observations may read deferred system variables, so
    blending state never has to round-trip through a body's `flux_<band>`
    variable (which would be the body→deferred-system cycle codegen
    rejects).

 3. **Scan time.** The same construction inside the epoch loop, for
    membership that genuinely varies per transit — a resolution taper in
    separation, or a scan-angle-dependent window. `WeightedPoint` is
    `isbits`, so one per epoch is free.

Tiers 2 and 3 are where instrument-specific blending behaviour belongs:
PlanetOrbits owns per-body states and the two generic linear reductions
(mass-weighted `barycentre`, flux-weighted `photocentre`), and everything an
instrument does on top of those is the observation's business.
"""
abstract type AbstractObs end
TypedTables.Table(obs::AbstractObs) = obs.table

"""
    ln_like(obs::AbstractObs, ctx::ObsContext)

Log-likelihood of `obs`'s data under the parameters and solved trajectory in
`ctx`.
"""
function ln_like end

"""
    likelihoodname(obs)

Name of an observation. It labels the observation's variables in the chain,
so it must be unique within a system.
"""
likelihoodname(obs::AbstractObs) = obs.name
export likelihoodname

"""
    refspecs(obs)

Tuple of reference specs (see [`Barycentre`](@ref)) this observation
resolves against the solved system. Used to validate names at model-build
time and to document the observation in `show`.
"""
refspecs(::AbstractObs) = ()

"""
    epochs(obs)

The epochs [MJD] this observation needs the system solved at, in table
order. Observations with no table return an empty vector.
"""
epochs(obs::AbstractObs) = hasproperty(obs, :table) && hasproperty(obs.table, :epoch) ?
                           collect(Float64, obs.table.epoch) : Float64[]

"""Whether this "observation" is really a prior term (excluded from epoch and likelihood counts)."""
_isprior(::AbstractObs) = false

"""
    check_siblings(obs, all_obs, ctx)

Validation hook called once per observation by [`System`](@ref), with the
model's whole observation tuple and a context

    (; name::Symbol, framevars::Vector{Symbol}, sysnames::Vector{Symbol})

— the system's name, the frame variables it defines, and every variable name
in its block. Default: nothing to check.

This exists because some settings are only wrong *in company*. An
observation's own constructor cannot see how many siblings it has, nor
whether the system it is about to join defines an absolute frame or a
variable it will fall through to — and the failure modes here (`G23HObs(...;
frame_shift=true)` on two Gaia sources; two observations with different
companion counts sharing one system-level flux-ratio vector) are a silently
wrong likelihood and an error naming no observation. Raise from here rather
than warning: a model that is wrong for a structural reason should not run.
"""
check_siblings(@nospecialize(obs::AbstractObs), @nospecialize(all_obs), ctx) = nothing

function likeobj_from_epoch_subset(obs::AbstractObs, obs_inds)
    error("""
    Data subsetting is not supported for observation type $(typeof(obs)).

    Implement `likeobj_from_epoch_subset` for it if you need cross-validation
    or PSIS-LOO.
    """)
end

# ---------------------------------------------------
# Evaluation context
# ---------------------------------------------------

"""
    ObsContext

What a likelihood is handed: the sample's parameters, its own variables, the
`PlanetOrbits.System` built from those parameters, and the `Trajectory`
solved once for the whole model at the union of every observation's epochs.

Use [`solutionat`](@ref) rather than indexing `traj` directly — the
trajectory is over the *deduplicated, sorted* epoch union, not this
observation's table order.
"""
struct ObsContext{Tθ,TO,TS,TT,TE,TB}
    θ_system::Tθ
    θ_obs::TO
    system::TS          # PlanetOrbits.System for this sample
    traj::TT            # solved PlanetOrbits.Trajectory
    epoch_index::TE     # row of this observation's table → column of `traj`
    # The per-sample scratch arena the trajectory was carved out of. A
    # likelihood needing its own temporaries should `@no_escape ctx.buf`
    # rather than reach for `Bumper.default_buffer()`, so that one arena —
    # sized to the model at build time — covers the whole evaluation. See
    # `_slab_size` in codegen.jl for why the size matters.
    buf::TB
    # Solve configuration a likelihood needs to see. `observing_geometry` and
    # `barycentric_lighttime` are consumed by `orbitsolve!` before any
    # likelihood runs, so they never appear here; `einstein_rv` is consumed by
    # the radial-velocity forward model itself, so it does. Interpolated as a
    # literal by codegen, exactly like the other two.
    einstein_rv::Bool
end
# The convenience forms default `einstein_rv` to `true`, the accurate setting.
# Anything that evaluates a likelihood on a *built model* must instead forward
# `sys.einstein_rv`, or a fit run with `einstein_rv=:off` would silently be
# scored against a different physics than it sampled — see `_obsctx`.
ObsContext(θ_system, θ_obs, system, traj, epoch_index) =
    ObsContext(θ_system, θ_obs, system, traj, epoch_index, Bumper.default_buffer(), true)
ObsContext(θ_system, θ_obs, system, traj, epoch_index, buf) =
    ObsContext(θ_system, θ_obs, system, traj, epoch_index, buf, true)

"""
    solutionat(ctx, i)

Per-epoch solution for row `i` of this observation's table.
"""
@inline solutionat(ctx::ObsContext, i::Integer) =
    @inbounds ctx.traj[ctx.epoch_index[i]]

"""
    ref(ctx, spec)

Resolve a reference spec against this sample's system. Constant-folds to a
`BodyRef` or a `WeightedPoint`; call it once outside the epoch loop.
"""
@inline ref(ctx::ObsContext, spec) = resolveref(ctx.system, spec)

"""
    resolverefs(ctx, specs::Tuple)

Resolve a *tuple* of reference specs in one go, returning a tuple of
`BodyRef`/`WeightedPoint` values in the same order. For observations that
name several targets (`ImageObs`, `InterferometryObs`, `G23HObs`), call it
once outside the epoch loop and index the result.

It accepts only a `Tuple`, deliberately. The obvious spelling —

    refs = [ref(ctx, s) for s in obs.targets]     # don't

builds a `Vector` whose element type is the *join* of the resolved types, so
a mixed `(BodyRef, WeightedPoint)` list widens to `Vector{Any}`; it
allocates once per likelihood evaluation, and the widening propagates into
every `raoff` call downstream. `map` over a tuple keeps each element's
concrete type, constant-folds (the specs carry their content in type
parameters), and allocates nothing. Restricting the signature to `Tuple`
makes the fast path the easy path rather than a thing to remember.
"""
@inline resolverefs(ctx::ObsContext, specs::Tuple) = map(s -> ref(ctx, s), specs)
@noinline resolverefs(::ObsContext, @nospecialize specs) = error(
    "`resolverefs` takes a `Tuple` of reference specs, not a $(typeof(specs)). " *
    "Store an observation's target list as a tuple — a vector of specs allocates " *
    "and loses inference on every likelihood evaluation.")

export ObsContext, solutionat, resolverefs

# ---------------------------------------------------
# Prior-shaped terms emitted by `@variables`
#
# These are not observations of anything; they exist so that `~` and `LL +=`
# inside a variables block can contribute to the log density. They read from
# `ctx.θ_obs`, which for a node-owned term is that node's namespace.
# ---------------------------------------------------

"""
    UserLikelihood

`lhs ~ Distribution` written inside a `@variables` block, where `lhs` is an
expression rather than a fresh variable name.
"""
struct UserLikelihood{TSym_LHS,TSym_RHS} <: AbstractObs
    priors::Priors
    derived::Derived
    name::String
end
UserLikelihood(lhs::Symbol, rhs::Symbol, name) =
    UserLikelihood{lhs,rhs}(Priors(), Derived(), String(name))

_isprior(::UserLikelihood) = true
likeobj_from_epoch_subset(l::UserLikelihood, _) = l
TypedTables.Table(::UserLikelihood) = nothing

function ln_like(::UserLikelihood{L,R}, ctx::ObsContext) where {L,R}
    lhs = getproperty(ctx.θ_obs, L)
    rhs = getproperty(ctx.θ_obs, R)
    if rhs isa NTuple{N,<:Number} where {N}
        rhs = SVector(rhs)
    end
    if lhs isa Distribution
        return logpdf(lhs, rhs)
    elseif rhs isa Distribution
        return logpdf(rhs, lhs)
    else
        error("neither side of the `~` expression evaluated to a distribution")
    end
end
Base.show(io::IO, ::MIME"text/plain", ::UserLikelihood{L,R}) where {L,R} =
    println(io, "UserLikelihood: $L ~ $R")

"""
    DirectLLObs

`LL += expr` written inside a `@variables` block: the expression's value is
added straight to the log density. An escape hatch for likelihoods that are
not `logpdf(dist, value)` — marginalized analytic integrals, for instance.
"""
struct DirectLLObs{TSym} <: AbstractObs
    priors::Priors
    derived::Derived
    name::String
end
DirectLLObs(sym::Symbol, name::String) = DirectLLObs{sym}(Priors(), Derived(), name)

_isprior(::DirectLLObs) = true
likeobj_from_epoch_subset(l::DirectLLObs, _) = l
TypedTables.Table(::DirectLLObs) = nothing
ln_like(::DirectLLObs{S}, ctx::ObsContext) where {S} = getproperty(ctx.θ_obs, S)
Base.show(io::IO, ::MIME"text/plain", ::DirectLLObs{S}) where {S} =
    println(io, "DirectLLObs: LL += $S")

"""
    UnitLengthPrior{X,Y}

Keeps the (x, y) pair behind a `UniformCircular` parametrization off the
origin, where the angle is undefined.
"""
struct UnitLengthPrior{X,Y} <: AbstractObs
    varx::Symbol
    vary::Symbol
end
likelihoodname(::UnitLengthPrior{X,Y}) where {X,Y} = "unitlengthprior_$(X)_$(Y)"
_isprior(::UnitLengthPrior) = true
likeobj_from_epoch_subset(l::UnitLengthPrior, _) = l
TypedTables.Table(::UnitLengthPrior) = nothing

function ln_like(::UnitLengthPrior{X,Y}, ctx::ObsContext) where {X,Y}
    x = getproperty(ctx.θ_obs, X)
    y = getproperty(ctx.θ_obs, Y)
    return logpdf(LogNormal(log(1.0), 0.1), sqrt(x^2 + y^2))
end
Base.show(io::IO, ::MIME"text/plain", ::UnitLengthPrior{X,Y}) where {X,Y} =
    println(io, "UnitLengthPrior: √($X^2+$Y^2) ~ LogNormal(log(1), 0.1)")

"""
    BlankLikelihood()

A no-op term, for models that need a variable namespace with nothing to
constrain it (prior-predictive checks, for instance).
"""
struct BlankLikelihood <: AbstractObs
    priors::Priors
    derived::Derived
    name::String
    BlankLikelihood(variables::Tuple{Priors,Derived}=(Priors(), Derived()), name="") =
        new(variables[1], variables[2], name)
end
_isprior(::BlankLikelihood) = true
ln_like(::BlankLikelihood, ctx::ObsContext) = zero(_system_number_type(ctx.θ_system))
TypedTables.Table(::BlankLikelihood) = nothing
likeobj_from_epoch_subset(l::BlankLikelihood, _) = l

# ---------------------------------------------------

"""
    _refdesc(obs) -> String
    _refdesc(targets::Tuple, reference) -> String

One-line description of what an observation observes, for `show`: its
target(s), then the reference they are measured against.

The default reads `refspecs(obs)` as "everything but the last entry is a
target". That is right for the pair types and for the variadic ones as they
are written, but it is a guess — a type whose refspec tuple is ordered
differently, or which wants to say something more specific about the roles
(`G23HObs` distinguishes its `host` from the companions blended into it),
should define its own method.

`show` used to `join(map(_refstr, refspecs(obs)), " vs ")`, which reads as a
chain of differences: `G23HObs` printed `A vs Ab vs Barycentre`, as though
three things were being subtracted from one another rather than one blended
source being measured against a reference point.
"""
_refdesc(@nospecialize obs::AbstractObs) = _refdesc_default(refspecs(obs))

_refdesc_default(::Tuple{}) = ""
_refdesc_default(specs::Tuple{Any}) = _refstr(specs[1])
_refdesc_default(specs::Tuple) = _refdesc(Base.front(specs), last(specs))

_refdesc(targets::Tuple{Any}, reference) = _refstr(targets[1]) * " vs " * _refstr(reference)
_refdesc(targets::Tuple, reference) =
    "(" * join(map(_refstr, targets), ", ") * ") vs " * _refstr(reference)

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize obs::AbstractObs)
    nm = string(typeof(obs).name.name)
    print(io, nm, " \"", likelihoodname(obs), "\"")
    desc = _refdesc(obs)
    if !isempty(desc)
        print(io, "  ", desc)
    end
    println(io)
    if hasproperty(obs, :table)
        show(io, mime, TypedTables.Table(obs))
        println(io)
    end
end
