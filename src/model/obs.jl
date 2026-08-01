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
end
ObsContext(θ_system, θ_obs, system, traj, epoch_index) =
    ObsContext(θ_system, θ_obs, system, traj, epoch_index, Bumper.default_buffer())

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

export ObsContext, solutionat

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

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize obs::AbstractObs)
    nm = string(typeof(obs).name.name)
    print(io, nm, " \"", likelihoodname(obs), "\"")
    if !isempty(refspecs(obs))
        print(io, "  ", join(map(_refstr, refspecs(obs)), " vs "))
    end
    println(io)
    if hasproperty(obs, :table)
        show(io, mime, TypedTables.Table(obs))
        println(io)
    end
end
