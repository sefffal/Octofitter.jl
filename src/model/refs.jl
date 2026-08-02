# ---------------------------------------------------
# Reference specs
#
# Every observation names *what it observes* and *what it is measured
# against*, using the same grammar orbits use for `about=`: a body, the
# barycentre of the system or a subsystem, or a flux-weighted photocentre.
#
# Specs carry their content in type parameters, so resolving one against a
# solved `PlanetOrbits.System` constant-folds to a `BodyRef` or a
# `WeightedPoint` — there is no runtime name lookup in the hot loop.
# ---------------------------------------------------

"""
    BodyRefSpec{Name}

A single body, by name. Written as the `Body` model node itself
(`target=b`) or as a `Symbol` (`target=:b`).
"""
struct BodyRefSpec{Name} end

"""
    Barycentre
    Barycentre(A, b, …)

The mass-weighted barycentre of the whole system, or of the subsystem
spanned by the given bodies. Use as an observation's `target` or `ref`:

    RadialVelocityObs(rvdata; target=A, ref=Barycentre, …)
"""
struct BarycentreSpec{Names} end
const Barycentre = BarycentreSpec{()}()

"""
    Photocentre
    Photocentre(:G)
    Photocentre(:G, (Aa, Ab))
    Photocentre((Aa, Ab))

The flux-weighted photocentre of the system — what a blended astrometric
measurement actually tracks — or, given a member list, of that subset only.
Bodies carry per-band fluxes (a `flux` or `flux_<band>` variable in their
block); pass the band when more than one is defined.

A member list says *this source blends exactly these bodies, in every
draw*. Use it when membership is structurally certain — two pairs several
arcseconds apart, say, where only intra-pair blending is ever possible, and
each catalog source is its own observation:

    GaiaDR4AstromObs(dataA; target=Photocentre(:G, (Aa, Ab)), ref=Barycentre, name="srcA")
    GaiaDR4AstromObs(dataB; target=Photocentre(:G, (Ba, Bb)), ref=Barycentre, name="srcB")

Member names are validated against the system's bodies at model-build time,
and the whole spec lives in type parameters, so it constant-folds. A subset
that is dark in the selected band is an error — a structural membership
declaration over bodies with no flux is not a point on the sky.

Membership that varies *per draw* (a sampled resolved-flag) or *per epoch*
(a scan-angle-dependent window) is not this: it belongs to the observation,
which reads `PlanetOrbits.fluxes(sys, band)` and builds its own
`WeightedPoint`. See the PlanetOrbits "Blended sources & photocentres" docs.
"""
struct PhotocentreSpec{Band,Names} end
const Photocentre = PhotocentreSpec{nothing,()}()

const AbstractRefSpec = Union{BodyRefSpec,BarycentreSpec,PhotocentreSpec}

export Barycentre, Photocentre

# --- construction from user-facing values ------------------------------------

(::BarycentreSpec)(members...) = BarycentreSpec{map(_refname, members)}()

# `Photocentre(band)`, `Photocentre(band, members)`, `Photocentre(members)`.
# The band always comes first, so a leading `Symbol` is never read as a
# member; members may be given as a tuple or as trailing arguments.
(::PhotocentreSpec)(band::Symbol) = PhotocentreSpec{band,()}()
(::PhotocentreSpec)(band::Symbol, members::Union{Tuple,AbstractVector}) =
    PhotocentreSpec{band,_membernames(members)}()
(::PhotocentreSpec)(band::Symbol, m1, members...) =
    PhotocentreSpec{band,_membernames((m1, members...))}()
(::PhotocentreSpec)(members::Union{Tuple,AbstractVector}) =
    PhotocentreSpec{nothing,_membernames(members)}()

function _membernames(members)
    ms = Tuple(members)
    isempty(ms) && error(
        "`Photocentre` was given an empty member list. Write `Photocentre` " *
        "(or `Photocentre(:band)`) for the whole system.")
    return map(_refname, ms)
end

_refname(s::Symbol) = s
_refname(b::BodyRefSpec{N}) where {N} = N
@noinline _refname(@nospecialize x) = error(
    "a reference member must be a `Body` model node or a `Symbol` naming one; " *
    "got a value of type $(typeof(x))")

"""
    refspec(x)

Normalize a user-supplied reference into a spec. Accepts a `Body`/`Orbit`
model node, a `Symbol`, or an already-built spec.
"""
refspec(s::AbstractRefSpec) = s
refspec(s::Symbol) = BodyRefSpec{s}()
@noinline refspec(@nospecialize x) = error(
    "an observation reference must be a `Body` model node, a `Symbol` naming one, " *
    "`Barycentre` (optionally `Barycentre(A, b)`), or `Photocentre` (optionally " *
    "`Photocentre(:band, (Aa, Ab))`); got a value of type $(typeof(x))")

# --- resolution against a solved system --------------------------------------
#
# Every spec's content is in its type parameters, so each of these folds to a
# `BodyRef` or a `WeightedPoint` with no runtime name lookup.

@inline resolveref(posys, ::BodyRefSpec{Name}) where {Name} =
    PlanetOrbits._resolve(PlanetOrbits._names(posys), Name)
@inline resolveref(posys, ::BarycentreSpec{()}) = PlanetOrbits.barycentre(posys)
@inline resolveref(posys, ::BarycentreSpec{Names}) where {Names} =
    PlanetOrbits.barycentre(posys, Names...)
@inline resolveref(posys, ::PhotocentreSpec{Band,()}) where {Band} =
    PlanetOrbits.photocentre(posys; band=Band)
@inline resolveref(posys, ::PhotocentreSpec{Band,Names}) where {Band,Names} =
    PlanetOrbits.photocentre(posys, Names...; band=Band)

# --- display ------------------------------------------------------------------

_refstr(::BodyRefSpec{N}) where {N} = string(N)
_refstr(::BarycentreSpec{()}) = "Barycentre"
_refstr(::BarycentreSpec{Names}) where {Names} = "Barycentre(" * join(Names, ", ") * ")"
_refstr(::PhotocentreSpec{nothing,()}) = "Photocentre"
_refstr(::PhotocentreSpec{B,()}) where {B} = "Photocentre(:$B)"
_refstr(::PhotocentreSpec{nothing,Names}) where {Names} = "Photocentre(" * _memberstr(Names) * ")"
_refstr(::PhotocentreSpec{B,Names}) where {B,Names} = "Photocentre(:$B, " * _memberstr(Names) * ")"
_memberstr(Names) = "(" * join(Names, ", ") * (length(Names) == 1 ? ",)" : ")")

# Names a spec depends on, for validation against the system's body list.
_refnames(::BodyRefSpec{N}) where {N} = (N,)
_refnames(::BarycentreSpec{Names}) where {Names} = Names
_refnames(::PhotocentreSpec{Band,Names}) where {Band,Names} = Names

# Bands a spec needs the bodies to declare, for the same validation. A
# photocentre over the whole system with no band named picks the sole band,
# which is a runtime question (there may be none, or several).
_refbands(::AbstractRefSpec) = ()
_refbands(::PhotocentreSpec{nothing}) = ()
_refbands(::PhotocentreSpec{B}) where {B} = (B,)
