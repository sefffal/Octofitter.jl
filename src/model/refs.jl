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

The flux-weighted photocentre of the system — what a blended astrometric
measurement actually tracks. Bodies carry per-band fluxes (a `flux` or
`flux_<band>` variable in their block); pass the band when more than one is
defined.
"""
struct PhotocentreSpec{Band} end
const Photocentre = PhotocentreSpec{nothing}()

const AbstractRefSpec = Union{BodyRefSpec,BarycentreSpec,PhotocentreSpec}

export Barycentre, Photocentre

# --- construction from user-facing values ------------------------------------

(::BarycentreSpec)(members...) = BarycentreSpec{map(_refname, members)}()
(::PhotocentreSpec)(band::Symbol) = PhotocentreSpec{band}()

_refname(s::Symbol) = s
_refname(b::BodyRefSpec{N}) where {N} = N

"""
    refspec(x)

Normalize a user-supplied reference into a spec. Accepts a `Body`/`Orbit`
model node, a `Symbol`, or an already-built spec.
"""
refspec(s::AbstractRefSpec) = s
refspec(s::Symbol) = BodyRefSpec{s}()
@noinline refspec(@nospecialize x) = error(
    "an observation reference must be a `Body` model node, a `Symbol` naming one, " *
    "`Barycentre` (optionally `Barycentre(A, b)`), or `Photocentre`; got a value of " *
    "type $(typeof(x))")

# --- resolution against a solved system --------------------------------------

@inline resolveref(posys, ::BodyRefSpec{Name}) where {Name} =
    PlanetOrbits._resolve(PlanetOrbits._names(posys), Name)
@inline resolveref(posys, ::BarycentreSpec{()}) = PlanetOrbits.barycentre(posys)
@inline resolveref(posys, ::BarycentreSpec{Names}) where {Names} =
    PlanetOrbits.barycentre(posys, Names...)
@inline resolveref(posys, ::PhotocentreSpec{Band}) where {Band} =
    PlanetOrbits.photocentre(posys; band=Band)

# --- display ------------------------------------------------------------------

_refstr(::BodyRefSpec{N}) where {N} = string(N)
_refstr(::BarycentreSpec{()}) = "Barycentre"
_refstr(::BarycentreSpec{Names}) where {Names} = "Barycentre(" * join(Names, ", ") * ")"
_refstr(::PhotocentreSpec{nothing}) = "Photocentre"
_refstr(::PhotocentreSpec{B}) where {B} = "Photocentre(:$B)"

# Names a spec depends on, for validation against the system's body list.
_refnames(::BodyRefSpec{N}) where {N} = (N,)
_refnames(::BarycentreSpec{Names}) where {Names} = Names
_refnames(::PhotocentreSpec) = ()
