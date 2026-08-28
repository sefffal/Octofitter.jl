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

struct BarycentreSpec{Names} end

"""
    Barycentre
    Barycentre(A, b, …)

The mass-weighted barycentre of the whole system, or of the subsystem
spanned by the given bodies. Use as an observation's `target` or `ref`:

    RadialVelocityObs(rvdata; target=A, ref=Barycentre, …)

# Why absolute astrometry also says `ref=Barycentre`

`ref` names the point an observation's modelled **offsets** are measured
from. It is not a claim that the point is at rest. The system barycentre
does move across the sky in general — and that motion is exactly what the
absolute frame in the **system** block describes: `ra`, `dec`, `plx`,
`pmra`, `pmdec` and `rv` at `ref_epoch` are the barycentre's own catalogue
quantities, which PlanetOrbits propagates rigorously in 3D, so its track
carries perspective acceleration and light-travel curvature rather than
being a straight line in (α, δ).

An absolute-astrometry prediction is the sum of the two:

    absolute position = the barycentre's propagated track   # from the frame
                      + (target − ref) at this epoch        # from `target`/`ref`

with the annual parallax supplied by the observation's own parallax factors.
`G23HObs` composes them literally — `frame_pmra(sol) + Δpmra`, where `Δpmra`
is the five-parameter refit of `raoff(sol, target, Barycentre)` over the
transits. So `ref=Barycentre` says "this source's excursion is measured from
the system's centre of mass"; the centre of mass moving is not an omission,
it is the frame's job.

Naming a **sub**-barycentre is legal and is usually *not* what you want for a
catalogue source in a multi-source system: `ref=Barycentre(Aa, Ab)` measures
the A pair's photocentre against its own centre of mass, which removes that
pair's motion about the system barycentre from the prediction — and for a
wide pair that motion is the signal. Both sources take the whole-system
`ref=Barycentre`, because both are predicted from the one shared frame; see
the [G23H tutorial](@ref fit-g23h).
"""
const Barycentre = BarycentreSpec{()}()

struct PhotocentreSpec{Band,Names} end

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

# ---------------------------------------------------
# Declaration specs are not queries
#
# `Barycentre` and `Photocentre` are *declaration* specs: singletons whose
# content lives in type parameters so that `resolveref` constant-folds. They
# are valid wherever a model says what an observation looks at (`target=`,
# `ref=`, `ObservableQuery`), and they mean nothing on their own — resolving
# one needs a solved system, because the weights come from the bodies' masses
# or fluxes.
#
# The natural thing to write after reading any tutorial is
# `radvel(sol, :A, Barycentre)`, which without these methods is a MethodError
# with a page of candidate signatures and no hint. Say what to do instead.
# ---------------------------------------------------
for fn in (:posx, :posy, :posz, :velx, :vely, :velz, :radvel,
           :raoff, :decoff, :pmra, :pmdec, :projectedseparation, :posangle)
    @eval begin
        @noinline PlanetOrbits.$fn(sol::PlanetOrbits.TrajectorySolution, t, r::Union{BarycentreSpec,PhotocentreSpec}) =
            _spec_is_not_a_ref($(QuoteNode(fn)), r)
        @noinline PlanetOrbits.$fn(sol::PlanetOrbits.TrajectorySolution, t::Union{BarycentreSpec,PhotocentreSpec}, r) =
            _spec_is_not_a_ref($(QuoteNode(fn)), t)
        @noinline PlanetOrbits.$fn(sol::PlanetOrbits.TrajectorySolution, t::Union{BarycentreSpec,PhotocentreSpec}, r::Union{BarycentreSpec,PhotocentreSpec}) =
            _spec_is_not_a_ref($(QuoteNode(fn)), t)
    end
end

@noinline function _spec_is_not_a_ref(fn::Symbol, spec)
    s = _refstr(spec)
    resolver = spec isa BarycentreSpec ? "PlanetOrbits.barycentre(posys)" :
                                         "PlanetOrbits.photocentre(posys; band=:H)"
    error("""
    `$s` is a declaration spec, not a resolved reference, so it cannot be passed
    directly to `$fn`. It says *what* an observation looks at (in `target=`, `ref=`
    or an `ObservableQuery`); turning it into a point needs the solved system,
    because its weights come from the bodies' masses or fluxes.

    Given a `PlanetOrbits.System` `posys` and a solution `sol`:

        $fn(sol, :A, Octofitter.resolveref(posys, $s))

    or equivalently `$fn(sol, :A, $resolver)`.
    """)
end
