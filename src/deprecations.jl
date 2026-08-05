# ---------------------------------------------------
# Names retired by the v2 model surface.
#
# Every entry here is a name that a working v1 script contains and that v2
# does not define. Left undefined, each of them fails as a bare
# `UndefVarError` from inside a `@variables` block or a `System(...)` call —
# accurate, and useless. A stub that names the replacement costs one method
# and saves the reader a search through the migration guide.
#
# These deliberately *error*: none of them can be aliased to a v2 equivalent,
# because in every case the replacement takes different arguments or lives
# somewhere else in the model. See `docs/src/v2-migration.md` for the full
# mapping table.
#
# (The HGCA family has its own stubs in `likelihoods/hgca-compat.jl`, next to
# the helper that replaces it.)
# ---------------------------------------------------

const _V2_GUIDE = "See `docs/src/v2-migration.md` (\"Migrating to v2\")."

@noinline θ_at_epoch_to_tperi(args...; kwargs...) = error("""
`θ_at_epoch_to_tperi` no longer exists, and nothing needs it: `θ` (position
angle at a reference epoch) and `epoch` are orbital-element keywords now, so
the conversion happens inside the orbit constructor.

    # v1
    tp = θ_at_epoch_to_tperi(θ, 58849; M=system.M, e, a, i, ω, Ω)

    # v2 — declare θ and epoch as elements and drop tp entirely
    θ ~ UniformCircular()
    epoch = 58849.0

$(_V2_GUIDE)""")
export θ_at_epoch_to_tperi

@noinline Planet(args...; kwargs...) = error("""
`Planet` no longer exists — v2 has one node type, `Body`, and observations do
not live on a companion.

    # v1
    b = Planet(name="b", basis=Visual{KepOrbit}, observations=[astrom],
               variables=@variables begin a ~ ...; M = system.M end)

    # v2
    b = Body(name="b", about=A, variables=@variables begin a ~ ...; end)
    sys = System(name="sys", bodies=[A, b], observations=[astrom], variables=...)

`about=` names the body (or bodies) the orbit is around, which is what the
v1 `basis=`/`M =` pair was standing in for; the frame is chosen by which
frame variables the `System` block declares. $(_V2_GUIDE)""")
export Planet

@noinline PlanetRelAstromObs(args...; kwargs...) = error("""
`PlanetRelAstromObs` is now `RelAstromObs`, and it names both ends of the
measurement explicitly instead of inheriting one from the planet it was
attached to:

    RelAstromObs(tab; target=b, ref=A, name="GPI")

It goes in the `System`'s `observations=` list. $(_V2_GUIDE)""")
export PlanetRelAstromObs

@noinline StarAbsoluteRVObs(args...; kwargs...) = error("""
`StarAbsoluteRVObs` is now `RadialVelocityObs` with the reflex pair named
explicitly, and it lives in core Octofitter rather than in
OctofitterRadialVelocity:

    RadialVelocityObs(tab; target=A, ref=Barycentre, name="HARPS", variables=...)

Note that v2 never invents a prior: v1 auto-injected `offset ~ Uniform(-1000,
1000)` and `jitter ~ LogUniform(0.001, 100)` when you gave no `variables=`
block, and v2 does not. Declare both explicitly or the fit has no zero point
and no jitter. $(_V2_GUIDE)""")
export StarAbsoluteRVObs

@noinline PlanetRelativeRVObs(args...; kwargs...) = error("""
`PlanetRelativeRVObs` is now `RadialVelocityObs` with the pair named
explicitly — the same type as the former `StarAbsoluteRVObs`, differing only
in its refs:

    RadialVelocityObs(tab; target=b, ref=A, name="relative RV", variables=...)

$(_V2_GUIDE)""")
export PlanetRelativeRVObs

@noinline MarginalizedStarAbsoluteRVObs(args...; kwargs...) = error("""
`MarginalizedStarAbsoluteRVObs` is now `OctofitterRadialVelocity.MarginalizedRVObs`
("Star" and "Absolute" were only ever ref choices), and it requires `target=`:

    using OctofitterRadialVelocity
    MarginalizedRVObs(tab; target=A, ref=Barycentre, name="HARPS", variables=...)

It also errors now if you declare an `offset` variable, which v1 silently
added on top of the marginalization. $(_V2_GUIDE)""")
export MarginalizedStarAbsoluteRVObs

@noinline masspostplot(args...; kwargs...) = error("""
`masspostplot` has been removed deliberately, with no replacement. It drew a
one-panel histogram of each body's mass posterior — which is one column of
`octocorner(model, chain)`, or one line of Julia:

    hist(vec(chain[:b_mass]))

$(_V2_GUIDE)""")
export masspostplot
