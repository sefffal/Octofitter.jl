# ---------------------------------------------------
# The two retired names worth catching.
#
# `Planet` and `θ_at_epoch_to_tperi` are the first things an old script hits,
# and left undefined each fails as a bare `UndefVarError` from inside a
# `@variables` block or a `System(...)` call — accurate, and useless. A stub
# that names the replacement costs one method.
#
# They deliberately *error* rather than aliasing: in both cases the
# replacement takes different arguments, so an alias would not have helped.
# Every other retired name is simply gone, so that the name itself is free
# again; `docs/src/v9-migration.md` carries the full mapping table.
#
# (The HGCA family has its own stubs in `likelihoods/hgca-compat.jl`, next to
# the helper that replaces it.)
# ---------------------------------------------------

const _V2_GUIDE = "See `docs/src/v9-migration.md` (\"Migrating to Octofitter v9\")."

@noinline θ_at_epoch_to_tperi(args...; kwargs...) = error("""
`θ_at_epoch_to_tperi` no longer exists, and nothing needs it: `θ` (position
angle at a reference epoch) and `epoch` are orbital-element keywords now, so
the conversion happens inside the orbit constructor.

    # v8
    tp = θ_at_epoch_to_tperi(θ, 58849; M=system.M, e, a, i, ω, Ω)

    # v9 — declare θ and epoch as elements and drop tp entirely
    θ ~ UniformCircular()
    epoch = 58849.0

$(_V2_GUIDE)""")
export θ_at_epoch_to_tperi

@noinline Planet(args...; kwargs...) = error("""
`Planet` no longer exists — v9 has one node type, `Body`, and observations do
not live on a companion.

    # v8
    b = Planet(name="b", basis=Visual{KepOrbit}, observations=[astrom],
               variables=@variables begin a ~ ...; M = system.M end)

    # v9
    b = Body(name="b", about=A, variables=@variables begin a ~ ...; end)
    sys = System(name="sys", bodies=[A, b], observations=[astrom], variables=...)

`about=` names the body (or bodies) the orbit is around, which is what the
v8 `basis=`/`M =` pair was standing in for; the frame is chosen by which
frame variables the `System` block declares. $(_V2_GUIDE)""")
export Planet


# `rvpostplot` is a working alias rather than an error: the figure still
# exists, it just has an honest name now. It never plotted a posterior — it
# plots one draw against every instrument's data, which is exactly why it can
# share a panel between instruments at all.
"""
    rvpostplot(model, chain, [sample_idx]; kwargs...)

Deprecated alias for [`rvplot`](@ref), which is the same figure under a name
that does not claim to show the posterior. $(_V2_GUIDE)
"""
function rvpostplot(args...; kwargs...)
    Base.depwarn("`rvpostplot` is now `rvplot`: the figure shows a single " *
                 "posterior draw, not the posterior.", :rvpostplot)
    return rvplot(args...; kwargs...)
end

"""
    rvpostplot_animated(model, chain; kwargs...)

Deprecated alias for [`rvplot_animated`](@ref). $(_V2_GUIDE)
"""
function rvpostplot_animated(args...; kwargs...)
    Base.depwarn("`rvpostplot_animated` is now `rvplot_animated`.", :rvpostplot_animated)
    return rvplot_animated(args...; kwargs...)
end

export rvpostplot, rvpostplot_animated
