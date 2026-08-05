# ---------------------------------------------------
# `HGCAObs` — a helper constructor over `G23HObs`   (agent E)
#
# The HGCA family is retired as a set of types: `HGCAObs`,
# `HGCAInstantaneousObs` and `GaiaCatalogFitObs` are all subsumed by
# `G23HObs` restricted to the channels the HGCA actually constrains
# (Hipparcos, Hipparcos–Gaia, DR3 — no DR2, no DR32, no RV), with
# `ueva_mode = :none`. See `design/observation-types-migration.md` §3.6.
#
# Release-note item: HGCA fits are NOT bit-identical to v1, because those
# channels are now modelled by G23H's treatment rather than by the HGCA's own
# cross-calibration.
# ---------------------------------------------------

"""
    HGCAObs(; gaia_id, host, companions=(), ref=Barycentre, kwargs...)

Proper-motion-anomaly astrometry in the style of the Hipparcos–Gaia Catalog
of Accelerations: the Hipparcos and Gaia DR3 catalog proper motions and the
long-baseline Hipparcos–Gaia scaled position difference.

This is a **helper constructor**, not a type. It builds a [`G23HObs`](@ref)
restricted to the six channels the HGCA constrains:

    HGCAObs(; gaia_id, host, companions=(), ref=Barycentre, kwargs...) =
        G23HObs(; gaia_id, host, companions, ref,
                  channels    = (:ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr3, :dec_dr3),
                  ueva_mode   = :none,
                  include_iad = false,
                  include_rv  = false,
                  kwargs...)

so everything `G23HObs` documents — source membership through
`host`/`companions`, flux ratios defaulting to the bodies' own `flux_G` /
`flux_Hp`, `variables=`, the offline `catalog=`/`forecast_table=`/
`hipparcos=` inputs — applies unchanged, and `likeobj_from_epoch_subset`,
`generate_from_params` and cross-validation come along for free.

Three parts of the mapping are load-bearing:

  - **`ueva_mode = :none`, not `false`.** It is a three-valued symbol
    (`:RUWE`/`:EAN`/`:none`), and `:none` drops both the `:ueva_dr3` datum
    *and* the UEVA-driven deflation of the DR3 covariance. The HGCA has no
    excess-noise channel, so both go.
  - **No DR2 channels.** `:ra_dr2`/`:dec_dr2` and the DR3−DR2 difference
    `:ra_dr32`/`:dec_dr32` are dropped along with `:rv_dr3`: the HGCA is
    Hipparcos + Hipparcos–Gaia + DR3.
  - **Not bit-identical to v1.** Those six channels are now modelled by
    G23H's treatment — its epoch-selection model, its five-parameter refits,
    its photocentre — rather than by the HGCA's own cross-calibration. The
    numbers move. This deserves one validation run against a published fit.

# Example

```julia
A = Body(name=:A, variables=@variables begin
    mass ~ truncated(Normal(1.0, 0.1), lower=0.1)
    flux_G = 1.0                 # the host defines the contrast scale
end)
b = Body(name=:b, about=A, variables=@variables begin
    mass ~ Uniform(0, 100) * mjup2msol
    flux_G = 0.0                 # dark companion
    flux_Hp = 0.0
    a ~ LogUniform(1, 100); e ~ Uniform(0, 0.99)
    i ~ Sine(); ω ~ UniformCircular(); Ω ~ UniformCircular(); tp ~ Uniform(50000, 60000)
end)

sys = System(name=:HD1234, bodies=(A, b), observations=(
    HGCAObs(; gaia_id=756291174721509376, host=A, companions=(b,), ref=Barycentre),
), variables=@variables begin
    plx ~ truncated(Normal(24.0, 0.1), lower=0)
end)
```

See also [`G23HObs`](@ref), which is what you want if you have the full G23H
catalog row and can afford its extra channels.
"""
HGCAObs(; host, companions=(), ref=Barycentre, name::AbstractString="HGCA", kwargs...) =
    G23HObs(; host, companions, ref, name,
        channels=(:ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr3, :dec_dr3),
        ueva_mode=:none,
        include_iad=false,
        include_rv=false,
        kwargs...)

export HGCAObs

# The two other HGCA-family types are gone rather than deprecated: their
# modelling code (`hgca-linfit.jl`, `hgca.jl`, and the `GaiaCatalogFitObs`
# half of the legacy fitting code) is not ported at all. A model that names
# them should say so clearly instead of failing on an undefined variable.
@noinline HGCAInstantaneousObs(args...; kwargs...) = error(
    "`HGCAInstantaneousObs` no longer exists. The instantaneous proper-motion " *
    "treatment it implemented is what `G23HObs` does for every channel, so use " *
    "`HGCAObs(; gaia_id, host, companions, ref)` — or `G23HObs` directly for the " *
    "full channel set. See the v2 migration guide.")

@noinline GaiaCatalogFitObs(args...; kwargs...) = error(
    "`GaiaCatalogFitObs` no longer exists. It was only ever the Gaia half of " *
    "`HGCAObs`, and refitting the five astrometric parameters to a modelled sky " *
    "path is now `G23HObs`'s `:ra_dr3`/`:dec_dr3` channels (or `GaiaDR4AstromObs` " *
    "for epoch astrometry). See the v2 migration guide.")

export HGCAInstantaneousObs, GaiaCatalogFitObs
