# ---------------------------------------------------
# The sky-plane offset front-end
#
# Five likelihoods ask the trajectory the same question — "where is this
# target on the sky relative to this reference, in mas, as *this instrument*
# would have reported it" — and then differ only in the noise model they wrap
# around the answer: `RelAstromObs` (2-D Gaussian), `ImageObs` (interpolate a
# reduced image), `LogLikelihoodMapObs` (interpolate a precomputed logL
# surface), and the interferometric types (Fourier-transform the flux
# distribution). See `design/observation-types-migration.md` §2.A.
#
# In v1 the rotate-and-scale block below was copy-pasted verbatim into four
# of them, each next to its own copy of the epicyclic superposition loop.
# Both are gone: the position is `raoff(sol, target, ref)`, and the
# calibration is applied here, once.
# ---------------------------------------------------

"""
    sky_calibration(ctx) -> (platescale, northangle)

The two instrument-calibration nuisance parameters, read from this
observation's own variables with their identity defaults (`platescale = 1`,
`northangle = 0`) and promoted to the sample's number type so they stay
ForwardDiff-clean when only one of the two is being sampled.

Every likelihood with a plate scale spells the lookup the same way; this is
that spelling, so a new one cannot quietly pick a different default.
"""
@inline function sky_calibration(ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    platescale = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : one(T)
    northangle = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : zero(T)
    return (platescale, northangle)
end

"""
    sky_offset(sol, target, reference; platescale=1, northangle=0) -> (Δα✱, Δδ)

Single-epoch sky-plane offset [mas] of `target` from `reference`, carried
onto the detector by the instrument's calibration. `sol` is a per-epoch
solution (`solutionat(ctx, i)`) and the two references are already resolved
— by `ref(ctx, spec)` or [`resolverefs`](@ref), outside the epoch loop.
Anything `raoff` accepts works, including a body name.

# Convention

`platescale` and `northangle` mean exactly what they mean in
[`RelAstromObs`](@ref): they describe the instrument's calibration *error*,
so that a **reported** measurement is corrected to the true sky by

    sep_true = sep_reported × platescale
    pa_true  = pa_reported + northangle        (position angle: N through E)

This function applies the **inverse** — it takes the model, which lives on
the true sky, and produces what the instrument would have reported:

    Δα✱ = ( Δα✱_model·cos(northangle) − Δδ_model·sin(northangle) ) / platescale
    Δδ  = ( Δα✱_model·sin(northangle) + Δδ_model·cos(northangle) ) / platescale

which is a rotation of the position angle by `−northangle` and a shrink by
`platescale`. So the values it returns are directly comparable with the
numbers in the data table — the pixel grid of an image, the phase reference
of a visibility — with no further correction.

That is the direction the four consumers want, and it is the opposite of
what `RelAstromObs` does internally: `RelAstromObs` moves the *data* onto
the true sky instead, because its residual and its σ live in the data's own
(sep, pa) or (ra, dec) frame. The two are not interchangeable — a residual
formed in the model frame is the data-frame residual rotated and scaled,
which reweights it against an untransformed σ — so do not "unify" them
without deciding which frame the σ belongs to.

!!! note "v1 discrepancy, deliberately not reproduced"
    v1's `OctofitterInterferometry` **multiplied** its offsets by
    `platescale` rather than dividing, i.e. it used the reciprocal of the
    convention above, while `OctofitterImages` (which divides the model
    position by the image's mas/pixel scale times the multiplier) and
    `PlanetRelAstromObs` both agree with it. The majority spelling is the
    one implemented here; a port of the interferometry likelihood therefore
    changes the sense of a non-unity `platescale`.
"""
@inline function sky_offset(sol, target, reference; platescale=1, northangle=0)
    sn, cs = sincos(northangle)
    return _rotscale(raoff(sol, target, reference), decoff(sol, target, reference),
                     sn, cs, inv(platescale))
end

"""
    sky_offset!(Δα✱, Δδ, ctx, target, reference; platescale=1, northangle=0)

Fill `Δα✱` and `Δδ` with the sky-plane offset [mas] of `target` from
`reference` at every epoch of this observation's table, in table order, with
the instrument calibration applied. Returns `(Δα✱, Δδ)`.

The buffers are **overwritten**, not accumulated into — unlike
[`accumulate_offsets!`](@ref), whose callers lay down a reference-point
linear motion first. One entry is written per row of the observation's
table; pass `@alloc`'d storage from `ctx.buf` rather than fresh vectors.

`target` and `reference` may be specs or already-resolved references; either
way they are resolved once, outside the loop, and `sincos(northangle)` is
computed once. Nothing is allocated and nothing is typed to `Float64`, so
this is safe in the hot loop and under ForwardDiff.

See [`sky_offset`](@ref) for the sign and scale convention — it is stated
there in full, because a sign error here propagates into four likelihoods.
"""
function sky_offset!(Δα, Δδ, ctx::ObsContext, target, reference;
                     platescale=1, northangle=0)
    axes(Δα) == axes(Δδ) || error(
        "sky_offset!: the two output buffers must have the same axes; got " *
        "$(axes(Δα)) and $(axes(Δδ)).")
    length(Δα) == length(ctx.epoch_index) || error(
        "sky_offset!: the output buffers are indexed by table row, so they must " *
        "have length $(length(ctx.epoch_index)); got $(length(Δα)).")
    t = _resolvedref(ctx, target)
    r = _resolvedref(ctx, reference)
    sn, cs = sincos(northangle)
    invps = inv(platescale)
    @inbounds for i in eachindex(Δα)
        sol = solutionat(ctx, i)
        Δα[i], Δδ[i] = _rotscale(raoff(sol, t, r), decoff(sol, t, r), sn, cs, invps)
    end
    return (Δα, Δδ)
end

# Exported because `OctofitterImages` and `OctofitterInterferometry` are
# separate packages and this is their shared front-end.
export sky_offset, sky_offset!, sky_calibration

# Rotate the position angle by −northangle and divide out the plate scale.
# Written over `sincos` products rather than over (sep, pa) so that the
# identity calibration is bit-exact: with northangle = 0 and platescale = 1
# this is `x*1 - y*0` and `x*1 + 0`, which is `x` for every finite input.
@inline _rotscale(x, y, sn, cs, invps) =
    ((x * cs - y * sn) * invps, (x * sn + y * cs) * invps)

# Specs resolve; already-resolved references pass through. Both branches
# constant-fold, so callers may pass whichever they have without paying for
# the choice.
@inline _resolvedref(ctx::ObsContext, spec::AbstractRefSpec) = ref(ctx, spec)
@inline _resolvedref(ctx::ObsContext, spec::Symbol) = ref(ctx, refspec(spec))
@inline _resolvedref(::ObsContext, r) = r
