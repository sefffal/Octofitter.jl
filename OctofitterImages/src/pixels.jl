# ---------------------------------------------------
# Where a source lands on the detector
#
# Both likelihoods in this package ask the trajectory the same question —
# "which pixel of *this* image is this body in, at this epoch" — and then
# differ only in what they read out of the image there. That question is
# Octofitter's shared `sky_offset` (design/observation-types-migration.md
# §2.A) followed by the RA sign flip.
#
# What is gone: v1 open-coded, separately in `images.jl` and
# `likelihood-maps.jl`, (a) its own copy of the rotate-by-northangle,
# divide-by-platescale block and (b) the epicyclic superposition loop, which
# rebuilt a companion's apparent position by summing the reflex motion of
# every body it decided at runtime was "inner" (`semimajoraxis(other) <
# semimajoraxis(this)`). That loop had no answer for crossing orbits, for
# equal semi-major axes, or for a moon, and it silently skipped companions
# with no `mass` variable. `raoff(sol, target, ref)` is the whole model now.
# ---------------------------------------------------

"""
    pixel_position(sol, target, reference, platescale, northangle) -> (x, y)

Position of `target` relative to `reference` at one epoch, in the image's own
pixel coordinates. `sol` is a per-epoch solution (`Octofitter.solutionat`),
`target` and `reference` are already-resolved references (see
`Octofitter.resolverefs`), and `platescale` is this image's mas/pixel scale
*including* the sampled multiplier.

`Octofitter.sky_offset` carries the model onto the detector — rotating the
position angle by `−northangle` and dividing by the plate scale — so its
output is already in pixels. The only thing left is the axis flip: right
ascension increases to the east, which is to the left in a north-up
east-left image, so the image's x axis runs opposite to Δα✱. (v1 wrote this
as `x = -ra_rotated`, `y = +dec_rotated`, and it means the same thing.)

The image's coordinates are pixel *offsets from the star* — every image is
expected to be recentred so index `[0,0]` is the reference position, which
is what makes `reference` the origin of this frame.
"""
@inline function pixel_position(sol, target, reference, platescale, northangle)
    Δα, Δδ = Octofitter.sky_offset(sol, target, reference; platescale, northangle)
    return (-Δα, Δδ)
end
