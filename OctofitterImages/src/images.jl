# ---------------------------------------------------
# Direct images
#
# One image contains every companion at once, so in v2 one `ImageObs`
# describes every companion in it (`targets=(b, c)`). v1 forced one
# `ImageObs` per planet and added the resulting likelihoods, which counts the
# same image's background once per companion and is only defensible when they
# are well separated.
#
# The per-target flux moved with it: v1 read a single `θ_obs.flux` from the
# observation's own variables, so two companions in one image needed two
# observations to have two fluxes. Each body now carries its own
# `flux_<band>` variable and this likelihood reads them through
# `PlanetOrbits.fluxes` — which is also what photocentres are built from, so
# a body's brightness means one thing across the whole model.
# ---------------------------------------------------

const images_cols = (:image, :epoch, :platescale,)

"""
    ImageObs(data; targets, ref, band=nothing, name="images", variables=@variables begin end)

A set of reduced images, and the companions modelled in them.

    ImageObs(image_dat; targets=(b, c), ref=A, band=:H, name="SPHERE")

`data` needs the columns $images_cols — one row per image — where `image` is
an `AstroImage` **recentred so that index `[0,0]` is `ref`**, `epoch` is MJD,
and `platescale` is that image's plate scale in mas/pixel.

# Sources
`targets` lists the bodies whose signal this likelihood models; `ref` is the
point the images are centred on (usually the host). `targets` may be a single
body or a tuple, and both take the usual grammar — a `Body` node, a `Symbol`,
or a `Photocentre` for a pair the imager does not resolve.

`targets` is a *structural* statement about which sources the forward model
includes, so it is deliberately not inferred from which bodies happen to have
a flux in this band: "fit the image term with `c` alone" and "include `b` but
let its flux go to zero" are different models, and only the first can be said
by leaving `b` out.

Each target's brightness is its own `flux_<band>` variable, in the image's
units (contrast, Jy, counts — whatever the pixels are, consistently). `band`
selects among them and may be omitted when the bodies declare exactly one.

# Contrast
Pass a `contrast` column (a callable mapping separation in *pixels* to the 1σ
flux uncertainty there, e.g. `contrast_interp`) or a `contrastmap`
column (an image of the same geometry). With neither, a contrast curve is
measured from each image by `contrast`.

# Variables
  - `platescale` — multiplier on the plate scale of every image [default 1].
  - `northangle` [rad] — rotation of the images relative to true north
    [default 0].

Neither is a flux: `flux` on an `ImageObs` is a v8 spelling and is rejected.

# Example
```julia
image_dat = Table(;
    epoch      = [mjd("2016-01-01"), mjd("2017-01-01")],
    image      = [AstroImages.recenter(img1), AstroImages.recenter(img2)],
    platescale = [19.4, 19.4],
)

ImageObs(image_dat; targets=(b, c), ref=A, band=:H, name="SPHERE",
    variables=@variables begin
        platescale = 1.0    # or: ~ truncated(Normal(1, 0.01), lower=0)
        northangle = 0.0    # or: ~ Normal(0, deg2rad(1))
    end)
```
with `flux_H` declared on `b` and on `c`.
"""
struct ImageObs{TTable<:Table,TTargets<:Tuple,TRef} <: Octofitter.AbstractObs
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    targets::TTargets
    ref::TRef
    band::Union{Nothing,Symbol}
    name::String
end

function ImageObs(
    observations;
    targets,
    ref,
    band::Union{Nothing,Symbol}=nothing,
    name::String="images",
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin end)
)
    (priors, derived) = variables
    _reject_obs_flux_variable(priors, derived, name)

    # Sources first: measuring contrast from a stack of images is not cheap,
    # and a mis-spelled target should not have to wait for it.
    tspecs = map(t -> _image_target(Octofitter.refspec(t), band), _astuple(targets))
    isempty(tspecs) && error(
        "`targets` is empty: an `ImageObs` with no source in it has no likelihood. " *
        "Name the bodies modelled in these images, e.g. targets=(b, c).")
    r = Octofitter.refspec(ref)

    table = Table(observations)
    if !issubset(images_cols, Tables.columnnames(table))
        error("Expected columns $images_cols")
    end
    if any(>=(mjd("2050")), table.epoch) || any(<=(mjd("1950")), table.epoch)
        @warn "Epochs fell outside 1950–2050; the expected format is MJD. Double-check your input."
    end

    # Fall back to measuring the contrast from the images themselves.
    if !in(:contrast, Tables.columnnames(table)) && !in(:contrastmap, Tables.columnnames(table))
        @info "Measuring contrast from image"
        table = Table(table, contrast=contrast_interp.(table.image))
    end

    # Linear interpolators over the images, and over the contrast maps if
    # given. Positions outside the data extrapolate to NaN, which `ln_like`
    # reads as "no information here" — see there.
    imageinterp = map(table.image) do img
        linear_interpolation(parent.(dims(img)), img, extrapolation_bc=convert(eltype(img), NaN))
    end
    table = Table(table; imageinterp)
    if hasproperty(table, :contrastmap)
        contrastmapinterp = map(table.contrastmap) do img
            linear_interpolation(parent.(dims(img)), img, extrapolation_bc=convert(eltype(img), NaN))
        end
        table = Table(table; contrastmapinterp)
    end

    return ImageObs{typeof(table),typeof(tspecs),typeof(r)}(
        table, priors, derived, tspecs, r, band, String(name))
end

export ImageObs

@inline _astuple(t::Tuple) = t
_astuple(v::AbstractVector) = Tuple(v)
@inline _astuple(x) = (x,)

# A target must be something that can *emit* in this band. A body can; an
# unresolved pair can (its flux is the sum of its members', and its position
# their photocentre); a barycentre is a dynamical point, not a light source.
@inline _image_target(spec::Octofitter.BodyRefSpec, band) = spec
# A photocentre target inherits the image's band unless it names one itself,
# so `targets=(Photocentre((Aa, Ab)),)` weights by the same fluxes the
# likelihood reads. A different band would weight the position by one band
# and scale the flux by another.
@inline _image_target(::Octofitter.PhotocentreSpec{nothing,Names}, band) where {Names} =
    Octofitter.PhotocentreSpec{band,Names}()
@inline function _image_target(spec::Octofitter.PhotocentreSpec{B,Names}, band) where {B,Names}
    B === band || error(
        "an `ImageObs` target photocentre is weighted by band :$B but the image's " *
        "band is $(band === nothing ? "unset" : ":$band"). Give both the same band.")
    return spec
end
@noinline _image_target(spec::Octofitter.BarycentreSpec, band) = error(
    "$(Octofitter._refstr(spec)) cannot be an `ImageObs` target: a barycentre has no " *
    "flux, so there is nothing for the image to have detected. Name the bodies " *
    "themselves, or a `Photocentre` if the imager does not resolve them.")

function _reject_obs_flux_variable(priors, derived, name)
    vars = vcat(collect(keys(priors.priors)), collect(keys(derived.variables)))
    bad = filter(v -> v === :flux || startswith(string(v), "flux_"), vars)
    isempty(bad) && return nothing
    error("""
    `$(join(bad, "`, `"))` was declared on ImageObs "$name", but an image's
    fluxes are per-body variables in v9, not per-observation ones: one image
    holds every companion, so one flux on the observation cannot describe it.

    Move each one to the body it belongs to and give it the image's band:

        b = Body(name="b", about=A, variables=@variables begin
                flux_H ~ Normal(3.8, 0.5)
                ...
            end)
        ImageObs(dat; targets=(b,), ref=A, band=:H, name="$name")
    """)
end

Octofitter.refspecs(obs::ImageObs) = (obs.targets..., obs.ref)

Octofitter.likeobj_from_epoch_subset(obs::ImageObs, inds) = ImageObs(
    obs.table[inds]; targets=obs.targets, ref=obs.ref, band=obs.band, name=obs.name,
    variables=(obs.priors, obs.derived))

# ---------------------------------------------------
# Fluxes
#
# `PlanetOrbits.fluxes(sys, band)` is a body-vector built from the bodies'
# own `flux_<band>` variables — the same numbers `photocentre` weights by.
# Reading them here rather than from `θ_obs` is what lets one image hold
# several companions, and what makes "b is 3× fainter than c in H" a
# statement the whole model shares.
# ---------------------------------------------------

@inline _bandfluxes(sys, band::Symbol) = PlanetOrbits.fluxes(sys, band)
@inline function _bandfluxes(sys, ::Nothing)
    fl = PlanetOrbits.fluxes(sys)
    length(fl) == 1 || error(
        isempty(fl) ?
        "no body in this system declares a flux, so an `ImageObs` has nothing to fit. " *
        "Give each target a `flux_<band>` variable." :
        "several bands are defined ($(keys(fl))); pass `band=` to `ImageObs` to say " *
        "which one these images are in.")
    return @inbounds values(fl)[1]
end

# Flux of one *source*, which is a body or an unresolved set of them. Written
# over the spec rather than the resolved reference because a `WeightedPoint`
# has already normalized its weights — the total flux is not recoverable from
# it.
@inline _sourceflux(sys, fl, ::Octofitter.BodyRefSpec{N}) where {N} =
    @inbounds fl[getproperty(PlanetOrbits.bodies(sys), N).idx]
@inline _sourceflux(sys, fl, ::Octofitter.PhotocentreSpec{B,()}) where {B} = sum(fl)
@inline _sourceflux(sys, fl, ::Octofitter.PhotocentreSpec{B,Names}) where {B,Names} =
    sum(map(n -> (@inbounds fl[getproperty(PlanetOrbits.bodies(sys), n).idx]), Names))

@inline function _sourcefluxes(ctx::Octofitter.ObsContext, obs::ImageObs)
    fl = _bandfluxes(ctx.system, obs.band)
    return map(spec -> _sourceflux(ctx.system, fl, spec), obs.targets)
end

# ---------------------------------------------------
# The likelihood
# ---------------------------------------------------

"""
    simulate(obs::ImageObs, ctx) -> (; x, y, flux, epochs)

Where this observation's targets are predicted to land, in each image's own
pixel coordinates, and how bright they are. `x` and `y` are
`(n_targets, n_epochs)` matrices in target-then-epoch order; `flux` is one
value per target.

The plate scale and north angle *are* applied here — unlike
`RelAstromObs.simulate`, which leaves its model on the true sky. An image is
never corrected onto the sky (that would resample the pixels and their
noise), so this observation carries its model onto the detector instead. See
`Octofitter.sky_offset`.
"""
function Octofitter.simulate(obs::ImageObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    platescale_mult, northangle = Octofitter.sky_calibration(ctx)
    targets = Octofitter.resolverefs(ctx, obs.targets)
    reference = Octofitter.ref(ctx, obs.ref)
    flux = _sourcefluxes(ctx, obs)
    L = length(obs.table.epoch)
    N = length(targets)
    x = Matrix{T}(undef, N, L)
    y = Matrix{T}(undef, N, L)
    for i in 1:L
        sol = Octofitter.solutionat(ctx, i)
        ps = obs.table.platescale[i] * platescale_mult
        for j in 1:N
            x[j, i], y[j, i] = pixel_position(sol, targets[j], reference, ps, northangle)
        end
    end
    return (; x, y, flux, epochs=obs.table.epoch)
end

"""
    ln_like(obs::ImageObs, ctx)

Matched-filter log-likelihood of the modelled companions in these images.

# The estimator, and why it is written as a quadratic form

Per image, with `N` modelled sources of flux `f`, this is

    ll = −½ (fᵀ A f − 2 fᵀ b)

where `b_j` is the matched-filter output at source `j`'s modelled position
and `A` is the Gram matrix of the source templates after whatever linear
operator the reduction applied. That is the general form of the profiled
linear estimator of Ruffio et al. 2017 (eq. 31) / Mawet et al. 2019 (eq. 8),
with the data-only term dropped as an additive constant.

The reduction assumed here — one image, per-position 1σ contrast `σₓ`, no
modelling of how one companion's PSF leaks into another's aperture — makes
`A` diagonal, `A_jj = 1/σₓ_j²` and `b_j = f̃ₓ_j/σₓ_j²`, and the expression
collapses to a sum of independent per-companion terms

    ll_i = −1/(2σₓ²) (f² − 2 f f̃ₓ)

which is exactly what v8 computed for its one companion. **The independence
comes from the estimator, not from the images**: a joint reduction (T-LOCI,
a forward-modelled matched filter) keeps `A` full, and its off-diagonal
entries are precisely the PSF-overlap coupling between companions. That case
is not implemented, but it is a different `A` rather than a different
likelihood — which is why the diagonal is written as a diagonal instead of
being folded away. The diagonal path never materializes `A` and allocates
nothing.

# Off the detector
A position that extrapolates outside the image gets `f̃ₓ = NaN` from the
interpolator, which is read as no flux there (`f̃ₓ = 0`), leaving the term
`−f²/2σₓ²`: a companion the images should have shown and did not is evidence
*against* a bright companion, not an absent constraint. A non-finite or zero
contrast is different — that is a σ the data cannot supply — and makes the
whole sample impossible, as in v1.
"""
function Octofitter.ln_like(obs::ImageObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    tbl = obs.table
    platescale_mult, northangle = Octofitter.sky_calibration(ctx)

    # Resolved once, outside the epoch loop: each of these is a name lookup
    # that constant-folds, and inside the loop it would not.
    targets = Octofitter.resolverefs(ctx, obs.targets)
    reference = Octofitter.ref(ctx, obs.ref)
    flux = _sourcefluxes(ctx, obs)

    ll = zero(T)
    for i in eachindex(tbl.epoch)
        sol = Octofitter.solutionat(ctx, i)
        ps = tbl.platescale[i] * platescale_mult
        # (f̃ₓ, σₓ) per source. A tuple, so `N` is a compile-time constant and
        # nothing here reaches the heap.
        readout = map(targets) do t
            x, y = pixel_position(sol, t, reference, ps, northangle)
            f̃ = tbl.imageinterp[i](x, y)
            σ = _contrast_at(tbl, i, x, y)
            # Outside the image the interpolator extrapolates to NaN. v1 read
            # that as zero flux measured there, and so do we.
            return (isfinite(f̃) ? f̃ : zero(f̃), σ)
        end
        # A contrast the data cannot supply is not a weak constraint, it is a
        # missing one; v1 rejected the whole sample and so do we.
        any(rd -> !isfinite(rd[2]) || iszero(rd[2]), readout) && return convert(T, -Inf)
        ll += _matched_filter_ll(flux, readout)
    end
    return ll
end

# −½ (fᵀ A f − 2 fᵀ b) for diagonal `A`. Kept as its own function so the
# full-`A` case has somewhere to go: the caller supplies `A` and `b`, not a
# per-companion log-likelihood.
@inline function _matched_filter_ll(f::NTuple{N,Any}, readout::NTuple{N,Any}) where {N}
    return sum(ntuple(Val(N)) do j
        f̃, σ = readout[j]
        A = inv(σ^2)          # the j-th diagonal entry of A
        b = f̃ * A             # the j-th matched-filter output
        return -(A * f[j]^2 - 2 * f[j] * b) / 2
    end)
end

# 2-D contrast map if there is one, 1-D contrast curve otherwise. The branch
# is on the table's column names, so it folds away.
@inline function _contrast_at(tbl, i, x, y)
    if hasproperty(tbl, :contrastmap)
        return tbl.contrastmapinterp[i](x, y)
    else
        return tbl.contrast[i](hypot(x, y))
    end
end

# ---------------------------------------------------
# Simulating images
# ---------------------------------------------------

"""
    generate_from_params(obs::ImageObs, ctx; add_noise)

A replicate image set drawn from the model — for posterior-predictive checks
and simulation-based calibration.

The replicate is generated in the space this likelihood actually reads: each
image is treated as the *matched-filter output map* the estimator above
assumes, in which the expected value at a position is the flux of a source
there and the uncertainty is that position's contrast. So each replicate is
a fresh map holding the modelled sources, and nothing else — the original
pixels are not reused.

Each source is injected as a point: the four pixels bracketing its modelled
position are raised so that interpolating there returns exactly its flux.
With `add_noise=false` the replicate therefore reproduces its own generating
parameters exactly. With `add_noise=true` every pixel gets an independent
draw of scale `σₓ`, the same contrast the likelihood reads.

!!! note
    This is a simulator for the likelihood, not for a coronagraph: a real
    matched-filter map has a PSF-shaped response and noise correlated on that
    scale, and neither is reproduced here. A source predicted outside the
    image is simply absent from the replicate.
"""
function Octofitter.generate_from_params(obs::ImageObs, ctx::Octofitter.ObsContext; add_noise)
    sim = Octofitter.simulate(obs, ctx)
    images = map(eachindex(obs.table.epoch)) do i
        img = obs.table.image[i]
        E = eltype(img)
        out = copy(img)
        xs, ys = parent.(dims(img))
        # Noise first, then the sources, so an injected source is not itself
        # perturbed away from the flux it was given.
        if add_noise
            for (jx, x) in enumerate(xs), (jy, y) in enumerate(ys)
                σ = _contrast_at(obs.table, i, x, y)
                out[jx, jy] = isfinite(σ) ? E(randn() * σ) : E(NaN)
            end
        else
            fill!(out, zero(E))
        end
        for j in eachindex(sim.flux)
            _inject_point!(out, xs, ys, sim.x[j, i], sim.y[j, i], E(sim.flux[j]))
        end
        return out
    end
    return ImageObs(Table(obs.table; image=images); targets=obs.targets, ref=obs.ref,
        band=obs.band, name=obs.name, variables=(obs.priors, obs.derived))
end

# Raise the four pixels bracketing (x, y) so that bilinear interpolation
# there returns `f`. The bilinear weights `w` sum to one, so adding `c·w`
# changes the interpolated value by `c·Σw²`; scaling by the inverse puts the
# requested flux exactly at the requested sub-pixel position.
function _inject_point!(out, xs, ys, x, y, f)
    (isfinite(x) && isfinite(y)) || return out
    (first(xs) <= x <= last(xs) && first(ys) <= y <= last(ys)) || return out
    ix = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    iy = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
    tx = (x - xs[ix]) / (xs[ix+1] - xs[ix])
    ty = (y - ys[iy]) / (ys[iy+1] - ys[iy])
    w = ((1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty)
    c = f / sum(wk -> wk^2, w)
    out[ix, iy] += c * w[1]
    out[ix+1, iy] += c * w[2]
    out[ix, iy+1] += c * w[3]
    out[ix+1, iy+1] += c * w[4]
    return out
end

# ---------------------------------------------------
# Contrast measurement (unchanged from v1)
# ---------------------------------------------------

"""
    contrast_interp(image; step=2)

Linear interpolation over the results of `contrast`, flat outside the
measured range. This is the callable an `ImageObs` `contrast` column holds:
separation in pixels → 1σ flux uncertainty.
"""
function contrast_interp(image::AstroImage; step=2)
    cont = contrast(image; step)
    mask = findfirst(isfinite, cont.contrast):findlast(isfinite, cont.contrast)
    return linear_interpolation(cont.separation[mask], cont.contrast[mask], extrapolation_bc=Flat())
end

"""
    contrast(image; step=2)

Measure the contrast of an image, in the sense of high contrast imaging.
That is, divide the image into annuli moving outwards from the centre
(index 0,0 if offset image) and calculate the standard deviation in
each.

Returns a vector of annulus locations in pixels and a vector of standard
deviations.

*NOTE* This is the 1σ contrast. Multiply by five to get the usual confidence
value.
"""
function contrast(image::AstroImage; step=2)
    dx = dims(image, X)
    dy = collect(dims(image, Y))
    dr = sqrt.(
        dx .^ 2 .+ (dy') .^ 2
    )

    c_img = collect(image)

    bins = 0:step:maximum(dr)
    contrast = zeros(size(bins))
    mask = falses(size(image))
    mask2 = isfinite.(c_img)
    for i in eachindex(bins)
        bin = bins[i]
        mask .= (bin .- step / 2) .< dr .< (bin .+ step / 2)
        mask .&= mask2
        c = std(view(c_img, mask))
        contrast[i] = c
    end

    return (; separation=bins, contrast)
end

function imgsep(image::AstroImage)
    dx = dims(image, X)
    dy = collect(dims(image, Y))
    dr = sqrt.(
        dx .^ 2 .+ (dy') .^ 2
    )
    return dr
end

# ---------------------------------------------------
# Plotting protocol
#
# One channel per modelled source: the matched-filter readout at that source's
# predicted position (the "data"), against the source's fitted flux (the
# "model"), per image. That is a genuine epoch series — one number per image
# per source — even though the underlying datum is a 2-D frame, and it is the
# quantity the estimator actually scores: the diagonal path's per-companion
# term is −(f² − 2 f f̃)/2σ², whose χ² is ((f̃ − f)/σ)² up to the additive
# constant f̃²/2σ² that `ln_like` drops.
#
# There is no smooth model curve: flux is epoch-independent and no
# `PlanetOrbits` observable produces it, so `query=nothing` and the panel shows
# model points, residuals and the marginal histogram (the same shape as the
# Gaia along-scan channel).
# ---------------------------------------------------

# A channel name per target. `_refstr` renders `Photocentre(:H, (Aa, Ab))` with
# punctuation, which is not a legal field name, so anything that is not a plain
# identifier falls back to its position in `targets`.
function _img_channelname(spec, k::Integer)
    s = Octofitter._refstr(spec)
    return Base.isidentifier(s) ? Symbol(:flux_, s) : Symbol(:flux_, k)
end

function Octofitter.plotchannels(obs::ImageObs)
    band = isnothing(obs.band) ? "" : " (" * string(obs.band) * ")"
    return ntuple(length(obs.targets)) do k
        Octofitter.PlotChannel(_img_channelname(obs.targets[k], k),
            "flux of " * Octofitter._refstr(obs.targets[k]) * band, "")
    end
end

function Octofitter.residuals(obs::ImageObs, ctx::Octofitter.ObsContext)
    sim = Octofitter.simulate(obs, ctx)
    tbl = obs.table
    ep = collect(Float64, tbl.epoch)
    L = length(ep)
    out = NamedTuple()
    for (k, spec) in enumerate(obs.targets)
        f = Float64(sim.flux[k])
        data = Vector{Float64}(undef, L)
        σ = Vector{Float64}(undef, L)
        use = trues(L)
        for i in 1:L
            x, y = Float64(sim.x[k, i]), Float64(sim.y[k, i])
            f̃ = tbl.imageinterp[i](x, y)
            # Off the detector the interpolator gives NaN, which `ln_like`
            # reads as no flux measured there — a real (and constraining)
            # zero, not a missing point, so `use` stays true.
            data[i] = isfinite(f̃) ? f̃ : 0.0
            s = Float64(_contrast_at(tbl, i, x, y))
            # A σ the data cannot supply makes the whole sample impossible in
            # `ln_like`; here it just means this image says nothing.
            σ[i] = s
            use[i] = isfinite(s) && !iszero(s)
        end
        σ[.!use] .= NaN
        out = merge(out, NamedTuple{(_img_channelname(spec, k),)}((
            (; epoch=ep, data, model=fill(f, L), resid=data .- f,
                σ, σ_eff=copy(σ), use),)))
    end
    return out
end
