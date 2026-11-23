# ---------------------------------------------------
# The observation: raw DR4 epoch products of one source
#
# One scene per transit (where each emitter is, how fast it moves along the
# line of sight, what its spectrum is), from which every block predicts its
# own data: the RVS epoch spectra and the epoch images.
# ---------------------------------------------------

using Octofitter: AbstractObs, ObsContext, Priors, Derived, Table, solutionat

const RVS_COLUMNS = (:epoch, :scan_pos_angle, :parallax_factor_al, :parallax_factor_ac,
    :ref_pos_al, :ref_pos_ac, :flux, :flux_error)

# The RVS epoch spectra of the source and what places the emitters in them.
struct RVSBlock{TTable<:Table,TR<:RVSEpochSpectra}
    table::TTable
    spectra::TR
    membership::Matrix{Float64}
    sinψ::Vector{Float64}
    cosψ::Vector{Float64}
end

# The epoch images of the source; `functionals[k]` turns emitter k's spectrum
# into its G photon rate and chromatic moments.
struct ImageBlock{TTable<:Table,TF<:Tuple,TD<:EpochImages}
    table::TTable
    data::TD
    functionals::TF
    sinψ::Vector{Float64}           # per window
    cosψ::Vector{Float64}
end

"""
    GaiaDR4RawObs(; emitters, rvs = nothing, images = nothing, name = "GaiaDR4Raw",
                  tangent_point = nothing, variables = @variables begin end,
                  rvs_instrument = RVSInstrument(), wavelength = DR4_RVS_WAVELENGTH, rvs_p_bad = 1e-3,
                  mismatch = rvs_instrument.mismatch_amplitude, membership = nothing,
                  image_instrument = nothing, image_p_bad = 1e-3)

The raw Gaia DR4 epoch products of one source, modelled from the scene of
the system's emitters: the RVS epoch spectra (`rvs`), the epoch images
(`images`), or both.

`emitters` names the bodies that emit light and what they are, for example
`(A = Star(evol, sed, rvs), b = Star(evol, sed, rvs))`; the keys are body
names. With RVS spectra, all must share one RVS wavelength grid.

The tangent point comes from a table's `ra0`/`dec0` columns or from
`tangent_point = (ra0, dec0)` [deg], as for `GaiaDR4AstromObs`.

The system must have an absolute frame (`plx`, `ra`, `dec`, `pmra`, `pmdec`,
`rv`, `ref_epoch`), which is Gaia DR4's. Every block places an emitter at
the same local-plane position, (w, z) of [`transit_position`](@ref).

# RVS epoch spectra

`rvs` is a table with one row per transit:

  - `epoch` [MJD], `scan_pos_angle` [**degrees**], `parallax_factor_al`,
    `parallax_factor_ac`: the transit, as in `epoch_astrometry`;
  - `ref_pos_al`, `ref_pos_ac` [mas]: where the pipeline took the source to
    be when it extracted the spectrum and set its wavelength scale, in the
    same local-plane coordinates as `centroid_pos_al` (see
    [`catalogue_track`](@ref)). The window is taken to be centred there;
  - `flux`, `flux_error`: the normalised spectrum of the transit, a vector
    per row on the bins `wavelength`, `NaN` where there is no data (see
    [`RVSEpochSpectra`](@ref), which also documents `rvs_instrument`,
    `rvs_p_bad` and `mismatch`). Rows must be in order of epoch.

`membership[emitter, transit]` is 1 where an emitter's light is in the
spectrum as the window caught it (the default everywhere) and 0 where the
pipeline subtracted it, the emitter having had its own window (a deblended
transit).

Per transit and emitter:

    velocity = frame rv + orbit + redshift + drv + κ (w − ref_pos_al) + κ_ac (z − ref_pos_ac)
    weight   = fraction of the light inside the window at z − ref_pos_ac, times membership

So the systemic velocity is the frame's `rv`, the one that drives
perspective acceleration, and a companion displaced along scan appears
Doppler-shifted by `κ` per mas.

# Epoch images

`images` is a table with one row per window (one CCD of one transit), as
[`EpochImages`](@ref) documents, and `image_instrument` the
[`ImageInstrument`](@ref) that holds the devices' image bases. Per window,
every point of every emitter is imaged at its own position, with its own G
photon rate times the window's exposure and its own chromatic moments, all
from the emitter's spectrum.

# Variables

Read from each emitter's body: whatever its emitter needs (`teff`, `logg`,
`mh`, `radius` for a [`Star`](@ref)); and for the RVS spectra `vsini`
[km/s], the rotational broadening (none if absent), and `drv` [**m/s**, as
every radial velocity in Octofitter], the zero point of that body's
template (zero if absent).

Of the observation itself, all optional:

  - `kappa_ac` [km s⁻¹ mas⁻¹], the across-scan term of the RVS wavelength
    scale, and `deblend_residual`, the fraction of a subtracted emitter's
    light that is left in a deblended spectrum;
  - `flux_scale`, a factor on every emitter's electrons in the images (the
    absolute flux calibration of the model against the windows);
  - `psf1`, `psf2`, …: coefficients of the devices' shape modes that are
    common to all the source's windows;
  - `jitter_al`, `jitter_ac` [mas]: the width of the per-transit shift of
    the scene, instead of the instrument's `transit_shift_sigma`.
"""
struct GaiaDR4RawObs{Names,TE<:Tuple,TB<:Tuple,TTable<:Table,TR<:Union{Nothing,RVSBlock},TI<:Union{Nothing,ImageBlock}} <: AbstractObs
    table::TTable                   # one row per row of the blocks' tables: `epoch`, `block`, `row`
    rvs::TR
    images::TI
    emitters::NamedTuple{Names,TE}
    bodies::TB
    tangent_point::NTuple{2,Float64}
    tangent_basis::NTuple{9,Float64}
    priors::Priors
    derived::Derived
    name::String
end

function _rvs_block(rvs, emitters, instrument, wavelength, p_bad, mismatch, membership)
    table = Table(rvs)
    issubset(RVS_COLUMNS, Octofitter.Tables.columnnames(table)) ||
        error("GaiaDR4RawObs: the `rvs` table needs the columns $RVS_COLUMNS")
    issorted(table.epoch) || error("GaiaDR4RawObs: the rows of `rvs` must be in order of epoch")
    grids = map(rvs_grid, values(emitters))
    all(==(first(grids)), grids) || error("GaiaDR4RawObs: the emitters' RVS spectra are on different wavelength grids")
    nt, ne = length(table.epoch), length(emitters)
    membership = membership === nothing ? ones(ne, nt) : Matrix{Float64}(membership)
    size(membership) == (ne, nt) || throw(DimensionMismatch("membership must be [emitter, transit] = ($ne, $nt)"))
    all(m -> 0 <= m <= 1, membership) || throw(ArgumentError("membership must be in [0, 1]"))
    all(>(0), sum(membership; dims=1)) || throw(ArgumentError("a transit has no emitter in it"))
    flux = reduce(hcat, table.flux)
    flux_error = reduce(hcat, table.flux_error)
    spectra = RVSEpochSpectra(flux, flux_error, first(grids); instrument, wavelength, p_bad, mismatch)
    return RVSBlock(table, spectra, membership, collect(Float64, sind.(table.scan_pos_angle)), collect(Float64, cosd.(table.scan_pos_angle)))
end

function _image_block(images, emitters, instrument, p_bad)
    instrument isa ImageInstrument || error("GaiaDR4RawObs: `images` needs `image_instrument`, an ImageInstrument")
    table = Table(images)
    data = EpochImages(table, instrument; p_bad)
    issorted(table.epoch) || error("GaiaDR4RawObs: the rows of `images` must be in order of epoch")
    degree = chromatic_degree(instrument)
    functionals = map(e -> SEDFunctionals(sed_bin_edges(e); degree), values(emitters))
    for f in functionals
        collect(f.wavenumber_range) ≈ collect(first(values(instrument.bases)).wavenumber_range) ||
            error("GaiaDR4RawObs: the image bases expand over a different range of wavenumber than the G passband")
    end
    return ImageBlock(table, data, functionals, collect(Float64, sind.(table.scan_pos_angle)), collect(Float64, cosd.(table.scan_pos_angle)))
end

function GaiaDR4RawObs(; emitters::NamedTuple, rvs=nothing, images=nothing, name="GaiaDR4Raw", tangent_point=nothing,
    variables::Tuple{Priors,Derived}=(Priors(), Derived()),
    rvs_instrument::RVSInstrument=RVSInstrument(), wavelength=DR4_RVS_WAVELENGTH, rvs_p_bad=1e-3,
    mismatch=rvs_instrument.mismatch_amplitude, membership=nothing,
    image_instrument=nothing, image_p_bad=1e-3)
    isempty(emitters) && error("GaiaDR4RawObs: no emitters")
    (rvs === nothing && images === nothing) && error("GaiaDR4RawObs: no data; pass `rvs`, `images` or both")
    rvs_block = rvs === nothing ? nothing : _rvs_block(rvs, emitters, rvs_instrument, wavelength, rvs_p_bad, mismatch, membership)
    image_block = images === nothing ? nothing : _image_block(images, emitters, image_instrument, image_p_bad)
    tps = [Octofitter._dr4_tangent_point(b.table, tangent_point) for b in (rvs_block, image_block) if b !== nothing]
    any(isnothing, tps) && error("GaiaDR4RawObs \"$name\": no tangent point. Keep the tables' `ra0`/`dec0` columns " *
                                 "or pass `tangent_point = (ra0, dec0)` [deg].")
    allequal(tps) || error("GaiaDR4RawObs \"$name\": the `rvs` and `images` tables have different tangent points")
    tp = first(tps)
    nr = rvs_block === nothing ? 0 : length(rvs_block.table.epoch)
    ni = image_block === nothing ? 0 : length(image_block.table.epoch)
    block_epochs(b) = b === nothing ? Float64[] : collect(Float64, b.table.epoch)
    table = Table(epoch=[block_epochs(rvs_block); block_epochs(image_block)], block=[fill(:rvs, nr); fill(:images, ni)], row=[1:nr; 1:ni])
    bodies = map(Octofitter.refspec, keys(emitters))
    priors, derived = variables
    return GaiaDR4RawObs(table, rvs_block, image_block, emitters, bodies, tp, tangent_basis(tp...), priors, derived, String(name))
end

Octofitter.refspecs(obs::GaiaDR4RawObs) = (obs.bodies..., Octofitter.Barycentre)

function Octofitter.check_siblings(obs::GaiaDR4RawObs, @nospecialize(all_obs), ctx)
    :ra in ctx.framevars || error("""
        System $(ctx.name): GaiaDR4RawObs "$(obs.name)" needs an absolute frame.

        The raw epoch products are placed in the sky and in velocity by the barycentre's \
        Gaia DR4 position, proper motion, parallax and radial velocity, which PlanetOrbits \
        propagates: define `plx`, `ra`, `dec`, `pmra`, `pmdec`, `rv` and `ref_epoch` in the \
        system block (see the "absolute frame" section of the Gaia DR4 tutorial for the \
        mas-offset form of `ra` and `dec`).
        """)
    return nothing
end

_nrvs(obs::GaiaDR4RawObs) = obs.rvs === nothing ? 0 : length(obs.rvs.table.epoch)

# The observation with its blocks' tables replaced by the given rows of the
# same tables (`nothing` drops a block), or by tables holding new data.
function _rebuild(obs::GaiaDR4RawObs, rvs_rows, image_rows; rvs_table=nothing, image_table=nothing)
    r, im = obs.rvs, obs.images
    keep_r = r !== nothing && rvs_rows !== nothing && !isempty(rvs_rows)
    keep_i = im !== nothing && image_rows !== nothing && !isempty(image_rows)
    (keep_r || keep_i) || error("GaiaDR4RawObs: the subset holds no data")
    kw = (; obs.emitters, obs.name, tangent_point=obs.tangent_point, variables=(obs.priors, obs.derived))
    if keep_r
        kw = (; kw..., rvs=something(rvs_table, r.table[rvs_rows]), rvs_instrument=r.spectra.instrument,
            wavelength=r.spectra.wavelength, membership=r.membership[:, rvs_rows],
            rvs_p_bad=r.spectra.p_bad[rvs_rows], mismatch=r.spectra.mismatch[:, rvs_rows])
    end
    if keep_i
        table = something(image_table, im.table[image_rows])
        # the transits of the rows that are kept, in order
        kept = unique(t for t in 1:ntransits(im.data) if any(in(image_rows), im.data.transit_windows[t]))
        kw = (; kw..., images=table, image_instrument=im.data.instrument, image_p_bad=im.data.p_bad[kept])
    end
    return GaiaDR4RawObs(; kw...)
end

function Octofitter.likeobj_from_epoch_subset(obs::GaiaDR4RawObs, inds)
    nr = _nrvs(obs)
    rows = collect(Int, (1:length(obs.table.epoch))[inds])
    return _rebuild(obs, [i for i in rows if i <= nr], [i - nr for i in rows if i > nr])
end

# ---------------------------------------------------
# The scene: RVS
# ---------------------------------------------------

"Velocity [km/s] an emitter's spectrum is offset by at its own surface; the gravitational redshift of a star."
rvs_velocity_offset(::AbstractEmitter, θ_body) = 0.0
rvs_velocity_offset(::Star, θ_body) = gravitational_redshift(θ_body.logg, θ_body.radius)

@inline _vsini(θ_body) = hasproperty(θ_body, :vsini) ? θ_body.vsini : 0.0
@inline _drv_kms(θ_body) = hasproperty(θ_body, :drv) ? θ_body.drv / 1000 : 0.0

# Prepare every emitter's spectrum once: `work[:, 4k-3:4k]` holds emitter k's
# unbroadened spectrum and the three arrays of its `RVSComponent`.
function _rvs_components(work::AbstractMatrix, obs::GaiaDR4RawObs{Names}, θ_system) where {Names}
    θ_bodies = values(NamedTuple{Names}(θ_system.bodies))
    spectra = obs.rvs.spectra
    grid = spectra.grid
    limb = spectra.instrument.limb_darkening
    ks = ntuple(identity, length(Names))
    return map(ks, values(obs.emitters), θ_bodies) do k, emitter, θ_body
        raw, flux, dflux, d2flux = (view(work, :, 4k - 3), view(work, :, 4k - 2), view(work, :, 4k - 1), view(work, :, 4k))
        rvs_spectrum!(raw, emitter, θ_body, θ_system)
        rvs_component!(flux, dflux, d2flux, raw, _vsini(θ_body), grid; limb_darkening=limb)
    end
end

# Velocities and weights of every emitter in transit `t`.
@inline function _rvs_scene(obs::GaiaDR4RawObs{Names}, ctx::ObsContext, t::Integer, refs::Tuple, barycentre) where {Names}
    θ_bodies = values(NamedTuple{Names}(ctx.θ_system.bodies))
    θ_obs = ctx.θ_obs
    blk = obs.rvs
    inst = blk.spectra.instrument
    tab = blk.table
    sol = solutionat(ctx, t)
    kappa_ac = hasproperty(θ_obs, :kappa_ac) ? θ_obs.kappa_ac : 0.0
    residual = hasproperty(θ_obs, :deblend_residual) ? θ_obs.deblend_residual : 0.0
    ks = ntuple(identity, length(Names))
    scene = map(ks, values(obs.emitters), θ_bodies, refs) do k, emitter, θ_body, body
        w, z = @inbounds transit_position(sol, body, barycentre, obs.tangent_basis, blk.sinψ[t], blk.cosψ[t],
            tab.parallax_factor_al[t], tab.parallax_factor_ac[t])
        Δη, Δζ = @inbounds (w - tab.ref_pos_al[t], z - tab.ref_pos_ac[t])
        v = transit_velocity(sol, body, barycentre, ctx.einstein_rv) + rvs_velocity_offset(emitter, θ_body) + _drv_kms(θ_body)
        m = @inbounds blk.membership[k, t]
        (rvs_velocity(inst, v, Δη; Δζ, kappa_ac), window_fraction(inst, Δζ) * (m + (1 - m) * residual))
    end
    return map(first, scene), map(last, scene)
end

function _rvs_loglike(obs::GaiaDR4RawObs{Names}, ctx::ObsContext, refs::Tuple, barycentre) where {Names}
    T = Octofitter._system_number_type(ctx.θ_system)
    spectra = obs.rvs.spectra
    @no_escape ctx.buf begin
        work = @alloc(T, length(spectra.grid), 4 * length(Names))
        components = _rvs_components(work, obs, ctx.θ_system)
        ll = zero(T)
        for t in 1:ntransits(spectra)
            velocities, weights = _rvs_scene(obs, ctx, t, refs, barycentre)
            ll += first(rvs_transit_loglike(spectra, t, components, velocities, weights; buf=ctx.buf))
        end
        ll
    end
end

# ---------------------------------------------------
# The scene: images
# ---------------------------------------------------

# Shape coefficients common to all windows of the source: `psf1`, `psf2`, … as far as they are defined.
@inline _shape(θ) = hasproperty(θ, :psf1) ? (θ.psf1, _shape2(θ)...) : ()
@inline _shape2(θ) = hasproperty(θ, :psf2) ? (θ.psf2, _shape3(θ)...) : ()
@inline _shape3(θ) = hasproperty(θ, :psf3) ? (θ.psf3, _shape4(θ)...) : ()
@inline _shape4(θ) = hasproperty(θ, :psf4) ? (θ.psf4,) : ()

@inline function _jitter(θ_obs, inst::ImageInstrument)
    (hasproperty(θ_obs, :jitter_al) || hasproperty(θ_obs, :jitter_ac)) || return nothing
    return (hasproperty(θ_obs, :jitter_al) ? θ_obs.jitter_al : inst.transit_shift_sigma[1],
        hasproperty(θ_obs, :jitter_ac) ? θ_obs.jitter_ac : inst.transit_shift_sigma[2])
end

# Number of points the emitters put in the scene.
_npoints(emitters::Tuple, θ_bodies::Tuple, θ_system, epoch) =
    sum(map((e, θ_body) -> length(scene_points(e, θ_body, θ_system, epoch)), emitters, θ_bodies))

# G photon rate and chromatic moments of every point: `rate[k]`, `moments[:, k]`, from point `k + 1` on.
@inline _fill_spectra!(rate, moments, k, ::Tuple{}, ::Tuple{}, ::Tuple{}, θ_system, epoch, buf) = k
@inline function _fill_spectra!(rate, moments, k, emitters::Tuple, θ_bodies::Tuple, functionals::Tuple, θ_system, epoch, buf)
    emitter, θ_body, f = first(emitters), first(θ_bodies), first(functionals)
    @no_escape buf begin
        sed = @alloc(eltype(rate), nbins(emitter))
        sed!(sed, emitter, θ_body, θ_system)
        k0 = k + 1
        total = chromatic_moments!(view(moments, :, k0), f, sed)
        for (_, _, weight) in scene_points(emitter, θ_body, θ_system, epoch)
            k += 1
            rate[k] = weight * total
            for n in axes(moments, 1)
                moments[n, k] = moments[n, k0]
            end
        end
        nothing
    end
    return _fill_spectra!(rate, moments, k, Base.tail(emitters), Base.tail(θ_bodies), Base.tail(functionals), θ_system, epoch, buf)
end

# Local-plane positions of every point in one window: `w[k]`, `z[k]`, from point `k + 1` on.
@inline _fill_positions!(w, z, k, ::Tuple{}, ::Tuple{}, ::Tuple{}, args...) = k
@inline function _fill_positions!(w, z, k, emitters::Tuple, θ_bodies::Tuple, refs::Tuple, θ_system, epoch, sol, barycentre,
    basis, sinψ, cosψ, f_al, f_ac)
    emitter, θ_body = first(emitters), first(θ_bodies)
    wb, zb = transit_position(sol, first(refs), barycentre, basis, sinψ, cosψ, f_al, f_ac)
    for (Δα, Δδ, _) in scene_points(emitter, θ_body, θ_system, epoch)
        k += 1
        dw, dz = scan_coordinates(Δα, Δδ, sinψ, cosψ)
        @inbounds w[k], z[k] = wb + dw, zb + dz
    end
    return _fill_positions!(w, z, k, Base.tail(emitters), Base.tail(θ_bodies), Base.tail(refs), θ_system, epoch, sol, barycentre,
        basis, sinψ, cosψ, f_al, f_ac)
end

# Run `f(t, w, z, rate, moments)` for every transit of the image block, with
# the scene's points placed in each of its windows.
function _image_transits(f::F, obs::GaiaDR4RawObs{Names}, ctx::ObsContext, refs::Tuple, barycentre, init) where {F,Names}
    T = Octofitter._system_number_type(ctx.θ_system)
    blk = obs.images
    d = blk.data
    θ_system = ctx.θ_system
    θ_bodies = values(NamedTuple{Names}(θ_system.bodies))
    emitters = values(obs.emitters)
    epoch = blk.table.epoch
    row0 = _nrvs(obs)
    K = _npoints(emitters, θ_bodies, θ_system, first(epoch))
    @no_escape ctx.buf begin
        rate = @alloc(T, K)
        moments = @alloc(T, chromatic_degree(d.instrument), K)
        _fill_spectra!(rate, moments, 0, emitters, θ_bodies, blk.functionals, θ_system, first(epoch), ctx.buf)
        nmax = maximum(length, d.transit_windows)
        w, z = @alloc(T, K, nmax), @alloc(T, K, nmax)
        acc = init
        for t in 1:ntransits(d)
            wins = d.transit_windows[t]
            for (q, iwin) in enumerate(wins)
                _fill_positions!(view(w, :, q), view(z, :, q), 0, emitters, θ_bodies, refs, θ_system, epoch[iwin],
                    solutionat(ctx, row0 + iwin), barycentre, obs.tangent_basis, blk.sinψ[iwin], blk.cosψ[iwin], d.f_al[t], d.f_ac[t])
            end
            nw = length(wins)
            acc = f(acc, t, view(w, :, 1:nw), view(z, :, 1:nw), rate, moments)
        end
        acc
    end
end

@inline _flux_scale(θ_obs) = hasproperty(θ_obs, :flux_scale) ? θ_obs.flux_scale : 1.0

function _image_loglike(obs::GaiaDR4RawObs, ctx::ObsContext, refs::Tuple, barycentre)
    T = Octofitter._system_number_type(ctx.θ_system)
    d = obs.images.data
    shape, flux_scale, jitter = _shape(ctx.θ_obs), _flux_scale(ctx.θ_obs), _jitter(ctx.θ_obs, d.instrument)
    return _image_transits(obs, ctx, refs, barycentre, zero(T)) do ll, t, w, z, rate, moments
        ll + _image_transit_ll(d, t, w, z, rate, moments, shape, flux_scale, jitter, ctx.buf)
    end
end

# With first-order dual numbers the gradient is carried by hand (`image/adjoint.jl`); any other number type goes through the generic code.
@inline _image_transit_ll(d, t, w, z, rate, moments, shape, flux_scale, jitter, buf) =
    first(image_transit_loglike(d, t, w, z, rate, moments; shape, flux_scale, jitter, buf))
@inline _image_transit_ll(d, t, w::AbstractMatrix{<:ForwardDiff.Dual{<:Any,Float64}}, z, rate, moments, shape, flux_scale, jitter, buf) =
    image_transit_loglike_dual(d, t, w, z, rate, moments; shape, flux_scale, jitter, buf)

# ---------------------------------------------------
# Likelihood, and data from the model
# ---------------------------------------------------

function Octofitter.ln_like(obs::GaiaDR4RawObs, ctx::ObsContext)
    refs = Octofitter.resolverefs(ctx, obs.bodies)
    barycentre = Octofitter.ref(ctx, Octofitter.Barycentre)
    ll = zero(Octofitter._system_number_type(ctx.θ_system))
    obs.rvs === nothing || (ll += _rvs_loglike(obs, ctx, refs, barycentre))
    obs.images === nothing || (ll += _image_loglike(obs, ctx, refs, barycentre))
    return ll
end

"""
    Octofitter.simulate(obs::GaiaDR4RawObs, ctx) -> (; rvs, images)

The model's data at the prior mean of the marginalised modes; a block the
observation does not have is `nothing`.

  - `rvs = (; flux, velocity, weight, epochs)`: the epoch spectra
    `flux[bin, transit]`, and each emitter's velocity [km/s] and weight in
    each transit, `[emitter, transit]`;
  - `images = (; samples, w, z, rate)`: the expected electrons of each
    window, a matrix per row of the image table; the local-plane position
    [mas] of each point of the scene in each window, `[point, window]`; and
    each point's G photon rate [e⁻ s⁻¹].
"""
function Octofitter.simulate(obs::GaiaDR4RawObs, ctx::ObsContext)
    return (; rvs=obs.rvs === nothing ? nothing : _simulate_rvs(obs, ctx, (t, c, v, w) -> rvs_transit_model(obs.rvs.spectra, t, c, v, w)),
        images=obs.images === nothing ? nothing : _simulate_images(obs, ctx, nothing))
end

function _simulate_rvs(obs::GaiaDR4RawObs{Names}, ctx::ObsContext, spectrum::F) where {Names,F}
    T = Octofitter._system_number_type(ctx.θ_system)
    refs = Octofitter.resolverefs(ctx, obs.bodies)
    barycentre = Octofitter.ref(ctx, Octofitter.Barycentre)
    spectra = obs.rvs.spectra
    nt, ne = ntransits(spectra), length(Names)
    flux = Matrix{Float64}(undef, length(spectra.wavelength), nt)
    velocity, weight = Matrix{T}(undef, ne, nt), Matrix{T}(undef, ne, nt)
    work = Matrix{T}(undef, length(spectra.grid), 4ne)
    components = _rvs_components(work, obs, ctx.θ_system)
    for t in 1:nt
        v, w = _rvs_scene(obs, ctx, t, refs, barycentre)
        velocity[:, t] .= v
        weight[:, t] .= w
        flux[:, t] .= spectrum(t, components, v, w)
    end
    return (; flux, velocity, weight, epochs=obs.rvs.table.epoch)
end

# `rng === nothing`: the mean; otherwise a draw from the good state.
function _simulate_images(obs::GaiaDR4RawObs, ctx::ObsContext, rng)
    refs = Octofitter.resolverefs(ctx, obs.bodies)
    barycentre = Octofitter.ref(ctx, Octofitter.Barycentre)
    d = obs.images.data
    shape, flux_scale, jitter = _shape(ctx.θ_obs), _flux_scale(ctx.θ_obs), _jitter(ctx.θ_obs, d.instrument)
    samples = Matrix{Float64}[]
    ws, zs = Vector{Float64}[], Vector{Float64}[]
    rates = Float64[]
    _image_transits(obs, ctx, refs, barycentre, nothing) do _, t, w, z, rate, moments
        append!(samples, rng === nothing ? image_transit_model(d, t, w, z, rate, moments; shape, flux_scale) :
                         image_transit_simulate(rng, d, t, w, z, rate, moments; shape, flux_scale, jitter))
        for q in axes(w, 2)
            push!(ws, collect(Float64, view(w, :, q)))
            push!(zs, collect(Float64, view(z, :, q)))
        end
        isempty(rates) && append!(rates, flux_scale .* rate)
        nothing
    end
    return (; samples, w=reduce(hcat, ws), z=reduce(hcat, zs), rate=rates)
end

function Octofitter.generate_from_params(obs::GaiaDR4RawObs, ctx::ObsContext; add_noise)
    rng = Octofitter.Random.default_rng()
    rvs_table = image_table = nothing
    if obs.rvs !== nothing
        spectra = obs.rvs.spectra
        sim = _simulate_rvs(obs, ctx, (t, c, v, w) -> rvs_transit_simulate(rng, spectra, t, c, v, w; add_noise))
        # bins without data stay without data
        flux = [ifelse.(isfinite.(spectra.flux[:, t]) .& isfinite.(spectra.flux_error[:, t]), sim.flux[:, t], NaN) for t in axes(sim.flux, 2)]
        rvs_table = Table(obs.rvs.table; flux)
    end
    if obs.images !== nothing
        image_table = Table(obs.images.table; samples=_simulate_images(obs, ctx, add_noise ? rng : nothing).samples)
    end
    return _rebuild(obs, 1:_nrvs(obs), obs.images === nothing ? nothing : (1:nwindows(obs.images.data)); rvs_table, image_table)
end
