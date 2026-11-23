# ---------------------------------------------------
# The observation: raw DR4 epoch products of one source
#
# One scene per transit (where each emitter is, how fast it moves along the
# line of sight, what its spectrum is), from which every block predicts its
# own data. So far the block is the RVS epoch spectra.
# ---------------------------------------------------

using Octofitter: AbstractObs, ObsContext, Priors, Derived, Table, solutionat

const RVS_COLUMNS = (:epoch, :scan_pos_angle, :parallax_factor_al, :parallax_factor_ac,
    :ref_pos_al, :ref_pos_ac, :flux, :flux_error)

"""
    GaiaDR4RawObs(; rvs, emitters, name = "GaiaDR4Raw", tangent_point = nothing,
                  instrument = RVSInstrument(), wavelength = DR4_RVS_WAVELENGTH, p_bad = 1e-3,
                  mismatch = instrument.mismatch_amplitude, membership = nothing,
                  variables = @variables begin end)

The raw Gaia DR4 epoch products of one source, modelled from the scene of
the system's emitters. So far: the RVS epoch spectra.

`emitters` names the bodies that emit light and what they are, for example
`(A = Star(evol, sed, rvs), b = Star(evol, sed, rvs))`; the keys are body
names. All must share one RVS wavelength grid.

`rvs` is a table with one row per transit:

  - `epoch` [MJD], `scan_pos_angle` [**degrees**], `parallax_factor_al`,
    `parallax_factor_ac`: the transit, as in `epoch_astrometry`;
  - `ref_pos_al`, `ref_pos_ac` [mas]: where the pipeline took the source to
    be when it extracted the spectrum and set its wavelength scale, in the
    same local-plane coordinates as `centroid_pos_al` (see
    [`catalogue_track`](@ref)). The window is taken to be centred there;
  - `flux`, `flux_error`: the normalised spectrum of the transit, a vector
    per row on the bins `wavelength`, `NaN` where there is no data (see
    [`RVSEpochSpectra`](@ref), which also documents `instrument`, `p_bad`
    and `mismatch`). Rows must be in order of epoch.

The tangent point comes from the table's `ra0`/`dec0` columns or from
`tangent_point = (ra0, dec0)` [deg], as for `GaiaDR4AstromObs`.

`membership[emitter, transit]` is 1 where an emitter's light is in the
spectrum as the window caught it (the default everywhere) and 0 where the
pipeline subtracted it, the emitter having had its own window (a deblended
transit).

# The model

The system must have an absolute frame (`plx`, `ra`, `dec`, `pmra`, `pmdec`,
`rv`, `ref_epoch`), which is Gaia DR4's. Per transit and emitter:

    velocity = frame rv + orbit + redshift + drv + κ (w − ref_pos_al) + κ_ac (z − ref_pos_ac)
    weight   = fraction of the light inside the window at z − ref_pos_ac, times membership

with (w, z) the emitter's local-plane position ([`transit_position`](@ref)).
So the systemic velocity is the frame's `rv`, the one that drives
perspective acceleration, and a companion displaced along scan appears
Doppler-shifted by `κ` per mas.

# Variables

Read from each emitter's body: `vsini` [km/s], the rotational broadening
(none if absent); `drv` [**m/s**, as every radial velocity in Octofitter],
the zero point of that body's template (zero if absent); and whatever its
emitter needs (`teff`, `logg`, `mh`, `radius` for a [`Star`](@ref)).

Of the observation itself, both optional: `kappa_ac` [km s⁻¹ mas⁻¹], the
across-scan term of the wavelength scale, and `deblend_residual`, the
fraction of a subtracted emitter's light that is left.
"""
struct GaiaDR4RawObs{Names,TE<:Tuple,TB<:Tuple,TTable<:Table,TR<:RVSEpochSpectra} <: AbstractObs
    table::TTable
    rvs::TR
    emitters::NamedTuple{Names,TE}
    bodies::TB
    tangent_point::NTuple{2,Float64}
    tangent_basis::NTuple{9,Float64}
    membership::Matrix{Float64}
    sinψ::Vector{Float64}
    cosψ::Vector{Float64}
    priors::Priors
    derived::Derived
    name::String
end

function GaiaDR4RawObs(; rvs, emitters::NamedTuple, name="GaiaDR4Raw", tangent_point=nothing,
    instrument::RVSInstrument=RVSInstrument(), wavelength=DR4_RVS_WAVELENGTH, p_bad=1e-3,
    mismatch=instrument.mismatch_amplitude, membership=nothing,
    variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    table = Table(rvs)
    issubset(RVS_COLUMNS, Octofitter.Tables.columnnames(table)) ||
        error("GaiaDR4RawObs: the `rvs` table needs the columns $RVS_COLUMNS")
    issorted(table.epoch) || error("GaiaDR4RawObs: the rows of `rvs` must be in order of epoch")
    isempty(emitters) && error("GaiaDR4RawObs: no emitters")
    grids = map(rvs_grid, values(emitters))
    all(==(first(grids)), grids) || error("GaiaDR4RawObs: the emitters' RVS spectra are on different wavelength grids")
    tp = Octofitter._dr4_tangent_point(table, tangent_point)
    tp === nothing && error("GaiaDR4RawObs \"$name\": no tangent point. Keep the table's `ra0`/`dec0` columns " *
                            "or pass `tangent_point = (ra0, dec0)` [deg].")
    nt, ne = length(table.epoch), length(emitters)
    membership = membership === nothing ? ones(ne, nt) : Matrix{Float64}(membership)
    size(membership) == (ne, nt) || throw(DimensionMismatch("membership must be [emitter, transit] = ($ne, $nt)"))
    all(m -> 0 <= m <= 1, membership) || throw(ArgumentError("membership must be in [0, 1]"))
    all(>(0), sum(membership; dims=1)) || throw(ArgumentError("a transit has no emitter in it"))
    flux = reduce(hcat, table.flux)
    flux_error = reduce(hcat, table.flux_error)
    spectra = RVSEpochSpectra(flux, flux_error, first(grids); instrument, wavelength, p_bad, mismatch)
    bodies = map(Octofitter.refspec, keys(emitters))
    priors, derived = variables
    return GaiaDR4RawObs(table, spectra, emitters, bodies, tp, tangent_basis(tp...), membership,
        collect(Float64, sind.(table.scan_pos_angle)), collect(Float64, cosd.(table.scan_pos_angle)),
        priors, derived, String(name))
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

function _rebuild(obs::GaiaDR4RawObs, table, membership)
    return GaiaDR4RawObs(; rvs=table, obs.emitters, obs.name, tangent_point=obs.tangent_point,
        instrument=obs.rvs.instrument, wavelength=obs.rvs.wavelength, membership, variables=(obs.priors, obs.derived),
        p_bad=obs.rvs.p_bad[_rows(obs, table)], mismatch=obs.rvs.mismatch[:, _rows(obs, table)])
end
# rows of `obs.table` that `table` holds, by epoch (a subset or the whole, in order)
_rows(obs::GaiaDR4RawObs, table) = [findfirst(==(e), obs.table.epoch) for e in table.epoch]

Octofitter.likeobj_from_epoch_subset(obs::GaiaDR4RawObs, inds) =
    _rebuild(obs, obs.table[inds], obs.membership[:, inds])

# ---------------------------------------------------
# The scene
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
    grid = obs.rvs.grid
    limb = obs.rvs.instrument.limb_darkening
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
    inst = obs.rvs.instrument
    tab = obs.table
    sol = solutionat(ctx, t)
    kappa_ac = hasproperty(θ_obs, :kappa_ac) ? θ_obs.kappa_ac : 0.0
    residual = hasproperty(θ_obs, :deblend_residual) ? θ_obs.deblend_residual : 0.0
    ks = ntuple(identity, length(Names))
    scene = map(ks, values(obs.emitters), θ_bodies, refs) do k, emitter, θ_body, body
        w, z = @inbounds transit_position(sol, body, barycentre, obs.tangent_basis, obs.sinψ[t], obs.cosψ[t],
            tab.parallax_factor_al[t], tab.parallax_factor_ac[t])
        Δη, Δζ = @inbounds (w - tab.ref_pos_al[t], z - tab.ref_pos_ac[t])
        v = transit_velocity(sol, body, barycentre, ctx.einstein_rv) + rvs_velocity_offset(emitter, θ_body) + _drv_kms(θ_body)
        m = @inbounds obs.membership[k, t]
        (rvs_velocity(inst, v, Δη; Δζ, kappa_ac), window_fraction(inst, Δζ) * (m + (1 - m) * residual))
    end
    return map(first, scene), map(last, scene)
end

function Octofitter.ln_like(obs::GaiaDR4RawObs{Names}, ctx::ObsContext) where {Names}
    T = Octofitter._system_number_type(ctx.θ_system)
    refs = Octofitter.resolverefs(ctx, obs.bodies)
    barycentre = Octofitter.ref(ctx, Octofitter.Barycentre)
    @no_escape ctx.buf begin
        work = @alloc(T, length(obs.rvs.grid), 4 * length(Names))
        components = _rvs_components(work, obs, ctx.θ_system)
        ll = zero(T)
        for t in 1:ntransits(obs.rvs)
            velocities, weights = _rvs_scene(obs, ctx, t, refs, barycentre)
            ll += first(rvs_transit_loglike(obs.rvs, t, components, velocities, weights; buf=ctx.buf))
        end
        ll
    end
end

"""
    Octofitter.simulate(obs::GaiaDR4RawObs, ctx) -> (; flux, velocity, weight, epochs)

The model's epoch spectra at the prior mean of the marginalised modes,
`flux[bin, transit]`, with each emitter's velocity [km/s] and weight in each
transit, `[emitter, transit]`.
"""
function Octofitter.simulate(obs::GaiaDR4RawObs{Names}, ctx::ObsContext) where {Names}
    return _simulate(obs, ctx, (t, c, v, w) -> rvs_transit_model(obs.rvs, t, c, v, w))
end

function _simulate(obs::GaiaDR4RawObs{Names}, ctx::ObsContext, spectrum::F) where {Names,F}
    T = Octofitter._system_number_type(ctx.θ_system)
    refs = Octofitter.resolverefs(ctx, obs.bodies)
    barycentre = Octofitter.ref(ctx, Octofitter.Barycentre)
    nt, ne = ntransits(obs.rvs), length(Names)
    flux = Matrix{Float64}(undef, length(obs.rvs.wavelength), nt)
    velocity, weight = Matrix{T}(undef, ne, nt), Matrix{T}(undef, ne, nt)
    work = Matrix{T}(undef, length(obs.rvs.grid), 4ne)
    components = _rvs_components(work, obs, ctx.θ_system)
    for t in 1:nt
        v, w = _rvs_scene(obs, ctx, t, refs, barycentre)
        velocity[:, t] .= v
        weight[:, t] .= w
        flux[:, t] .= spectrum(t, components, v, w)
    end
    return (; flux, velocity, weight, epochs=obs.table.epoch)
end

function Octofitter.generate_from_params(obs::GaiaDR4RawObs, ctx::ObsContext; add_noise)
    rng = Octofitter.Random.default_rng()
    sim = _simulate(obs, ctx, (t, c, v, w) -> rvs_transit_simulate(rng, obs.rvs, t, c, v, w; add_noise))
    # bins without data stay without data
    flux = [ifelse.(isfinite.(obs.rvs.flux[:, t]) .& isfinite.(obs.rvs.flux_error[:, t]), sim.flux[:, t], NaN) for t in axes(sim.flux, 2)]
    return _rebuild(obs, Table(obs.table; flux), obs.membership)
end
