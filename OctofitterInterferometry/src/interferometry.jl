# ---------------------------------------------------
# Interferometric observables
#
# One observation type, replacing v1's `InterferometryObs` +
# `GRAVITYWideKPObs`. The two shared the OI-FITS reader, the u/v tables and
# `cvis_bin!`, and — once the primary stops being privileged — the entire
# position/flux front-end as well. What actually differed was one stage: the
# observable projection and its covariance (plain closure phases with diagonal
# σ, versus the kernel-phase basis with the `CKP` correlation model and a
# block Cholesky). That difference now lives in the `phases` field and two
# methods of `_epoch_ln_like`; everything up to `cps_model` is shared.
#
# Three separate places in v1 baked in a single privileged companion, and all
# three are gone:
#
#  1. `cvis_model .+= 1.0` — the host's flux was literally hard-coded to 1,
#     and the sum normalized by `1/(1 + Σ contrast)`. The host is now an
#     ordinary entry in `targets` reading an ordinary `flux_<band>`.
#  2. Companion positions came from `raoff(sol)` minus a hand-rolled
#     photocentre perturbation summed over "inner" planets, selected at
#     runtime by `semimajoraxis(other) < semimajoraxis(this)`. That is
#     `raoff(sol, target, ref)` now — via the shared `sky_offset` front-end.
#  3. The fibre pointing was `f·ρ/(1+f)`, a closed-form *two-body*
#     photocentre with no meaning for three bodies. It is a reference now.
#
# Setting the host's `flux_<band> = 1.0` reproduces (1) and (2) bit-for-bit;
# `test/runtests.jl` asserts it. (3) genuinely changes numbers — see the
# docstring.
# ---------------------------------------------------

"""
    ClosurePhases()

Model closure phases directly, with the per-observable uncertainties from the
data and an optional `σ_cp_jitter` added in quadrature. Optionally also models
squared visibilities, for rows whose `use_vis2` is true.

Works for any array geometry.
"""
struct ClosurePhases end

"""
    KernelPhases(; correlated=true)

Project the closure-phase residuals onto the kernel-phase basis (the orthonormal
basis of the closure design's row space built by the internal
`OctofitterInterferometry.kernel_phase_basis`) and evaluate them under the block-structured
correlation model of `CKP`. `correlated=false` keeps the projection but
assumes independent spectral channels.

Requires GRAVITY's array geometry (6 baselines, 4 closure triangles).
"""
struct KernelPhases
    correlated::Bool
end
KernelPhases(; correlated::Bool=true) = KernelPhases(correlated)

const AbstractPhaseModel = Union{ClosurePhases,KernelPhases}

export ClosurePhases, KernelPhases

const interferometry_cols = (:epoch, :u, :v, :cps_data, :dcps,
                             :index_cps1, :index_cps2, :index_cps3)
const interferometry_cols_vis2 = (:vis2_data, :dvis2)

"""
    InterferometryObs(data...; targets, ref=Barycentre, band, name, …)

Interferometric observables — closure phases, optionally squared
visibilities, optionally projected onto a kernel-phase basis — of a
collection of point sources.

    vis = InterferometryObs(tab;
        targets = (A, b, c),
        ref     = A,
        band    = :K,
        name    = "GRAVITY",
        variables = @variables begin
            σ_cp_jitter ~ LogUniform(0.1, 100)
            platescale = 1.0
            northangle = 0.0
        end)

# The forward model

    V(u,v) = Σ_j f_j exp(−2πi (u Δα✱_j + v Δδ_j)) / Σ_j f_j

`targets` names the bodies in the sum, `ref` is the phase centre the offsets
are measured from, and `f_j` is each body's `flux_<band>` variable — the
**host included**, so a model reproducing v8's contrast-ratio convention gives
the host `flux_K = 1.0` and the companions their ratios.

Shifting the phase centre multiplies `V` by a global phase, which every
supported observable (|V|², closure phases, kernel phases) is invariant to, so
`ref` is a free choice; `Barycentre` is the default and `A` matches v8's
spelling exactly.

!!! warning "…invariant modulo 360°"
    `closurephase!` folds each *baseline* phase into (−180°, 180°] and
    then sums the triangle, without folding the sum. So a phase centre far from
    the sources — far enough that the individual baseline phases wrap — shifts a
    modelled closure phase by a multiple of 360° relative to a nearby one, and
    the residual against data that live in (−180°, 180°] changes with it. This
    is v8's behaviour, kept deliberately (v8 has the folding line commented out
    in the source), but it means `ref` should be somewhere near the flux
    centroid in practice. `Barycentre` is that for a faint companion, and
    `Photocentre(band)` is that in general.

`targets` is a **structural declaration** of which sources this observation's
forward model contains, and is deliberately not defaulted to "every body with
a flux in this band": "fit with `c` only" and "include `b`, but let its flux
go to zero" are different models, and you must be able to say the first.

# Data

Each row is one exposure. Give either a `filename` column pointing at an
OI-FITS file (read by `_prepare_input_row`, with optional
`wavelength_min_meters` / `wavelength_max_meters` cuts) or the prepared
columns directly: `epoch` [MJD], `u`, `v` [inverse wavelengths, baseline ×
channel], `cps_data`, `dcps` [degrees, triangle × channel], `index_cps1`,
`index_cps2`, `index_cps3`, and `eff_wave` [m]. Add `vis2_data`, `dvis2` and a
true `use_vis2` to include squared visibilities. Rows are sorted by epoch.

# Options

  - `kernel_phases=true` selects [`KernelPhases`](@ref) over the default
    [`ClosurePhases`](@ref); `kp_correlation=false` drops the spectral
    correlation within a kernel phase. In kernel-phase mode a `jitter` column
    and an optional `kp_Cy` column name, per exposure, the variables holding
    the kernel-phase jitter [degrees] and the spectral correlation
    coefficient; they may live in this observation's variables or the
    system's.
  - `fiber_coupling=true` models single-mode fibre injection losses.
    `fiber_pointing` is where the fibre is pointed, as a reference —
    `Photocentre(band)` by default, or a body, or `Photocentre(band, (A, b))`
    for an explicit subset — and each source's throughput is a function of its
    own offset from that point.

# Variables

  - `σ_cp_jitter` [degrees] — added in quadrature to the closure-phase
    uncertainties (closure-phase mode only; kernel-phase mode uses the
    per-exposure `jitter` variable named by the table).
  - `platescale`, `northangle` [rad] — the instrument calibration, applied as
    in [`Octofitter.sky_offset`](@ref).

!!! warning "Two v8 behaviours this does not reproduce"
    **`platescale` is now a divisor.** v8's interferometry likelihood
    *multiplied* the modelled offsets by `platescale`, the reciprocal of the
    convention used by relative astrometry and by images. The shared
    `sky_offset` front-end uses the majority convention, so a non-unity
    `platescale` now means the opposite of what it did. `platescale = 1`
    (the default) is unaffected.

    **Fibre coupling is computed differently.** v8 evaluated the throughput of
    *every* companion at the host-to-photocentre distance `f·ρ/(1+f)` and left
    the host at 1.0. Each source's throughput is now evaluated at its own
    offset from `fiber_pointing`, which is the quantity the injection
    efficiency actually depends on, and which is defined for any number of
    bodies. For a faint companion the two differ by roughly the full coupling
    loss at the companion's separation.
"""
struct InterferometryObs{TTable<:Table,TT<:Tuple,TR,TP<:AbstractPhaseModel,TF,TI} <: AbstractInterferometryObs
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    targets::TT                 # tuple of BodyRefSpec — the sources in the sum
    ref::TR                     # phase centre
    band::Union{Symbol,Nothing}
    name::String
    phases::TP
    fiber_pointing::TF          # `nothing`, or a reference spec
    fiber_coupling::TI          # `nothing`, or (sep_mas, λ_m) -> throughput
    # Buffer sizes, so the hot loop does not recompute them per evaluation.
    n_baselines_max::Int
    n_cps_max::Int
    n_resid_max::Int
    n_kp_max::Int
end

export InterferometryObs

function InterferometryObs(
    observations...;
    targets,
    ref=Barycentre,
    band::Union{Symbol,Nothing}=nothing,
    name,
    kernel_phases::Bool=false,
    kp_correlation::Bool=true,
    fiber_coupling::Bool=false,
    fiber_pointing=nothing,
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.Priors(), Octofitter.Derived()),
)
    (priors, derived) = variables

    input_table = Table(observations...)
    table = if :filename ∈ Tables.columnnames(input_table)
        Table(map(_prepare_input_row, eachrow(input_table)))
    else
        input_table
    end

    issubset(interferometry_cols, Tables.columnnames(table)) || error(
        "Expected columns $interferometry_cols, got $(Tables.columnnames(table))")
    if !hasproperty(table, :use_vis2)
        table = Table(table; use_vis2=fill(false, length(table.epoch)))
    end
    if any(r -> r === true, table.use_vis2)
        issubset(interferometry_cols_vis2, Tables.columnnames(table)) || error(
            "`use_vis2` is set on at least one row, so columns " *
            "$interferometry_cols_vis2 are required; got $(Tables.columnnames(table))")
    end
    if !hasproperty(table, :eff_wave)
        # Only needed for the wavelength axis of the fibre-coupling grid and
        # for the kernel-phase block sizes; the visibility model itself works
        # off `u`/`v`, which are already in inverse wavelengths.
        (kernel_phases || fiber_coupling) && error(
            "`kernel_phases` and `fiber_coupling` need an `eff_wave` column [m] " *
            "giving each exposure's spectral channels.")
    end
    table = table[sortperm(vec(table.epoch))]

    if kernel_phases
        table = Table(map(r -> _prepare_kernel_phases(only(r)), eachrow(table)))
    end

    targetspecs = _target_specs(targets)
    phasespec = Octofitter.refspec(ref)
    phases = kernel_phases ? KernelPhases(kp_correlation) : ClosurePhases()

    coupling, pointing = if fiber_coupling
        @info "Pre-calculating fiber coupling efficiency over grid"
        interp = fiber_coupling_interpolator(reduce(vcat, map(vec, table.eff_wave)))
        pt = isnothing(fiber_pointing) ?
             (isnothing(band) ? Photocentre : Photocentre(band)) :
             Octofitter.refspec(fiber_pointing)
        (interp, pt)
    else
        isnothing(fiber_pointing) || error(
            "`fiber_pointing` was given but `fiber_coupling=false`, so it would " *
            "have no effect. Pass `fiber_coupling=true` as well.")
        (nothing, nothing)
    end

    return _build(table, priors, derived, targetspecs, phasespec, band, String(name),
                  phases, pointing, coupling)
end
InterferometryObs(observations::NamedTuple...; kwargs...) =
    InterferometryObs(observations; kwargs...)

# Assemble the struct from an already-prepared table, recomputing the cached
# buffer sizes. `likeobj_from_epoch_subset` and `generate_from_params` go
# through here rather than through the constructor: the prepared table still
# carries its `filename` column, so re-running the constructor would re-read
# every OI-FITS file (and, in kernel-phase mode, refactorize every design
# matrix) to arrive back where it started.
function _build(table, priors, derived, targets, ref, band, name, phases, pointing, coupling)
    n_bl = maximum(m -> size(m, 1), table.u; init=0)
    n_cp = maximum(m -> size(m, 1), table.cps_data; init=0)
    n_resid, n_kp = if phases isa KernelPhases
        (maximum(i -> size(table.cps_data[i], 1) * length(table.eff_wave[i]),
                 eachindex(table.epoch); init=0),
         maximum(m -> size(m, 1), table.P₁; init=0))
    else
        (0, 0)
    end
    return InterferometryObs{typeof(table),typeof(targets),typeof(ref),typeof(phases),
                             typeof(pointing),typeof(coupling)}(
        table, priors, derived, targets, ref, band, name, phases, pointing, coupling,
        n_bl, n_cp, n_resid, n_kp)
end

"""
    GRAVITYWideKPObs(data...; targets, ref, band, name="GRAVITY-WIDE", …)

GRAVITY-WIDE preset: [`InterferometryObs`](@ref) with `kernel_phases=true` and
`fiber_coupling=true`.

This is a *preset*, not the v8 type — v8's `GRAVITYWideKPObs` took its
companion fluxes from a positionally-indexed `flux` vector on the observation
and had no notion of which bodies it was modelling. `targets` and `band` are
required, and the fluxes come from the bodies' `flux_<band>` variables.
"""
GRAVITYWideKPObs(observations...; name="GRAVITY-WIDE", kernel_phases=true,
                 fiber_coupling=true, kwargs...) =
    InterferometryObs(observations...; name, kernel_phases, fiber_coupling, kwargs...)

export GRAVITYWideKPObs

# --- the reference declarations -----------------------------------------------

# Targets must be bodies. A `Barycentre`/`Photocentre` spec resolves to a
# `WeightedPoint`, whose weights are normalized — the *total* band flux of the
# set it stands for is not recoverable from it, and this model needs exactly
# that. An unresolved pair standing in as one source is spelled as its two
# bodies, which is also the more accurate model.
function _target_specs(targets)
    ts = targets isa Union{Tuple,AbstractVector} ? Tuple(targets) : (targets,)
    isempty(ts) && error(
        "`targets` is empty: an interferometric observation needs at least one source.")
    specs = map(Octofitter.refspec, ts)
    all(s -> s isa Octofitter.BodyRefSpec, specs) || error("""
        `targets` must name bodies — each source in the visibility sum carries its
        own `flux_<band>` and its own position. `Barycentre`/`Photocentre` specs
        resolve to normalized weighted points, which do not carry a total flux, so
        they cannot stand in for a source here. (They are fine as `ref` or
        `fiber_pointing`, which are positions only.)
        """)
    return specs
end

Octofitter.refspecs(obs::InterferometryObs) = isnothing(obs.fiber_pointing) ?
    (obs.targets..., obs.ref) : (obs.targets..., obs.ref, obs.fiber_pointing)

# The generic reading of `refspecs` would take the fibre pointing for the
# reference. Say what each entry is instead.
function Octofitter._refdesc(obs::InterferometryObs)
    s = Octofitter._refdesc(obs.targets, obs.ref)
    isnothing(obs.fiber_pointing) && return s
    return s * " [fiber at " * Octofitter._refstr(obs.fiber_pointing) * "]"
end

function Octofitter.likeobj_from_epoch_subset(obs::InterferometryObs, obs_inds)
    return _build(obs.table[obs_inds], obs.priors, obs.derived, obs.targets, obs.ref,
                  obs.band, obs.name, obs.phases, obs.fiber_pointing, obs.fiber_coupling)
end

# --- the front end ------------------------------------------------------------

# Per-target band fluxes. v1 read `θ_obs.flux[i_planet]`: a vector on the
# *observation*, indexed by the order the planets happened to be declared in,
# with the host's value implied to be 1. Fluxes are body variables now, so
# there is no index to get wrong and the host is not special.
@inline function _target_fluxes(ctx, obs::InterferometryObs, targets)
    fl = _bandfluxes(ctx.system, obs.band)
    return map(t -> (@inbounds fl[t.idx]), targets)
end

@inline _bandfluxes(sys, band::Symbol) = PlanetOrbits.fluxes(sys, band)
@inline function _bandfluxes(sys, ::Nothing)
    allbands = PlanetOrbits.fluxes(sys)
    length(allbands) == 1 || _err_band(keys(allbands))
    return values(allbands)[1]
end
@noinline function _err_band(bands)
    isempty(bands) && error(
        "InterferometryObs: no body in this system declares a flux, so there is " *
        "nothing to build a visibility out of. Give each source a `flux_<band>` " *
        "variable in its own block (the host included).")
    error("InterferometryObs: this system declares the flux bands " *
          join(bands, ", ") * ", so the observation must say which one it " *
          "observes: `band=:" * string(first(bands)) * "`.")
end

@noinline _err_legacy_flux() = error("""
    This model gives the interferometry observation a `flux` variable. That was
    v8's spelling: a vector of companion contrast ratios, indexed by planet
    order, against a host whose flux was hard-coded to 1.

    In v9 each body carries its own `flux_<band>` variable and the observation
    names the bodies it models:

        A = Body(name="A", variables=@variables begin
            mass = 1.2
            flux_K = 1.0            # the host is an ordinary source now
        end)
        b = Body(name="b", about=A, variables=@variables begin
            mass = 10mjup
            flux_K ~ Uniform(0, 1)  # what `flux[1]` used to be
            …
        end)
        InterferometryObs(dat; targets=(A, b), ref=A, band=:K, name="…")
    """)

@noinline _err_dark() = error(
    "every body in `targets` has zero flux in this band, so the modelled " *
    "visibility is 0/0. Give at least one source a non-zero `flux_<band>`.")

@inline function _fiber_ref(ctx, obs::InterferometryObs)
    isnothing(obs.fiber_pointing) && return nothing
    return Octofitter.ref(ctx, obs.fiber_pointing)
end

# Sky offsets of every source (and of the fibre) from the phase centre, at one
# epoch. `map` over the target tuple, so the per-source types stay concrete and
# the loop unrolls.
@inline function _epoch_offsets(ctx, i_epoch, front)
    sol = Octofitter.solutionat(ctx, i_epoch)
    ps, na = front.platescale, front.northangle
    offsets = map(front.targets) do t
        Octofitter.sky_offset(sol, t, front.phasecentre; platescale=ps, northangle=na)
    end
    fiberoff = isnothing(front.fiber) ? nothing :
               Octofitter.sky_offset(sol, front.fiber, front.phasecentre;
                                     platescale=ps, northangle=na)
    return offsets, fiberoff
end

# Effective brightness of each source at one wavelength: its band flux, times
# its own single-mode-fibre throughput when one is being modelled.
@inline _effective_weights(::InterferometryObs, fluxes, offsets, ::Nothing, λ) = fluxes
@inline function _effective_weights(obs::InterferometryObs, fluxes, offsets, fiberoff, λ)
    return map(fluxes, offsets) do f, o
        sep = hypot(o[1] - fiberoff[1], o[2] - fiberoff[2])
        f * obs.fiber_coupling(sep, λ)
    end
end

# A per-exposure variable named by a table column (kernel-phase jitter, and the
# spectral correlation coefficient). It may live on this observation or on the
# system; derived values are already materialized into both by the time a
# likelihood runs, so v1's third `obs.derived[name](θ_system)` branch is gone.
@inline function _named_var(ctx, name::Symbol, ::Type{T}) where {T}
    if hasproperty(ctx.θ_obs, name)
        return convert(T, getproperty(ctx.θ_obs, name))
    elseif hasproperty(ctx.θ_system, name)
        return convert(T, getproperty(ctx.θ_system, name))
    end
    _err_named_var(name)
end
@noinline _err_named_var(name) = error(
    "InterferometryObs: the data table names `:$name` as a variable, but neither " *
    "the observation's nor the system's variables block defines it.")

# --- the likelihood -----------------------------------------------------------

function Octofitter.ln_like(obs::InterferometryObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    hasproperty(ctx.θ_obs, :flux) && _err_legacy_flux()

    front = _front_end(obs, ctx, T)
    iszero(sum(front.fluxes)) && _err_dark()

    ll = zero(T)
    @no_escape ctx.buf begin
        cvis = @alloc(Complex{T}, obs.n_baselines_max)
        cps = @alloc(T, obs.n_cps_max)
        resid = @alloc(T, obs.n_resid_max)
        kp_resid = @alloc(T, obs.n_kp_max)
        bufs = (; cvis, cps, resid, kp_resid)
        # The kernel-phase covariance is factorized in place and handed to
        # `PDMat`, so it has to be an ordinary owned matrix; it is allocated
        # once here at the largest size any exposure needs and viewed per
        # exposure rather than reallocated per exposure as v1 did.
        Σ, C = _kp_workspace(obs.phases, obs.n_kp_max, T)
        # `ll` is threaded through rather than summed per epoch: v1 accumulated
        # one running total over every (epoch, wavelength) term, and
        # floating-point addition is not associative, so summing per-epoch
        # subtotals would break the bit-for-bit agreement the tests pin.
        for i_epoch in eachindex(obs.table.epoch)
            ll = _epoch_ln_like(obs.phases, obs, ctx, i_epoch, front, bufs, Σ, C, ll)
        end
    end
    return ll
end

@inline function _front_end(obs::InterferometryObs, ctx, ::Type{T}) where {T}
    targets = Octofitter.resolverefs(ctx, obs.targets)
    platescale, northangle = Octofitter.sky_calibration(ctx)
    return (;
        targets,
        phasecentre=Octofitter.ref(ctx, obs.ref),
        fiber=_fiber_ref(ctx, obs),
        fluxes=_target_fluxes(ctx, obs, targets),
        platescale, northangle,
        σ_cp_jitter=hasproperty(ctx.θ_obs, :σ_cp_jitter) ? ctx.θ_obs.σ_cp_jitter : zero(T),
    )
end

_kp_workspace(::ClosurePhases, _, ::Type{T}) where {T} = (nothing, nothing)
_kp_workspace(::KernelPhases, n, ::Type{T}) where {T} =
    (Matrix{T}(undef, n, n), Matrix{T}(undef, n, n))

# --- closure phases -----------------------------------------------------------

function _epoch_ln_like(::ClosurePhases, obs::InterferometryObs, ctx, i, front, bufs,
                        _Σ, _C, ll::T) where {T}
    tab = obs.table
    u_all, v_all = tab.u[i], tab.v[i]
    cps_data, dcps = tab.cps_data[i], tab.dcps[i]
    index_cps1, index_cps2, index_cps3 = tab.index_cps1[i], tab.index_cps2[i], tab.index_cps3[i]
    cvis = @view bufs.cvis[1:size(u_all, 1)]
    cps = @view bufs.cps[1:size(cps_data, 1)]
    jitter = front.σ_cp_jitter
    use_vis2 = tab.use_vis2[i] === true

    offsets, fiberoff = _epoch_offsets(ctx, i, front)

    for i_wave in axes(u_all, 2)
        u = @view u_all[:, i_wave]
        v = @view v_all[:, i_wave]
        w = _effective_weights(obs, front.fluxes, offsets, fiberoff, _λ(tab, i, i_wave))
        _cvis_epoch!(cvis, u, v, offsets, w)
        closurephase!(cps; vis=cvis, index_cps1, index_cps2, index_cps3)

        # Written in v1's two-pass order — all the normalizations, then all the
        # residuals — because the bit-for-bit check runs through here and
        # floating-point addition is not associative.
        σ_cp = @view dcps[:, i_wave]
        const_cp = zero(T)
        for I in eachindex(σ_cp)
            const_cp -= log(2π * (σ_cp[I]^2 + jitter^2))
        end
        const_cp /= 2
        lnlike_cp = zero(T)
        data = @view cps_data[:, i_wave]
        for I in eachindex(σ_cp, data, cps)
            σ² = jitter^2 + σ_cp[I]^2
            lnlike_cp -= 0.5 * (data[I] - cps[I])^2 / σ²
        end
        ll += lnlike_cp + const_cp

        if use_vis2
            # v1's squared-visibility branch read an uninitialized
            # `lnlike_v2` and so could only ever throw; this is the term it
            # was reaching for.
            vis2_data = @view tab.vis2_data[i][:, i_wave]
            dvis2 = @view tab.dvis2[i][:, i_wave]
            lnlike_v2 = zero(T)
            for I in eachindex(vis2_data, dvis2, cvis)
                σ² = dvis2[I]^2
                lnlike_v2 -= 0.5 * ((vis2_data[I] - abs2(cvis[I]))^2 / σ² + log(2π * σ²))
            end
            ll += lnlike_v2
        end
    end
    return ll
end

# `eff_wave` is only needed when something is wavelength-dependent beyond
# u/v — currently just the fibre throughput — so a plain closure-phase table
# is not required to carry it.
@inline _λ(tab, i, i_wave) = hasproperty(tab, :eff_wave) ? tab.eff_wave[i][i_wave] : NaN

# --- kernel phases ------------------------------------------------------------

function _epoch_ln_like(kp::KernelPhases, obs::InterferometryObs, ctx, i, front, bufs,
                        Σ_buf, C_buf, ll::T) where {T}
    tab = obs.table
    u_all, v_all = tab.u[i], tab.v[i]
    cps_data = tab.cps_data[i]
    index_cps1, index_cps2, index_cps3 = tab.index_cps1[i], tab.index_cps2[i], tab.index_cps3[i]
    Λ = length(tab.eff_wave[i])
    n_T3 = size(cps_data, 1)
    cvis = @view bufs.cvis[1:size(u_all, 1)]
    cps = @view bufs.cps[1:n_T3]
    resid = @view bufs.resid[1:(n_T3*Λ)]

    offsets, fiberoff = _epoch_offsets(ctx, i, front)

    for i_wave in axes(u_all, 2)
        u = @view u_all[:, i_wave]
        v = @view v_all[:, i_wave]
        w = _effective_weights(obs, front.fluxes, offsets, fiberoff, tab.eff_wave[i][i_wave])
        _cvis_epoch!(cvis, u, v, offsets, w)
        closurephase!(cps; vis=cvis, index_cps1, index_cps2, index_cps3)
        # Grouped by triangle, channels within — the layout `P₁` was built for.
        for i_T3 in 1:n_T3
            resid[(i_T3-1)*Λ+i_wave] = cps_data[i_T3, i_wave] - cps[i_T3]
        end
    end

    P₁ = tab.P₁[i]
    σ_kp = tab.σ_kp[i]
    L = size(P₁, 1)
    kp_resid = @view bufs.kp_resid[1:L]
    mul!(kp_resid, P₁, resid)

    jitter = hasproperty(tab, :jitter) ? max(eps(), _named_var(ctx, tab.jitter[i], T)) : zero(T)
    Cy = kp.correlated && hasproperty(tab, :kp_Cy) ?
         max(eps(), _named_var(ctx, tab.kp_Cy[i], T)) : zero(T)

    C = @view C_buf[1:L, 1:L]
    Σ = @view Σ_buf[1:L, 1:L]
    CKP!(C, tab[i], Cy)
    _kp_covariance!(Σ, C, σ_kp, jitter)
    # Three kernel phases per spectral channel with GRAVITY, each block
    # spanning all Λ channels — the partition `CKP!` just populated.
    return ll + _kp_logpdf(Σ, kp_resid, 3)
end

# --- simulation ---------------------------------------------------------------

"""
    simulate(obs::InterferometryObs, ctx) -> Vector of per-exposure NamedTuples

Model observables at this observation's epochs: `cps_model` [degrees]
(triangle × channel), `vis2_model` (baseline × channel) and the complex
`cvis_model` they were formed from.

The values are what the instrument would have *reported* — `platescale` and
`northangle` are already applied, in the direction described by
[`Octofitter.sky_offset`](@ref) — so they are directly comparable with
`obs.table.cps_data`.
"""
function Octofitter.simulate(obs::InterferometryObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    front = _front_end(obs, ctx, T)
    tab = obs.table
    return map(eachindex(tab.epoch)) do i
        u_all, v_all = tab.u[i], tab.v[i]
        n_bl, n_wave = size(u_all)
        n_T3 = size(tab.cps_data[i], 1)
        cvis_model = zeros(Complex{T}, n_bl, n_wave)
        cps_model = zeros(T, n_T3, n_wave)
        offsets, fiberoff = _epoch_offsets(ctx, i, front)
        for i_wave in 1:n_wave
            cvis = @view cvis_model[:, i_wave]
            w = _effective_weights(obs, front.fluxes, offsets, fiberoff, _λ(tab, i, i_wave))
            _cvis_epoch!(cvis, (@view u_all[:, i_wave]), (@view v_all[:, i_wave]), offsets, w)
            closurephase!((@view cps_model[:, i_wave]); vis=cvis,
                index_cps1=tab.index_cps1[i], index_cps2=tab.index_cps2[i],
                index_cps3=tab.index_cps3[i])
        end
        return (; epoch=tab.epoch[i], cps_model, vis2_model=abs2.(cvis_model), cvis_model)
    end
end

"""
    generate_from_params(obs::InterferometryObs, ctx; add_noise)

Replace this observation's closure phases (and squared visibilities, where
present) with the model's, optionally perturbed by their own uncertainties.

v8 accepted `add_noise` and ignored it, which made simulation-based
calibration on this type meaningless; the noise is drawn here.
"""
function Octofitter.generate_from_params(obs::InterferometryObs, ctx::Octofitter.ObsContext;
                                         add_noise)
    sims = Octofitter.simulate(obs, ctx)
    tab = obs.table
    cps_data = map(eachindex(sims)) do i
        cp = collect(sims[i].cps_model)
        add_noise ? cp .+ randn.() .* tab.dcps[i] : cp
    end
    table = Table(tab; cps_data)
    if hasproperty(tab, :vis2_data)
        vis2_data = map(eachindex(sims)) do i
            v2 = collect(sims[i].vis2_model)
            add_noise ? v2 .+ randn.() .* tab.dvis2[i] : v2
        end
        table = Table(table; vis2_data)
    end
    return _build(table, obs.priors, obs.derived, obs.targets, obs.ref, obs.band,
                  obs.name, obs.phases, obs.fiber_pointing, obs.fiber_coupling)
end

# ---------------------------------------------------
# Plotting protocol
#
# One exposure is not one number: it is an `(n_triangles × n_channels)` block
# of closure phases (and, optionally, an `(n_baselines × n_channels)` block of
# squared visibilities). The channels below flatten those blocks and repeat the
# epoch, so every measured closure phase is one point in the generic panel —
# which is the honest data-space view, and the one that makes a systematic
# offset in a single triangle visible.
#
# There is no `PlanetOrbits` observable behind a closure phase, so
# `query=nothing`: the panel draws model points, residuals and the marginal
# histogram, not a smooth curve.
#
# Under `KernelPhases`, the quantity `ln_like` actually scores is the *kernel*
# projection `P₁·resid` against a dense within-block covariance, which has no
# per-exposure data counterpart to plot against. The closure phases are still
# the measurement, so they are what is shown; the residual strip is then a
# data-space diagnostic rather than a whitened likelihood residual, and the
# reported σ is the per-point `dcps` rather than the projected covariance.
# ---------------------------------------------------

function Octofitter.plotchannels(obs::InterferometryObs)
    chs = Any[Octofitter.PlotChannel(:closure_phase, "closure phase", "°")]
    any(obs.table.use_vis2 .=== true) &&
        push!(chs, Octofitter.PlotChannel(:vis2, "squared visibility", ""))
    return Tuple(chs)
end

function Octofitter.residuals(obs::InterferometryObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :σ_cp_jitter) ? Float64(ctx.θ_obs.σ_cp_jitter) : 0.0
    sims = Octofitter.simulate(obs, ctx)
    tab = obs.table

    ep = Float64[]; data = Float64[]; model = Float64[]; σ = Float64[]
    ep2 = Float64[]; data2 = Float64[]; model2 = Float64[]; σ2 = Float64[]
    for i in eachindex(tab.epoch)
        cp_d = vec(collect(Float64, tab.cps_data[i]))
        cp_m = vec(collect(Float64, sims[i].cps_model))
        append!(ep, fill(Float64(tab.epoch[i]), length(cp_d)))
        append!(data, cp_d); append!(model, cp_m)
        append!(σ, vec(collect(Float64, tab.dcps[i])))
        if tab.use_vis2[i] === true
            v_d = vec(collect(Float64, tab.vis2_data[i]))
            v_m = vec(collect(Float64, sims[i].vis2_model))
            append!(ep2, fill(Float64(tab.epoch[i]), length(v_d)))
            append!(data2, v_d); append!(model2, v_m)
            append!(σ2, vec(collect(Float64, tab.dvis2[i])))
        end
    end

    out = (; closure_phase=(; epoch=ep, data, model, resid=data .- model,
        σ, σ_eff=hypot.(σ, jitter), use=trues(length(ep))))
    isempty(ep2) && return out
    # Squared visibilities carry no jitter term in `ln_like`.
    return merge(out, (; vis2=(; epoch=ep2, data=data2, model=model2,
        resid=data2 .- model2, σ=σ2, σ_eff=copy(σ2), use=trues(length(ep2)))))
end
