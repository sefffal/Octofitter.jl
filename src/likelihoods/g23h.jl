# ---------------------------------------------------
# G23H — joint Gaia (DR2/DR3) + Hipparcos catalog astrometry
#
# One observation, several channels: the Hipparcos and Gaia catalog proper
# motions, the Hipparcos–Gaia and DR3–DR2 scaled position differences, the
# Gaia DR3 unit-weight-error variance (UEVA/RUWE), the Gaia RV variability,
# and optionally the Hipparcos per-transit abscissa residuals. They stay in
# one observation because the published catalog covariance couples them: the
# DR2/DR3 proper motions are correlated through `rho_dr2_dr3`, and the UEVA
# channel *deflates* the DR3 block.
#
# What the v2 port changed, and why
# ---------------------------------
#
# 1. `_simulate_skypath_perturbations!`'s per-companion linear photocentre
#    coefficient `(−m_k + f_k·m_host)/(M_k(1 + f_k))` is gone. It is
#    algebraically exact for **one** luminous companion, but for several it
#    normalizes each companion's term by its own `(1 + f_k)` and superposes —
#    and a photocentre is not a superposition (the point BINARYS makes).
#    v2 builds one `WeightedPoint` per draw from the effective flux ratios
#    and asks for `raoff(sol, photocentre, barycentre)`, which is the exact
#    flux-weighted mean of apparent positions for any number of companions
#    and any hierarchy. Single-companion configurations are unchanged to
#    roundoff; multi-companion ones differ, and the v2 answer is the correct
#    one. There is deliberately no legacy-compatibility mode.
#
# 2. `propagate_astrom` and its differential-light-travel-time term are gone.
#    That term was flagged `# This isn't right!` in the source (it double
#    counted the proper motion) and the absolute-frame branch now reads the
#    frame observables `frame_ra`/`frame_dec`/`frame_pmra`/`frame_pmdec` at
#    the catalog reference epochs, which PlanetOrbits propagates rigorously —
#    including the barycentric light-travel sign that v1 had inverted.
#
# 3. The ~250-line per-planet `planet_sols_cache` is gone. `epochs(obs)`
#    declares the Hipparcos transits, the Gaia forecast pool and the six
#    catalog reference epochs, and the model solves the union once.
#
# 4. `fluxratio`/`fluxratio_hip` remain *observation-consumed values* (they
#    are instrument-band contrast ratios against this source's host, not
#    intrinsic body fluxes), but the constructor now declares `host=` and
#    `companions=(…)` explicitly, so tuple order ↔ body mapping is written
#    down rather than implied by companion declaration order.
# ---------------------------------------------------

using Arrow

# Diagnostic hook (off by default; no hot-path cost when `nothing`). Set to a
# `Vector` and `ln_like` pushes a NamedTuple of the per-channel catalog-vs-model
# residuals and their marginal σ for the coupled PM/UEVA block — used to check
# that the simulator and the likelihood agree (pulls ~ N(0,1) at truth).
const _G23H_DEBUG_PULLS = Ref{Any}(nothing)
const _G23H_DEBUG_PULLS_LOCK = ReentrantLock()

"""
    G23HObs(; gaia_id, host, companions, …)

Joint Gaia DR2/DR3 + Hipparcos catalog astrometry from the G23H catalog:
calibrated proper motions at five epochs, the Gaia DR3 astrometric
excess-noise (UEVA/RUWE) channel, the Gaia RV variability channel, and
optionally the Hipparcos per-transit abscissae.

# Source membership
    host=A, companions=(b, c, d)

`host` is the body the catalog source is centred on; `companions` are the
bodies that may blend into it, **in the order the flux-ratio variables are
indexed**. Both take `Body` model nodes or `Symbol`s.

`host`/`companions` is not a single `Photocentre` spec, and cannot be: the
Hipparcos branch is a grating response rather than a linear reduction (see
`_hippacentre!`), and *which* bodies could ever blend is a structural,
build-time statement while *how much* each one blends is a per-draw one.

# Flux ratios

Each companion contributes light in proportion to its flux ratio against the
host, in two bands: `fluxratio` (G, for the Gaia DR2/DR3 photocentre) and
`fluxratio_hip` (Hp, for the Hipparcos abscissa). They are looked up in this
order, and the first hit wins:

 1. **this observation's own variables** — `fluxratio = (f_b, f_c)`. Use this
    when the ratios are *dynamic*: derived from system variables, including
    deferred ones, so a resolved-flag latent at system level can gate a
    companion out of the photocentre for that draw. That is the escape hatch
    the three-tier design needs, and it is why blending state never has to
    round-trip through a body's own flux variable (which would be the
    body→deferred-system cycle codegen rejects).
 2. **a system-level vector of the same name** — `system.fluxratio`.
 3. **the bodies' own fluxes** — `flux_G` / `flux_Hp` declared in each body's
    variables block, as `f_k = flux_<band>(companion_k) / flux_<band>(host)`.
    This is the common, *static* case, and it needs no vector at all: the
    ratios stop being positional and the same flux variable feeds every
    observation in that band.

A vector given under (1) or (2) must be a length-`length(companions)`
container (a tuple, `SVector`, or vector) matching `companions` element for
element; a bare scalar is accepted only when there is exactly one companion.
A model that declares no fluxes at all leaves every companion dark in every
band, which is the right answer for a fit with no photometry; a model that
declares *some* band but not the one asked for is a name mismatch and errors.

# Data
`catalog` is the G23H catalog: a path to the Arrow file, or an already-loaded
table (row selected by `gaia_id`), or a single row as a NamedTuple. Give
either `gaia_id` or `hip_id`; `hip_id` is resolved to a Gaia source id
through the catalog.

`dr2_transits_catalog` supplies the Gaia DR2 matched-transit count
(`astrometric_matched_observations_dr2`), which sizes the DR2 epoch
selection and is not carried by the published catalog. There is no fallback.

The Gaia scan geometry comes from one of:

  - `forecast_table` — a prepared table with `epoch` [MJD],
    `scanAngle_rad` and `parallaxFactorAlongScan`. Nothing is fetched and no
    ephemeris is needed; this is also how the test suite runs offline.
  - `scanlaw_table` — the `scanninglaw` Python package's output (`times`,
    `angles`); parallax factors are computed from the Earth ephemeris.
  - neither — the GOST web service is queried for `ra`/`dec`.

`hipparcos` supplies the Hipparcos IAD as `(; table, hip_sol)` (see
`hipparcos_iad`); by default it is loaded from the `Hipparcos_IAD`
data dependency using the catalog's `hip_id`. A catalog row with
`hip_id = NaN` drops every Hipparcos channel.

# Channels
`table.kind` names the channels, one row each: `:iad_hip`, `:ra_hip`,
`:dec_hip`, `:ra_hg`, `:dec_hg`, `:ra_dr2`, `:dec_dr2`, `:ra_dr32`,
`:dec_dr32`, `:ra_dr3`, `:dec_dr3`, `:ueva_dr3`, `:rv_dr3`.
`likeobj_from_epoch_subset` selects on those rows, so a caller can exclude a
channel by dropping its row.

`channels=` restricts the set at construction, and does it by filtering the
same `table.kind` rows — one code path, so the two spellings cannot diverge:

    channels = (:ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr3, :dec_dr3)

is what [`HGCAObs`](@ref) is. Dropping the last Hipparcos channel drops the
Hipparcos catalog distributions with it, and dropping `:iad_hip` drops the
six abscissa nuisance parameters, exactly as a post-hoc subset would.

`ueva_mode` is `:RUWE` (default), `:EAN`, or `:none`. `:none` drops the
`:ueva_dr3` datum and the UEVA-driven DR3 covariance deflation — for stars
whose `sig_AL`/`sig_att_radec`/`sig_cal` calibration is absent and which
therefore cannot form the priors that channel needs — and leaves every other
channel untouched.

# `frame_shift`
Default `true`, and single-source behaviour is unchanged. Every channel —
DR3's own included — is expressed relative to the model's DR3-epoch proper
motion of the *host body* rather than of the barycentre, so the frame `pmra`
stops meaning "the barycentre's" and starts meaning "the host's, as DR3
measured it". The Hipparcos–Gaia channel keeps its acceleration content
because the shift is common-mode across epochs, and the reparameterization
pays for itself: the frame proper motion is otherwise degenerate with the
reflex, and this points the sampled direction at what the data pin hardest.

It is a *reparameterization* only under a wide prior on `pmra`/`pmdec`. Put
an informative prior on either and `frame_shift=true` silently changes the
posterior, because the prior is then a prior on a different quantity than it
reads as. Nothing checks that for you.

**With more than one `G23HObs` on the same system frame it is an error**, not
a warning. Each observation would subtract its own `Δpm` while redefining the
one shared `pmra`, so every source's DR3 channel would predict the same
number — for a wide pair like ups And, tens of σ from two catalogue values
that genuinely differ. The relative proper motion is not merely uninformative
under that shift; it is removed from the likelihood, and it is the entire
wide-orbit signal. Use [`AnchoredFrame`](@ref), which reconditions the same
degeneracy in the model rather than inside one observation, and composes.

# Variables
Defaults are derived from the catalog row unless `variables` is given:
`σ_AL`, `σ_att`, `σ_calib` (mas, from the catalog's calibration), the epoch
selections `transit_priorities`/`transits`/`transits_dr2`/`transits_rv`, the
Hipparcos abscissa nuisance block (`hip_iad_jitter`, `iad_Δra`, `iad_Δdec`,
`iad_Δplx`, `iad_Δpmra`, `iad_Δpmdec`), and `σ_rv_per_transit`. Optional
`σ_hip_pmra`, `σ_hg_pmra`, `σ_dr2_pmra`, `σ_dr32_pmra`, `σ_dr3_pmra` (and
the `pmdec` counterparts) add per-channel proper-motion jitter in quadrature.

The system's `pmra`/`pmdec` (or this observation's, if the system has no
absolute frame) define the reference-point proper motion every modelled
channel is expressed against.

# Several sources in one system
Two Gaia sources belonging to one physical system are two `G23HObs` on one
`System`, sharing one absolute frame. That shared frame is the whole point —
it is what makes the wide pair's relative astrometry constrain the wide orbit
— so give each observation its own `host=`/`companions=` and let the frame do
the binding. Set `frame_shift=false` on all of them and parameterize the frame
with [`AnchoredFrame`](@ref).

**Every nuisance parameter here stays per source, deliberately.**
`transit_priorities`, `σ_AL`, `σ_att`, `σ_calib` and `u_dup_dr2` are declared
in each observation's own namespace and are not shared, even for two sources a
few arcseconds apart whose forecast pools are near-identical:

  - the losses genuinely common to both (dead time, decontamination, the OBMT
    gap lists) are *already* removed by `dr2_ok_mask`/`dr3_ok_mask` before
    priorities are applied. What `transit_priorities` marginalizes is the
    residual, source-specific loss: gating, window class, saturation,
    per-source AGIS outlier rejection. For a pair 8.5 magnitudes apart these
    differ a lot, and not monotonically — the *brighter* star can have fewer
    usable transits, because bright-star handling discards them.
  - sharing one priorities vector would force the smaller selection to be the
    top-k subset of the larger, which asserts a nesting the two selection
    mechanisms do not produce.
  - `σ_AL`, `σ_att` and `σ_calib` are magnitude-dependent calibration terms
    and must not be shared under any scheme.

If you are tempted to share anyway, note the silent failure mode:
`transit_priorities` indexes each observation's *own* forecast pool, built by
a per-source GOST query. Two queries at different sky positions need not
return the same number of rows, and if they differ by one, index `i` denotes a
different transit in each — with no error, no warning, and a likelihood that
is merely wrong.
"""
struct G23HObs{TTable,TTableH,TTableG,TCat,THip,THost,TComp,TRef} <: AbstractObs
    table::TTable                       # one row per likelihood channel
    priors::Priors
    derived::Derived
    hip_table::TTableH
    gaia_table::TTableG                 # the fixed Gaia forecast pool
    catalog::TCat
    hip_sol::THip
    A_prepared_5_hip::Matrix{Float64}
    A_prepared_5_dr2::Matrix{Float64}
    A_prepared_5_dr3::Matrix{Float64}
    # Cached weighted pseudo-inverse of the Hipparcos design matrix,
    # Q = pinv(A ./ σ) ./ σ', so the per-sample weighted five-parameter LSQ is
    # one 5×N matrix-vector product. A and σ are both fixed at construction.
    pinv_5_hip::Matrix{Float64}
    # Q_hip · (constant catalog residuals) — the all-inactive fast path's
    # Hipparcos five-parameter answer, in (Δα, Δδ, Δpmra, Δpmdec) order.
    hip_x_const::NTuple{4,Float64}
    include_iad::Bool
    ueva_mode::Symbol
    # Whether every channel is expressed relative to the *host body's* DR3
    # proper motion rather than the barycentre's. See the `frame_shift`
    # keyword on the constructor.
    frame_shift::Bool
    name::String
    host::THost
    companions::TComp
    ref::TRef
    # Epoch layout of `epochs(obs)`: Hipparcos transits, then the Gaia
    # forecast pool, then the catalog reference epochs (0 = absent).
    n_hip::Int
    n_gaia::Int
    i_ra_hip::Int
    i_dec_hip::Int
    i_ra_dr2::Int
    i_dec_dr2::Int
    i_ra_dr3::Int
    i_dec_dr3::Int
    all_epochs::Vector{Float64}
end

export G23HObs

likelihoodname(obs::G23HObs) = obs.name
epochs(obs::G23HObs) = obs.all_epochs
refspecs(obs::G23HObs) = (obs.host, obs.companions..., obs.ref)

# The generic `targets… vs ref` reading would print `(A, Ab) vs Barycentre`,
# which loses the asymmetry that matters here: the catalog source is centred
# on `host`, and the companions are bodies that *blend into* it (with their
# own per-band flux ratios) rather than co-equal targets. Shared with
# `HipparcosIADObs`, which names its refs the same way.
_blend_refdesc(host, companions, reference) = isempty(companions) ?
    _refstr(host) * " vs " * _refstr(reference) :
    _refstr(host) * " (blended with " *
    join(map(_refstr, companions), ", ") * ") vs " * _refstr(reference)

_refdesc(obs::G23HObs) = _blend_refdesc(obs.host, obs.companions, obs.ref)

# ──────────────────────────────────────────────────────────────────────
# Layout of the coupled catalog block
#
# The block is a list of *sub-blocks*, each naming the channel kinds it
# contributes in component order. Its width, each sub-block's offset, the
# `channels=` mask and the debug labels are all derived from this list, so
# adding a channel is an edit here plus a value in the `blockvals` tuple in
# `_g23h_ln_like` — never an index-arithmetic change to the covariance
# assembly, which is the most numerically delicate code in this file. (It
# used to be literal `1:2`/`3:4`/`Val(11)` arithmetic in five places, and a
# new channel meant editing all of them in step.)
#
# A sub-block is the unit *within which* the catalog quotes correlations: the
# five-parameter fits correlate α with δ, so those go together. Correlations
# *between* sub-blocks are declared separately (see `_g23h_place_cross!`) —
# today just DR2↔DR3 via `rho_dr2_dr3`.
const _g23h_blocks = (
    (name=:hip,  kinds=(:ra_hip, :dec_hip)),
    (name=:hg,   kinds=(:ra_hg, :dec_hg)),
    (name=:dr2,  kinds=(:ra_dr2, :dec_dr2)),
    (name=:dr32, kinds=(:ra_dr32, :dec_dr32)),
    (name=:dr3,  kinds=(:ra_dr3, :dec_dr3)),
    (name=:ueva, kinds=(:ueva_dr3,)),
)

# The channels that enter the coupled catalog MvNormal, in its component
# order. `:iad_hip` and `:rv_dr3` are deliberately absent: they are separate,
# uncorrelated terms.
const _g23h_channel_kinds =
    reduce((acc, b) -> (acc..., b.kinds...), _g23h_blocks; init=())
const _g23h_nchan = length(_g23h_channel_kinds)

# Zero-based offset of each sub-block into the assembled vector/matrix.
# Written as an independent sum per entry rather than a running accumulator:
# it does not depend on the iteration order of whatever builds it.
const _g23h_block_offset = ntuple(length(_g23h_blocks)) do bi
    sum(length(_g23h_blocks[k].kinds) for k in 1:(bi-1); init=0)
end

# Position of a sub-block in `_g23h_blocks`, by name — so a cross-covariance
# is declared as `_g23h_offsetof(:dr2)` rather than a literal `5`.
@inline _g23h_offsetof(name::Symbol) =
    _g23h_block_offset[findfirst(b -> b.name === name, _g23h_blocks)]

# Place a sub-block's covariance on the diagonal at the offset its position in
# `_g23h_blocks` implies. Unrolled through `ntuple(Val(…))` so that `bi` is a
# literal at each step: `blockvals` is a heterogeneous tuple (the UEVA
# sub-block is 1×1 where the rest are 2×2), and a runtime loop over it would
# be type-unstable and send the whole assembly to the heap.
@inline function _g23h_place_diag!(Σ_full, blockvals::Tuple)
    ntuple(Val(length(_g23h_blocks))) do bi
        o = _g23h_block_offset[bi]
        Σb = blockvals[bi].Σ
        @inbounds for j in axes(Σb, 2), i in axes(Σb, 1)
            Σ_full[o+i, o+j] = Σb[i, j]
        end
        nothing
    end
    return nothing
end

# Cross-covariance between two sub-blocks, written symmetrically.
@inline function _g23h_place_cross!(Σ_full, oa::Int, ob::Int, K)
    @inbounds for j in axes(K, 2), i in axes(K, 1)
        Σ_full[oa+i, ob+j] = K[i, j]
        Σ_full[ob+j, oa+i] = K[i, j]
    end
    return nothing
end

# Every channel a G23HObs can carry — what `channels=` is validated against.
const _g23h_all_kinds = (
    :iad_hip, :ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr2, :dec_dr2,
    :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3, :ueva_dr3, :rv_dr3)

"""
    _g23h_channel_rows(table, channels) -> row indices

Rows of the channel table the constructor's `channels=` keyword keeps.
`nothing` keeps everything.

A requested channel this source does not have is a warning rather than an
error: a `channels=` list is written for a *family* of sources (every HGCA
fit asks for the Hipparcos channels), and one member of that family having no
Hipparcos entry is data, not a mistake. A channel name that is not a channel
at all is a typo, and errors.
"""
function _g23h_channel_rows(table, channels)
    isnothing(channels) && return 1:length(table.kind)
    chans = Tuple(channels)
    isempty(chans) && error(
        "`channels=()` selects no channels at all. Pass `channels=nothing` (the " *
        "default) to keep every channel this source has.")
    for k in chans
        k ∈ _g23h_all_kinds || error(
            "`:$k` is not a G23HObs channel. The channels are $(_g23h_all_kinds).")
    end
    absent = filter(k -> k ∉ table.kind, collect(chans))
    isempty(absent) || @warn(
        "`channels=` asked for channels this source does not have; they are skipped.",
        absent, available = Tuple(table.kind))
    keep = findall(k -> k ∈ chans, table.kind)
    isempty(keep) && error(
        "`channels=$(chans)` selected none of this source's channels " *
        "($(Tuple(table.kind))); the observation would constrain nothing.")
    return keep
end

# The two consequences of a channel table losing rows, in one place so that
# `channels=` (construction) and `likeobj_from_epoch_subset` (cross-validation)
# cannot drift apart. Both are idempotent, so applying this to an unrestricted
# table is a no-op.
function _g23h_restrict(table, catalog, priors, derived)
    # The Hipparcos catalog proper-motion distributions exist only to be
    # compared against; with every Hipparcos channel gone there is nothing to
    # compare, and `_g23h_simulate!` reads `isnothing(dist_hip)` as "this
    # source has no Hipparcos".
    if all(k -> k ∉ table.kind, (:iad_hip, :ra_hip, :dec_hip, :ra_hg, :dec_hg))
        catalog = (; catalog..., dist_hip=nothing, dist_hg=nothing)
    end
    if :iad_hip ∉ table.kind
        # Six nuisance parameters and two derived quantities that no longer
        # touch the likelihood anywhere; sampled, they would inflate the model
        # dimension by 6 for nothing.
        iad_keys = (:hip_iad_jitter, :iad_Δra, :iad_Δdec, :iad_Δplx,
                    :iad_Δpmra, :iad_Δpmdec, :iad_pmra, :iad_pmdec)
        if any(k -> k ∈ iad_keys, keys(priors.priors)) ||
           any(k -> k ∈ iad_keys, keys(derived.variables))
            priors = Priors(OrderedDict(k => v for (k, v) in priors.priors if k ∉ iad_keys))
            derived = Derived(
                OrderedDict(k => v for (k, v) in derived.variables if k ∉ iad_keys),
                derived.captured_names, derived.captured_vals)
        end
    end
    return (table, catalog, priors, derived)
end

# ──────────────────────────────────────────────────────────────────────
# DR2/DR3 matched-transit sidecar plumbing
#
# The Gaia DR2 astrometric solution shares its input-span start with DR3 but
# stops earlier, and the two transit sets overlap without nesting: DR3's
# rebuilt crossmatch commonly used more transits inside the DR2 window than
# DR2 itself matched, and DR2 occasionally used transits DR3 later dropped.
# The epoch simulator therefore selects the two sets separately — DR3's from
# the whole AGIS span, DR2's from the DR2 window — with the DR2 selection
# sized by the DR2 matched-transit count. For bright stars that count also
# includes doubly-downlinked transits (two windows per FoV crossing, both
# tallied by DR2's crossmatch, resolved to one in (E)DR3); those are modelled
# as repeated epochs with the distinct-crossing count marginalized.
#
# The DR2/DR3 catalog-covariance correlation is NOT derived from these
# counts; it is the published `rho_dr2_dr3`, against which the catalog
# uncertainties were calibrated.
# ──────────────────────────────────────────────────────────────────────

const _G23H_DR2_SIDECAR = Ref{Any}(nothing)
const _G23H_DR2_SIDECAR_LOCK = ReentrantLock()

function _g23h_dr2_sidecar_datadep()
    lock(_G23H_DR2_SIDECAR_LOCK) do
        if isnothing(_G23H_DR2_SIDECAR[])
            dir = datadep"G23H_DR2Transits"
            feathers = filter(f -> endswith(f, ".feather"), readdir(dir; join=true))
            isempty(feathers) && error("No .feather file in the G23H_DR2Transits DataDep at $dir")
            _G23H_DR2_SIDECAR[] = Arrow.Table(first(feathers))
        end
        return _G23H_DR2_SIDECAR[]
    end
end

# Merge the DR2 matched-transit count for `gaia_id` into the catalog row. No-op
# if the catalog already carries the column; errors (no fallback) if it cannot
# be resolved — the count is what sizes the DR2 epoch selection.
function _g23h_merge_dr2_sidecar(catalog, gaia_id, dr2_transits_catalog)
    hasproperty(catalog, :astrometric_matched_observations_dr2) && return catalog
    tbl = if !isnothing(dr2_transits_catalog)
        Tables.istable(dr2_transits_catalog) ? dr2_transits_catalog :
        Arrow.Table(dr2_transits_catalog)
    else
        _g23h_dr2_sidecar_datadep()
    end
    idx = findfirst(==(gaia_id), Tables.getcolumn(tbl, :gaia_source_id))
    isnothing(idx) && error(
        "Gaia source $gaia_id was not found in the G23H_DR2Transits sidecar. The DR2 " *
        "matched-transit count (`astrometric_matched_observations_dr2`) sizes the DR2 " *
        "epoch selection; there is no fallback.")
    srow = NamedTuple(Table(tbl)[idx])
    # The published sidecar carries the count under Gaia DR2's native column
    # name (`astrometric_matched_observations`, renamed
    # `astrometric_matched_transits` in (E)DR3); accept either spelling.
    n_dr2 = if haskey(srow, :astrometric_matched_observations_dr2)
        srow.astrometric_matched_observations_dr2
    elseif haskey(srow, :astrometric_matched_observations)
        srow.astrometric_matched_observations
    else
        error("The G23H_DR2Transits sidecar is missing the DR2 matched-observation count.")
    end
    return merge(catalog, (; astrometric_matched_observations_dr2=n_dr2))
end

# Total size of the DR2 epoch selection. Includes DR2's doubly-downlinked
# bright-star transits, so it can exceed the number of distinct crossings.
function _g23h_dr2_target_transits(catalog, n_dr3::Integer)
    hasproperty(catalog, :astrometric_matched_observations_dr2) || error(
        "G23HObs requires the Gaia DR2 matched-transit count " *
        "(`astrometric_matched_observations_dr2`), from the G23H_DR2Transits sidecar " *
        "or the `dr2_transits_catalog` keyword; it was not found for this source.")
    v = catalog.astrometric_matched_observations_dr2
    (ismissing(v) || !isfinite(v)) && error(
        "Gaia DR2 `astrometric_matched_observations_dr2` is missing/non-finite for this source.")
    n_dr2 = round(Int, v)
    if n_dr2 > 3 * n_dr3
        # Modest excess over n_dr3 is genuine; a several-fold excess suggests
        # the column holds CCD-level AL observations (~9 per FoV transit).
        @warn "DR2 matched transits ($n_dr2) far exceed DR3 matched transits ($n_dr3); the " *
              "sidecar column may be CCD-level AL observations rather than FoV transits."
    end
    return max(n_dr2, 0)
end

"""
    _g23h_select_dr2_epochs(priorities, pool, n_distinct, n_total)

Materialize the DR2 epoch set from the sampled transit priorities:
`n_distinct` distinct crossings (the top-priority epochs of the DR2-window
pool) padded with repeats up to `n_total` entries. A repeated index puts that
epoch into the DR2 least-squares fit twice — the model of a doubly-downlinked
bright-star transit that DR2's matched-observation count tallied twice. Which
crossings carry the repeats is taken from the top of the priority ordering;
by exchangeability of the iid priorities that is distributed as a uniform
choice among the selected crossings. The returned length is always `n_total`,
so chain storage and likelihood buffers stay rectangular.
"""
function _g23h_select_dr2_epochs(priorities::AbstractVector, pool::AbstractVector{<:Integer},
                                 n_distinct::Integer, n_total::Integer)
    (n_distinct <= 0 || n_total <= 0) && return Int[]
    sel = pool[partialsortperm(priorities[pool], 1:n_distinct, rev=true)]  # priority order
    n_rep = n_total - n_distinct
    # mod1 wrap covers the extreme case n_rep > n_distinct (multiplicity ≥ 3).
    reps = n_rep > 0 ? sel[mod1.(1:n_rep, n_distinct)] : Int[]
    return sort!(vcat(sel, reps))
end

# ──────────────────────────────────────────────────────────────────────
# Construction
# ──────────────────────────────────────────────────────────────────────

function G23HObs(;
        host,
        companions=(),
        ref=Barycentre,
        gaia_id=nothing,
        hip_id=nothing,
        catalog=nothing,
        forecast_table=nothing,
        scanlaw_table=nothing,
        hipparcos=nothing,
        variables::Union{Nothing,Tuple{Priors,Derived}}=nothing,
        channels=nothing,
        include_rv=true,
        include_iad::Bool=false,
        rv_ln_uncert_err_floor::Union{Nothing,Real}=0.30,
        ueva_mode::Symbol=:RUWE,
        frame_shift::Bool=true,
        freeze_epochs::Bool=false,
        dr2_transits_catalog=nothing,
        # G-band threshold below which the DR2 duplicate-transit count is
        # marginalized (doubly-downlinked bright-star transits). -Inf disables.
        dr2_dup_gmag_threshold::Real=6.5,
        name::AbstractString="G23H",
    )
    ueva_mode in (:RUWE, :EAN, :none) ||
        error("`ueva_mode` must be :RUWE, :EAN or :none; got :$ueva_mode")
    isnothing(gaia_id) && isnothing(hip_id) &&
        error("Either `gaia_id` or `hip_id` must be specified")
    !isnothing(gaia_id) && !isnothing(hip_id) &&
        error("Specify either `gaia_id` or `hip_id`, not both")

    hostspec = refspec(host)
    compspecs = map(refspec, Tuple(companions))
    refspec_ = refspec(ref)

    catalog = _g23h_catalog_row(catalog, gaia_id, hip_id)
    if isnothing(gaia_id) && hasproperty(catalog, :gaia_source_id)
        gaia_id = catalog.gaia_source_id
    end

    # Catalog epochs are Julian years, not decimal years (T. Brandt, priv. comm.).
    J2000_mjd = 51544.5
    _yr2mjd(y) = (y - 2000) * julian_year + J2000_mjd
    catalog = (; catalog...,
        epoch_ra_hip_mjd=_yr2mjd(catalog.epoch_ra_hip),
        epoch_dec_hip_mjd=_yr2mjd(catalog.epoch_dec_hip),
        epoch_ra_dr2_mjd=_yr2mjd(catalog.epoch_ra_dr2),
        epoch_dec_dr2_mjd=_yr2mjd(catalog.epoch_dec_dr2),
        epoch_ra_dr3_mjd=_yr2mjd(catalog.epoch_ra_dr3),
        epoch_dec_dr3_mjd=_yr2mjd(catalog.epoch_dec_dr3),
    )

    if !hasproperty(catalog, :astrometric_chi2_al_dr3) || !hasproperty(catalog, :rv_nb_transits)
        @warn "Column missing from catalog; querying the Gaia DR3 TAP server (or its cache)"
        dr3 = _query_gaia_dr3(; gaia_id)
        catalog = (; catalog...,
            astrometric_chi2_al_dr3=dr3.astrometric_chi2_al,
            parallax_error=dr3.parallax_error,
            rv_nb_transits=dr3.rv_nb_transits,
            radial_velocity_error=dr3.radial_velocity_error)
    end

    if !isnothing(gaia_id)
        catalog = _g23h_merge_dr2_sidecar(catalog, gaia_id, dr2_transits_catalog)
    end

    # Floor the uncertainty on the per-transit RV ln-σ calibration. The GP
    # behind `rv_ln_uncert_err_dr3` is Malmquist-biased for cool dwarfs —
    # their (colour, magnitude) calibration bins are dominated by distant
    # giants with sharper lines — so the quoted uncertainty on ln σ_rv is
    # unrealistically small and over-constrains the σ_rv_per_transit prior.
    if !isnothing(rv_ln_uncert_err_floor) && hasproperty(catalog, :rv_ln_uncert_err_dr3) &&
       !ismissing(catalog.rv_ln_uncert_err_dr3) && isfinite(catalog.rv_ln_uncert_err_dr3)
        catalog = merge(catalog, (;
            rv_ln_uncert_err_dr3=max(catalog.rv_ln_uncert_err_dr3, rv_ln_uncert_err_floor)))
    end

    # --- Hipparcos -----------------------------------------------------------
    has_hip = !(catalog.hip_id isa Number && isnan(catalog.hip_id))
    if !has_hip
        @info "No Hipparcos data for this source; the Hipparcos channels are dropped."
        hip_table = Table(_g23h_empty_hip_table())
        hip_sol = nothing
        dist_hip = nothing
        dist_hg = nothing
        A_prepared_5_hip = zeros(0, 0)
    else
        loaded = isnothing(hipparcos) ?
                 hipparcos_iad(; hip_id=Int(catalog.hip_id)) : hipparcos
        # Copy the columns: `FlexTable(::Table)` wraps the *same* vectors, and
        # the recalibration below mutates them. Without this, a caller who
        # passes a preloaded `hipparcos=` table gets it shifted underneath
        # them, and a second observation built from it is shifted twice.
        hip_table = FlexTable(map(copy, Tables.columntable(loaded.table)))
        hip_sol = loaded.hip_sol
        # Brandt, Michalik & Brandt: +0.140 mas on the residuals and 2.25 mas
        # of extra dispersion, which mitigates overfitting. This also
        # recomputes `proj_meas_alongscan`, which is derived from `res`.
        hipparcos_recalibrate!(hip_table)
        hip_table = Table(hip_table)
        A_prepared_5_hip = prepare_A_5param(hip_table.cosϕ, hip_table.sinϕ, hip_table.epoch,
            hip_table.parallaxFactorAlongScan,
            catalog.epoch_ra_hip_mjd, catalog.epoch_dec_hip_mjd)
        dist_hip = _g23h_pm_dist(catalog.pmra_hip, catalog.pmdec_hip,
            catalog.pmra_hip_error[1], catalog.pmdec_hip_error[1], catalog.pmra_pmdec_hip[1])
        dist_hg = _g23h_pm_dist(catalog.pmra_hg, catalog.pmdec_hg,
            catalog.pmra_hg_error[1], catalog.pmdec_hg_error[1], catalog.pmra_pmdec_hg[1])
    end

    dist_dr2 = _g23h_pm_dist(catalog.pmra_dr2, catalog.pmdec_dr2,
        catalog.pmra_dr2_error[1], catalog.pmdec_dr2_error[1], catalog.pmra_pmdec_dr2[1])
    dist_dr32 = _g23h_pm_dist(catalog.pmra_dr32, catalog.pmdec_dr32,
        catalog.pmra_dr32_error[1], catalog.pmdec_dr32_error[1], catalog.pmra_pmdec_dr32[1])
    dist_dr3 = _g23h_pm_dist(catalog.pmra_dr3, catalog.pmdec_dr3,
        catalog.pmra_dr3_error[1], catalog.pmdec_dr3_error[1], catalog.pmra_pmdec_dr3[1])
    catalog = (; catalog..., dist_hip, dist_hg, dist_dr2, dist_dr32, dist_dr3)

    # --- the Gaia forecast pool ---------------------------------------------
    gaia_table, dr2_ok_mask, dr3_ok_mask =
        _g23h_forecast_pool(catalog, forecast_table, scanlaw_table)

    A_prepared_5_dr2 = prepare_A_5param(gaia_table.cosϕ, gaia_table.sinϕ, gaia_table.epoch,
        gaia_table.parallaxFactorAlongScan, catalog.epoch_ra_dr2_mjd, catalog.epoch_dec_dr2_mjd)
    A_prepared_5_dr3 = prepare_A_5param(gaia_table.cosϕ, gaia_table.sinϕ, gaia_table.epoch,
        gaia_table.parallaxFactorAlongScan, catalog.epoch_ra_dr3_mjd, catalog.epoch_dec_dr3_mjd)

    # --- the channel table ---------------------------------------------------
    has_rv = include_rv && hasproperty(catalog, :rv_ln_uncert_dr3) &&
             !ismissing(catalog.rv_ln_uncert_dr3) && !ismissing(catalog.rv_ln_uncert_err_dr3) &&
             isfinite(catalog.rv_ln_uncert_dr3) && isfinite(catalog.rv_ln_uncert_err_dr3)
    table = _g23h_channel_table(catalog, hip_table, gaia_table, has_hip, has_rv, ueva_mode)
    # `channels=` and `likeobj_from_epoch_subset` are the same operation on the
    # same rows, so they share the row selection *and* its two consequences.
    table = table[_g23h_channel_rows(table, channels)]

    # --- default variables ---------------------------------------------------
    if isnothing(variables)
        # Defaults follow the *restricted* table: a dropped `:rv_dr3` must not
        # leave `σ_rv_per_transit` sampled but unused.
        variables = _g23h_default_variables(catalog, gaia_table, dr2_ok_mask, dr3_ok_mask,
            has_hip, :rv_dr3 ∈ table.kind, hip_sol, table, freeze_epochs,
            dr2_dup_gmag_threshold, ueva_mode)
    end
    (priors, derived) = variables
    table, catalog, priors, derived = _g23h_restrict(table, catalog, priors, derived)

    # §11: a row with a Hipparcos entry but null perspective-acceleration terms
    # would have silently NaN'd the likelihood under an absolute frame. Caught
    # here, where the catalogue row can be named, rather than as a NaN log
    # density inside the sampler.
    _nl(k) = hasproperty(catalog, k) ? getproperty(catalog, k) : missing
    if has_hip && any(k -> k in (:ra_hg, :dec_hg, :ra_hip, :dec_hip), table.kind) &&
       !all(v -> v isa Real && isfinite(v), (_nl(:nonlinear_dpmra), _nl(:nonlinear_dpmdec)))
        error("G23HObs \"$name\": Gaia DR3 $(gaia_id) has Hipparcos data, so its " *
              "Hipparcos and Hipparcos–Gaia proper-motion channels carry the HGCA " *
              "perspective-acceleration correction — but the catalog row's " *
              "`nonlinear_dpmra`/`nonlinear_dpmdec` are " *
              "$(_nl(:nonlinear_dpmra))/$(_nl(:nonlinear_dpmdec)), not finite " *
              "numbers. Under an absolute frame that correction would be added to the " *
              "model proper motions. Supply the terms in the catalog row, or drop the " *
              "affected channels with " *
              "`channels=(:ra_dr2, :dec_dr2, :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3, " *
              ":ueva_dr3)`.")
    end

    pinv_5_hip = _g23h_pinv_hip(A_prepared_5_hip, hip_table)
    hip_x_const = _g23h_hip_x_const(pinv_5_hip, hip_table, include_iad, has_hip)

    n_hip = size(hip_table, 1)
    n_gaia = length(gaia_table.epoch)
    all_ep = Float64[]
    append!(all_ep, hip_table.epoch)
    append!(all_ep, gaia_table.epoch)
    k = n_hip + n_gaia
    i_ra_hip = i_dec_hip = 0
    if has_hip
        push!(all_ep, catalog.epoch_ra_hip_mjd); i_ra_hip = (k += 1)
        push!(all_ep, catalog.epoch_dec_hip_mjd); i_dec_hip = (k += 1)
    end
    push!(all_ep, catalog.epoch_ra_dr2_mjd); i_ra_dr2 = (k += 1)
    push!(all_ep, catalog.epoch_dec_dr2_mjd); i_dec_dr2 = (k += 1)
    push!(all_ep, catalog.epoch_ra_dr3_mjd); i_ra_dr3 = (k += 1)
    push!(all_ep, catalog.epoch_dec_dr3_mjd); i_dec_dr3 = (k += 1)

    return G23HObs{typeof(table),typeof(hip_table),typeof(gaia_table),typeof(catalog),
        typeof(hip_sol),typeof(hostspec),typeof(compspecs),typeof(refspec_)}(
        table, priors, derived, hip_table, gaia_table, catalog, hip_sol,
        A_prepared_5_hip, A_prepared_5_dr2, A_prepared_5_dr3,
        pinv_5_hip, hip_x_const, include_iad, ueva_mode, frame_shift, String(name),
        hostspec, compspecs, refspec_,
        n_hip, n_gaia, i_ra_hip, i_dec_hip, i_ra_dr2, i_dec_dr2, i_ra_dr3, i_dec_dr3,
        all_ep)
end

# The multi-source guard. `frame_shift` redefines the system's `pmra`/`pmdec`
# to mean "this observation's host body, at DR3", which two observations
# cannot both do — they would each subtract their own Δpm from the one shared
# parameter and predict identical DR3 proper motions for two sources the
# catalogue separates by many σ. That is not a degradation, it is the deletion
# of the wide-orbit signal, so it errors.
#
# Keyed on `isempty(framevars)`, not on the observation count alone: with no
# system frame each `G23HObs` supplies its own `pmra`/`pmdec` in its own
# namespace, nothing is shared, and the shift is per-observation in fact as
# well as in name.
function check_siblings(obs::G23HObs, @nospecialize(all_obs), ctx)
    others = G23HObs[o for o in all_obs if o isa G23HObs && o !== obs]
    isempty(others) && return nothing
    _g23h_check_frame_shift(obs, others, ctx)
    _g23h_check_shared_fluxratios(obs, others, ctx)
    return nothing
end

function _g23h_check_frame_shift(obs::G23HObs, others::Vector{G23HObs}, ctx)
    obs.frame_shift || return nothing
    isempty(ctx.framevars) && return nothing
    error("""
    System $(ctx.name): "$(obs.name)" has `frame_shift=true`, but it shares the \
    system's frame with $(length(others)) other G23HObs \
    ($(join(('"' * o.name * '"' for o in others), ", "))).

    `frame_shift` expresses every channel relative to *this* observation's host body \
    at the DR3 epoch, which redefines the system-level `pmra`/`pmdec`. Two \
    observations cannot both redefine the same parameter: each would subtract its own \
    Δpm, so both sources' DR3 channels would predict one and the same proper motion, \
    and the relative proper motion — the entire wide-orbit signal — would be removed \
    from the likelihood rather than merely weakened.

    Pass `frame_shift=false` to every G23HObs in this system, and recover the sampling \
    efficiency it was buying with `AnchoredFrame`, which reconditions the same \
    degeneracy in the model instead of inside one observation.
    """)
end

# §10: the tier-2 flux-ratio fallthrough is keyed by the bare name, so a
# system-level `fluxratio` is read by *every* G23HObs, each validating it
# against its own companion count. Two observations with different counts can
# therefore never both be right — and the runtime failure ("received 3 flux
# ratios but the observation declares 1 companion") names neither of them.
#
# Statically decidable: which tier an observation lands on is a question about
# variable *names*, and the counts are structural.
function _g23h_check_shared_fluxratios(obs::G23HObs, others::Vector{G23HObs}, ctx)
    for key in (:fluxratio, :fluxratio_hip)
        _g23h_tier2(obs, key, ctx) || continue
        clash = [o for o in others if _g23h_tier2(o, key, ctx) &&
                                      length(o.companions) != length(obs.companions)]
        isempty(clash) && continue
        error("""
        System $(ctx.name): "$(obs.name)" declares $(length(obs.companions)) \
        companion(s) and "$(clash[1].name)" declares $(length(clash[1].companions)) \
        companion(s), but both fall through to the system-level `$key` vector — which is keyed by name \
        alone, so both read the same values and each checks them against its own count. \
        At most one can be right.

        Give each observation its own `$key` in its `variables=` block, or drop the \
        system-level one and let the flux ratios come from the bodies' own \
        `flux_<band>` variables (tier 3), which is per-body and so cannot be \
        ambiguous.
        """)
    end
    return nothing
end

# Tier 1 is the observation's own namespace; tier 2 is the system's. An
# observation that declares the key itself never reaches tier 2.
_g23h_tier2(obs::G23HObs, key::Symbol, ctx) =
    !(key in keys(obs.priors.priors)) && !(key in keys(obs.derived.variables)) &&
    key in ctx.sysnames

# ──────────────────────────────────────────────────────────────────────
# The full 5×5, fetched per source (§7)
# ──────────────────────────────────────────────────────────────────────

"""
    _g23h_solution_5x5(catalog, gaia_id; release) -> SMatrix{5,5} or nothing

That source's Gaia astrometric covariance, in (α*, δ, ϖ, μα*, μδ) order at the
release's reference epoch, fetched from the archive and cached per source.
Returns `nothing` if the fetch is unavailable — position/parallax/mean-RV
channels then refuse to build, with a message saying so; the proper-motion
channels need none of this and are unaffected.

**The fetch is checked against the catalog, not merely against HTTP 200.**
The catalog's own columns are reconstructible from this matrix:

  - `epoch_ra_dr3` is the central epoch, `ref − cov(α*,μα*)/var(μα*)`;
  - `ra_error_central_dr3` is σ(α*) propagated to that epoch, times
    [`G23H_ERROR_INFLATION`](@ref);
  - `pmra_dr3_error` is σ(μα*) times the same factor;
  - `ra_dec_corr_central_dr3` and `pmra_pmdec_dr3` are the raw correlations.

so we rebuild them and compare. A mismatch means the catalog row and the
archive row describe *different solutions* — a wrong source id, a re-reduction,
DR3 vs DR3-lite — and mixing the two would produce a covariance that is
internally inconsistent in a way nothing downstream could detect. Verified on
both ups And sources: central epochs to 2e-8 yr, every error ratio exactly
1.370000.
"""
function _g23h_solution_5x5(catalog, gaia_id; release::Symbol=:dr3, verbose::Bool=true)
    isnothing(gaia_id) && return nothing
    sol = try
        release === :dr3 ? _query_gaia_dr3(; gaia_id) : _query_gaia_dr2(; gaia_id)
    catch err
        err isa InterruptException && rethrow()
        verbose && @warn "Could not fetch the Gaia $(uppercase(String(release))) astrometric solution for $gaia_id; position, parallax and mean-RV channels are unavailable for this source." exception = err
        return nothing
    end
    C = try
        _gaia_5x5(sol)
    catch err
        err isa InterruptException && rethrow()
        verbose && @warn "Gaia $(uppercase(String(release))) solution for $gaia_id has no usable 5×5 covariance." exception = err
        return nothing
    end
    isposdef(C) || error(
        "the Gaia $(uppercase(String(release))) 5×5 covariance for source $gaia_id is not " *
        "positive definite (minimum eigenvalue $(minimum(eigvals(C)))). The archive's " *
        "errors and correlations disagree; refusing to build a likelihood on it.")
    release === :dr3 && _g23h_check_5x5(C, catalog, gaia_id)
    return C
end

# The reconstruction gate described above. Tolerances are loose next to the
# agreement actually measured (1e-8 yr, exactly 1.370000) and tight next to any
# real mismatch, which would be tens of percent.
function _g23h_check_5x5(C, catalog, gaia_id)
    ref_yr = 2000 + (meta_gaia_DR3.ref_epoch_mjd - 51544.5) / julian_year
    for (coord, epochkey, σkey, pmσkey) in (
            (1, :epoch_ra_dr3, :ra_error_central_dr3, :pmra_dr3_error),
            (2, :epoch_dec_dr3, :dec_error_central_dr3, :pmdec_dr3_error))
        hasproperty(catalog, epochkey) || continue
        t_c = _central_epoch(C, coord)
        # `epoch_*_dr3` is stored as Float32 in the catalog, so ~1e-5 yr of
        # representation error is expected and is not a mismatch.
        Δ = ref_yr + t_c - Float64(getproperty(catalog, epochkey))
        abs(Δ) < 1e-3 || _err_5x5_mismatch(gaia_id, "central epoch $epochkey",
            ref_yr + t_c, Float64(getproperty(catalog, epochkey)))
        Cc = _propagate_5x5(C, t_c)
        _check_ratio(gaia_id, σkey, sqrt(Cc[coord, coord]),
            Float64(getproperty(catalog, σkey)))
        _check_ratio(gaia_id, pmσkey, sqrt(C[coord+3, coord+3]),
            Float64(getproperty(catalog, pmσkey)))
    end
    # Parallax carries no inflation — the asymmetry this whole check exists to
    # pin down, since it is the one a position/parallax channel gets wrong by
    # analogy.
    hasproperty(catalog, :parallax_error) &&
        _check_ratio(gaia_id, :parallax_error, sqrt(C[3, 3]),
            Float64(catalog.parallax_error), 1.0)
    return nothing
end

function _check_ratio(gaia_id, key, raw, stored, expect=G23H_ERROR_INFLATION)
    (isfinite(raw) && isfinite(stored) && raw > 0) || return nothing
    r = stored / raw
    isapprox(r, expect; rtol=1e-3) ||
        _err_5x5_mismatch(gaia_id, "$key (expected $(expect)× the archive value, got $(r)×)",
            raw * expect, stored)
    return nothing
end

@noinline _err_5x5_mismatch(gaia_id, what, rebuilt, stored) = error("""
G23HObs: the Gaia archive solution fetched for source $gaia_id does not reproduce \
the G23H catalog row's $what — rebuilt $rebuilt, catalog has $stored.

The catalog's columns are derived from the archive's five-parameter solution \
(central epochs from the position/proper-motion covariance, errors inflated by \
$(G23H_ERROR_INFLATION)), so they must agree. That they do not means the two rows \
describe different solutions: a wrong source id, a re-reduction, or a DR3 vs \
DR3-lite mismatch. Combining them would give a covariance that is internally \
inconsistent in a way no downstream check could catch, so this refuses rather \
than proceeding.

If you do not need the position, parallax or mean-RV channels, drop them from \
`channels=` — the proper-motion channels use only the catalog row and are \
unaffected by this.
""")

# HGCA's nonlinear (perspective-acceleration) correction, as the pair of
# additive proper-motion offsets (HG epoch, Hipparcos epoch). The factor of
# two on the Hipparcos epoch is because `nonlinear_dpm*` is defined to the HG
# epoch.
#
# Catalog provenance, so nobody "fixes" this layer: the G23H catalog's DR2,
# DR32 and DR3 columns carry NO baked-in perspective correction, while the
# `_hip` and `_hg` columns come verbatim from HGCA v2, which DOES apply one —
# and records how much in `nonlinear_dpmra`/`nonlinear_dpmdec`. Applying
# these offsets undoes HGCA's correction, so all five channels flow through
# the same uncorrected code path (the model propagates perspective itself).
# This is NOT a double count, and the signs are right — verified numerically
# at Barnard's star (audit 2026-08-14, FINDINGS.md layer 2): the model spread
# (μ_h − μ_dr3) = (+2.529, −32.058) mas/yr plus the applied 2·dpm =
# (−2.531, +32.079) cancels to 0.02 mas/yr, and the HG version to 5e-4.
# Removing or re-signing this layer breaks that cancellation by ~32 mas/yr.
# Side effect worth knowing: through this term the channel data pin the
# system rv at ~0.29 mas/yr per km/s in pmdec.
#
# The `isfinite` test is the point of the helper (§11). The guard used to be
# `absolute && !isnothing(dist_hip)`, and the two catalogue rows anyone had
# checked happened to have finite `nonlinear_dpm*` exactly where they had
# Hipparcos — so a source with a Hipparcos entry but null nonlinear terms
# would have added `NaN` to `μ_hg` and NaN'd the whole log-likelihood, with
# nothing in the message to say why. `G23HObs`'s constructor rejects that row
# outright; this stays as the defence for a struct assembled by other means.
@inline function _g23h_nonlinear(cat, absolute::Bool, has_hip::Bool)
    z = zero(SVector{2,Float64})
    (absolute && has_hip) || return (z, z)
    (isfinite(cat.nonlinear_dpmra) && isfinite(cat.nonlinear_dpmdec)) || return (z, z)
    d = SVector{2,Float64}(cat.nonlinear_dpmra, cat.nonlinear_dpmdec)
    return (d, 2d)
end

function _g23h_pm_dist(μ_ra, μ_dec, σ_ra, σ_dec, ρ)
    c = ρ * σ_ra * σ_dec
    return MvNormal(@SVector([μ_ra, μ_dec]), @SArray [σ_ra^2 c; c σ_dec^2])
end

function _g23h_catalog_row(catalog, gaia_id, hip_id)
    isnothing(catalog) && (catalog = joinpath(datadep"G23H_Catalog", "G23H-v1.0.feather"))
    if catalog isa NamedTuple
        return catalog
    elseif Tables.istable(catalog)
        return _g23h_row_from_table(catalog, gaia_id, hip_id)
    else
        return _g23h_row_from_table(Arrow.Table(catalog), gaia_id, hip_id)
    end
end

function _g23h_row_from_table(t, gaia_id, hip_id)
    if !isnothing(hip_id)
        matches = findall(==(hip_id), Tables.getcolumn(t, :hip_id))
        isempty(matches) && error("Hipparcos ID $hip_id was not found in the catalog.")
        gaia_id = Tables.getcolumn(t, :gaia_source_id)[matches[1]]
        @info "Resolved HIP $hip_id to Gaia DR3 source ID $gaia_id"
    end
    idx = findfirst(==(gaia_id), Tables.getcolumn(t, :gaia_source_id))
    isnothing(idx) && error("Gaia source ID $gaia_id was not found in the catalog.")
    return NamedTuple(Table(t)[idx])
end

_g23h_empty_hip_table() = (; iorb=Int[], epoch_yrs=Float64[], parf=Float64[],
    cosϕ=Float64[], sinϕ=Float64[], res=Float64[], sres=Float64[], reject=Bool[],
    sres_renorm=Float64[], epoch=Float64[], x=Float64[], y=Float64[], z=Float64[],
    rv_kms=Float64[], plx_vs_time=Float64[], Δα✱=Float64[], Δδ=Float64[],
    scanAngle_rad=Float64[], parallaxFactorAlongScan=Float64[],
    proj_meas_alongscan=Float64[])

"""
Build the fixed Gaia forecast pool and the per-release usability masks.

The pool is trimmed to the (E)DR3 AGIS input span — GOST forecasts start at
the beginning of science operations and can run well past the DR3 data span,
and out-of-span epochs can never contribute to any modelled channel but would
consume `transit_priorities` selection slots.

Known data gaps are then applied **per release**: the DR2 and (E)DR3 gap
lists differ inside the DR2 window by ~13 days of DR2-valid time that only
the (E)DR3 processing excluded, so removing the union from one shared pool
undercounts the epochs DR2 could use. Rows dead for both releases are
dropped. Non-persistent DR2 entries (lunar eclipses, say) stay
DR2-selectable: empirically DR2's AGIS did use those transits.
"""
function _g23h_forecast_pool(catalog, forecast_table, scanlaw_table)
    ft = if !isnothing(forecast_table)
        t = FlexTable(forecast_table)
        for col in (:epoch, :scanAngle_rad, :parallaxFactorAlongScan)
            hasproperty(t, col) || error(
                "`forecast_table` needs an `$col` column (it has " *
                "$(Tables.columnnames(t)))")
        end
        t
    elseif !isnothing(scanlaw_table)
        @info "Using the supplied `scanninglaw` table; GOST will not be queried."
        t = FlexTable(scanlaw_table)
        t.epoch = tcb_at_gaia_2mjd.(t.times)
        t.scanAngle_rad = deg2rad.(t.angles)
        e = FlexTable(geocentre_position_query.(t.epoch))
        f = @. e.x * sind(catalog.ra) - e.y * cosd(catalog.ra)
        g = @. e.x * cosd(catalog.ra) * sind(catalog.dec) +
               e.y * sind(catalog.ra) * sind(catalog.dec) - e.z * cosd(catalog.dec)
        t.parallaxFactorAlongScan = @. f * sin(t.scanAngle_rad) + g * cos(t.scanAngle_rad)
        t
    else
        t = FlexTable(GOST_forecast(catalog.ra, catalog.dec))
        t.epoch = jd2mjd.(t.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        t.scanAngle_rad = t.scanAngle_rad_
        t
    end
    # Hipparcos' scan-angle convention: ϕ = π/2 + scanAngle.
    ft.cosϕ = cos.(π / 2 .+ ft.scanAngle_rad)
    ft.sinϕ = sin.(π / 2 .+ ft.scanAngle_rad)

    # `FlexTable` row masking returns a 2-D table; flatten back to one row per
    # transit so every column is a plain vector.
    _flatten(t) = Table(map(vec, Tables.columntable(t)))
    tbl = _flatten(ft[gaia_agis_span_dr3.start_mjd .<= vec(ft.epoch) .<= gaia_agis_span_dr3.stop_mjd, :])

    gaps_dr2 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiadr2_08252020.csv"), FlexTable)
    gaps_edr23 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiaedr3_12232020.csv"), FlexTable)
    _persistent(x) = x === true || (x isa AbstractString && uppercase(strip(x)) == "TRUE")
    dr2_hard = _persistent.(gaps_dr2.persistent)
    gap_starts_dr2 = obmt2mjd.(gaps_dr2.start[dr2_hard])
    gap_stops_dr2 = obmt2mjd.(gaps_dr2.end[dr2_hard])
    gap_starts_dr3 = obmt2mjd.(gaps_edr23.start)
    gap_stops_dr3 = obmt2mjd.(gaps_edr23.end)
    _in_gaps(e, starts, stops) = any(s <= e <= p for (s, p) in zip(starts, stops))

    ep = vec(tbl.epoch)
    dr2_ok = (ep .<= gaia_agis_span_dr2.stop_mjd) .&
             .!_in_gaps.(ep, Ref(gap_starts_dr2), Ref(gap_stops_dr2))
    dr3_ok = .!_in_gaps.(ep, Ref(gap_starts_dr3), Ref(gap_stops_dr3))
    keep = dr2_ok .| dr3_ok
    all(keep) || @info "Removed forecast transits in data gaps of every applicable release." n_removed = count(.!keep)
    n_dr2_only = count(dr2_ok .& .!dr3_ok)
    n_dr2_only > 0 && @info "Keeping transits usable by DR2 but excluded by the EDR3 gap list." n = n_dr2_only
    return _flatten(tbl[keep, :]), dr2_ok[keep], dr3_ok[keep]
end

# One row per likelihood channel. This table is what `likeobj_from_epoch_subset`
# selects on, and what the plotting layer reads; the actual reductions run
# against the prepared design matrices.
function _g23h_channel_table(catalog, hip_table, gaia_table, has_hip, has_rv, ueva_mode)
    # The Hipparcos and Hipparcos–Gaia rows are built unconditionally and
    # dropped at the end (`has_hip || append!(drop, 1:5)`), so for a Gaia-only
    # source these conversions run on the catalog's `NaN` epoch fields — and
    # `years2mjd(NaN)` throws `InexactError: Int64(NaN)` from `Dates.Date`,
    # long before the drop. v8 guarded each conversion individually; folding
    # all eleven into one `_e` helper during the v9 port dropped the guard,
    # which made `G23HObs` unconstructible for *any* source with no Hipparcos
    # entry (ups And B, and every faint secondary like it).
    #
    # Keyed on the value rather than on `has_hip`: that also covers a source
    # that has a `hip_id` but non-finite Hipparcos–Gaia epochs, where the row
    # is kept and a NaN epoch is a loud downstream failure rather than a
    # construction-time crash naming the wrong thing.
    _e(y) = isfinite(y) ? years2mjd(y) : oftype(float(y), NaN)
    hipmean = has_hip ? mean(hip_table.epoch) : NaN
    hipmin = has_hip ? minimum(hip_table.epoch) : 0.0
    hipmax = has_hip ? maximum(hip_table.epoch) : 0.0
    ep_dr2_first = first(gaia_table.epoch[gaia_table.epoch .>= gaia_agis_span_dr2.start_mjd])
    ep_dr2_last = last(gaia_table.epoch[gaia_table.epoch .<= gaia_agis_span_dr2.stop_mjd])
    ep_dr3_first = first(gaia_table.epoch[gaia_table.epoch .>= gaia_agis_span_dr3.start_mjd])
    ep_dr3_last = last(gaia_table.epoch[gaia_table.epoch .<= gaia_agis_span_dr3.stop_mjd])

    table = FlexTable(
        epoch=[hipmean, _e(catalog.epoch_ra_hip), _e(catalog.epoch_dec_hip),
            _e(catalog.epoch_ra_hg), _e(catalog.epoch_dec_hg),
            _e(catalog.epoch_ra_dr2), _e(catalog.epoch_dec_dr2),
            _e(catalog.epoch_ra_dr32), _e(catalog.epoch_dec_dr32),
            _e(catalog.epoch_ra_dr3), _e(catalog.epoch_dec_dr3),
            _e((catalog.epoch_dec_dr3 + catalog.epoch_dec_dr2) / 2)],
        start_epoch=[hipmin, hipmin, hipmin, _e(catalog.epoch_ra_hip), _e(catalog.epoch_dec_hip),
            ep_dr2_first, ep_dr2_first, _e(catalog.epoch_ra_dr2), _e(catalog.epoch_dec_dr2),
            ep_dr3_first, ep_dr3_first, ep_dr3_first],
        stop_epoch=[hipmax, hipmax, hipmax, _e(catalog.epoch_ra_dr3), _e(catalog.epoch_dec_dr3),
            ep_dr2_last, ep_dr2_last, _e(catalog.epoch_ra_dr3), _e(catalog.epoch_dec_dr3),
            ep_dr3_last, ep_dr3_last, ep_dr3_last],
        pm=[NaN, catalog.pmra_hip, catalog.pmdec_hip, catalog.pmra_hg, catalog.pmdec_hg,
            catalog.pmra_dr2, catalog.pmdec_dr2, catalog.pmra_dr32, catalog.pmdec_dr32,
            catalog.pmra_dr3, catalog.pmdec_dr3,
            ueva_mode == :RUWE ? catalog.ruwe_dr3 : catalog.astrometric_excess_noise_dr3],
        σ_pm=[NaN, catalog.pmra_hip_error, catalog.pmdec_hip_error,
            catalog.pmra_hg_error, catalog.pmdec_hg_error,
            catalog.pmra_dr2_error, catalog.pmdec_dr2_error,
            catalog.pmra_dr32_error, catalog.pmdec_dr32_error,
            catalog.pmra_dr3_error, catalog.pmdec_dr3_error, NaN],
        kind=[:iad_hip, :ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr2, :dec_dr2,
            :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3, :ueva_dr3],
    )
    if has_rv
        push!(table.epoch, mean(gaia_table.epoch))
        push!(table.start_epoch, first(gaia_table.epoch))
        push!(table.stop_epoch, last(gaia_table.epoch))
        push!(table.pm, NaN)
        push!(table.σ_pm, NaN)
        push!(table.kind, :rv_dr3)
    end
    tab = Table(table)
    drop = Int[]
    has_hip || append!(drop, 1:5)
    # :none is not a UEVA model: the datum and the DR3 covariance deflation it
    # drives are both dropped, and every other channel is untouched.
    ueva_mode === :none && push!(drop, 12)
    isempty(drop) && return tab
    return tab[setdiff(1:length(tab.kind), drop)]
end

function _g23h_pinv_hip(A, hip_table)
    isempty(A) && return zeros(0, 0)
    # `sres ≤ 0` marks scans the Hipparcos reduction rejected; give them zero
    # weight (σ = Inf) so this catalog-reproduction LSQ matches what the
    # catalog actually fit. A raw negative sres would sneak the scan back in
    # at weight 1/|sres|, and sres == 0 would produce an Inf-weighted row.
    σ = map(s -> s > 0 ? float(s) : Inf, hip_table.sres)
    return prepare_pinv_5param(A, σ)
end

# Q_hip · b at zero perturbation: b = residuals (constant), so this is the
# whole Hipparcos five-parameter answer on the all-inactive fast path.
function _g23h_hip_x_const(Q, hip_table, include_iad, has_hip)
    (isempty(Q) || !has_hip || !include_iad) && return (0.0, 0.0, 0.0, 0.0)
    x = Q * collect(hip_table.res)
    # `fit_5param_pinv` reorders to (Δα, Δδ, μα*, μδ, ϖ).
    return (x[1], x[2], x[4], x[5])
end

# ──────────────────────────────────────────────────────────────────────
# Default variables
# ──────────────────────────────────────────────────────────────────────

function _g23h_default_variables(catalog, gaia_table, dr2_ok_mask, dr3_ok_mask,
                                 has_hip, has_rv, hip_sol, table, freeze_epochs,
                                 dr2_dup_gmag_threshold, ueva_mode)
    len_epochs = length(gaia_table.epoch)   # union pool; transit_priorities spans it
    astrometric_matched_transits_dr3 = catalog.astrometric_matched_transits_dr3

    # No `fluxratio`/`fluxratio_hip` here, deliberately. They used to be
    # emitted as derived variables defaulting to a vector of zeros, which meant
    # the *default* variable set always shadowed the bodies' own `flux_G` /
    # `flux_Hp` and silently made every companion dark. The lookup now falls
    # through to the system-level vector and then to the body fluxes (see
    # `_g23h_fluxratios`), which is only possible if nothing is declared here.
    # `:none` fixes the three σ nuisances rather than sampling them. They enter
    # the likelihood only through the UEVA channel, which `:none` switches off
    # entirely, so sampling them would add three dimensions the data cannot
    # constrain. More importantly it would make the observation *unconstructable*
    # for exactly the sources `:none` exists to serve: the G23H
    # sig_AL/sig_att_radec/sig_cal calibration is absent (NaN) for the very
    # brightest stars, where it was not extrapolated, and
    # `truncated(Normal(NaN, NaN), …)` cannot be built. The values below are
    # representative catalog medians, present only so that
    # `σ_formal = √(σ_att² + σ_AL²)` stays finite. (Ported from main's 531db2c,
    # which the v2 port had implemented only halfway — the deflation half.)
    variables = if ueva_mode === :none
        @variables begin
            σ_AL = 0.132
            σ_att = 0.0779
            σ_calib = 0.0795
        end
    else
        @variables begin
            σ_AL ~ truncated(Normal(catalog.sig_AL, catalog.sig_AL_sigma), lower=eps(), upper=10.0)
            σ_att ~ truncated(Normal(catalog.sig_att_radec, catalog.sig_att_radec_sigma), lower=eps(), upper=10.0)
            σ_calib ~ truncated(Normal(catalog.sig_cal, catalog.sig_cal_sigma), lower=eps(), upper=10.0)
        end
    end

    # Per-release selection pools. The table is already trimmed to the common
    # AGIS start (DR2 and DR3 both begin at OBMT 1192.13 rev), so the DR3-only
    # "tail" pool is exactly the DR3-usable epochs after the DR2 stop.
    dr2_stop_mjd = gaia_agis_span_dr2.stop_mjd
    epochs_mjd = vec(gaia_table.epoch)
    dr2_pool = findall(dr2_ok_mask)
    dr3_win = findall(dr3_ok_mask .& (epochs_mjd .<= dr2_stop_mjd))
    dr3_tail = findall(dr3_ok_mask .& (epochs_mjd .> dr2_stop_mjd))
    n_dr3_pool = length(dr3_win) + length(dr3_tail)
    @info "Count of missed or rejected transits:" dr3 = max(0, n_dr3_pool - astrometric_matched_transits_dr3)

    # DR3 epoch selection: sample the *indices* of the forecast epochs Gaia
    # actually used, via continuous `transit_priorities` (highest wins).
    # `astrometric_matched_transits_dr3` epochs, split between the DR2 window
    # and the DR3-only tail in proportion to the pool sizes, clamped into the
    # hard feasibility bounds.
    if n_dr3_pool < astrometric_matched_transits_dr3
        @warn "Fewer usable epochs in the forecast than `astrometric_matched_transits`; results may be inaccurate."
        n2_win = length(dr3_win)
        n_tail = length(dr3_tail)
    else
        n2_win = clamp(
            round(Int, astrometric_matched_transits_dr3 * length(dr3_win) / n_dr3_pool),
            max(astrometric_matched_transits_dr3 - length(dr3_tail), 0),
            min(length(dr3_win), astrometric_matched_transits_dr3))
        n_tail = astrometric_matched_transits_dr3 - n2_win
    end

    # DR2 epoch selection, made separately (the two sets are not nested).
    # `n_dr2_total` counts DR2's matched observations *including*
    # doubly-downlinked bright-star transits, so the DR2 set is `n_dr2_total`
    # entries drawn as `n_dr2_distinct` distinct crossings plus repeats. For
    # bright stars the distinct count is a latent marginalized between "every
    # crossing doubled" and "no duplicates"; for fainter stars duplicates are
    # empirically rare, so only geometrically-forced repeats occur.
    n_dr2_total = _g23h_dr2_target_transits(catalog, astrometric_matched_transits_dr3)
    n_dr2_hi = min(n_dr2_total, length(dr2_pool))
    length(dr2_pool) < n_dr2_total && @warn(
        "Gaia DR2 matched-transit count exceeds the geometric DR2-window pool; the excess " *
        "must be doubly-downlinked transits and is modelled as repeated epochs.",
        n_pool = length(dr2_pool), n_dr2_total)
    gmag = hasproperty(catalog, :phot_g_mean_mag_dr3) ? catalog.phot_g_mean_mag_dr3 : NaN
    bright = !ismissing(gmag) && isfinite(gmag) && gmag < dr2_dup_gmag_threshold
    n_dr2_lo = bright ? clamp(cld(n_dr2_total, 2), min(1, n_dr2_hi), n_dr2_hi) : n_dr2_hi
    marginalize_dup = n_dr2_lo < n_dr2_hi
    @info "DR2/DR3 epoch selection" n2_win n_tail n_dr2_total n_dr2_distinct_range = (n_dr2_lo, n_dr2_hi)

    # Sort: the window logic downstream needs chronological order, while
    # `partialsortperm` returns priority order. Both selections read the SAME
    # priorities vector, so within the DR2 window the smaller selection is
    # automatically the top-k subset of the larger — maximal epoch overlap, on
    # the grounds that a transit usable by DR2's pipeline was almost certainly
    # reused by DR3.
    if freeze_epochs
        transit_priorities = (randn(len_epochs)...,)
        n_dr2_distinct = rand(n_dr2_lo:max(n_dr2_lo, n_dr2_hi))
        transits = sort(vcat(
            dr3_win[partialsortperm(SVector(transit_priorities)[dr3_win], 1:n2_win, rev=true)],
            dr3_tail[partialsortperm(SVector(transit_priorities)[dr3_tail], 1:n_tail, rev=true)]))
        transits_dr2 = _g23h_select_dr2_epochs(SVector(transit_priorities), dr2_pool,
            n_dr2_distinct, n_dr2_total)
        variables = vcat(variables, @variables begin
            transit_priorities = $transit_priorities
            transits = $transits
            transits_dr2 = $transits_dr2
        end)
    else
        variables = vcat(variables, @variables begin
            transit_priorities ~ MvNormal(zeros(len_epochs), I)
            transits = sort(vcat(
                $(dr3_win)[partialsortperm(SVector(transit_priorities)[$(dr3_win)], 1:$(n2_win), rev=true)],
                $(dr3_tail)[partialsortperm(SVector(transit_priorities)[$(dr3_tail)], 1:$(n_tail), rev=true)]))
        end)
        if marginalize_dup
            # u_dup_dr2 maps uniformly onto the integer range of distinct
            # crossing counts [n_dr2_lo, n_dr2_hi].
            variables = vcat(variables, @variables begin
                u_dup_dr2 ~ Uniform(0, 1)
                transits_dr2 = _g23h_select_dr2_epochs(SVector(transit_priorities), $(dr2_pool),
                    min($(n_dr2_lo) + floor(Int, u_dup_dr2 * $(n_dr2_hi - n_dr2_lo + 1)), $(n_dr2_hi)),
                    $(n_dr2_total))
            end)
        else
            variables = vcat(variables, @variables begin
                transits_dr2 = _g23h_select_dr2_epochs(SVector(transit_priorities), $(dr2_pool),
                    $(n_dr2_hi), $(n_dr2_total))
            end)
        end
    end

    # The IAD nuisance parameters are only meaningful while the `:iad_hip`
    # channel exists (`likeobj_from_epoch_subset` strips them if a subset
    # later removes that row).
    if has_hip && :iad_hip ∈ table.kind
        variables = vcat(variables, @variables begin
            hip_iad_jitter ~ LogUniform(0.001, 100)
            iad_Δra ~ Uniform(-1000, 1000)
            iad_Δdec ~ Uniform(-1000, 1000)
            iad_Δplx ~ Uniform(-10, 10)
            iad_Δpmra ~ Uniform(-1000, 1000)
            iad_Δpmdec ~ Uniform(-1000, 1000)
            iad_pmra = $(hip_sol.pm_ra) + iad_Δpmra
            iad_pmdec = $(hip_sol.pm_de) + iad_Δpmdec
        end)
    end

    if has_rv
        # The paired GP calibration reports the per-transit RV uncertainty in
        # log space: `rv_ln_uncert_dr3` is the posterior mean of ln σ and
        # `rv_ln_uncert_err_dr3` its posterior standard deviation, so σ itself
        # is LogNormal and is sampled directly.
        variables = vcat(variables, @variables begin
            σ_rv_per_transit ~ LogNormal(catalog.rv_ln_uncert_dr3, catalog.rv_ln_uncert_err_dr3)
        end)
        n_rv = Int(catalog.rv_nb_transits)
        n_astro_sel = min(Int(astrometric_matched_transits_dr3), n_dr3_pool)
        @info "Count of RV transits:" n_rv n_astro_sel
        # RV transits are modelled as a subset of the astrometric-used
        # transits: an entirely-missed visit is missed for both. The
        # astrometric set is not a plain global top-k, so select the top-n_rv
        # priorities from *within* `transits`.
        if 0 < n_rv < n_astro_sel
            variables = vcat(variables, @variables begin
                transits_rv = sort(transits[partialsortperm(SVector(transit_priorities)[transits], 1:$n_rv, rev=true)])
            end)
        elseif n_rv >= n_astro_sel > 0
            n_rv > n_astro_sel && @warn(
                "More Gaia RV transits than modelled astrometric transits; using all astrometric transits for RV.",
                n_rv, n_astro_sel)
            variables = vcat(variables, @variables begin
                transits_rv = transits
            end)
        end
    end
    return variables
end

# ──────────────────────────────────────────────────────────────────────
# Cross-validation subsets
# ──────────────────────────────────────────────────────────────────────

function likeobj_from_epoch_subset(obs::G23HObs, obs_inds)
    # Restrict on the POST-subset table, not `obs.table`: the whole point is to
    # react to the rows the subset removed. Same call the `channels=` keyword
    # makes at construction.
    # `table[inds]`, not `table[inds, :]`: the latter returns a 2-D table, and
    # the constructor's `channels=` path produces a 1-D one. Both spellings
    # work downstream, but "the two paths cannot diverge" should be literal.
    table, catalog, priors, derived =
        _g23h_restrict(obs.table[obs_inds], obs.catalog, obs.priors, obs.derived)
    return G23HObs{typeof(table),typeof(obs.hip_table),typeof(obs.gaia_table),typeof(catalog),
        typeof(obs.hip_sol),typeof(obs.host),typeof(obs.companions),typeof(obs.ref)}(
        table, priors, derived, obs.hip_table, obs.gaia_table, catalog, obs.hip_sol,
        obs.A_prepared_5_hip, obs.A_prepared_5_dr2, obs.A_prepared_5_dr3,
        obs.pinv_5_hip, obs.hip_x_const, obs.include_iad, obs.ueva_mode, obs.frame_shift,
        obs.name,
        obs.host, obs.companions, obs.ref,
        obs.n_hip, obs.n_gaia, obs.i_ra_hip, obs.i_dec_hip,
        obs.i_ra_dr2, obs.i_dec_dr2, obs.i_ra_dr3, obs.i_dec_dr3, obs.all_epochs)
end

# ──────────────────────────────────────────────────────────────────────
# Per-draw source membership
#
# Tier 2 (per draw): the Gaia DR2/DR3 photocentre. Weights are the host (1.0)
# and each companion's effective G-band flux ratio, normalized — one
# `WeightedPoint`, one `raoff` per epoch, exact for any number of luminous
# companions.
#
# Tier 3 (per epoch): the Hipparcos abscissa. The grating response makes
# membership a function of the per-transit projected separation, so the
# weighting happens inside the scan loop and is not a photocentre at all —
# see `_hippacentre!`.
# ──────────────────────────────────────────────────────────────────────

# Effective per-companion flux ratios, zeroed wherever the companion is
# massless (an absent companion contributes no reflex, so it must contribute
# no light either). Tuple recursion, never a loop — a loop with a growing
# accumulator infers as `Tuple` and that instability propagates into
# everything downstream.
#
# Three sources, first hit wins:
#
#   1. `θ_obs.<key>` — the per-draw override. Only this tier can express
#      marginalized resolvedness, because it may read *deferred* system
#      variables, and a resolved-flag latent gating a companion's light to
#      zero for one draw is exactly that.
#   2. `θ_system.<key>` — a system-level vector of the same name. This is what
#      the old default-variables block forwarded, kept so models written
#      against it keep working.
#   3. the bodies' own `flux_<band>` variables, as F_k / F_host. The static
#      case, and the one that makes the common model non-positional.
#
# Every branch is resolved statically: `hasproperty` on a NamedTuple with a
# constant-propagated key folds away, so the hot loop sees one of the three.
@inline function _g23h_fluxratios(θ_obs, θ_system, key::Symbol, ::Val{Band}, sys,
                                  hostidx::Int, cidx::NTuple{N,Int},
                                  active::NTuple{N,Bool}, ::Type{T}) where {Band,N,T}
    f = if hasproperty(θ_obs, key)
        _g23h_asratios(getproperty(θ_obs, key), Val(N), T)
    elseif hasproperty(θ_system, key)
        _g23h_asratios(getproperty(θ_system, key), Val(N), T)
    else
        _g23h_bandratios(sys, Val(Band), hostidx, cidx, T)
    end
    return ntuple(k -> active[k] ? f[k] : zero(T), Val(N))
end

# Flux ratios read straight off the bodies: f_k = flux_<Band>(companion_k) /
# flux_<Band>(host). The ratio, not the flux, because that is what both
# instrument responses consume — the Hipparcos grating phase is defined
# against the host's own signal (`Re` starts at 1), and the Gaia photocentre
# weights are normalized anyway.
@inline function _g23h_bandratios(sys, ::Val{Band}, hostidx::Int,
                                  cidx::NTuple{N,Int}, ::Type{T}) where {Band,N,T}
    bands = keys(PlanetOrbits.fluxes(sys))
    if !(Band in bands)
        # A model with no photometry at all: every companion dark, which is
        # what v1 did when the flux-ratio vector was omitted. A model that
        # declares *some* band but not this one is a name mismatch
        # (`flux_H` for `flux_Hp`, say) and is worth failing on rather than
        # silently modelling a dark companion.
        isempty(bands) || _g23h_err_band(Band, bands)
        return ntuple(_ -> zero(T), Val(N))
    end
    fl = PlanetOrbits.fluxes(sys, Band)
    f_host = @inbounds fl[hostidx]
    # Test the primal: a differentiated zero flux is a Dual whose value is
    # zero but whose partials are not, and `iszero` on that is false.
    iszero(PlanetOrbits._primal(f_host)) && _g23h_err_darkhost(Band)
    return ntuple(k -> T((@inbounds fl[cidx[k]]) / f_host), Val(N))
end

@noinline _g23h_err_band(band, bands) = error(
    "G23HObs needs each body's flux in band :$band to form its flux ratios, but " *
    "this system's bodies declare $(bands). Either rename the body variable to " *
    "`flux_$band`, or give the observation an explicit `$(band === :Hp ? "fluxratio_hip" : "fluxratio")` vector.")
@noinline _g23h_err_darkhost(band) = error(
    "G23HObs's host body has zero flux in band :$band, so the companions' flux " *
    "ratios against it are undefined. Give the host a flux (1.0 makes every " *
    "other body's flux a contrast ratio), or pass an explicit flux-ratio vector.")

@inline function _g23h_asratios(v::Union{Tuple,AbstractVector,StaticVector}, ::Val{N}, ::Type{T}) where {N,T}
    length(v) == N || _g23h_err_ratios(length(v), N)
    return ntuple(k -> T(v[k]), Val(N))
end
# A bare scalar is unambiguous only for a single companion. Legacy accepted it
# for any count and silently applied it to every companion, which made a
# mis-shaped flux vector invisible.
@inline _g23h_asratios(v::Number, ::Val{1}, ::Type{T}) where {T} = (T(v),)
@noinline _g23h_asratios(v::Number, ::Val{N}, ::Type{T}) where {N,T} = error(
    "G23HObs received a scalar flux ratio but the observation declares $N companions. " *
    "Give one value per companion, in `companions=` order, e.g. " *
    "`fluxratio = (f_b, f_c, f_d)`.")
@noinline _g23h_err_ratios(got, want) = error(
    "G23HObs received $got flux ratios but the observation declares $want companions. " *
    "The flux-ratio vector is indexed by `companions=` order, so its length must match.")

# w[host] = 1, w[companion k] = f_k, everything else 0 — then normalized by
# `photocentre`, which is the tier-2 entry point.
@inline function _g23h_photocentre(sys, hostidx::Int, cidx::NTuple{N,Int},
                                   f::NTuple{N,T}) where {N,T}
    w = SVector(ntuple(j -> _g23h_weight(j, hostidx, cidx, f, T),
        Val(PlanetOrbits.nbodies(sys))))
    return PlanetOrbits.photocentre(w)
end
@inline _g23h_weight(j, hostidx, ::Tuple{}, ::Tuple{}, ::Type{T}) where {T} =
    j == hostidx ? one(T) : zero(T)
@inline _g23h_weight(j, hostidx, cidx::Tuple, f::Tuple, ::Type{T}) where {T} =
    j == first(cidx) ? T(first(f)) : _g23h_weight(j, hostidx, Base.tail(cidx), Base.tail(f), T)

"""
    _hippacentre!(Δα, Δδ, σ_inflation, ctx, rows, cosϕ, sinϕ,
                  host, reference, comps, f_hip, active, s)

The Hipparcos instrument response: the BINARYS "Hippacentre" along-scan
offset from the *combined* multi-companion modulated signal (Leclerc et al.
2023, A&A 672 A82, Eqs. 13 and 15), accumulated into `(Δα, Δδ)` so that the
downstream scan projection `b = Δα·cosϕ + Δδ·sinϕ` recovers Δν_B exactly.
Cross-scan components are zero, being unobservable from an abscissa.

This is **not** a photocentre, which is why it lives here rather than in
PlanetOrbits: the modulating grid makes the response a periodic function of
the projected separation. For N companions at scan-projected separations
ρ_p^(k) from the host, with Hp-band flux ratios f_k gated by a per-transit
resolution taper α_k = `α_resolve_hip(|ρ^(k)|)`, the combined phase in the
host frame is

    φ = atan2( Σ_k f_k α_k sin ζ_k ,  1 + Σ_k f_k α_k cos ζ_k ),   ζ_k = 2π ρ_p^(k)/s

and the offset from the system barycentre, projected on scan, is

    Δν_B = (s/2π)·φ + host_along

where `host_along` is the host's reflex about the reference point. That
reflex is **one** query — `raoff(sol, host, reference)` — where v8 summed
`raoff(sol_k, m_k)` over companions; the single query is exact for any
hierarchy, and the sum was only ever right for a flat astrocentric set whose
per-row gravitating masses summed to the system's.

The σ-inflation factor is the combined first-harmonic amplitude reduction
(Eq. 15, generalized to N companions):

    f_σ = (1 + Σ_k f_k α_k) / √( (1 + Σ_k f_k α_k cos ζ_k)² + (Σ_k f_k α_k sin ζ_k)² )

`σ_inflation` is multiplied in place, so pass a buffer initialised to 1. It
is the **noise model** for comparing the predicted abscissa to the observed
residuals, and must not be folded into the weighting of the LSQ that
reproduces the published catalog five-parameter solution — that fit was
performed by the Hipparcos pipeline with point-source σ.

With every companion dark, or every one resolved (α → 0), this reduces to
Δν_B = host_along and f_σ = 1: the "resolved binary, primary alone" answer.
"""
function _hippacentre!(Δα, Δδ, σ_inflation, ctx::ObsContext, rows,
                       cosϕ, sinϕ, host, reference,
                       comps::NTuple{N,PlanetOrbits.BodyRef},
                       f_hip::NTuple{N,T}, active::NTuple{N,Bool},
                       s::Float64) where {N,T}
    any(active) || return
    s_over_2π = s / (2π)
    two_π_over_s = (2π) / s
    # Squared resolution scale in mas⁻², saving a sqrt per transit per
    # companion: α(ρ) = exp(−(ρ_arcsec/RES)²) and ρ_arcsec = √(ra²+dec²)/1000.
    inv_res_mas2 = 1 / (1000 * HIPPARCOS_RESOLUTION_ARCSEC)^2
    @inbounds for i in eachindex(rows)
        sol = solutionat(ctx, rows[i])
        c = cosϕ[i]
        sn = sinϕ[i]
        # The host's reflex about the reference point — one query, all levels
        # of the hierarchy included, not gated by α (the host's barycentric
        # motion is physical whether or not Hipparcos resolved the pair).
        host_along = raoff(sol, host, reference) * c + decoff(sol, host, reference) * sn
        Re = one(T)
        Im = zero(T)
        f_total = zero(T)
        for k in 1:N
            active[k] || continue
            ck = comps[k]
            ra_p = raoff(sol, ck, host)
            dec_p = decoff(sol, ck, host)
            ρ_pk = ra_p * c + dec_p * sn
            α_k = exp(-(ra_p * ra_p + dec_p * dec_p) * inv_res_mas2)
            ζ_k = two_π_over_s * ρ_pk
            f_k = f_hip[k] * α_k
            # A degenerate orbit proposal can make raoff/decoff — and hence
            # ζ_k — non-finite, and `sincos` throws a DomainError on ±Inf/NaN,
            # which would take down the whole evaluation. The host-reflex term
            # is already non-finite for such a proposal, so propagate NaN and
            # let the sample be cleanly rejected.
            if isfinite(ζ_k)
                sin_ζ, cos_ζ = sincos(ζ_k)
            else
                sin_ζ = cos_ζ = oftype(ζ_k, NaN)
            end
            Re += f_k * cos_ζ
            Im += f_k * sin_ζ
            f_total += f_k
        end
        Δν = s_over_2π * atan(Im, Re) + host_along
        Δα[i] += Δν * c
        Δδ[i] += Δν * sn
        if σ_inflation !== nothing
            σ_inflation[i] *= (1 + f_total) / sqrt(Re * Re + Im * Im)
        end
    end
    return
end

# ──────────────────────────────────────────────────────────────────────
# Reference-point astrometry
#
# v1's `propagate_astrom`, minus the differential-light-travel term its own
# source flagged as double counting the proper motion. On an absolute frame
# the answer is the propagated frame at the catalog reference epochs, which
# `observe_pass!` already produced; otherwise it is the passthrough of the
# model's own proper motion, exactly as v1.
# ──────────────────────────────────────────────────────────────────────

@inline _g23h_frame_astrom(sra::PlanetOrbits._AbsSol, sdec::PlanetOrbits._AbsSol, pmra, pmdec) =
    (frame_ra(sra), frame_dec(sdec), frame_pmra(sra), frame_pmdec(sdec))
@inline _g23h_frame_astrom(sra, sdec, pmra, pmdec) =
    (zero(pmra), zero(pmdec), pmra, pmdec)

@inline function _g23h_astrom_at(ctx::ObsContext, i_ra::Int, i_dec::Int, pmra, pmdec)
    i_ra == 0 && return (zero(pmra), zero(pmdec), pmra, pmdec)
    return _g23h_frame_astrom(solutionat(ctx, i_ra), solutionat(ctx, i_dec), pmra, pmdec)
end

@inline _g23h_isabs(::PlanetOrbits.AbsoluteFrame) = true
@inline _g23h_isabs(::PlanetOrbits.AbstractFrame) = false

# The reference point's proper motion. The system block owns it whenever the
# model has an absolute frame; a parallax-only model declares it on this
# observation instead (v2 reserves the system-level `pmra`/`pmdec` names for
# the frame, and a partial frame is rejected at model-build time).
@inline function _g23h_pm(θ_system, θ_obs, ::Type{T}) where {T}
    if hasproperty(θ_system, :pmra)
        return T(θ_system.pmra), T(θ_system.pmdec)
    elseif hasproperty(θ_obs, :pmra)
        return T(θ_obs.pmra), T(θ_obs.pmdec)
    else
        _g23h_err_pm()
    end
end
@noinline _g23h_err_iad_pm() = error(
    "G23HObs's Hipparcos abscissa channel needs `iad_pmra` and `iad_pmdec` — the " *
    "proper motion of the Hipparcos frame, which the measured abscissae carry. " *
    "They are generated automatically unless custom `variables` are supplied; a " *
    "custom set must define them (e.g. `iad_pmra = hip_sol.pm_ra + iad_Δpmra`), or " *
    "drop the `:iad_hip` channel.")

@noinline _g23h_err_pm() = error("""
    G23HObs needs the reference point's proper motion, and this system does not define it.

    Declare the full absolute frame in the **system** block — it is all-or-nothing:

        variables = @variables begin
            plx  ~ truncated(Normal(cat.plx, cat.plx_error), lower=0)
            ra   = cat.ra;   dec  = cat.dec           # degrees, at ref_epoch
            pmra ~ Normal(cat.pmra, 10); pmdec ~ Normal(cat.pmdec, 10)   # mas/yr
            rv   = 0.0                                # km/s
            ref_epoch = 57388.0                       # MJD
        end

    (`pmra`/`pmdec` could also be defined in this observation's own `variables=`
    block, but note that supplying `variables=` **replaces** G23HObs's default block
    outright — you would lose `σ_AL`, `σ_att`, `σ_calib`, the `transit_priorities` /
    `transits` / `transits_dr2` machinery, the `iad_*` frame-offset nuisances and
    `σ_rv_per_transit`, and would have to write all of them yourself. The system
    block is almost always what you want.)
    """)

# ──────────────────────────────────────────────────────────────────────
# Epoch selection
# ──────────────────────────────────────────────────────────────────────

# Row of `epochs(obs)` for pool index `j`.
@inline _g23h_gaia_row(obs::G23HObs, j::Integer) = obs.n_hip + Int(j)

"""
Materialize this draw's epoch selections into bump-allocated index buffers:
`ii` the DR3-used pool rows, `kk` the DR2-used pool rows (repeats legal), and
`jj` the RV-used pool rows. Returns `ok=false` when a selection contains a
duplicate that it must not (`transits` and `transits_rv` index distinct
crossings; only DR2 may repeat).
"""
function _g23h_selection(obs::G23HObs, ctx::ObsContext)
    θ_obs = ctx.θ_obs
    buf = ctx.buf
    n_pool = obs.n_gaia
    seen = Bumper.alloc!(buf, Bool, n_pool)

    has_tr = hasproperty(θ_obs, :transits)
    n_sel = has_tr ? length(θ_obs.transits) : n_pool
    ii = Bumper.alloc!(buf, Int, n_sel)
    if has_tr
        tr = θ_obs.transits
        @inbounds for i in 1:n_sel
            ii[i] = Int(tr[i])
        end
        fill!(seen, false)
        @inbounds for i in 1:n_sel
            t = ii[i]
            (1 <= t <= n_pool) || return (; ok=false, ii, kk=ii, jj=ii, istart=1, iend=0)
            seen[t] && return (; ok=false, ii, kk=ii, jj=ii, istart=1, iend=0)
            seen[t] = true
        end
    else
        @inbounds for i in 1:n_sel
            ii[i] = i
        end
    end

    hasproperty(θ_obs, :transits_dr2) || _g23h_err_dr2()
    tr2 = θ_obs.transits_dr2
    n_dr2 = length(tr2)
    kk = Bumper.alloc!(buf, Int, n_dr2)
    @inbounds for i in 1:n_dr2
        kk[i] = Int(tr2[i])
    end

    if hasproperty(θ_obs, :transits_rv)
        trv = θ_obs.transits_rv
        n_rv = length(trv)
        jj = Bumper.alloc!(buf, Int, n_rv)
        @inbounds for i in 1:n_rv
            jj[i] = Int(trv[i])
        end
        fill!(seen, false)
        @inbounds for i in 1:n_rv
            t = jj[i]
            (1 <= t <= n_pool) || return (; ok=false, ii, kk, jj, istart=1, iend=0)
            seen[t] && return (; ok=false, ii, kk, jj, istart=1, iend=0)
            seen[t] = true
        end
    else
        jj = ii
    end

    # The DR3 window within the selection. The pool is already trimmed to the
    # DR3 AGIS span, so in practice this is the whole selection; kept because
    # a caller-supplied `forecast_table` need not be.
    ep = obs.gaia_table.epoch
    istart = 0
    iend = 0
    @inbounds for i in 1:n_sel
        e = ep[ii[i]]
        if istart == 0 && e >= gaia_agis_span_dr3.start_mjd
            istart = i
        end
        if e <= gaia_agis_span_dr3.stop_mjd
            iend = i
        end
    end
    istart == 0 && (istart = 1)
    iend == 0 && (iend = n_sel)
    return (; ok=true, ii, kk, jj, istart, iend)
end

@noinline _g23h_err_dr2() = error(
    "G23HObs requires a `transits_dr2` observation variable (the DR2-used epoch " *
    "selection). It is generated automatically unless custom `variables` are supplied — " *
    "custom variable sets must define it.")

# ──────────────────────────────────────────────────────────────────────
# The forward model
# ──────────────────────────────────────────────────────────────────────

"""
    _g23h_simulate!(bufs, sel, obs, ctx)

Model values for every channel: the five proper-motion pairs, the cube-root
UEVA statistic and its uncertainty, the Gaia RV sample variance, and the
Hipparcos abscissa residuals.
"""
function _g23h_simulate!(bufs, sel, obs::G23HObs, ctx::ObsContext)
    θ_system = ctx.θ_system
    θ_obs = ctx.θ_obs
    T = _system_number_type(θ_system)
    (; Δα_hip, Δδ_hip, σ_infl_hip, iad_resid, Δα_dr2, Δδ_dr2, Δα_dr3, Δδ_dr3) = bufs
    (; ii, kk, jj, istart, iend) = sel

    (; σ_att, σ_AL, σ_calib) = θ_obs
    σ_formal = sqrt(σ_att^2 + σ_AL^2)
    pmra_sys, pmdec_sys = _g23h_pm(θ_system, θ_obs, T)

    sys = ctx.system
    hostref = ref(ctx, obs.host)
    reference = ref(ctx, obs.ref)
    comps = resolverefs(ctx, obs.companions)
    cidx = map(c -> c.idx, comps)
    masses = sys.masses
    # `iszero` on a Dual is false when the value is zero but the partials are
    # not, so a differentiated zero mass would slip past this gate and land on
    # the very code paths it exists to skip. Test the primal.
    active = map(c -> !iszero(PlanetOrbits._primal(masses[c.idx])), comps)
    any_active = any(active)

    has_hip = !isnothing(obs.catalog.dist_hip)
    has_iad = :iad_hip ∈ obs.table.kind
    has_rv = :rv_dr3 ∈ obs.table.kind

    # ---- Hipparcos ---------------------------------------------------------
    hip_bias_pm_sq = zero(T)
    Δα_h = Δδ_h = Δpmra_h = Δpmdec_h = zero(T)
    if has_hip
        if any_active
            f_hip = _g23h_fluxratios(θ_obs, θ_system, :fluxratio_hip, Val(:Hp), sys,
                hostref.idx, cidx, active, T)
            _hippacentre!(Δα_hip, Δδ_hip, σ_infl_hip, ctx, 1:obs.n_hip,
                obs.hip_table.cosϕ, obs.hip_table.sinϕ, hostref, reference,
                comps, f_hip, active, HIPPARCOS_GRID_STEP_ARCSEC)
            # Extract the catalog five-parameter bias with the *uninflated* σ:
            # the pipeline that produced the catalog used point-source σ, so
            # to reproduce the bias it absorbed we must weight the LSQ the
            # same way. σ_infl_hip feeds only the abscissa residual noise and
            # the catalog-point covariance.
            resid = obs.include_iad ? obs.hip_table.res : 0.0
            out = fit_5param_pinv(obs.pinv_5_hip, obs.hip_table.cosϕ, obs.hip_table.sinϕ,
                Δα_hip, Δδ_hip, resid; buf=ctx.buf)
            Δα_h, Δδ_h, Δpmra_h, Δpmdec_h = out.parameters
        else
            # All companions inactive: every perturbation is zero and the LSQ
            # collapses to the cached projection of the constant residuals.
            Δα_h = T(obs.hip_x_const[1])
            Δδ_h = T(obs.hip_x_const[2])
            Δpmra_h = T(obs.hip_x_const[3])
            Δpmdec_h = T(obs.hip_x_const[4])
        end
        hip_bias_pm_sq = Δpmra_h^2 + Δpmdec_h^2
    end

    α_h₀, δ_h₀, pmra_h₀, pmdec_h₀ =
        _g23h_astrom_at(ctx, obs.i_ra_hip, obs.i_dec_hip, pmra_sys, pmdec_sys)
    μ_h = has_hip ? @SVector([pmra_h₀ + Δpmra_h, pmdec_h₀ + Δpmdec_h]) :
          @SVector([zero(T), zero(T)])

    # ---- Hipparcos abscissa residuals --------------------------------------
    if has_hip && has_iad
        # The five-parameter Hipparcos frame offset — the general facility of
        # `skypath.jl`, not something Hipparcos-specific. The parallax anchors
        # on the *Hipparcos catalog* value rather than the system's: this
        # channel is here for the companion curvature, and the frame is pure
        # nuisance. (`HipparcosIADObs`, where the abscissae are the only data,
        # anchors on the system `plx` instead so that they constrain it.)
        # `frame_offset` defaults a missing component to zero, which is right
        # for a position or parallax offset and catastrophic for a proper
        # motion: the measured abscissa carries the catalog's full sky path,
        # so a zero `iad_pmra` silently models a star that does not move.
        # Keep the loud failure the destructuring form used to give.
        (hasproperty(θ_obs, :iad_pmra) && hasproperty(θ_obs, :iad_pmdec)) ||
            _g23h_err_iad_pm()
        off0 = frame_offset(θ_obs, obs.hip_sol.plx, T)
        # The published five-parameter solution already absorbed the
        # companion's apparent position and proper motion, so subtract that
        # bias back out before comparing against the residuals it left.
        off = FrameOffset(off0.Δra - Δα_h, off0.Δdec - Δδ_h, off0.plx,
            off0.pmra - Δpmra_h, off0.pmdec - Δpmdec_h)
        inv_jy = inv(julian_year)
        ep = obs.hip_table.epoch
        cϕ = obs.hip_table.cosϕ
        sϕ = obs.hip_table.sinϕ
        pf = obs.hip_table.parallaxFactorAlongScan
        meas = obs.hip_table.proj_meas_alongscan
        # The residual is the perpendicular distance from the model point to
        # the measured abscissa line; because the line's direction is
        # (−sinϕ, cosϕ) that distance reduces exactly to the scalar
        # scan projection, which is what `proj_meas_alongscan` precomputes.
        #
        # Stored *signed* (`measured − model`). The only consumer,
        # `_g23h_ln_like`, squares it, so the likelihood is bit-identical to
        # the `abs` this used to hold; the sign is what a residual plot needs
        # (an unsigned residual has a half-normal marginal, which reads as a
        # bias rather than as noise). `simulate` still reports the magnitude
        # under the `iad_resid` name and the signed vector as
        # `iad_resid_signed`, so nothing outside this file changes.
        @inbounds @simd for i in eachindex(ep, iad_resid)
            Δt = (ep[i] - hipparcos_catalog_epoch_mjd) * inv_jy
            proj = frame_offset_alongscan(off, Δt, cϕ[i], sϕ[i], pf[i],
                Δα_hip[i], Δδ_hip[i])
            iad_resid[i] = meas[i] - proj
        end
    end

    # ---- Gaia DR2 and DR3 --------------------------------------------------
    #
    # Tier 2: one `WeightedPoint` per draw over the host and its companions,
    # weighted by their effective G-band flux ratios, and one `raoff` per
    # transit. For several luminous companions this is the exact flux-weighted
    # mean of apparent positions, which is *not* the superposition of
    # per-companion photocentres v1 computed.
    if any_active
        f_G = _g23h_fluxratios(θ_obs, θ_system, :fluxratio, Val(:G), sys,
            hostref.idx, cidx, active, T)
        photo = _g23h_photocentre(sys, hostref.idx, cidx, f_G)
        @inbounds for i in istart:iend
            sol = solutionat(ctx, _g23h_gaia_row(obs, ii[i]))
            Δα_dr3[i-istart+1] = raoff(sol, photo, reference)
            Δδ_dr3[i-istart+1] = decoff(sol, photo, reference)
        end
        @inbounds for i in eachindex(kk)
            sol = solutionat(ctx, _g23h_gaia_row(obs, kk[i]))
            Δα_dr2[i] = raoff(sol, photo, reference)
            Δδ_dr2[i] = decoff(sol, photo, reference)
        end
        A_dr3 = @view obs.A_prepared_5_dr3[view(ii, istart:iend), :]
        cϕ_dr3 = @view obs.gaia_table.cosϕ[view(ii, istart:iend)]
        sϕ_dr3 = @view obs.gaia_table.sinϕ[view(ii, istart:iend)]
        out_dr3 = fit_5param_prepared(A_dr3, cϕ_dr3, sϕ_dr3, Δα_dr3, Δδ_dr3, 0.0, σ_formal;
            include_chi2=Val(true), buf=ctx.buf)
        Δα_dr3_p, Δδ_dr3_p, Δpmra_dr3, Δpmdec_dr3 = out_dr3.parameters
        chi2_astro = out_dr3.chi_squared_astro

        A_dr2 = @view obs.A_prepared_5_dr2[kk, :]
        cϕ_dr2 = @view obs.gaia_table.cosϕ[kk]
        sϕ_dr2 = @view obs.gaia_table.sinϕ[kk]
        out_dr2 = fit_5param_prepared(A_dr2, cϕ_dr2, sϕ_dr2, Δα_dr2, Δδ_dr2; buf=ctx.buf)
        Δα_dr2_p, Δδ_dr2_p, Δpmra_dr2, Δpmdec_dr2 = out_dr2.parameters
    else
        # Every per-transit perturbation is exactly zero, and a five-parameter
        # fit to zero data returns zero parameters and zero χ². Skipping both
        # solves is what makes the all-inactive draws cheap — and in a model
        # that marginalizes over the companion count they are roughly half of
        # them.
        z = zero(T)
        Δα_dr3_p = Δδ_dr3_p = Δpmra_dr3 = Δpmdec_dr3 = z
        Δα_dr2_p = Δδ_dr2_p = Δpmra_dr2 = Δpmdec_dr2 = z
        chi2_astro = z
    end
    α_dr3₀, δ_dr3₀, pmra_dr3₀, pmdec_dr3₀ =
        _g23h_astrom_at(ctx, obs.i_ra_dr3, obs.i_dec_dr3, pmra_sys, pmdec_sys)
    μ_dr3 = @SVector [pmra_dr3₀ + Δpmra_dr3, pmdec_dr3₀ + Δpmdec_dr3]
    α_dr2₀, δ_dr2₀, pmra_dr2₀, pmdec_dr2₀ =
        _g23h_astrom_at(ctx, obs.i_ra_dr2, obs.i_dec_dr2, pmra_sys, pmdec_sys)
    μ_dr2 = @SVector [pmra_dr2₀ + Δpmra_dr2, pmdec_dr2₀ + Δpmdec_dr2]

    # ---- H-G and DR3-DR2 scaled position differences -----------------------
    Δt_hg_ra = obs.catalog.epoch_ra_dr3_mjd - obs.catalog.epoch_ra_hip_mjd
    Δt_hg_dec = obs.catalog.epoch_dec_dr3_mjd - obs.catalog.epoch_dec_hip_mjd
    Δt_32_ra = obs.catalog.epoch_ra_dr3_mjd - obs.catalog.epoch_ra_dr2_mjd
    Δt_32_dec = obs.catalog.epoch_dec_dr3_mjd - obs.catalog.epoch_dec_dr2_mjd
    if _g23h_isabs(sys.frame)
        if has_hip
            Δα_hg_prop = (α_dr3₀ - α_h₀) * 60 * 60 * 1000 * cosd((δ_dr3₀ + δ_h₀) / 2)
            Δδ_hg_prop = (δ_dr3₀ - δ_h₀) * 60 * 60 * 1000
            pmra_hg = (Δα_dr3_p - Δα_h + Δα_hg_prop) / Δt_hg_ra * julian_year
            pmdec_hg = (Δδ_dr3_p - Δδ_h + Δδ_hg_prop) / Δt_hg_dec * julian_year
        else
            pmra_hg = zero(T)
            pmdec_hg = zero(T)
        end
        Δα_32_prop = (α_dr3₀ - α_dr2₀) * 60 * 60 * 1000 * cosd((δ_dr3₀ + δ_dr2₀) / 2)
        Δδ_32_prop = (δ_dr3₀ - δ_dr2₀) * 60 * 60 * 1000
        pmra_dr32 = (Δα_dr3_p - Δα_dr2_p + Δα_32_prop) / Δt_32_ra * julian_year
        pmdec_dr32 = (Δδ_dr3_p - Δδ_dr2_p + Δδ_32_prop) / Δt_32_dec * julian_year
    else
        if has_hip
            pmra_hg = (Δα_dr3_p - Δα_h) / Δt_hg_ra * julian_year + pmra_sys
            pmdec_hg = (Δδ_dr3_p - Δδ_h) / Δt_hg_dec * julian_year + pmdec_sys
        else
            pmra_hg = zero(T)
            pmdec_hg = zero(T)
        end
        pmra_dr32 = (Δα_dr3_p - Δα_dr2_p) / Δt_32_ra * julian_year + pmra_sys
        pmdec_dr32 = (Δδ_dr3_p - Δδ_dr2_p) / Δt_32_dec * julian_year + pmdec_sys
    end
    μ_hg = @SVector [pmra_hg, pmdec_hg]
    μ_dr32 = @SVector [pmra_dr32, pmdec_dr32]

    # ---- UEVA and the DR3 uncertainty deflation ----------------------------
    (; astrometric_chi2_al_dr3, astrometric_n_good_obs_al_dr3,
       astrometric_matched_transits_dr3, astrometric_excess_noise_dr3, ruwe_dr3) = obs.catalog
    gaia_n_dof = obs.catalog.astrometric_params_solved_dr3 == 31 ? 5 : 6
    N = astrometric_n_good_obs_al_dr3
    N_FoV = astrometric_matched_transits_dr3
    N_AL = N / N_FoV
    # Expected UEVA for a single star, and its standard deviation
    # (Eqs. D.8 and D.9).
    μ_UEVA_single = (N_AL / (N - gaia_n_dof)) *
                    ((N_FoV - gaia_n_dof) * σ_calib^2 + N_FoV * σ_AL^2)
    σ_UEVA_single = sqrt(2 * N_AL / (N - gaia_n_dof)^2 * (
        N_AL * (N_FoV - gaia_n_dof) * σ_calib^4 + N_FoV * σ_AL^4 +
        2 * N_FoV * σ_AL^2 * σ_calib^2))
    if obs.ueva_mode === :none
        μ_1_3 = zero(T)
        UEVA_unc = zero(T)
        UEVA_model = zero(T)
        deflation_factor_dr3 = one(T)
    else
        UEVA_Gaia = if obs.ueva_mode === :EAN
            astrometric_excess_noise_dr3^2 + σ_att^2 + σ_AL^2
        else
            u0 = 1 / ruwe_dr3 * sqrt(astrometric_chi2_al_dr3 / (N - gaia_n_dof))
            (ruwe_dr3 * u0)^2 * σ_formal^2
        end
        μ_1_3 = UEVA_Gaia^(1 / 3)
        # /3 from the cube-root transformation.
        UEVA_unc = σ_UEVA_single * μ_UEVA_single^(-2 / 3) / 3
        # chi_squared_astro sums over the transits actually modelled (the DR3
        # window slice) while the normalizations below assume N_FoV of them.
        # The counts agree except in the forecast-shortfall case; rescale so
        # the predicted excess per companion stays consistent with the N_FoV
        # normalization.
        n_dr3_modeled = iend - istart + 1
        chi2_scaled = chi2_astro * N_AL * (N_FoV / n_dr3_modeled)
        UEVA_model_raw = (chi2_scaled * σ_formal^2) / (N - gaia_n_dof)
        # Cube-root transform (Eq. 27, Sect. 5.1.1).
        UEVA_model = cbrt((chi2_scaled * σ_formal^2) / (N_AL * N_FoV - gaia_n_dof) + μ_UEVA_single)
        d_raw = sqrt(μ_UEVA_single / UEVA_Gaia)
        deflation_factor_dr3 = d_raw > 1 ? one(d_raw) : d_raw
    end

    # ---- Gaia RV variability -----------------------------------------------
    if has_rv
        n_rv_ep = length(jj)
        rv_model = Bumper.alloc!(ctx.buf, T, n_rv_ep)
        fill!(rv_model, zero(T))
        if any_active
            @inbounds for i in 1:n_rv_ep
                sol = solutionat(ctx, _g23h_gaia_row(obs, jj[i]))
                rv_model[i] = radvel(sol, hostref, reference) / 1e3   # km/s
            end
        end
        rv_sum = zero(T)
        @inbounds for i in 1:n_rv_ep
            rv_sum += rv_model[i]
        end
        rv_mean = rv_sum / n_rv_ep
        rv_sumsq = zero(T)
        @inbounds for i in 1:n_rv_ep
            d = rv_model[i] - rv_mean
            rv_sumsq += d * d
        end
        sample_variance = rv_sumsq / (n_rv_ep - 1)
        N_rv = obs.catalog.rv_nb_transits
        rv_dof = T(N_rv - 1)
        ε_catalog = obs.catalog.radial_velocity_error
        # Convert the catalog error back to a sample variance (Eq. A4).
        s_catalog_squared = T((2 * N_rv / π) * (ε_catalog^2 - 0.113^2))
    else
        rv_dof = convert(T, NaN)
        rv_mean = convert(T, NaN)
        sample_variance = convert(T, NaN)
        s_catalog_squared = convert(T, NaN)
    end

    # Shift the whole reference frame so the model proper motion refers to the
    # primary rather than the barycentre. This vastly improves sampling
    # efficiency, and is why Δpmra_dr3 is subtracted from *every* channel
    # including DR3's own.
    #
    # It is also a redefinition of a *system*-level parameter made inside one
    # observation, which is why it is switchable and why more than one
    # `G23HObs` on a shared frame with it on is an error — see the keyword's
    # docstring, and `AnchoredFrame` for the parameterization that gets the
    # same conditioning without the redefinition.
    shift = obs.frame_shift ? SVector{2,T}(Δpmra_dr3, Δpmdec_dr3) : zero(SVector{2,T})
    return (;
        UEVA_model, UEVA_unc, μ_1_3,
        μ_h=μ_h .- shift,
        μ_hg=μ_hg .- shift,
        μ_dr2=μ_dr2 .- shift,
        μ_dr32=μ_dr32 .- shift,
        μ_dr3=μ_dr3 .- shift,
        hip_bias_pm_sq,
        n_dr3=iend - istart + 1,
        n_dr2=length(kk),
        rv_dof, rv_mean, sample_variance, s_catalog_squared,
        deflation_factor_dr3,
        Δα_dr3=Δα_dr3_p, Δδ_dr3=Δδ_dr3_p, Δpmra_dr3, Δpmdec_dr3,
    )
end

"""
    simulate(obs::G23HObs, ctx) -> NamedTuple

Model values for every channel of this observation at the sample in `ctx`,
plus the per-transit Hipparcos abscissa residuals and σ inflation. Allocates
ordinary arrays; `ln_like` runs the same code against bump-allocated ones.

`iad_resid` is the unsigned perpendicular distance from the model to each
measured abscissa line; `iad_resid_signed` is the same residual as
`measured − model`.
"""
function simulate(obs::G23HObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    @no_escape ctx.buf begin
        sel = _g23h_selection(obs, ctx)
        sel.ok || error("G23HObs: the sampled transit selection contains duplicate indices.")
        n_hip = obs.n_hip
        n_dr3 = sel.iend - sel.istart + 1
        bufs = (; Δα_hip=zeros(T, n_hip), Δδ_hip=zeros(T, n_hip),
            σ_infl_hip=ones(T, n_hip), iad_resid=zeros(T, n_hip),
            Δα_dr2=zeros(T, length(sel.kk)), Δδ_dr2=zeros(T, length(sel.kk)),
            Δα_dr3=zeros(T, n_dr3), Δδ_dr3=zeros(T, n_dr3))
        sim = _g23h_simulate!(bufs, sel, obs, ctx)
        out = merge(sim, (; μ=@SVector([sim.μ_h[1], sim.μ_h[2], sim.μ_hg[1], sim.μ_hg[2],
                sim.μ_dr2[1], sim.μ_dr2[2], sim.μ_dr32[1], sim.μ_dr32[2],
                sim.μ_dr3[1], sim.μ_dr3[2], sim.UEVA_model, sim.sample_variance]),
            # `iad_resid` keeps its published meaning (the perpendicular
            # distance, unsigned); `iad_resid_signed` is the same quantity with
            # the sign the plotting protocol needs.
            iad_resid=abs.(bufs.iad_resid), iad_resid_signed=bufs.iad_resid,
            σ_inflation_hip=bufs.σ_infl_hip,
            Δα_mas_hip=bufs.Δα_hip, Δδ_mas_hip=bufs.Δδ_hip))
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────
# Correction flags
# ──────────────────────────────────────────────────────────────────────

# G23HObs deliberately does NOT declare `reduced_lighttime_free`, though its
# channels are built from catalog quantities that were reduced light-time-free.
#
# The trait says "this observation IS a reduction in that convention, so
# propagate it in that convention". That reading is right for a single
# reduction window: within one, the light-time-free model is the
# parameterization the solution was expressed in. It is the wrong reading for
# G23H's dominant channels, which are INTER-window — the Hipparcos–Gaia
# position difference and DR3−DR2 connect two solutions decades apart, and no
# reduction convention governs the path between them. Hipparcos reports the
# apparent direction at J1991.25 and Gaia at J2016.0; both are statements about
# the real apparent path, locally linearized. Connecting them is a propagation
# question, and the rigorous apparent path is the right answer to it.
#
# The pin was added (2026-08-14) when it WAS load-bearing: PlanetOrbits then
# double-counted the catalog Doppler factor, so `barycentric_lighttime=true`
# contradicted the data by μ·|v_r|/c = 3.8 mas/yr at Barnard and manufactured a
# 0.74 M_jup companion against a true ~0.13. That bug is fixed — `_dedoppler`
# now propagates the true worldline and reads out apparent quantities — so the
# flag no longer separates "true vs apparent" but "rigorous apparent vs the
# catalogs' linear-in-3D approximation to it". Both are self-consistent, the
# difference is a genuine second-order term, and `:auto` can therefore measure
# it like any other correction rather than being told the answer.
#
# Measured, for scale: on Barnard's FGS product the two settings differ by
# 0.082 mas in predicted position over a 23-year lever arm (the FGS reference
# track is the pipeline's light-time-free propagation, so it isolates exactly
# this term). That is ~8% of the FGS−Gaia proper-motion-anomaly σ — small, but
# a systematic, and worth measuring rather than assuming.
#
# `reduced_lighttime_free` itself stays. It is still the right declaration for
# an observation that is purely a single reduction; see its docstring.

# What a correction does to this observation: the five proper-motion channel
# pairs of the coupled catalog block, judged against the tightest catalog σ
# among the channels present. UEVA and the RV-variance channel are deliberately
# excluded — different units, and both derive from the same per-transit
# Δα/Δδ the PM channels already probe.
has_correction_impact(::Type{<:G23HObs}) = true
function correction_impact(obs::G23HObs, a::ObsContext, b::ObsContext)
    cat = obs.catalog
    kinds = obs.table.kind
    μa = simulate(obs, a).μ
    μb = simulate(obs, b).μ
    m = 0.0
    n = 0
    σs = Float64[]
    for (bi, blk) in enumerate(_g23h_blocks)
        blk.name === :ueva && continue
        all(k -> k ∈ kinds, blk.kinds) || continue
        off = _g23h_block_offset[bi]
        for j in 1:2
            d = abs(float(μa[off+j] - μb[off+j]))
            isfinite(d) || return (; delta=NaN, sigma=_tightest(σs), n=0)
            m = max(m, d)
            n += 1
        end
        push!(σs, float(getproperty(cat, Symbol(:pmra_, blk.name, :_error))),
              float(getproperty(cat, Symbol(:pmdec_, blk.name, :_error))))
    end
    return (; delta=m, sigma=_tightest(σs), n)
end

# ──────────────────────────────────────────────────────────────────────
# Likelihood
# ──────────────────────────────────────────────────────────────────────

function ln_like(obs::G23HObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    ll = zero(T)
    @no_escape ctx.buf begin
        sel = _g23h_selection(obs, ctx)
        if !sel.ok
            ll = convert(T, -Inf)
        else
            n_hip = obs.n_hip
            n_dr3 = sel.iend - sel.istart + 1
            n_dr2 = length(sel.kk)
            Δα_hip = @alloc(T, n_hip); fill!(Δα_hip, zero(T))
            Δδ_hip = @alloc(T, n_hip); fill!(Δδ_hip, zero(T))
            iad_resid = @alloc(T, n_hip); fill!(iad_resid, zero(T))
            # σ-inflation buffer for the BINARYS first-harmonic correction.
            # Initialised to 1 (no inflation); populated per transit from the
            # *combined* multi-companion modulated signal. Consumed in two
            # places: the abscissa residual likelihood, and the catalog
            # covariance Σ_h below. Never in the catalog-bias LSQ.
            σ_infl_hip = @alloc(T, n_hip); fill!(σ_infl_hip, one(T))
            Δα_dr2 = @alloc(T, n_dr2); fill!(Δα_dr2, zero(T))
            Δδ_dr2 = @alloc(T, n_dr2); fill!(Δδ_dr2, zero(T))
            Δα_dr3 = @alloc(T, n_dr3); fill!(Δα_dr3, zero(T))
            Δδ_dr3 = @alloc(T, n_dr3); fill!(Δδ_dr3, zero(T))
            bufs = (; Δα_hip, Δδ_hip, σ_infl_hip, iad_resid, Δα_dr2, Δδ_dr2, Δα_dr3, Δδ_dr3)
            sim = _g23h_simulate!(bufs, sel, obs, ctx)
            ll = _g23h_ln_like(obs, ctx, sim, iad_resid, σ_infl_hip, T)
        end
    end
    return isnan(ll) ? convert(T, -Inf) : ll
end

function _g23h_ln_like(obs::G23HObs, ctx::ObsContext, sim, iad_resid, σ_infl_hip, ::Type{T0}) where {T0}
    θ_obs = ctx.θ_obs
    ll = zero(T0)
    cat = obs.catalog
    kinds = obs.table.kind

    # Optional per-channel proper-motion jitter, added in quadrature.
    dist_hip = _g23h_jittered(cat.dist_hip, θ_obs, :σ_hip_pmra, :σ_hip_pmdec,
        cat.pmra_hip_error, cat.pmdec_hip_error, cat.pmra_pmdec_hip, T0)
    dist_hg = _g23h_jittered(cat.dist_hg, θ_obs, :σ_hg_pmra, :σ_hg_pmdec,
        cat.pmra_hg_error, cat.pmdec_hg_error, cat.pmra_pmdec_hg, T0)
    dist_dr2 = _g23h_jittered(cat.dist_dr2, θ_obs, :σ_dr2_pmra, :σ_dr2_pmdec,
        cat.pmra_dr2_error, cat.pmdec_dr2_error, cat.pmra_pmdec_dr2, T0)
    dist_dr32 = _g23h_jittered(cat.dist_dr32, θ_obs, :σ_dr32_pmra, :σ_dr32_pmdec,
        cat.pmra_dr32_error, cat.pmdec_dr32_error, cat.pmra_pmdec_dr32, T0)
    dist_dr3 = _g23h_jittered(cat.dist_dr3, θ_obs, :σ_dr3_pmra, :σ_dr3_pmdec,
        cat.pmra_dr3_error, cat.pmdec_dr3_error, cat.pmra_pmdec_dr3, T0)

    absolute = _g23h_isabs(ctx.system.frame)
    Δ_hg, Δ_h = _g23h_nonlinear(cat, absolute, !isnothing(dist_hip))
    μ_h = sim.μ_h + Δ_h
    μ_hg = sim.μ_hg + Δ_hg

    mask = ntuple(i -> _g23h_channel_kinds[i] ∈ kinds, Val(_g23h_nchan))
    n_components = count(mask)

    # Change-of-variables Jacobian for the UEVA component. The MvNormal below
    # treats t = UEVA^(1/3) as the datum, but t is a parameter-dependent
    # transform of the raw catalog datum:
    #   :EAN  t = (EAN² + σ_att² + σ_AL²)^(1/3),          datum EAN
    #   :RUWE t = (χ²_AL·(σ_att² + σ_AL²)/(N−dof))^(1/3), datum χ²_AL
    # A proper density in the raw datum needs log|dt/d(datum)|; the
    # parameter-dependent part is what appears here. (Same class of bug as the
    # rv_dr3 ξ² term below; caught by SBC 2026-07-04.)
    if :ueva_dr3 ∈ kinds
        if obs.ueva_mode === :EAN
            # EAN == 0 is a boundary atom of the Gaia reduction (the
            # excess-noise fit pinned at zero); the continuous
            # change-of-variables does not apply there.
            if cat.astrometric_excess_noise_dr3 > 0
                ll += -2 * log(sim.μ_1_3)
            end
        elseif obs.ueva_mode === :RUWE
            ll += (T0(1) / 3) * log(θ_obs.σ_att^2 + θ_obs.σ_AL^2)
        end
    end

    # The Hipparcos abscissa channel is uncorrelated with the catalog block,
    # so it is handled separately rather than factorizing a much bigger matrix.
    if :iad_hip ∈ kinds
        (; hip_iad_jitter) = θ_obs
        half_log2π = T0(0.5) * log(T0(2π))
        jitter² = hip_iad_jitter * hip_iad_jitter
        rej = obs.hip_table.reject
        sres = obs.hip_table.sres_renorm
        @inbounds for i in eachindex(iad_resid)
            rej[i] && continue
            # Inflate the per-transit residual σ by the BINARYS first-harmonic
            # factor: it is 1 in the unresolved limit and grows where the
            # binary modulation reduces the signal amplitude, and the residual
            # noise scales the same way.
            s = sres[i] * σ_infl_hip[i]
            σ² = s * s + jitter²
            r = iad_resid[i]
            ll += -T0(0.5) * (r * r / σ² + log(σ²)) - half_log2π
        end
    end

    # Gaia RV variability (Chance et al. 2022).
    if :rv_dr3 ∈ kinds
        ε_catalog = cat.radial_velocity_error
        N_rv = cat.rv_nb_transits
        σ_rv = θ_obs.σ_rv_per_transit                      # km/s
        s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)   # Eq. A4
        # Non-centrality parameter (Eq. C2). Since sample_variance is computed
        # from the modelled RV at N_rv epochs, (N_rv − 1)·s² = Σ(μ_n − μ̄)²
        # reproduces λ = Σ((μ_n − μ̄)/σ)² exactly.
        ncp = (N_rv - 1) * sim.sample_variance / σ_rv^2
        # The paper's Eq. C2 states N_k degrees of freedom while its null (Eq.
        # A6) uses N_k − 1. The standard sampling-theory result for
        # Σ(v_n − v̄)²/σ² is non-central χ²_{N−1}(λ), which reduces to the
        # central χ²_{N−1} under the null; we use N−1 to stay self-consistent.
        try
            d = NoncentralChisq(N_rv - 1, ncp)
            ξ² = (N_rv - 1) * s_catalog_squared / σ_rv^2   # Eq. A5
            # Change of variables: the raw datum is ε, but the density is
            # evaluated in ξ², a transform that depends on σ_rv. The
            # parameter-dependent part of log|dξ²/dε| is −2·log σ; without it
            # the posterior is biased high by 2·sd(ln σ)² in log space.
            ll += logpdf(d, ξ²) - 2 * log(σ_rv)
        catch err
            # Inherited from v1, where this branch was silent. Two very
            # different things land in it and only one of them is a
            # rejection:
            #   * a genuinely out-of-domain ncp or ξ² — a bad proposal, which
            #     −Inf correctly rejects;
            #   * `Distributions.NoncentralChisq`'s log-density not accepting
            #     ForwardDiff `Dual`s at all, which makes the *whole model*
            #     −Inf under any gradient-based sampler while the primal path
            #     stays finite. That is a limitation, not a rejection, and it
            #     is invisible unless something says so.
            # The numbers are unchanged; the second case is now diagnosable.
            err isa MethodError && @warn(
                "The Gaia RV variability channel (:rv_dr3) could not be evaluated: " *
                "`NoncentralChisq`'s log-density does not support this number type. " *
                "Under automatic differentiation this makes the entire log-density " *
                "-Inf. Use a derivative-free sampler, or drop the :rv_dr3 row with " *
                "`likeobj_from_epoch_subset` / `include_rv=false`.",
                exception = err, maxlog = 1)
            ll += -Inf
        end
    end

    # ---- the coupled catalog block -----------------------------------------
    μ_dr3_cat, Σ_dr3 = params(dist_dr3)
    μ_dr2_cat, Σ_dr2 = params(dist_dr2)
    μ_dr32_cat, Σ_dr32 = params(dist_dr32)
    μ_h_cat, Σ_h = isnothing(dist_hip) ? (@SVector[0.0, 0.0], @SMatrix zeros(2, 2)) : params(dist_hip)
    μ_hg_cat, Σ_hg = isnothing(dist_hg) ? (@SVector[0.0, 0.0], @SMatrix zeros(2, 2)) : params(dist_hg)

    # DR3 (and DR2) position covariance at the central epoch, already inflated
    # by Gaia. The UEVA channel says how much of that inflation our companion
    # model explains, and the DR3 block is deflated by it — which also feeds
    # back into the DR3−DR2 difference through the cross term.
    σ_ra3, σ_dec3, ρ3 = cat.ra_error_central_dr3, cat.dec_error_central_dr3, cat.ra_dec_corr_central_dr3
    σ_ra2, σ_dec2, ρ2 = cat.ra_error_central_dr2, cat.dec_error_central_dr2, cat.ra_dec_corr_central_dr2
    Σ_pos_dr3 = @SMatrix [σ_ra3^2 ρ3*σ_ra3*σ_dec3; ρ3*σ_ra3*σ_dec3 σ_dec3^2]
    ρ_dr3_dr2 = cat.rho_dr2_dr3
    Σ_cross = @SMatrix [
        ρ_dr3_dr2*σ_ra3*σ_ra2 ρ_dr3_dr2*ρ3*σ_ra3*σ_dec2
        ρ_dr3_dr2*ρ2*σ_dec3*σ_ra2 ρ_dr3_dr2*σ_dec3*σ_dec2
    ]
    Δt_ra = (cat.epoch_ra_dr3_mjd - cat.epoch_ra_dr2_mjd) / julian_year
    Δt_dec = (cat.epoch_dec_dr3_mjd - cat.epoch_dec_dr2_mjd) / julian_year
    d = sim.deflation_factor_dr3
    ΔΣ_pos = (d^2 - 1) * Σ_pos_dr3 - (d - 1) * (Σ_cross + Σ_cross')
    Tr = @SMatrix [1/Δt_ra 0.0; 0.0 1/Δt_dec]
    ΔΣ_dr32 = Tr * ΔΣ_pos * Tr'

    T = promote_type(eltype(Σ_h), eltype(Σ_hg), eltype(Σ_dr2), eltype(Σ_dr32),
        eltype(ΔΣ_dr32), eltype(Σ_dr3), typeof(d), T0)
    Σ_h = SMatrix{2,2,T,4}(Σ_h)
    Σ_hg = SMatrix{2,2,T,4}(Σ_hg)
    Σ_dr2 = SMatrix{2,2,T,4}(Σ_dr2)

    # BINARYS f_σ inflation of the Hipparcos catalog covariance. The catalog's
    # own five-parameter fit used point-source σ, so the LSQ reproducing its
    # bias must too — but the *uncertainty* on the catalog point should
    # reflect the binary-induced increase in per-transit noise. Multiply Σ_h
    # by the transit-averaged f_σ². In the dark-companion / single-star /
    # wide-resolved limit f_σ → 1 and this is a no-op. Σ_hg is left
    # uninflated: its long baseline has a Gaia endpoint as well, so the right
    # multiplier would be < f_σ², and HGCA's renormalisation self-calibrates.
    if !isnothing(dist_hip) && obs.n_hip > 0
        n_used = 0
        sumsq = zero(T)
        rej = obs.hip_table.reject
        @inbounds for i in eachindex(σ_infl_hip)
            rej[i] && continue
            n_used += 1
            sumsq += σ_infl_hip[i]^2
        end
        n_used > 0 && (Σ_h = Σ_h * (sumsq / n_used))
    end

    # BINARYS epistemic uncertainty on the catalog-bias correction. The
    # model's predicted bias is the photocentre modulation absorbed by the
    # published Hipparcos five-parameter fit, and it is correct only to the
    # extent that this likelihood matches what the Hipparcos pipeline did.
    # Known approximations: the H1+H2 composite catalog point modelled with a
    # single-reduction basis matrix; Hp-band mass–luminosity systematics in
    # the per-companion flux ratios; a resolution gate anchored to the grating
    # step rather than an empirical detection efficiency. Absorb the residual
    # model error as ε²·|Δpm_h|² added isotropically — zero at the
    # dark-companion limit, growing where the bias does.
    ε_binarys = T(0.3)
    if !isnothing(dist_hip) && sim.hip_bias_pm_sq > zero(T)
        Σ_h = Σ_h + (ε_binarys^2 * sim.hip_bias_pm_sq) * SMatrix{2,2,T,4}(I)
    end

    Σ_dr32 = SMatrix{2,2,T,4}(Σ_dr32 .+ ΔΣ_dr32)
    Σ_dr3 = SMatrix{2,2,T,4}(Σ_dr3 .* d^2)

    # The one place the sub-block list is paired with values, in its order.
    # Widths and offsets come from the list itself, so a new channel is added
    # here and in `_g23h_blocks` and nowhere else.
    blockvals = (
        (μc=μ_h_cat,    μm=μ_h,        Σ=Σ_h),
        (μc=μ_hg_cat,   μm=μ_hg,       Σ=Σ_hg),
        (μc=μ_dr2_cat,  μm=sim.μ_dr2,  Σ=Σ_dr2),
        (μc=μ_dr32_cat, μm=sim.μ_dr32, Σ=Σ_dr32),
        (μc=μ_dr3_cat,  μm=sim.μ_dr3,  Σ=Σ_dr3),
        (μc=SVector(sim.μ_1_3), μm=SVector(sim.UEVA_model),
         Σ=SMatrix{1,1,T,1}(sim.UEVA_unc^2)),
    )
    μ_catalog_full = vcat(map(b -> b.μc, blockvals)...)
    μ_model_full = vcat(map(b -> b.μm, blockvals)...)

    @no_escape ctx.buf begin
        Σ_full = @alloc(T, _g23h_nchan, _g23h_nchan); fill!(Σ_full, zero(T))
        _g23h_place_diag!(Σ_full, blockvals)
        # Cross-epoch correlation between the DR2 and DR3 proper motions.
        K = ρ_dr3_dr2 * sqrt(Σ_dr2) * sqrt(Σ_dr3)'
        _g23h_place_cross!(Σ_full, _g23h_offsetof(:dr2), _g23h_offsetof(:dr3), K)

        idx = @alloc(Int, n_components)
        let k = 0
            @inbounds for i in 1:_g23h_nchan
                mask[i] && (idx[k+=1] = i)
            end
        end
        μ_cat = @alloc(T, n_components)
        μ_mod = @alloc(T, n_components)
        @inbounds for k in 1:n_components
            μ_cat[k] = μ_catalog_full[idx[k]]
            μ_mod[k] = μ_model_full[idx[k]]
        end
        Σ_sel = @alloc(T, n_components, n_components)
        @inbounds for kj in 1:n_components, ki in 1:n_components
            Σ_sel[ki, kj] = Σ_full[idx[ki], idx[kj]]
        end

        if _G23H_DEBUG_PULLS[] !== nothing
            lbls = [_g23h_channel_kinds[idx[k]] for k in 1:n_components]
            c = [Float64(μ_cat[k]) for k in 1:n_components]
            m = [Float64(μ_mod[k]) for k in 1:n_components]
            s = [sqrt(Float64(Σ_sel[k, k])) for k in 1:n_components]
            Σcap = [Float64(Σ_sel[ki, kj]) for ki in 1:n_components, kj in 1:n_components]
            push!(_G23H_DEBUG_PULLS[],
                (; labels=lbls, catalog=c, model=m, sigma=s, pull=(c .- m) ./ s, Σ=Σcap))
        end

        if n_components == 1
            ll += logpdf(Normal(μ_cat[1], sqrt(Σ_sel[1, 1])), μ_mod[1])
        elseif n_components > 1
            # MvNormal log-density by an in-place Cholesky on a bump buffer —
            # `Distributions.MvNormal` would re-wrap and re-factorize Σ on
            # every call through a heap-allocating PDMat.
            #   log p = −½(δ'Σ⁻¹δ + n·log 2π + log|Σ|),
            #   δ'Σ⁻¹δ = ‖L⁻¹δ‖²,  log|Σ| = 2·Σ log L_ii.
            L = @alloc(T, n_components, n_components)
            @inbounds for kj in 1:n_components, ki in 1:n_components
                L[ki, kj] = Σ_sel[ki, kj]
            end
            δ = @alloc(T, n_components)
            @inbounds for k in 1:n_components
                δ[k] = μ_mod[k] - μ_cat[k]
            end
            try
                F = cholesky!(Hermitian(L, :L))
                ldiv!(F.L, δ)
                quad = zero(T)
                logdetΣ = zero(T)
                @inbounds for k in 1:n_components
                    quad += δ[k] * δ[k]
                    logdetΣ += log(L[k, k])
                end
                ll += -T(0.5) * (quad + n_components * log(T(2π)) + 2logdetΣ)
            catch err
                err isa PosDefException ? (ll = convert(T0, -Inf)) : rethrow(err)
            end
        end
    end
    return ll
end

"""
    _g23h_catalog_moments(obs, ctx) -> (; labels, catalog, model, sigma, pull, Σ)

The coupled catalog block as the likelihood actually assembles it, at one
sample: the selected channel names, the catalog values, the model values, and
the marginal σ — after the per-channel proper-motion jitter, the UEVA-driven
DR3 deflation and its DR3−DR2 cross term, the BINARYS f_σ inflation of Σ_h and
its epistemic term.

It gets them by running `ln_like` with the diagnostic hook above armed, rather
than by re-deriving any of it. Reassembling ~60 lines of covariance algebra
plot-side is the drift the v9 plotting protocol exists to prevent — v8's
`absastromplot` did exactly that, with line-number citations into this file —
and the hook was already here to ask whether the simulator and the likelihood
agree, which is the same question a residual plot asks.
"""
function _g23h_catalog_moments(obs::G23HObs, ctx::ObsContext)
    return lock(_G23H_DEBUG_PULLS_LOCK) do
        sink = Any[]
        old = _G23H_DEBUG_PULLS[]
        _G23H_DEBUG_PULLS[] = sink
        try
            ln_like(obs, ctx)
        finally
            _G23H_DEBUG_PULLS[] = old
        end
        isempty(sink) && error(
            "G23HObs: the catalog block was not evaluated at this sample — the sampled " *
            "transit selection was rejected — so there are no residuals to report.")
        return last(sink)
    end
end

# Rebuild a channel's catalog MvNormal with extra proper-motion jitter, when
# the model samples any. Returns the original distribution untouched otherwise
# (and `nothing` stays `nothing`).
@inline function _g23h_jittered(dist, θ_obs, kra::Symbol, kdec::Symbol,
                                σ_ra, σ_dec, ρ, ::Type{T}) where {T}
    (isnothing(dist) || !(hasproperty(θ_obs, kra) || hasproperty(θ_obs, kdec))) && return dist
    jra = hasproperty(θ_obs, kra) ? getproperty(θ_obs, kra) : zero(T)
    jdec = hasproperty(θ_obs, kdec) ? getproperty(θ_obs, kdec) : zero(T)
    c = ρ[1] * σ_ra[1] * σ_dec[1]
    μ = mean(dist)
    return MvNormal(SVector(μ[1], μ[2]),
        @SArray [σ_ra[1]^2+jra^2 c; c σ_dec[1]^2+jdec^2])
end

# ──────────────────────────────────────────────────────────────────────
# Forward simulation of a synthetic catalog
# ──────────────────────────────────────────────────────────────────────

"""
    generate_from_params(obs::G23HObs, ctx; add_noise)

A new `G23HObs` whose catalog values are those this sample predicts. Every
covariance the likelihood will assemble at these parameters is mirrored here,
including the DR2/DR3 cross-correlation and the UEVA-driven deflation —
independent draws would leave the two conditional DR2/DR3 directions
over-dispersed by 1/(1−ρ²) in whitened space, and a fit reads that excess as
astrometric acceleration, i.e. spurious decades-period companions.
"""
function generate_from_params(obs::G23HObs, ctx::ObsContext; add_noise)
    sim = simulate(obs, ctx)
    (; μ_h, μ_hg, μ_dr2, μ_dr32, μ_dr3, UEVA_model, UEVA_unc, μ_1_3, sample_variance) = sim
    cat = obs.catalog
    θ_obs = ctx.θ_obs
    has_hip = !isnothing(cat.dist_hip)

    # `ln_like` adds the HGCA nonlinear correction to the model μ for an
    # absolute frame; without mirroring it here the synthetic data would sit
    # on a different convention than the likelihood's model — a phantom
    # proper-motion anomaly at the true parameters. Same helper, so the two
    # cannot come to disagree about when it applies.
    let (Δ_hg, Δ_h) = _g23h_nonlinear(cat, _g23h_isabs(ctx.system.frame), has_hip)
        μ_hg = μ_hg + Δ_hg
        μ_h = μ_h + Δ_h
    end

    (; σ_att, σ_AL, σ_calib) = θ_obs
    σ_formal = sqrt(σ_att^2 + σ_AL^2)
    N = cat.astrometric_n_good_obs_al_dr3
    N_FoV = cat.astrometric_matched_transits_dr3
    N_AL = N / N_FoV
    gaia_n_dof = cat.astrometric_params_solved_dr3 == 31 ? 5 : 6
    μ_UEVA_single = (N_AL / (N - gaia_n_dof)) *
                    ((N_FoV - gaia_n_dof) * σ_calib^2 + N_FoV * σ_AL^2)

    if obs.ueva_mode === :none
        # No UEVA datum, so no inflation/deflation to simulate.
        inflation_dr3 = 1.0
        deflation_truth = 1.0
        new_chi2_al = cat.astrometric_chi2_al_dr3
        new_ruwe = cat.ruwe_dr3
        new_ean = cat.astrometric_excess_noise_dr3
    else
        # A companion degrades Gaia's five-parameter fit, raising χ² and
        # RUWE/EAN; Gaia then inflates every formal uncertainty by that
        # excess, and the likelihood deflates them back when the companion
        # model explains it. The synthetic catalog must therefore be inflated
        # by the companion we inject, or the deflation has nothing to bite on.
        new_UEVA = max(UEVA_model + (add_noise ? randn() * UEVA_unc : 0.0), 0.0)^3
        UEVA_original = μ_1_3^3
        inflation_dr3 = sqrt(max(1.0, new_UEVA / max(eps(), UEVA_original)))
        new_chi2_al = max(Float64(N - gaia_n_dof), new_UEVA * (N - gaia_n_dof) / σ_formal^2)
        old_chi2 = cat.astrometric_chi2_al_dr3
        new_ruwe = old_chi2 > 0 ? cat.ruwe_dr3 * sqrt(new_chi2_al / old_chi2) : cat.ruwe_dr3
        new_ean = sqrt(max(0.0, new_UEVA - σ_att^2 - σ_AL^2))
        # The deflation the fit will apply at truth: UEVA reconstructed from
        # the NEW catalog is max(σ_formal², new_UEVA) in both modes (the
        # clamps above make the two expressions agree).
        deflation_truth = min(1.0, sqrt(μ_UEVA_single / max(σ_formal^2, new_UEVA)))
    end

    new_pmra_dr3_error = cat.pmra_dr3_error[1] * inflation_dr3
    new_pmdec_dr3_error = cat.pmdec_dr3_error[1] * inflation_dr3
    new_ra_err_c3 = cat.ra_error_central_dr3 * inflation_dr3
    new_dec_err_c3 = cat.dec_error_central_dr3 * inflation_dr3
    # DR32 errors inflate too (they depend on DR3 positions). This slightly
    # overestimates, since DR2 contributes as well, but the ΔΣ_dr32 correction
    # in `ln_like` handles the fine structure.
    new_pmra_dr32_error = cat.pmra_dr32_error[1] * inflation_dr3
    new_pmdec_dr32_error = cat.pmdec_dr32_error[1] * inflation_dr3

    function _draw_pm(μ, σa, σb, ρ)
        add_noise || return SVector{2,Float64}(μ)
        c = ρ * σa * σb
        Σ = @SArray [σa^2 c; c σb^2]
        return μ .+ cholesky(Hermitian(Σ)).L * @SVector [randn(), randn()]
    end

    if has_hip
        new_pm_hip = _draw_pm(μ_h, cat.pmra_hip_error[1], cat.pmdec_hip_error[1], cat.pmra_pmdec_hip[1])
        new_pm_hg = _draw_pm(μ_hg, cat.pmra_hg_error[1], cat.pmdec_hg_error[1], cat.pmra_pmdec_hg[1])
    end

    ρ_dr3_dr2 = cat.rho_dr2_dr3
    c2 = cat.pmra_pmdec_dr2[1] * cat.pmra_dr2_error[1] * cat.pmdec_dr2_error[1]
    Σ_dr2_gen = @SArray [cat.pmra_dr2_error[1]^2 c2; c2 cat.pmdec_dr2_error[1]^2]
    c3 = cat.pmra_pmdec_dr3[1] * new_pmra_dr3_error * new_pmdec_dr3_error
    Σ_dr3_gen = (@SArray [new_pmra_dr3_error^2 c3; c3 new_pmdec_dr3_error^2]) .* deflation_truth^2
    if add_noise
        K = ρ_dr3_dr2 * sqrt(Σ_dr2_gen) * sqrt(Σ_dr3_gen)'
        Σ_23 = SMatrix{4,4,Float64,16}([Σ_dr2_gen K; K' Σ_dr3_gen])
        n23 = cholesky(Hermitian(Σ_23)).L * @SVector [randn(), randn(), randn(), randn()]
        new_pm_dr2 = μ_dr2 .+ @SVector [n23[1], n23[2]]
        new_pm_dr3 = μ_dr3 .+ @SVector [n23[3], n23[4]]
    else
        new_pm_dr2 = SVector{2,Float64}(μ_dr2)
        new_pm_dr3 = SVector{2,Float64}(μ_dr3)
    end

    Σ_pos3 = @SMatrix [new_ra_err_c3^2 cat.ra_dec_corr_central_dr3*new_ra_err_c3*new_dec_err_c3
        cat.ra_dec_corr_central_dr3*new_ra_err_c3*new_dec_err_c3 new_dec_err_c3^2]
    Σ_cross = @SMatrix [
        ρ_dr3_dr2*new_ra_err_c3*cat.ra_error_central_dr2 ρ_dr3_dr2*cat.ra_dec_corr_central_dr3*new_ra_err_c3*cat.dec_error_central_dr2
        ρ_dr3_dr2*cat.ra_dec_corr_central_dr2*new_dec_err_c3*cat.ra_error_central_dr2 ρ_dr3_dr2*new_dec_err_c3*cat.dec_error_central_dr2
    ]
    ΔΣ_pos = (deflation_truth^2 - 1) * Σ_pos3 - (deflation_truth - 1) * (Σ_cross + Σ_cross')
    Δt_ra = (cat.epoch_ra_dr3_mjd - cat.epoch_ra_dr2_mjd) / julian_year
    Δt_dec = (cat.epoch_dec_dr3_mjd - cat.epoch_dec_dr2_mjd) / julian_year
    Tr = @SMatrix [1/Δt_ra 0.0; 0.0 1/Δt_dec]
    c32 = cat.pmra_pmdec_dr32[1] * new_pmra_dr32_error * new_pmdec_dr32_error
    Σ_dr32_gen = (@SArray [new_pmra_dr32_error^2 c32; c32 new_pmdec_dr32_error^2]) +
                 Tr * ΔΣ_pos * Tr'
    new_pm_dr32 = if add_noise
        μ_dr32 .+ cholesky(Hermitian(Σ_dr32_gen)).L * @SVector [randn(), randn()]
    else
        SVector{2,Float64}(μ_dr32)
    end

    # ---- Hipparcos abscissa residuals --------------------------------------
    n_hip = obs.n_hip
    new_hip_res = zeros(n_hip)
    if has_hip && n_hip > 0
        b = zeros(n_hip)
        @inbounds for i in 1:n_hip
            b[i] = sim.Δα_mas_hip[i] * obs.hip_table.cosϕ[i] +
                   sim.Δδ_mas_hip[i] * obs.hip_table.sinϕ[i]
        end
        # The five-parameter fit absorbs position, proper motion and parallax;
        # the residual is the curvature the companion leaves behind.
        new_hip_res .= b .- obs.A_prepared_5_hip * (obs.A_prepared_5_hip \ b)
        if add_noise
            # `hip_iad_jitter` only exists while the :iad_hip channel is
            # active; without it, simulate with zero excess jitter — the
            # residuals feed no active likelihood term then.
            jit = :iad_hip ∈ obs.table.kind ? θ_obs.hip_iad_jitter : 0.0
            @inbounds for i in 1:n_hip
                new_hip_res[i] += randn() *
                                  hypot(obs.hip_table.sres_renorm[i] * sim.σ_inflation_hip[i], jit)
            end
        end
    end

    # ---- Gaia RV -----------------------------------------------------------
    has_rv = :rv_dr3 ∈ obs.table.kind
    new_rv_error = hasproperty(cat, :radial_velocity_error) ? cat.radial_velocity_error : NaN
    if has_rv && isfinite(sample_variance)
        σ_rv = θ_obs.σ_rv_per_transit
        N_rv = cat.rv_nb_transits
        ncp = max(0.0, (N_rv - 1) * sample_variance / σ_rv^2)
        ξ² = add_noise ? rand(NoncentralChisq(max(1, N_rv - 1), ncp)) : ncp + (N_rv - 1)
        S² = ξ² * σ_rv^2 / max(1, N_rv - 1)
        # s² = (2N/π)(ε² − 0.113²)  ⟹  ε = √(s²π/(2N) + 0.113²)
        new_rv_error = sqrt(max(0.0, S² * π / (2 * N_rv) + 0.113^2))
    end

    # ---- rebuild ------------------------------------------------------------
    new_cat = (; cat...)
    if has_hip
        new_cat = (; new_cat..., pmra_hip=new_pm_hip[1], pmdec_hip=new_pm_hip[2],
            pmra_hg=new_pm_hg[1], pmdec_hg=new_pm_hg[2])
    end
    new_cat = (; new_cat...,
        pmra_dr2=new_pm_dr2[1], pmdec_dr2=new_pm_dr2[2],
        pmra_dr32=new_pm_dr32[1], pmdec_dr32=new_pm_dr32[2],
        pmra_dr3=new_pm_dr3[1], pmdec_dr3=new_pm_dr3[2],
        astrometric_chi2_al_dr3=new_chi2_al, ruwe_dr3=new_ruwe,
        astrometric_excess_noise_dr3=new_ean,
        pmra_dr3_error=new_pmra_dr3_error, pmdec_dr3_error=new_pmdec_dr3_error,
        ra_error_central_dr3=new_ra_err_c3, dec_error_central_dr3=new_dec_err_c3,
        pmra_dr32_error=new_pmra_dr32_error, pmdec_dr32_error=new_pmdec_dr32_error)
    has_rv && (new_cat = (; new_cat..., radial_velocity_error=new_rv_error))

    dist_hip = has_hip ? _g23h_pm_dist(new_cat.pmra_hip, new_cat.pmdec_hip,
        cat.pmra_hip_error[1], cat.pmdec_hip_error[1], cat.pmra_pmdec_hip[1]) : nothing
    dist_hg = has_hip ? _g23h_pm_dist(new_cat.pmra_hg, new_cat.pmdec_hg,
        cat.pmra_hg_error[1], cat.pmdec_hg_error[1], cat.pmra_pmdec_hg[1]) : nothing
    dist_dr2 = _g23h_pm_dist(new_cat.pmra_dr2, new_cat.pmdec_dr2,
        cat.pmra_dr2_error[1], cat.pmdec_dr2_error[1], cat.pmra_pmdec_dr2[1])
    dist_dr32 = _g23h_pm_dist(new_cat.pmra_dr32, new_cat.pmdec_dr32,
        new_pmra_dr32_error, new_pmdec_dr32_error, cat.pmra_pmdec_dr32[1])
    dist_dr3 = _g23h_pm_dist(new_cat.pmra_dr3, new_cat.pmdec_dr3,
        new_pmra_dr3_error, new_pmdec_dr3_error, cat.pmra_pmdec_dr3[1])
    new_cat = (; new_cat..., dist_hip, dist_hg, dist_dr2, dist_dr32, dist_dr3)

    new_hip_table = if n_hip > 0
        cols = Tables.columntable(obs.hip_table)
        Table(merge(cols, (; res=new_hip_res)))
    else
        obs.hip_table
    end

    new_pm = copy(obs.table.pm)
    new_σ_pm = copy(obs.table.σ_pm)
    for (i, kind) in enumerate(obs.table.kind)
        if kind === :ra_hip && has_hip
            new_pm[i] = new_cat.pmra_hip
        elseif kind === :dec_hip && has_hip
            new_pm[i] = new_cat.pmdec_hip
        elseif kind === :ra_hg && has_hip
            new_pm[i] = new_cat.pmra_hg
        elseif kind === :dec_hg && has_hip
            new_pm[i] = new_cat.pmdec_hg
        elseif kind === :ra_dr2
            new_pm[i] = new_cat.pmra_dr2
        elseif kind === :dec_dr2
            new_pm[i] = new_cat.pmdec_dr2
        elseif kind === :ra_dr32
            new_pm[i] = new_cat.pmra_dr32
            new_σ_pm[i] = new_pmra_dr32_error
        elseif kind === :dec_dr32
            new_pm[i] = new_cat.pmdec_dr32
            new_σ_pm[i] = new_pmdec_dr32_error
        elseif kind === :ra_dr3
            new_pm[i] = new_cat.pmra_dr3
            new_σ_pm[i] = new_pmra_dr3_error
        elseif kind === :dec_dr3
            new_pm[i] = new_cat.pmdec_dr3
            new_σ_pm[i] = new_pmdec_dr3_error
        elseif kind === :ueva_dr3
            new_pm[i] = obs.ueva_mode === :RUWE ? new_cat.ruwe_dr3 :
                        new_cat.astrometric_excess_noise_dr3
        end
    end
    new_table = Table(merge(Tables.columntable(obs.table), (; pm=new_pm, σ_pm=new_σ_pm)))

    # The caches derive from hip_table.{sres,res}, so recompute them.
    pinv_5_hip = obs.pinv_5_hip
    hip_x_const = _g23h_hip_x_const(pinv_5_hip, new_hip_table, obs.include_iad, has_hip)

    return G23HObs{typeof(new_table),typeof(new_hip_table),typeof(obs.gaia_table),
        typeof(new_cat),typeof(obs.hip_sol),typeof(obs.host),typeof(obs.companions),
        typeof(obs.ref)}(
        new_table, obs.priors, obs.derived, new_hip_table, obs.gaia_table, new_cat, obs.hip_sol,
        obs.A_prepared_5_hip, obs.A_prepared_5_dr2, obs.A_prepared_5_dr3,
        pinv_5_hip, hip_x_const, obs.include_iad, obs.ueva_mode, obs.frame_shift, obs.name,
        obs.host, obs.companions, obs.ref,
        obs.n_hip, obs.n_gaia, obs.i_ra_hip, obs.i_dec_hip,
        obs.i_ra_dr2, obs.i_dec_dr2, obs.i_ra_dr3, obs.i_dec_dr3, obs.all_epochs)
end
