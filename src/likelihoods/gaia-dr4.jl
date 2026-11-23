# ---------------------------------------------------
# Gaia DR4 epoch astrometry (along-scan)
#
# The v1 version carried three things that no longer exist here:
#   - a hand-summed loop over companions with `raoff(sol, mass)` scaling,
#     which only worked because every companion orbited the same star;
#   - a branch on the orbit *type* (`OrbitSolutionAbsoluteVisual`) to pick
#     between two barycentre-position models;
#   - a `fluxratio` vector indexed positionally by companion.
#
# All three are now one line: the photocentre's offset from whatever
# reference the observation declares. Blended photocentres of several
# luminous bodies are correct by construction, because the weights are
# recomputed from the bodies' fluxes every sample.
#
# The reference point — where the excursion is measured from — has two
# sources, chosen by the system's frame:
#
#   - parallax-only (`plx` alone): a straight line in the local plane built
#     from this observation's own `ra_offset_mas`/`dec_offset_mas`/`pmra`/
#     `pmdec`/`ref_epoch` variables, and the constant `plx`. Unchanged since
#     the v2 port, bit for bit.
#   - absolute frame (`plx, ra, dec, pmra, pmdec, rv, ref_epoch`): the
#     barycentre's *propagated* direction from PlanetOrbits (`frame_ra`/
#     `frame_dec`), projected gnomonically about the tangent point the data
#     are counted from (`ra0`, `dec0`), with the parallax at the epoch. Before
#     this, the frame's propagation — perspective acceleration, the
#     epoch-dependent parallax, the light-time terms — never reached the DR4
#     likelihood: the reference point was the straight line even when the
#     system carried the full frame.
#
# Under the absolute frame the excursion (gnomonic about the barycentre's
# direction, as PlanetOrbits reports it) is not *added* to the projected
# barycentre: the two are coordinates in different tangent planes, and the
# sum is wrong by ~ |excursion| · |offset from the tangent point| · tan δ —
# mas for an arcsecond-scale photocentre–barycentre separation. The source's
# direction is rebuilt on the sphere and projected about the tangent point
# instead (`_dr4_project`), and its parallax is its own (`_dr4_source_plx`).
# ---------------------------------------------------

"""
    GaiaDR4AstromObs(data; target=Photocentre, ref=Barycentre, name="GaiaDR4",
                     detrend=false, tangent_point=nothing, variables=…)

Gaia DR4 individual (epoch) astrometry: one along-scan abscissa per
transit.

`data` needs `:epoch` [MJD], `:scan_pos_angle` ψ [**degrees**],
`:parallax_factor_al` [dimensionless], `:centroid_pos_al` [mas] and
`:centroid_pos_error_al` [mas]; an `:outlier_flag` column is honoured if
present. The DR4 epoch-astrometry table's `:ra0`/`:dec0` columns [deg] —
the tangent point the local-plane coordinates are counted from — are read
if present, and are **required** when the system has an absolute frame (see
"Reference point" below); `tangent_point=(ra0, dec0)` [deg] supplies them
for a table that lacks the columns.

# Units
`scan_pos_angle` is the only angle-like column, and it is taken in
**degrees** — the unit the Gaia archive publishes it in (the DR4 VOTABLE
declares `unit="deg"`), so a table read straight out of the archive needs no
conversion. The conversion to radians happens once here, at construction:
`obs.table` keeps the data verbatim in degrees, while the `sinψ`/`cosψ` the
likelihood projects with are radian-based derived quantities. Everything
else is unchanged: epochs in MJD, along-scan positions and uncertainties in
mas, parallax factors dimensionless.

The model for each transit is

    w = Δα* sin ψ + Δδ cos ψ + ϖ · f_al

where `Δα*, Δδ` are the local-plane coordinates of the reference point (see
below) plus the offset of `target` from `ref` — normally the photocentre
relative to the system barycentre.

# Reference point: parallax-only vs absolute frame

Which reference-point model applies is decided by the frame the **system**
block defines, never by a keyword here.

**Parallax-only system** (`plx` alone). The reference point moves in a
straight line in the local plane, from this observation's own variables:

    Δα* = ra_offset_mas + pmra · (t − ref_epoch)
    Δδ  = dec_offset_mas + pmdec · (t − ref_epoch)
    ϖ   = plx                                      # constant

so `ra_offset_mas`/`dec_offset_mas` absorb the offset between the data's
tangent point and the barycentre, and no tangent point is needed.

**Absolute frame** (`plx, ra, dec, pmra, pmdec, rv, ref_epoch` in the system
block). The system block owns the barycentre's position and motion, and
PlanetOrbits propagates it rigorously in 3D. At each transit the modelled
point is the `target`'s apparent direction — the barycentre's propagated
direction `(frame_ra(sol), frame_dec(sol))` displaced on the sphere by the
target's excursion from it — projected **gnomonically** about the tangent
point `(ra0, dec0)` that DR4's `centroid_pos_al` is counted from. The
parallax term uses the target's own parallax *at that epoch*,
`1000 / (d(t) + Δz)` [mas] for the barycentre's propagated distance `d(t)`
and the target's line-of-sight offset `Δz` from it:

    (Δα*, Δδ) = gnomonic(direction of target; about = (ra0, dec0))   [mas]
    ϖ         = ϖ_target(t)

So perspective acceleration (≈ 5 mas over DR4 for Barnard's star), the
change of parallax with distance (≈ 0.18 mas across the DR4 window for
Barnard's star, ±0.09 mas about the reference epoch), and — when
`barycentric_lighttime` is on — the light-time terms all reach the
abscissae. The system's `plx` is the barycentre's parallax at the frame's
`ref_epoch`.

**The absolute frame is Gaia DR4's.** The system's `ra`, `dec`, `pmra`,
`pmdec` are the barycentre's position and proper motion *in the DR4 frame*
at `ref_epoch`, so this observation has no zero-point offset of its own.
Other absolute datasets in the same model (Hipparcos, earlier Gaia releases,
HST/FGS) carry their offsets *to* DR4 — the same arrangement as a single
reference timescale with every other clock converted into it. (`G23HObs`
does not yet carry a DR3→DR4 offset.) Several DR4 sources in one system —
the two components of a resolved wide pair, say — therefore share the frame
with nothing to tie, and their relative astrometry is kept whole.

Under an absolute frame the observation's reference-point variables
(`ra_offset_mas`, `dec_offset_mas`, `pmra`, `pmdec`, `ref_epoch`) must **not**
be declared, and the model refuses to build if any is: a second position or
proper motion here would either double count the motion or silently override
the frame, the same rule `G23HObs` follows for proper motion. `ref` must be
the whole-system `Barycentre`, because the frame describes the barycentre and
the excursion has to be measured from it.

The DR4 position is known only to a fraction of a mas, while the system's
`ra`/`dec` are in degrees. Write them as mas offsets from a catalogue
position (the tangent point is a convenient one) and derive the angles:

    Δra_mas  ~ Normal(0, 1000)
    Δdec_mas ~ Normal(0, 1000)
    ra  = \$ra0 + Δra_mas / (3.6e6 * cosd(\$dec0))
    dec = \$dec0 + Δdec_mas / 3.6e6

The same goes for the anchor variables of an [`AnchoredFrame`](@ref).

The tangent point is a property of the data, not a parameter: it comes from
the table's `ra0`/`dec0` columns (which DR4 repeats on every row of a source,
and which must therefore be constant) or from `tangent_point=(ra0, dec0)`.
It is required under an absolute frame and ignored otherwise; giving both
the columns and the keyword with different values is an error.

The excursion PlanetOrbits reports (`raoff`/`decoff` from the barycentre)
is a pair of gnomonic coordinates about the barycentre's direction. It is
turned back into a direction on the sphere and *that* is projected about the
tangent point — exactly, with no small-angle step. Simply adding it to the
projected barycentre would mix two tangent planes, an error of order
`|excursion| · |barycentre − tangent point| · tan δ`: sub-µas for a compact
orbit, but 7 µas for a 1″ excursion at δ = 41° and ~2 mas for 20″, i.e. for
the resolved components of a wide pair. `detrend` acts on the excursion
before that step, exactly as under a parallax-only frame.

# Parallax factors
`f_al` is Gaia's own `parallax_factor_al`, taken at face value: this
observation consumes SSB observables plus explicit parallax factors, and
never observer-aware observables. Every absolute-astrometry observation type
does one or the other, never both, and which one is literal code in the type
so that it can be reviewed per type.

What Gaia's factors omit is the annual–orbital (Kopeikin) coupling — the
dependence of the parallax factor on the *companion's* line-of-sight depth,
which is `≈ 4.85 · z[AU] / d[pc]²` µas per AU of observer displacement. For
any DR4 target that is sub-µas, well below the per-transit precision, so the
face-value factors are exact enough. PlanetOrbits' observer-aware
observables exist for the cases where it is not.

What *is* kept under an absolute frame is the target's own distance: the
parallax multiplying `f_al` is the target's, `1000 / (d + Δz)`, not the
barycentre's. The difference is `≈ ϖ² Δz / 2.06×10⁸` mas — 1e-8 of the
parallax for a planet an AU out, ~1e-4 for a component of a few-hundred-au
pair at 13 pc. Under a parallax-only frame the system's constant `plx` is
used, as before.

# `detrend`
With `detrend=true` the linear (constant + slope) part of the
target-versus-reference excursion is removed before it enters the model, so
only the curvature does. That breaks the degeneracy between the fitted
proper motion and a wide companion, and makes the reference-point position
and proper motion — this observation's variables, or the system's frame
under an absolute frame — describe the *photocentre* rather than the
barycentre.

# Variables
`astrometric_jitter` [mas] adds in quadrature to the per-transit formal
error. Under a parallax-only frame `ra_offset_mas`, `dec_offset_mas`,
`pmra`, `pmdec`, `ref_epoch` define the reference-point motion; under an
absolute frame they are the system's job and are rejected here (see above).
The parallax comes from the system block in both cases.

# Blended sources

`target` is what the *catalog source* is, and a catalog source is not
generally a body: it is whatever flux the pipeline blended into one
centroid. `Photocentre` (the default) is the whole system's flux-weighted
point; `Photocentre(:G, (Aa, Ab))` is the point over a named subset.

Two sources in a 2+2 quadruple — two tight pairs several arcseconds apart,
so that only intra-pair blending is possible — are two instances of this
observation, each with its own nuisance parameters, sharing the system's
`plx` and frame:

    System(name=:quad, bodies=(Aa, Ab, Ba, Bb, wide), observations=(
        GaiaDR4AstromObs(scans_A; target=Photocentre(:G, (Aa, Ab)),
                         ref=Barycentre, name="srcA", variables=@variables begin
                             ra_offset_mas ~ Normal(0, 100); dec_offset_mas ~ Normal(0, 100)
                             pmra ~ Normal(0, 100); pmdec ~ Normal(0, 100); ref_epoch = 57388.5
                         end),
        GaiaDR4AstromObs(scans_B; target=Photocentre(:G, (Ba, Bb)),
                         ref=Barycentre, name="srcB", variables=…),
    ), variables=…)

Each source's modelled signal then carries *both* its pair's wide-orbit
motion and the intra-pair photocentric wobble, because a photocentre is one
dot product over absolute body states — there is no per-level bookkeeping to
get wrong. Bodies declare `flux_G` in their own blocks; setting the host's
to 1.0 makes the others contrast ratios.

Membership that is not structurally fixed — a sampled resolved-flag, a
scan-angle-dependent window — is not expressible as a static spec, and is
not meant to be: an observation of that kind reads
`PlanetOrbits.fluxes(sys, band)` and builds its own `WeightedPoint` per draw
or per epoch.
"""
struct GaiaDR4AstromObs{TTable<:Table,TT,TR} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    target::TT
    ref::TR
    name::String
    detrend::Bool
    # Precomputed detrending coefficients: with the epochs fixed, the
    # least-squares constant and slope of any excursion are two dot products.
    detrend_Δt::Vector{Float64}
    detrend_inv_N::Float64
    detrend_inv_sum_Δt²::Float64
    # Per-transit constants of the data, hoisted out of the per-sample loops.
    # The scan geometry and formal errors are fixed at construction, so the
    # projection's sincos and the Gaussian normalization's inputs are data,
    # not per-sample work (the sincos alone was 10% of the DR4 primal).
    #
    # These are also the *only* place ψ is ever in radians: `table` keeps the
    # user's degrees verbatim, so reconstructing an observation from its own
    # table (`likeobj_from_epoch_subset`, `generate_from_params`) can never
    # double-convert. Anything downstream that needs the projection — the
    # plotting extension included — must read these, not `table.scan_pos_angle`.
    sinψ::Vector{Float64}
    cosψ::Vector{Float64}
    # Outlier handling as arithmetic rather than a branch: `w` is 1 for kept
    # transits and 0 for flagged ones, so the residual loop stays branch-free
    # and vectorizes. Flagged rows may carry NaN centroids/errors in real DR4
    # tables, and NaN·0 = NaN — so `b_al` and `σ²` are *sanitized* copies
    # (flagged rows read 0 and 1), while `table` keeps the data verbatim.
    w::Vector{Float64}
    b_al::Vector{Float64}
    σ²::Vector{Float64}
    invσ²::Vector{Float64}    # w/σ², i.e. already zero on flagged rows
    n_eff::Float64            # Σw
    # -½ Σ w·log(2π σ²): the whole normalization term when no jitter is
    # sampled, which is the production fast path (per-transit noise estimated
    # robustly up front rather than fit).
    lognorm0::Float64
    # The tangent point (ra0, dec0) [deg] the local-plane coordinates
    # `centroid_pos_al` are counted from, or `nothing` when neither the table
    # nor the caller says. Read only under an absolute frame, where the
    # reference point is the frame's direction projected about it; a
    # parallax-only model never touches it. See `_dr4_tangent_point`.
    tangent_point::Union{Nothing,NTuple{2,Float64}}
    # The tangent point's unit vector and local east/north axes, (u0, e_α0,
    # e_δ0) flattened, hoisted for the same reason as `sinψ`/`cosψ`: data,
    # not per-sample work. NaN when there is no tangent point.
    tangent_basis::NTuple{9,Float64}
end

const dr4_cols = (:epoch, :scan_pos_angle, :parallax_factor_al,
                  :centroid_pos_al, :centroid_pos_error_al)

# The observation variables that define the reference point under a
# parallax-only frame. Under an absolute frame the system block owns all of
# them, and declaring any here is a build-time error (`check_siblings`).
const DR4_REFPOINT_VARS = (:ra_offset_mas, :dec_offset_mas, :pmra, :pmdec, :ref_epoch)

"""
The tangent point `(ra0, dec0)` [deg] of a DR4 table, from its `ra0`/`dec0`
columns or the `tangent_point` keyword, or `nothing` if neither gives one.

DR4 repeats the per-source tangent point on every row ("ra0 and dec0 are
always the same for all transits of a given source_id"), so a column that
varies means the table mixes sources or has been edited, and that is an
error rather than something to average over. Both sources at once must
agree exactly: the keyword exists for tables that lack the columns, not to
override them.
"""
function _dr4_tangent_point(table, tangent_point)
    has_ra0 = hasproperty(table, :ra0)
    has_dec0 = hasproperty(table, :dec0)
    has_ra0 == has_dec0 || error(
        "GaiaDR4AstromObs: the table has a `$(has_ra0 ? :ra0 : :dec0)` column but no " *
        "`$(has_ra0 ? :dec0 : :ra0)` column. The tangent point is a pair; supply both " *
        "columns, or drop both and pass `tangent_point=(ra0, dec0)` [deg].")
    from_cols = if has_ra0
        ra0s = table.ra0
        dec0s = table.dec0
        (allequal(ra0s) && allequal(dec0s)) || error(
            "GaiaDR4AstromObs: the `ra0`/`dec0` columns are not constant. DR4 gives one " *
            "tangent point per source and repeats it on every transit, so a varying " *
            "column means this table mixes sources (or NaNs have crept in). One " *
            "observation is one source: split the table by `source_id`.")
        (Float64(first(ra0s)), Float64(first(dec0s)))
    else
        nothing
    end
    from_kw = if tangent_point === nothing
        nothing
    else
        length(tangent_point) == 2 || error(
            "GaiaDR4AstromObs: `tangent_point` takes `(ra0, dec0)` in degrees; got " *
            "$(repr(tangent_point)).")
        (Float64(tangent_point[1]), Float64(tangent_point[2]))
    end
    if from_cols !== nothing && from_kw !== nothing && from_cols != from_kw
        error("GaiaDR4AstromObs: `tangent_point=$(from_kw)` disagrees with the table's " *
              "own `ra0`/`dec0` = $(from_cols). The tangent point is where the data's " *
              "local-plane coordinates are counted from, so it is fixed by the table; " *
              "drop the keyword.")
    end
    tp = from_cols !== nothing ? from_cols : from_kw
    if tp !== nothing
        ra0, dec0 = tp
        (isfinite(ra0) && isfinite(dec0) && -90 <= dec0 <= 90) || error(
            "GaiaDR4AstromObs: the tangent point (ra0, dec0) = $(tp) is not a valid " *
            "direction in degrees.")
    end
    return tp
end

function GaiaDR4AstromObs(observations;
                          target=Photocentre, ref=Barycentre,
                          name="GaiaDR4",
                          detrend::Bool=false,
                          tangent_point=nothing,
                          variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    # Collect every column up front: a multithreaded `CSV.read` hands back
    # `ChainedVector`s, and the per-transit indexing below asserts on them.
    # See `materialize_cols`.
    table = materialize_cols(Table(observations))
    if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
        table = Table(table; epoch=jd2mjd.(table.obs_time_tcb))
    end
    equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    issubset(dr4_cols, Tables.columnnames(table)) ||
        error("Expected columns $dr4_cols")
    table = table[sortperm(vec(table.epoch))]

    ep = table.epoch
    Δt = collect((ep .- sum(ep) / length(ep)) ./ PlanetOrbits.year2day_julian)
    t, r = refspec(target), refspec(ref)

    # ψ is ingested in degrees (the archive's own unit) and converted here,
    # once. `sind`/`cosd` rather than `sin(deg2rad(·))`: the argument
    # reduction is exact on the degree scale, so a scan angle that is a whole
    # number of degrees projects exactly.
    sinψ = collect(Float64, sind.(table.scan_pos_angle))
    cosψ = collect(Float64, cosd.(table.scan_pos_angle))
    w = hasproperty(table, :outlier_flag) ?
        collect(Float64, table.outlier_flag .== 0) : ones(Float64, length(ep))
    b_al = collect(Float64, table.centroid_pos_al)
    σ² = collect(Float64, table.centroid_pos_error_al) .^ 2
    # Sanitize flagged rows (see the field comment): their values are never
    # *used* — every read is multiplied by w = 0 — but they must be finite.
    @. b_al = ifelse(w == 0, 0.0, b_al)
    @. σ² = ifelse(w == 0, 1.0, σ²)
    invσ² = w ./ σ²
    lognorm0 = -sum(w .* log.(2π .* σ²)) / 2
    tp = _dr4_tangent_point(table, tangent_point)
    tangent_basis = tp === nothing ? ntuple(_ -> NaN, 9) : _dr4_tangent_basis(tp...)
    return GaiaDR4AstromObs{typeof(table),typeof(t),typeof(r)}(
        table, priors, derived, t, r, String(name), detrend,
        Δt, 1.0 / length(ep), 1.0 / sum(Δt .^ 2),
        sinψ, cosψ, w, b_al, σ², invσ², sum(w), lognorm0,
        tp, tangent_basis)
end

export GaiaDR4AstromObs

refspecs(obs::GaiaDR4AstromObs) = (obs.target, obs.ref)

# The absolute-frame contract, checked once when the observation joins a
# system — the constructor cannot know which frame that system will have.
# Nothing to check for a parallax-only or frame-free system: that path is the
# one this observation always had.
function check_siblings(obs::GaiaDR4AstromObs, @nospecialize(all_obs), ctx)
    # A partial frame never gets this far (the `System` constructor rejects
    # it), so `ra` present means the full absolute frame.
    :ra in ctx.framevars || return nothing
    declared = Symbol[collect(keys(obs.priors.priors))...;
                      collect(keys(obs.derived.variables))...]
    clash = [v for v in DR4_REFPOINT_VARS if v in declared]
    isempty(clash) || error("""
        System $(ctx.name): GaiaDR4AstromObs "$(obs.name)" declares \
        $(join(("`" .* string.(clash) .* "`"), ", ")), but the system has an absolute \
        frame.

        Under an absolute frame the system block owns the barycentre's position and \
        proper motion (`ra`, `dec`, `pmra`, `pmdec` at `ref_epoch`), and this \
        observation's reference point is that frame, propagated by PlanetOrbits and \
        projected about the data's tangent point. A second position or proper motion \
        here would either double count the motion or override the frame, so it is \
        refused rather than silently resolved one way.

        Remove $(join(("`" .* string.(clash) .* "`"), ", ")) from this observation's \
        `variables`. The absolute frame is DR4's own: the system's `ra`/`dec` are the \
        barycentre's DR4 position, so there is no separate DR4 zero point to fit. \
        Write them as mas offsets from a catalogue position, e.g. \
        `Δra_mas ~ Normal(0, 1000)`, `ra = \$ra0 + Δra_mas / (3.6e6 * cosd(\$dec0))`. \
        To keep the linear reference-point model instead, define `plx` alone in the \
        system block (a parallax-only frame).
        """)
    obs.tangent_point === nothing && error("""
        System $(ctx.name): GaiaDR4AstromObs "$(obs.name)" has no tangent point, and the \
        system has an absolute frame.

        DR4's `centroid_pos_al` is a local-plane coordinate counted from a fixed \
        tangent point (`ra0`, `dec0` in the epoch-astrometry table). Under an absolute \
        frame the model predicts the barycentre's propagated direction and must project \
        it about that same point, so the point has to come from the data. Keep the \
        table's `ra0`/`dec0` columns, or pass `tangent_point=(ra0, dec0)` [deg] to \
        `GaiaDR4AstromObs`.
        """)
    obs.ref isa BarycentreSpec{()} || error("""
        System $(ctx.name): GaiaDR4AstromObs "$(obs.name)" has `ref=$(obs.ref)`, but the \
        system has an absolute frame.

        The absolute frame describes the whole system's barycentre, and this \
        observation adds the excursion of `target` *from `ref`* to the frame's \
        propagated track — so `ref` must be that barycentre. Use `ref=Barycentre`.
        """)
    return nothing
end

function likeobj_from_epoch_subset(obs::GaiaDR4AstromObs, inds)
    return GaiaDR4AstromObs(obs.table[inds, :, 1];
        target=obs.target, ref=obs.ref, obs.name, obs.detrend,
        tangent_point=obs.tangent_point,
        variables=(obs.priors, obs.derived))
end

# ---------------------------------------------------
# The absolute-frame reference point
# ---------------------------------------------------

@inline _dr4_isabs(::PlanetOrbits.AbsoluteFrame) = true
@inline _dr4_isabs(::PlanetOrbits.AbstractFrame) = false

const _DR4_RAD2MAS = rad2deg(1.0) * 3_600_000

"""
The unit vector `u0` of the tangent point `(α0, δ0)` [deg] and its local
east and north axes `e_α0`, `e_δ0`, flattened as a 9-tuple — the basis
[`_dr4_project`](@ref) projects onto.
"""
function _dr4_tangent_basis(α0, δ0)
    sα, cα = sincosd(α0)
    sδ, cδ = sincosd(δ0)
    return (cδ * cα, cδ * sα, sδ,          # u0
            -sα, cα, 0.0,                  # e_α0, toward increasing α
            -sδ * cα, -sδ * sα, cδ)        # e_δ0, toward increasing δ
end

"""
    _dr4_project(α, δ, x, y, basis)

Gnomonic coordinates `(ξ, η)` [mas] about the tangent point of the direction
whose gnomonic coordinates about `(α, δ)` [deg] are `(x, y)` [mas]. `basis`
is the tangent point's [`_dr4_tangent_basis`](@ref).

Gnomonic coordinates are central projection: `(x, y)` about a direction `u`
with local axes `e_α`, `e_δ` name the ray through `u + x e_α + y e_δ` (x, y in
radians). Projecting that ray about the tangent point is

    ξ = (r ⋅ e_α0) / (r ⋅ u0),   η = (r ⋅ e_δ0) / (r ⋅ u0)

and the normalization of `r` cancels. Both steps are exact, so the result is
the target's gnomonic position about the tangent point whatever the size of
the excursion or of the barycentre's offset from the tangent point — the
first-order alternative, `gnomonic(α, δ) + (x, y)`, mixes two tangent planes
and is wrong by `~ |(x, y)| ⋅ |(α, δ) − (α0, δ0)| ⋅ tan δ` (see the
"Reference point" section of [`GaiaDR4AstromObs`](@ref)). With `x = y = 0` it
is the gnomonic projection of `(α, δ)` itself.

ξ points toward increasing α and η toward increasing δ — the axes ψ is
measured in. The DR4 data model defines `ra0`/`dec0` as the tangent point
"from which the local plane coordinates centroid_pos_al and centroid_pos_ac
are counted". `r ⋅ e_α0` and `r ⋅ e_δ0` are differences of O(1) products
whose result is the arcsecond-scale offset; the cancellation costs ~1e-16
rad, i.e. ~2e-8 mas.
"""
@inline function _dr4_project(α, δ, x, y, basis)
    sα, cα = sincosd(α)
    sδ, cδ = sincosd(δ)
    xr = x / _DR4_RAD2MAS
    yr = y / _DR4_RAD2MAS
    # r = u + x e_α + y e_δ, unnormalized.
    rx = cδ * cα - xr * sα - yr * sδ * cα
    ry = cδ * sα + xr * cα - yr * sδ * sα
    rz = sδ + yr * cδ
    u0x, u0y, u0z, ax, ay, az, dx, dy, dz = basis
    scale = _DR4_RAD2MAS / (rx * u0x + ry * u0y + rz * u0z)
    return ((rx * ax + ry * ay + rz * az) * scale,
            (rx * dx + ry * dy + rz * dz) * scale)
end

"""
The system barycentre's propagated distance at this solution's epoch [pc].

PlanetOrbits has public accessors for every other propagated frame quantity
(`frame_ra`, `frame_dec`, `frame_pmra`, `frame_pmdec`, `frame_rv`) but not
for the distance, so this reads two of its internals: `frame(sol)` (the
effective frame of the solve — de-Dopplered under the light-time solve, see
its docstring) and `_propagate_dist`, at the solved emission epoch
`traj.t_em` — the same frame, function and epoch `frame_ra`/`frame_dec` are
computed from, so position and parallax describe one worldline. Replace with
a public accessor when PlanetOrbits gains one.
"""
@inline function _dr4_frame_dist(sol::PlanetOrbits._AbsSol)
    _, _, _, d_pc = PlanetOrbits._propagate_dist(PlanetOrbits.frame(sol),
                                                @inbounds sol.traj.t_em[sol.k])
    return d_pc
end

"""The system barycentre's parallax at this solution's epoch [mas], `1000 / d(t)`."""
@inline _dr4_frame_plx(sol::PlanetOrbits._AbsSol) = 1000 / _dr4_frame_dist(sol)

"""
The parallax [mas] of `target` at this solution's epoch: `1000 / (d + Δz)`
for the barycentre's propagated distance `d` [pc] and the target's
line-of-sight offset `Δz = posz(sol, target, reference)` [AU, +receding] from
it — the inverse of [`barycentre_parallax`](@ref). For a blended photocentre
`Δz` is the flux-weighted mean, which is its parallax to first order in
`Δz/d` (the next term is ~1e-10 of the parallax for a 100 au pair at 10 pc).
"""
@inline _dr4_source_plx(sol::PlanetOrbits._AbsSol, target, reference) =
    1000 / (_dr4_frame_dist(sol) + PlanetOrbits.posz(sol, target, reference) / PlanetOrbits.pc2au)

# The θ-level half of the `check_siblings` contract: `hasproperty` on the θ
# NamedTuple is a compile-time constant, so this folds to nothing on a valid
# model. It exists for contexts built by hand, which skip `System`'s checks.
@inline function _dr4_check_absolute_θ(θ_obs, name)
    if hasproperty(θ_obs, :ra_offset_mas) || hasproperty(θ_obs, :dec_offset_mas) ||
       hasproperty(θ_obs, :pmra) || hasproperty(θ_obs, :pmdec) ||
       hasproperty(θ_obs, :ref_epoch)
        _dr4_err_absolute_θ(name)
    end
    return nothing
end
@noinline _dr4_err_absolute_θ(name) = error(
    "GaiaDR4AstromObs \"$name\": this observation's variables define a reference " *
    "point ($(join(DR4_REFPOINT_VARS, ", "))), but the system has an absolute frame, " *
    "which owns it. Remove them, or define `plx` alone in the system block.")
@noinline _dr4_err_no_tangent(name) = error(
    "GaiaDR4AstromObs \"$name\": the system has an absolute frame but this observation " *
    "has no tangent point. Keep the table's `ra0`/`dec0` columns or pass " *
    "`tangent_point=(ra0, dec0)` [deg].")

"""
    simulate!(along_scan, ra_offset, dec_offset, obs, ctx)

Fill the modelled along-scan abscissae (and the RA/Dec offsets they were
projected from) for every transit.

`ra_offset`/`dec_offset` are local-plane coordinates [mas] of the modelled
point: about the tangent point `(ra0, dec0)` under an absolute frame, and
relative to the observation's own (arbitrary) origin under a parallax-only
one — see "Reference point" in the [`GaiaDR4AstromObs`](@ref) docstring.
"""
function simulate!(along_scan, ra_offset, dec_offset, obs::GaiaDR4AstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    tab = obs.table
    # A property of the frame's *type*, so every branch on it below folds at
    # compile time: each model compiles only its own reference-point path.
    absolute = _dr4_isabs(ctx.system.frame)

    if absolute
        # The excursion is accumulated from zero first, in the barycentre's
        # tangent plane; the loop after it moves the target onto the sphere
        # and projects it about the data's tangent point.
        _dr4_check_absolute_θ(θ_obs, obs.name)
        obs.tangent_point === nothing && _dr4_err_no_tangent(obs.name)
        fill!(ra_offset, zero(T))
        fill!(dec_offset, zero(T))
    else
        # Parallax-only: the reference point's own linear motion, from this
        # observation's variables. Unchanged since the v2 port.
        ref_epoch = hasproperty(θ_obs, :ref_epoch) ? θ_obs.ref_epoch : zero(T)
        ra0 = hasproperty(θ_obs, :ra_offset_mas) ? θ_obs.ra_offset_mas : zero(T)
        dec0 = hasproperty(θ_obs, :dec_offset_mas) ? θ_obs.dec_offset_mas : zero(T)
        pmra_ = hasproperty(θ_obs, :pmra) ? θ_obs.pmra : zero(T)
        pmdec_ = hasproperty(θ_obs, :pmdec) ? θ_obs.pmdec : zero(T)

        @inbounds for i in eachindex(tab.epoch)
            Δt = (tab.epoch[i] - ref_epoch) / PlanetOrbits.year2day_julian
            ra_offset[i] = ra0 + pmra_ * Δt
            dec_offset[i] = dec0 + pmdec_ * Δt
        end
    end

    # The excursion of the observed point about that reference. One call —
    # for a blended photocentre this is the flux-weighted mean of *apparent*
    # positions, not a superposition of per-companion terms, which is the
    # case v1's positional `fluxratio` vector could not represent.
    target = ref(ctx, obs.target)
    reference = ref(ctx, obs.ref)
    if obs.detrend
        sum_ra = zero(T); dot_ra = zero(T)
        sum_dec = zero(T); dot_dec = zero(T)
        @no_escape ctx.buf begin
            pra = @alloc(T, length(tab.epoch))
            pdec = @alloc(T, length(tab.epoch))
            @inbounds for i in eachindex(tab.epoch)
                sol = solutionat(ctx, i)
                pra[i] = raoff(sol, target, reference)
                pdec[i] = decoff(sol, target, reference)
                sum_ra += pra[i]; dot_ra += obs.detrend_Δt[i] * pra[i]
                sum_dec += pdec[i]; dot_dec += obs.detrend_Δt[i] * pdec[i]
            end
            mean_ra = sum_ra * obs.detrend_inv_N
            slope_ra = dot_ra * obs.detrend_inv_sum_Δt²
            mean_dec = sum_dec * obs.detrend_inv_N
            slope_dec = dot_dec * obs.detrend_inv_sum_Δt²
            @inbounds for i in eachindex(tab.epoch)
                ra_offset[i] += pra[i] - mean_ra - slope_ra * obs.detrend_Δt[i]
                dec_offset[i] += pdec[i] - mean_dec - slope_dec * obs.detrend_Δt[i]
            end
        end
    else
        accumulate_offsets!(ra_offset, dec_offset, ctx, target, reference)
    end

    if absolute
        # `ra_offset`/`dec_offset` hold the (possibly detrended) excursion,
        # gnomonic about the barycentre's apparent direction. Compose it with
        # that direction on the sphere and project about the tangent point
        # (`_dr4_project`); the parallax is the target's own at the epoch.
        basis = obs.tangent_basis
        @inbounds for i in eachindex(tab.epoch)
            sol = solutionat(ctx, i)
            ξ, η = _dr4_project(frame_ra(sol), frame_dec(sol),
                                ra_offset[i], dec_offset[i], basis)
            ra_offset[i] = ξ
            dec_offset[i] = η
            along_scan[i] = ξ * obs.sinψ[i] + η * obs.cosψ[i] +
                            _dr4_source_plx(sol, target, reference) * tab.parallax_factor_al[i]
        end
    else
        plx = ctx.θ_system.plx
        @inbounds @simd for i in eachindex(tab.epoch)
            along_scan[i] = ra_offset[i] * obs.sinψ[i] + dec_offset[i] * obs.cosψ[i] +
                            plx * tab.parallax_factor_al[i]
        end
    end
    return (; along_scan, ra_offset, dec_offset, epochs=tab.epoch)
end

function simulate(obs::GaiaDR4AstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    L = length(obs.table.epoch)
    return simulate!(Vector{T}(undef, L), Vector{T}(undef, L), Vector{T}(undef, L), obs, ctx)
end

function ln_like(obs::GaiaDR4AstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    ll = zero(T)
    L = length(obs.table.epoch)
    @no_escape ctx.buf begin
        along_scan = @alloc(T, L)
        ra_offset = @alloc(T, L)
        dec_offset = @alloc(T, L)
        simulate!(along_scan, ra_offset, dec_offset, obs, ctx)
        # `hasproperty` on the θ NamedTuple is a compile-time constant, so
        # this branch folds and each model compiles only its own loop.
        if hasproperty(ctx.θ_obs, :astrometric_jitter)
            jitter² = ctx.θ_obs.astrometric_jitter^2
            # log(2πσ²) = log2π + log(σ²), with the constant hoisted across
            # the loop and the remaining log through the branch-free kernel
            # (Float64 only — Duals fall back to Base.log inside `vlog`).
            ll -= T(obs.n_eff * log(2π) / 2)
            @inbounds @simd for i in 1:L
                σ² = jitter² + obs.σ²[i]
                resid = obs.b_al[i] - along_scan[i]
                ll -= obs.w[i] * (resid^2 / σ² + PlanetOrbits.vlog(σ²)) / 2
            end
        else
            # No jitter sampled: the whole normalization is data, computed at
            # construction, and the loop is two flops and a muladd per
            # transit. This is the production fast path for pipelines that
            # estimate per-transit noise robustly instead of fitting it.
            ll += T(obs.lognorm0)
            @inbounds @simd for i in 1:L
                resid = obs.b_al[i] - along_scan[i]
                ll -= resid^2 * obs.invσ²[i] / 2
            end
        end
    end
    return ll
end

function generate_from_params(obs::GaiaDR4AstromObs, ctx::ObsContext; add_noise)
    sim = simulate(obs, ctx)
    tab = Table(obs.table)
    al = copy(sim.along_scan)
    if add_noise
        al .+= randn.() .* tab.centroid_pos_error_al
    end
    newtab = Table(tab; centroid_pos_al=al)
    return GaiaDR4AstromObs(newtab; target=obs.target, ref=obs.ref, obs.name, obs.detrend,
        tangent_point=obs.tangent_point,
        variables=(obs.priors, obs.derived))
end

# Gaia measures one number per transit — the along-scan abscissa — so that is
# what the `:auto` correction test compares.
has_correction_impact(::Type{<:GaiaDR4AstromObs}) = true
correction_impact(obs::GaiaDR4AstromObs, a::ObsContext, b::ObsContext) =
    _simulate_impact(simulate(obs, a), simulate(obs, b), (:along_scan,),
                     _tightest(obs.table.centroid_pos_error_al))
