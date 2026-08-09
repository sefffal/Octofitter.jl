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
# ---------------------------------------------------

"""
    GaiaDR4AstromObs(data; target=Photocentre, ref=Barycentre, name="GaiaDR4",
                     detrend=false, variables=…)

Gaia DR4 individual (epoch) astrometry: one along-scan abscissa per
transit.

`data` needs `:epoch` [MJD], `:scan_pos_angle` ψ [rad],
`:parallax_factor_al`, `:centroid_pos_al` [mas] and
`:centroid_pos_error_al` [mas]; an `:outlier_flag` column is honoured if
present.

The model for each transit is

    η = Δα* sin ψ + Δδ cos ψ + ϖ · f_al

where `Δα*, Δδ` are the barycentre's own linear motion (from this
observation's `ra_offset_mas`, `dec_offset_mas`, `pmra`, `pmdec`,
`ref_epoch` variables) plus the offset of `target` from `ref` — normally the
photocentre relative to the system barycentre.

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

# `detrend`
With `detrend=true` the linear (constant + slope) part of the
target-versus-reference excursion is removed before it enters the model, so
only the curvature does. That breaks the degeneracy between the fitted
proper motion and a wide companion, and makes the observation's position and
proper-motion variables describe the *photocentre* rather than the
barycentre.

# Variables
`astrometric_jitter` [mas] adds in quadrature to the per-transit formal
error; `ra_offset_mas`, `dec_offset_mas`, `pmra`, `pmdec`, `ref_epoch`
define the reference-point motion. The parallax comes from the system block.

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
end

const dr4_cols = (:epoch, :scan_pos_angle, :parallax_factor_al,
                  :centroid_pos_al, :centroid_pos_error_al)

function GaiaDR4AstromObs(observations;
                          target=Photocentre, ref=Barycentre,
                          name="GaiaDR4",
                          detrend::Bool=false,
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

    sinψ = collect(Float64, sin.(table.scan_pos_angle))
    cosψ = collect(Float64, cos.(table.scan_pos_angle))
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
    return GaiaDR4AstromObs{typeof(table),typeof(t),typeof(r)}(
        table, priors, derived, t, r, String(name), detrend,
        Δt, 1.0 / length(ep), 1.0 / sum(Δt .^ 2),
        sinψ, cosψ, w, b_al, σ², invσ², sum(w), lognorm0)
end

export GaiaDR4AstromObs

refspecs(obs::GaiaDR4AstromObs) = (obs.target, obs.ref)

function likeobj_from_epoch_subset(obs::GaiaDR4AstromObs, inds)
    return GaiaDR4AstromObs(obs.table[inds, :, 1];
        target=obs.target, ref=obs.ref, obs.name, obs.detrend,
        variables=(obs.priors, obs.derived))
end

"""
    simulate!(along_scan, ra_offset, dec_offset, obs, ctx)

Fill the modelled along-scan abscissae (and the RA/Dec offsets they were
projected from) for every transit.
"""
function simulate!(along_scan, ra_offset, dec_offset, obs::GaiaDR4AstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    tab = obs.table

    ref_epoch = hasproperty(θ_obs, :ref_epoch) ? θ_obs.ref_epoch : zero(T)
    ra0 = hasproperty(θ_obs, :ra_offset_mas) ? θ_obs.ra_offset_mas : zero(T)
    dec0 = hasproperty(θ_obs, :dec_offset_mas) ? θ_obs.dec_offset_mas : zero(T)
    pmra_ = hasproperty(θ_obs, :pmra) ? θ_obs.pmra : zero(T)
    pmdec_ = hasproperty(θ_obs, :pmdec) ? θ_obs.pmdec : zero(T)

    # The reference point's own linear motion.
    @inbounds for i in eachindex(tab.epoch)
        Δt = (tab.epoch[i] - ref_epoch) / PlanetOrbits.year2day_julian
        ra_offset[i] = ra0 + pmra_ * Δt
        dec_offset[i] = dec0 + pmdec_ * Δt
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

    plx = ctx.θ_system.plx
    @inbounds @simd for i in eachindex(tab.epoch)
        along_scan[i] = ra_offset[i] * obs.sinψ[i] + dec_offset[i] * obs.cosψ[i] +
                        plx * tab.parallax_factor_al[i]
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
        variables=(obs.priors, obs.derived))
end

# Gaia measures one number per transit — the along-scan abscissa — so that is
# what the `:auto` correction test compares.
has_correction_impact(::Type{<:GaiaDR4AstromObs}) = true
correction_impact(obs::GaiaDR4AstromObs, a::ObsContext, b::ObsContext) =
    _simulate_impact(simulate(obs, a), simulate(obs, b), (:along_scan,),
                     _tightest(obs.table.centroid_pos_error_al))
