using FITSIO
using LinearAlgebra

# HST FGS epoch astrometry from the pickle-fgs reduction pipeline.
#
# Consumes the `PICKLEFGS-EPOCHASTROM` FITS product (pickle-fgs
# spec/09_octofitter_export.md): per-HST-orbit normal points as tangent-plane
# offsets (Δα*, Δδ) [mas, ICRS] from a reference 5-parameter solution
# (`REFSOL`, conventionally the target's Gaia DR3 row — a coordinate origin,
# not a prior), a dense 2m×2m covariance in factorized form, the F5ND
# cross-filter wedge as an explicit 2-parameter nuisance with a Gaussian
# prior, and the reference-frame offset nuisances in re-instantiable form.
#
# The reference implementation of the consumer algebra is
# `picklefgs.export.refit_from_product`; `Octofitter.fgs_gls_fit` ports it
# verbatim and must agree with the product's embedded validation numbers
# (`PLXCHK ± PLXCHKS`).

const FGS_PICKLE_FORMAT = "PICKLEFGS-EPOCHASTROM"
const FGS_PICKLE_FVERS = 4  # the format version this reader implements

"""
    FGSEpochAstromObs(filename; host, companions=(), ref=Barycentre,
                      name="FGS", marginalize=true, variables=nothing)

Likelihood over HST Fine Guidance Sensor epoch astrometry produced by the
pickle-fgs pipeline (`PICKLEFGS-EPOCHASTROM` format; see pickle-fgs
`spec/09_octofitter_export.md`).

# Source membership

    host=A, companions=(b, c, d)

`host` is the star the FGS normal points are centred on, `companions` are the
bodies whose light may blend into that source, and `ref` is what the host's
reflex is measured against — `Barycentre` by default.

Companion flux ratios come from the same three-tier lookup [`G23HObs`](@ref)
documents, in the `:FGS` band: this observation's own `fluxratio` vector,
then a system-level vector of that name, then the bodies' own `flux_FGS`
variables as `flux_FGS(companion) / flux_FGS(host)`. With none of those the
companions are dark and the source tracks the host's reflex alone, which is
the right answer for FGS planet astrometry and the case this type was
written for.

Note this is a genuine flux-weighted photocentre over the whole source, not
v8's superposition of per-companion photocentres — the two differ once more
than one companion is luminous.

# Forward model

Per epoch i (mas):

    model_x = [frame_ra(tᵢ) − RA_REF_DEGᵢ]·cos δ_ref + Δπ·PLXFAC_Xᵢ
              + Σ_planets pert_x + [Mᵢ w]_x·F5NDᵢ
    model_y = [frame_dec(tᵢ) − DEC_REF_DEGᵢ]         + Δπ·PLXFAC_Yᵢ
              + Σ_planets pert_y + [Mᵢ w]_y·F5NDᵢ

`frame_ra`/`frame_dec` are the barycentre's rigorous apparent direction, and
`RA_REF_DEG`/`DEC_REF_DEG` are the exported reference track — the pipeline's
light-time-free propagation of `REFSOL`, i.e. the coordinate origin the
exported offsets were differenced against. The mixture of conventions is
deliberate: the track defines where zero is, so re-deriving it rigorously
would move the origin rather than improve it. Solved at `JYEAR_TDB +
ROEMER_YR`, the epoch the track itself was propagated to.

Reading the frame rather than expanding about `REFSOL` is what makes the
perspective acceleration a modelled quantity instead of one frozen at
`REFSOL`'s `RV_KMS`. On a product whose epochs sit two decades from its
reference epoch that is not a refinement: for Barnard it is worth 3.05 mas per
km/s of `rv`, 0.62 mas per mas of Δπ and 0.03 mas per mas/yr of Δμ. It is also
why FVERS 4 is required rather than merely preferred — see the loader's error.

`Δπ` is the sampled parallax minus the file's `REFSOL` value, and stays
differential because the reference track already carries `REFSOL`'s own
parallax ellipse; `PLXFAC` are the exported parallax-factor columns, which
embed the producing pipeline's ephemeris conventions. Planet photocentre
perturbations use the standard Octofitter machinery, and `w = (wedge_x,
wedge_y)` is the F5ND cross-filter wedge in the ideal (detector) frame with its
per-epoch ideal→sky map `Mᵢ`.

The system must carry `ra`, `dec` [deg], `plx` [mas], `pmra`, `pmdec`
[mas/yr] variables, with `ref_epoch` equal to the REFSOL epoch (Gaia DR3:
J2016.0). The FGS epochs are absolute in the Gaia DR3 frame (anchored
through the field's reference stars) and statistically independent of the
target's own DR3 row — but see spec/09 §6 before combining with that row.

Observation variables (defaults built from the file if `variables` is not
given — pass your own to override, keeping the same names):

* `wedge_x`, `wedge_y` — the wedge nuisances, prior `Normal(WPRI, WSIG)`
  from the frozen cross-filter calibration embedded in the file. These MUST
  be fitted with this prior (never fixed to zero): corr(π, wedge_x) ≈ 0.99.
* `offsets_u` (only `marginalize=false`) — the k reference-frame offset
  nuisances, whitened: `offsets_u ~ Product(fill(Normal(0,1), k))` and
  `Δ = Δ̂ + L_Δ u` inside the likelihood.
* `locks` (when the product carries a `LOCKSTATE` block) — one free
  FineLock wrong-lock crossing offset [mas, ideal frame] per LOCKSTATE
  row, shared by the epochs with `LOCK_CLU == CLUSTER` and applied through
  the per-epoch ideal→sky map on the displaced axis only. Prior
  `Product(fill(Normal(0, LOCKPRI), n_lock))` with the file's broad
  conditioning prior — the physical slab prior arbitrated the discrete
  lock *states* in the producing pipeline (`picklefgs.lockstate`); the
  exported hard assignment had p ≈ 1 everywhere. These MUST stay free
  (never fixed): the affected epochs are only anchored through them.
* `drift` (when the product carries a `DRIFT` block) — one free
  ideal-frame linear drift rate [mas/yr] per DRIFT row, applied to EVERY
  epoch with the time weight `(JYEAR_TDB − DRFEPOCH)` through the per-epoch
  ideal→sky map on that axis. Prior `Product(fill(Normal(0, DRFPRI),
  n_drift))` with the file's broad conditioning prior. This is a per-field
  detector-frame drift identified and measured by the producing pipeline;
  the producer's fitted rates ship as diagnostics and MUST be refitted
  here (never fixed) — an unmodelled drift aliases into parallax at
  −0.42 mas per mas/yr, and a fixed point estimate would understate the
  π/PM covariance.

`marginalize=true` (default) uses the pre-marginalized covariance
`COV_MARG` — the reference-frame nuisances are integrated out analytically.
`marginalize=false` instantiates them as k extra parameters on the
conditional covariance `COV_COND`; the science posterior is identical
(linear-Gaussian; the producing pipeline locks the two modes to 1e-6 mas)
but the offset posteriors become inspectable — does the orbit fit drag any
reference star's correction? (cf. `MarginalizedRVObs` for the same pattern
in RV.)
"""
struct FGSEpochAstromObs{TTable<:Table,THost,TComp,TRef} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    name::String
    marginalize::Bool
    host::THost
    companions::TComp
    ref::TRef
    # REFSOL: the exported offsets' coordinate origin
    refsol::@NamedTuple{source_id::Int64, ra::Float64, dec::Float64, parallax::Float64,
        pmra::Float64, pmdec::Float64, rv_kms::Float64, epoch_jyear::Float64}
    ref_epoch_mjd::Float64
    # lower Cholesky factor of the active noise covariance [mas²]
    # (COV_MARG when marginalize, else COV_COND), plus its log-det half
    L_noise::Matrix{Float64}
    ldet_half::Float64
    # reference-frame offset nuisances (k = 2 × n_ref_stars; empty when
    # marginalize=true or the product carries none)
    offset_design::Matrix{Float64}   # (2m, k) B_eff [mas/mas]
    offset_Ld::Matrix{Float64}       # (k, k) lower Cholesky of Σ_Δ [mas]
    offset_mean::Vector{Float64}     # (k,) Δ̂ [mas] (already applied to the data)
    offset_source_id::Vector{Int64}
    # wedge prior (ideal frame, mas) — kept for reporting; the sampled prior
    # lives in `priors`
    wedge_prior_mean::Vector{Float64}
    wedge_prior_sigma::Vector{Float64}
    # LOCKSTATE block (docstring): one free crossing offset per
    # entry, axis ∈ (1, 2) = ideal (x, y), members = epoch indices with
    # LOCK_CLU == cluster. delta/sigma are the producer's fitted values
    # (diagnostic / initialization only). Empty when the product has none.
    lock_cols::Vector{@NamedTuple{cluster::Int, axis::Int, delta_mas::Float64, sigma_mas::Float64}}
    lock_members::Vector{Vector{Int}}
    lock_prior_mas::Float64
    # DRIFT block (docstring): one free ideal-frame linear drift
    # rate per entry, axis ∈ (1, 2) = ideal (x, y), applied to EVERY epoch
    # with time weight tdrift[i] = JYEAR_TDB[i] − DRFEPOCH. rate/sigma are
    # the producer's fitted values (diagnostic only). Empty when absent.
    drift_cols::Vector{@NamedTuple{axis::Int, rate_mas_yr::Float64, sigma_mas_yr::Float64}}
    tdrift::Vector{Float64}
    drift_prior_mas_yr::Float64
    # embedded validation gate from the producing pipeline
    meta::@NamedTuple{target::String, plx_check_mas::Float64, plx_check_sigma::Float64,
        chi2red_5p::Float64, chi2_ok::Bool}
end

likelihoodname(obs::FGSEpochAstromObs) = obs.name
refspecs(obs::FGSEpochAstromObs) = (obs.host, obs.companions..., obs.ref)
_refdesc(obs::FGSEpochAstromObs) = _blend_refdesc(obs.host, obs.companions, obs.ref)

function likeobj_from_epoch_subset(obs::FGSEpochAstromObs, obs_inds)
    error("""
    Data subsetting is not supported for FGSEpochAstromObs.

    The FGS normal points are correlated across epochs through the shared
    reference frame (dense covariance), so pointwise likelihoods are not
    independent and cross-validation subsetting is not defined for this
    observation type.
    """)
end

function FGSEpochAstromObs(
    filename::AbstractString;
    host,
    companions=(),
    ref=Barycentre,
    name::String="FGS",
    marginalize::Bool=true,
    variables::Union{Nothing,Tuple{Priors,Derived}}=nothing,
)
    hostspec = refspec(host)
    compspecs = map(refspec, Tuple(companions))
    refspec_ = refspec(ref)
    FITS(filename, "r") do f
        hdr = read_header(f[1])
        fmt = get(hdr, "FORMAT", nothing)
        fmt == FGS_PICKLE_FORMAT || error("$filename: FORMAT=$fmt, not a $FGS_PICKLE_FORMAT file")
        fvers = get(hdr, "FVERS", nothing)
        fvers == FGS_PICKLE_FVERS ||
            error("""
            $filename: FVERS=$fvers (this reader implements $FGS_PICKLE_FORMAT FVERS $FGS_PICKLE_FVERS).

            FVERS 4 added the reference track — RA_REF_DEG, DEC_REF_DEG, PLX_REF_MAS,
            ROEMER_YR — which this likelihood's forward model differences against. An
            earlier product cannot be read: without the track the only available model
            is a first-order expansion about REFSOL, which freezes the perspective
            acceleration at REFSOL's rv (worth 3 mas per km/s on a product whose epochs
            are two decades from its reference epoch).

            Regenerate with pickle-fgs `scripts/export_epoch_astrometry.py`. The columns
            are a function of REFSOL and the epoch times alone, so the reduction is
            unchanged — DX/DY and every covariance come back bit-identical.""")

        ep = f["EPOCHS"]
        col(nm) = read(ep, nm)
        epoch_utc = Vector{Float64}(col("TIME_MJD_UTC"))
        jyear_tdb = Vector{Float64}(col("JYEAR_TDB"))
        m = length(epoch_utc)

        rs_hdu = f["REFSOL"]
        refsol = (;
            source_id=Int64(read(rs_hdu, "SOURCE_ID")[1]),
            ra=Float64(read(rs_hdu, "RA_DEG")[1]),
            dec=Float64(read(rs_hdu, "DEC_DEG")[1]),
            parallax=Float64(read(rs_hdu, "PARALLAX_MAS")[1]),
            pmra=Float64(read(rs_hdu, "PMRA_MASYR")[1]),
            pmdec=Float64(read(rs_hdu, "PMDEC_MASYR")[1]),
            rv_kms=Float64(read(rs_hdu, "RV_KMS")[1]),
            epoch_jyear=Float64(read(rs_hdu, "EPOCH_JYEAR")[1]),
        )
        ref_epoch_mjd = 51544.5 + (refsol.epoch_jyear - 2000.0) * julian_year

        # The reference track (spec/09 §2): the BARYCENTRIC half of the origin
        # the exported offsets were differenced against.
        ra_ref = Vector{Float64}(col("RA_REF_DEG"))
        dec_ref = Vector{Float64}(col("DEC_REF_DEG"))
        roemer_yr = Vector{Float64}(col("ROEMER_YR"))

        # The epoch the trajectory is solved at, which is NOT the epoch the
        # product reports. The model differences against a track the pipeline
        # propagated to `JYEAR_TDB + ROEMER_YR`, so it has to be evaluated at
        # that same instant or the two sides are read at different times. Both
        # terms are small and both are real: the UTC/TDB gap that `TIME_MJD_UTC`
        # would carry is 0.023 mas at Barnard's µ, and the Roemer shift (≤ 499 s)
        # is 0.17 mas. `epoch_utc` is kept for display.
        epoch = @. 51544.5 + (jyear_tdb + roemer_yr - 2000.0) * julian_year

        table = Table(;
            epoch,                                  # MJD, TDB + Roemer: the solve epoch
            epoch_utc,                              # MJD UTC, as exported
            tau=jyear_tdb .- refsol.epoch_jyear,    # Julian years from REFSOL epoch (TDB)
            ra_ref, dec_ref, roemer_yr,             # reference track [deg, deg, jyear]
            dx=Vector{Float64}(col("DX_MAS")),      # mas, tangent plane vs REFSOL
            dy=Vector{Float64}(col("DY_MAS")),
            σ_x=Vector{Float64}(col("SIG_X_MAS")),  # √diag(COV_MARG); display only
            σ_y=Vector{Float64}(col("SIG_Y_MAS")),
            plxfac_x=Vector{Float64}(col("PLXFAC_X")),
            plxfac_y=Vector{Float64}(col("PLXFAC_Y")),
            M11=Vector{Float64}(col("M11")),        # ideal→sky per epoch
            M12=Vector{Float64}(col("M12")),
            M21=Vector{Float64}(col("M21")),
            M22=Vector{Float64}(col("M22")),
            is_f5nd=Vector{Bool}(col("IS_F5ND")),
            roll_v3_deg=Vector{Float64}(col("ROLL_V3_DEG")),
            orbit_id=Vector{String}(col("ORBIT_ID")),
            n_exp=Vector{Int}(col("N_EXP")),
        )

        # FITS images written from numpy C-order arrays come back transposed
        # in FITSIO's column-major read (NAXIS1 = numpy's last axis) —
        # immaterial for the symmetric covariances, essential for the
        # (2m, k) design matrix.
        read_image(nm) = Matrix{Float64}(transpose(read(f[nm])))

        cov_marg = read_image("COV_MARG")
        @assert size(cov_marg) == (2m, 2m)
        has_offsets = get(hdr, "HASOFFS", false) == true

        if marginalize || !has_offsets
            noise_cov = cov_marg
            offset_design = zeros(2m, 0)
            offset_Ld = zeros(0, 0)
            offset_mean = zeros(0)
            offset_source_id = Int64[]
        else
            noise_cov = read_image("COV_COND")
            offset_design = read_image("OFFSET_DESIGN")
            offset_cov = read_image("OFFSET_COV")
            op = f["OFFSET_POST"]
            offset_source_id = Vector{Int64}(read(op, "SOURCE_ID"))
            k = 2 * length(offset_source_id)
            @assert size(offset_design) == (2m, k) && size(offset_cov) == (k, k)
            offset_mean = Vector{Float64}(undef, k)
            offset_mean[1:2:end] .= read(op, "DRA_MAS")
            offset_mean[2:2:end] .= read(op, "DDEC_MAS")
            offset_Ld = Matrix(cholesky(Symmetric(offset_cov)).L)
        end
        k = size(offset_design, 2)

        L_noise = Matrix(cholesky(Symmetric(noise_cov)).L)
        ldet_half = sum(log, diag(L_noise))

        whdr = read_header(f["WEDGE"])
        wedge_prior_mean = [Float64(whdr["WPRI_X"]), Float64(whdr["WPRI_Y"])]
        wedge_prior_sigma = [Float64(whdr["WSIG_X"]), Float64(whdr["WSIG_Y"])]

        # LOCKSTATE block (docstring)
        lock_cols = @NamedTuple{cluster::Int, axis::Int, delta_mas::Float64, sigma_mas::Float64}[]
        lock_members = Vector{Int}[]
        lock_prior_mas = Inf
        if get(hdr, "HASLOCK", false) == true && get(hdr, "NLOCKC", 0) > 0
            lock_cluster = Vector{Int}(read(ep, "LOCK_CLU"))
            ls = f["LOCKSTATE"]
            lock_prior_mas = Float64(read_header(ls)["LOCKPRI"])
            clusters = Vector{Int}(read(ls, "CLUSTER"))
            axes = Vector{String}(read(ls, "AXIS"))
            deltas = Vector{Float64}(read(ls, "DELTA_MAS"))
            sigmas = Vector{Float64}(read(ls, "SIGMA_MAS"))
            for j in eachindex(clusters)
                ax = axes[j] == "x" ? 1 : axes[j] == "y" ? 2 :
                    error("$filename: LOCKSTATE AXIS '$(axes[j])' not in (x, y)")
                push!(lock_cols, (; cluster=clusters[j], axis=ax,
                    delta_mas=deltas[j], sigma_mas=sigmas[j]))
                push!(lock_members, findall(==(clusters[j]), lock_cluster))
            end
        end
        n_lock = length(lock_cols)

        # DRIFT block (docstring)
        drift_cols = @NamedTuple{axis::Int, rate_mas_yr::Float64, sigma_mas_yr::Float64}[]
        tdrift = zeros(m)
        drift_prior_mas_yr = Inf
        if get(hdr, "HASDRIFT", false) == true
            dr = f["DRIFT"]
            dhdr = read_header(dr)
            drift_epoch = Float64(dhdr["DRFEPOCH"])
            drift_prior_mas_yr = Float64(dhdr["DRFPRI"])
            tdrift = jyear_tdb .- drift_epoch
            daxes = Vector{String}(read(dr, "AXIS"))
            drates = Vector{Float64}(read(dr, "RATE_MASYR"))
            dsigmas = Vector{Float64}(read(dr, "SIGMA_MASYR"))
            for j in eachindex(daxes)
                ax = daxes[j] == "x" ? 1 : daxes[j] == "y" ? 2 :
                    error("$filename: DRIFT AXIS '$(daxes[j])' not in (x, y)")
                push!(drift_cols, (; axis=ax, rate_mas_yr=drates[j], sigma_mas_yr=dsigmas[j]))
            end
        end
        n_drift = length(drift_cols)

        meta = (;
            target=String(hdr["TARGET"]),
            plx_check_mas=Float64(hdr["PLXCHK"]),
            plx_check_sigma=Float64(hdr["PLXCHKS"]),
            chi2red_5p=Float64(hdr["CHI2R5P"]),
            chi2_ok=hdr["CHI2OK"] == true,
        )
        if !meta.chi2_ok
            @warn "This FGS product shipped with its per-field validation gate FLAGGED (χ²red of the embedded 5-parameter check outside [0.5, 2.0]). The covariance may not price all structure on this field; interpret posteriors with care. See pickle-fgs spec/09 §6." meta.chi2red_5p
        end

        if isnothing(variables)
            priors_dict = OrderedDict{Symbol,Distribution}(
                :wedge_x => Normal(wedge_prior_mean[1], wedge_prior_sigma[1]),
                :wedge_y => Normal(wedge_prior_mean[2], wedge_prior_sigma[2]),
            )
            if !marginalize && k > 0
                priors_dict[:offsets_u] = Product(fill(Normal(0.0, 1.0), k))
            end
            if n_lock > 0
                priors_dict[:locks] = Product(fill(Normal(0.0, lock_prior_mas), n_lock))
            end
            if n_drift > 0
                priors_dict[:drift] = Product(fill(Normal(0.0, drift_prior_mas_yr), n_drift))
            end
            priors = Priors(priors_dict)
            derived = Derived(OrderedDict{Symbol,Any}(), (), ())
        else
            (priors, derived) = variables
            for key in (:wedge_x, :wedge_y)
                haskey(priors.priors, key) || haskey(derived.variables, key) ||
                    error("observation variables must define `$key` (the F5ND cross-filter wedge nuisance; the file's prior is Normal($(wedge_prior_mean), $(wedge_prior_sigma)))")
            end
            if !marginalize && k > 0 &&
               !(haskey(priors.priors, :offsets_u) || haskey(derived.variables, :offsets_u))
                error("marginalize=false requires the `offsets_u` observation variable (whitened reference-frame offsets, e.g. `offsets_u ~ Product(fill(Normal(0,1), $k))`)")
            end
            if n_lock > 0 &&
               !(haskey(priors.priors, :locks) || haskey(derived.variables, :locks))
                error("this product carries a lock-state block: the observation variables must define `locks` (e.g. `locks ~ Product(fill(Normal(0, $lock_prior_mas), $n_lock))`; see the FGSEpochAstromObs docstring)")
            end
            if n_drift > 0 &&
               !(haskey(priors.priors, :drift) || haskey(derived.variables, :drift))
                error("this product carries a drift block: the observation variables must define `drift` (e.g. `drift ~ Product(fill(Normal(0, $drift_prior_mas_yr), $n_drift))`; see the FGSEpochAstromObs docstring)")
            end
        end

        obs = FGSEpochAstromObs{typeof(table),typeof(hostspec),typeof(compspecs),
            typeof(refspec_)}(
            table, priors, derived, name, marginalize || !has_offsets,
            hostspec, compspecs, refspec_,
            refsol, ref_epoch_mjd,
            L_noise, ldet_half,
            offset_design, offset_Ld, offset_mean, offset_source_id,
            wedge_prior_mean, wedge_prior_sigma,
            lock_cols, lock_members, lock_prior_mas,
            drift_cols, tdrift, drift_prior_mas_yr,
            meta,
        )

        # Self-check: the file-only GLS fit must reproduce the embedded
        # validation parallax (the producing pipeline guarantees 0.002 mas
        # closure; we allow 0.05 for numeric slack).
        gls = fgs_gls_fit(obs)
        dplx = gls.parallax_mas - meta.plx_check_mas
        if abs(dplx) > 0.05
            @warn "FGS product self-check: file-only GLS parallax differs from the embedded PLXCHK — loader/format mismatch?" gls.parallax_mas meta.plx_check_mas dplx
        else
            @info "FGS product loaded: $(meta.target), $m epochs, k=$k offset nuisances $(marginalize ? "(marginalized)" : "(explicit)"). GLS self-check π = $(round(gls.parallax_mas, digits=3)) ± $(round(gls.sigma[3], digits=3)) mas (file PLXCHK $(meta.plx_check_mas) ± $(meta.plx_check_sigma); Δ = $(round(dplx, sigdigits=2)) mas)."
        end

        return obs
    end
end
export FGSEpochAstromObs

# --- forward model ----------------------------------------------------------------

# --- the frame path ----------------------------------------------------------------
#
# Rather than expand `Δα*₀ + Δμ_α* τ` about REFSOL, this reads straight off the
# frame PlanetOrbits already propagated:
#
#     model = [frame_ra, frame_dec](t) − [RA_REF_DEG, DEC_REF_DEG](t)
#
# `frame_ra`/`frame_dec` are the system barycentre's rigorous **apparent**
# direction, so perspective acceleration, the frame's own light-time terms and
# every cross-derivative the Taylor expansion froze are modelled rather than
# approximated. `rv` becomes a real lever: 3.05 mas per km/s at Barnard.
#
# The subtraction is deliberately MIXED-CONVENTION and must stay that way. The
# exported track is the pipeline's light-time-free B&L14 propagation of REFSOL,
# and it is what `DX/DY` were differenced against — a coordinate origin, not a
# claim about the path. Re-deriving it rigorously here would silently move the
# origin. The right form is `rigorous_apparent(t) − REFSOL_lighttime_free(t)`,
# which is what this is.
#
# Two conventions come along with the track and are handled at load time: the
# solve epoch is TDB rather than the product's native UTC, and it carries the
# exported Roemer shift, because the track was propagated to
# `JYEAR_TDB + ROEMER_YR`.
#
# Precision note: `frame_ra` and `ra_ref` are both ≈ 269.4° for Barnard and
# differ by a few mas. The subtraction itself is exact (Sterbenz); the error is
# one ulp of the larger, ≈ 2e-4 mas from `frame_ra`'s own `atan`, four orders
# below a per-epoch σ. Under `Dual` the partials do not cancel at all.
@inline _fgs_wrap180(Δ) = Δ - 360 * round(Δ / 360)

function _fgs_frame!(model_x, model_y, obs::FGSEpochAstromObs,
                     ctx::ObsContext, Δplx)
    tbl = obs.table
    @inbounds for i in eachindex(tbl.epoch)
        sol = solutionat(ctx, i)
        dec_ref = tbl.dec_ref[i]
        model_x[i] += _fgs_wrap180(_fgs_frame_ra(sol) - tbl.ra_ref[i]) *
                      3.6e6 * cosd(dec_ref) + Δplx * tbl.plxfac_x[i]
        model_y[i] += (_fgs_frame_dec(sol) - dec_ref) * 3.6e6 +
                      Δplx * tbl.plxfac_y[i]
    end
    return nothing
end

# Reading a direction needs an absolute frame. FGS always has one in practice —
# the docstring requires ra/dec/plx/pmra/pmdec — but say so here rather than
# surfacing a `MethodError` from inside the epoch loop.
@inline _fgs_frame_ra(sol::PlanetOrbits._AbsSol) = frame_ra(sol)
@inline _fgs_frame_dec(sol::PlanetOrbits._AbsSol) = frame_dec(sol)
@noinline _fgs_frame_ra(@nospecialize sol) = _err_fgs_needs_absolute()
@noinline _fgs_frame_dec(@nospecialize sol) = _err_fgs_needs_absolute()
@noinline _err_fgs_needs_absolute() = error(
    "FGSEpochAstromObs reads the system's apparent sky path directly, which " *
    "needs an absolute frame: give the system " *
    "`ra`, `dec`, `plx`, `pmra`, `pmdec` (and `rv`) variables with " *
    "`ref_epoch` at the product's REFSOL epoch.")

# --- correction flags --------------------------------------------------------------
#
# What a correction does to this observation: the largest per-epoch shift in the
# predicted offset, against the tightest single-epoch σ in the product.
#
# The reference track is the pipeline's LIGHT-TIME-FREE propagation, so
# differencing a rigorously propagated path against it isolates exactly the
# barycentric light-time term — 0.082 mas over Barnard's 23-year lever arm,
# against a 1.6 mas tightest per-epoch σ. That makes this a real measurement
# rather than the "cannot say" that would default the flag on.
has_correction_impact(::Type{<:FGSEpochAstromObs}) = true
function correction_impact(obs::FGSEpochAstromObs, a::ObsContext, b::ObsContext)
    σ = min(minimum(obs.table.σ_x), minimum(obs.table.σ_y))
    sa = simulate(obs, a)
    sb = simulate(obs, b)
    d = 0.0
    for (xa, xb) in ((sa.model_x, sb.model_x), (sa.model_y, sb.model_y))
        for i in eachindex(xa)
            δ = abs(float(xa[i] - xb[i]))
            isfinite(δ) || return (; delta=NaN, sigma=σ, n=0)
            d = max(d, δ)
        end
    end
    return (; delta=d, sigma=σ, n=2 * length(obs.table.epoch))
end

# Fill `model_x`, `model_y` (length m, mas) with the model offsets vs REFSOL.
# Shared by ln_like and simulate. Buffers must be zeroed by the caller before
# entry — the source's sky excursion is accumulated into them.
function _fgs_model!(
    model_x::AbstractVector, model_y::AbstractVector,
    obs::FGSEpochAstromObs, ctx::ObsContext, ::Type{T},
) where {T}
    tbl = obs.table
    rs = obs.refsol
    θ_system = ctx.θ_system
    θ_obs = ctx.θ_obs
    sys = ctx.system

    # The source's own excursion: one flux-weighted photocentre over the host
    # and its companions, one `raoff`/`decoff` per epoch. Resolve every
    # reference once, outside the epoch loop.
    hostref = ref(ctx, obs.host)
    reference = ref(ctx, obs.ref)
    comps = resolverefs(ctx, obs.companions)
    cidx = map(c -> c.idx, comps)
    masses = sys.masses
    # Test the primal: a differentiated zero mass is a Dual whose value is
    # zero but whose partials are not, so `iszero` on it is false.
    active = map(c -> !iszero(PlanetOrbits._primal(masses[c.idx])), comps)
    if any(active)
        f = _g23h_fluxratios(θ_obs, θ_system, :fluxratio, Val(:FGS), sys,
            hostref.idx, cidx, active, T)
        photo = _g23h_photocentre(sys, hostref.idx, cidx, f)
        accumulate_offsets!(model_x, model_y, ctx, photo, reference)
    end

    # The parallax stays differential. `xi_ref` already carries REFSOL's full
    # parallactic displacement and `PLXFAC` is the pipeline's own
    # d(sky)/d(π at the REFSOL epoch) through the whole chain, so `Δπ · PLXFAC`
    # is exactly the residual ellipse — and both Δπ and PLXFAC are referred to
    # the same epoch, which is what makes that true.
    Δplx = θ_system.plx - rs.parallax
    wx = θ_obs.wedge_x
    wy = θ_obs.wedge_y

    _fgs_frame!(model_x, model_y, obs, ctx, Δplx)

    @inbounds for i in eachindex(tbl.epoch)
        if tbl.is_f5nd[i]
            model_x[i] += tbl.M11[i] * wx + tbl.M12[i] * wy
            model_y[i] += tbl.M21[i] * wx + tbl.M22[i] * wy
        end
    end

    # lock-state crossings: shared ideal-frame offset per LOCKSTATE
    # row on its member epochs, displaced axis only (docstring)
    if !isempty(obs.lock_cols)
        locks = θ_obs.locks
        @inbounds for j in eachindex(obs.lock_cols)
            δ = locks[j]
            ax = obs.lock_cols[j].axis
            for i in obs.lock_members[j]
                if ax == 1
                    model_x[i] += tbl.M11[i] * δ
                    model_y[i] += tbl.M21[i] * δ
                else
                    model_x[i] += tbl.M12[i] * δ
                    model_y[i] += tbl.M22[i] * δ
                end
            end
        end
    end

    # drift: free ideal-frame linear rate per DRIFT row, every
    # epoch, time weight tdrift = JYEAR_TDB − DRFEPOCH (docstring)
    if !isempty(obs.drift_cols)
        drift = θ_obs.drift
        @inbounds for j in eachindex(obs.drift_cols)
            d = drift[j]
            ax = obs.drift_cols[j].axis
            for i in eachindex(tbl.epoch)
                w = d * obs.tdrift[i]
                if ax == 1
                    model_x[i] += tbl.M11[i] * w
                    model_y[i] += tbl.M21[i] * w
                else
                    model_x[i] += tbl.M12[i] * w
                    model_y[i] += tbl.M22[i] * w
                end
            end
        end
    end
    return nothing
end

function ln_like(obs::FGSEpochAstromObs, ctx::ObsContext)
    θ_system = ctx.θ_system
    θ_obs = ctx.θ_obs
    T = _system_number_type(θ_system)
    tbl = obs.table
    m = length(tbl.epoch)
    n = 2m
    k = size(obs.offset_design, 2)

    # The Δ 5-parameter linearization is taken about REFSOL at its epoch; the
    # sampled (ra, dec, pmra, pmdec) must be defined at the same epoch.
    if hasproperty(θ_system, :ref_epoch) && abs(θ_system.ref_epoch - obs.ref_epoch_mjd) > 1.0
        error("FGSEpochAstromObs: system ref_epoch ($(θ_system.ref_epoch) MJD) must equal the product's REFSOL epoch ($(obs.ref_epoch_mjd) MJD = J$(obs.refsol.epoch_jyear)).")
    end

    ll = zero(T)
    @no_escape ctx.buf begin
        model_x = @alloc(T, m)
        model_y = @alloc(T, m)
        fill!(model_x, zero(T))
        fill!(model_y, zero(T))
        _fgs_model!(model_x, model_y, obs, ctx, T)

        resid = @alloc(T, n)
        @inbounds for i in 1:m
            resid[2i-1] = tbl.dx[i] - model_x[i]
            resid[2i] = tbl.dy[i] - model_y[i]
        end

        # explicit reference-frame nuisances: Δ = Δ̂ + L_Δ u, data already at Δ̂
        if !obs.marginalize && k > 0
            u = θ_obs.offsets_u
            q = @alloc(T, k)
            @inbounds for a in 1:k
                acc = zero(T)
                for b in 1:a  # L_Δ lower triangular
                    acc += obs.offset_Ld[a, b] * u[b]
                end
                q[a] = acc
            end
            D = obs.offset_design
            @inbounds for r in 1:n
                acc = zero(T)
                for a in 1:k
                    acc += D[r, a] * q[a]
                end
                resid[r] -= acc
            end
        end

        # whiten against the (constant) noise Cholesky: z = L⁻¹ resid
        L = obs.L_noise
        z = @alloc(T, n)
        @inbounds for i in 1:n
            acc = resid[i]
            for j in 1:(i-1)
                acc -= L[i, j] * z[j]
            end
            z[i] = acc / L[i, i]
        end
        zz = zero(T)
        @inbounds for i in 1:n
            zz += z[i]^2
        end
        ll = -zz / 2 - obs.ldet_half - n * T(log(2π)) / 2
    end
    return ll
end

# Model offsets at the data epochs for a given parameter draw (plotting /
# diagnostics; standard Octofitter simulate signature).
function simulate(obs::FGSEpochAstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    m = length(obs.table.epoch)
    model_x = zeros(T, m)
    model_y = zeros(T, m)
    _fgs_model!(model_x, model_y, obs, ctx, T)
    # `epoch_utc`, not `epoch`: the latter is the retarded
    # TDB instant the trajectory is solved at, which is right for the physics and
    # wrong for an axis — a plot's x should be the epoch the product reports.
    return (; model_x, model_y, epoch=obs.table.epoch_utc,
        resid_x=obs.table.dx .- model_x, resid_y=obs.table.dy .- model_y)
end

# --- file-only GLS reference fit --------------------------------------------------

"""
    fgs_gls_fit(obs::FGSEpochAstromObs) -> NamedTuple

Linear generalized-least-squares astrometric fit (5 parameters + the 2 wedge
nuisances under their Gaussian prior, plus the k reference-frame offsets when
`marginalize=false`) from the product contents only, with **no companion** —
a verbatim port of `picklefgs.export.refit_from_product`, the reference
consumer algebra of pickle-fgs spec/09 §3.

This is the acceptance check for the loader: `parallax_mas ± sigma[3]` must
reproduce the file's embedded `PLXCHK ± PLXCHKS` (the producing pipeline
guarantees 0.002 mas closure against its rigorous fit), identically in both
`marginalize` modes. It is run automatically at load time.

Returns absolute (REFSOL + fitted offset) astrometry: `ra_deg`, `dec_deg`,
`parallax_mas`, `pmra_mas_yr`, `pmdec_mas_yr`, `wedge_mas`, `sigma` (the 7
marginal σs on the mas grid), `cov` (7×7), `chi2red`, `offsets_mas`.
"""
function fgs_gls_fit(obs::FGSEpochAstromObs)
    tbl = obs.table
    m = length(tbl.epoch)
    n = 2m
    n_lock = length(obs.lock_cols)
    n_drift = length(obs.drift_cols)
    n_ast = 7 + n_lock + n_drift  # [5 astrometric | 2 wedge | locks | drift]
    k = obs.marginalize ? 0 : size(obs.offset_design, 2)
    n_par = n_ast + k

    y = Vector{Float64}(undef, n)
    A = zeros(n, n_par)
    for i in 1:m
        r = 2i - 1
        y[r] = tbl.dx[i]
        y[r+1] = tbl.dy[i]
        A[r, 1] = 1.0
        A[r+1, 2] = 1.0
        A[r, 3] = tbl.plxfac_x[i]
        A[r+1, 3] = tbl.plxfac_y[i]
        A[r, 4] = tbl.tau[i]
        A[r+1, 5] = tbl.tau[i]
        if tbl.is_f5nd[i]
            A[r, 6] = tbl.M11[i]
            A[r, 7] = tbl.M12[i]
            A[r+1, 6] = tbl.M21[i]
            A[r+1, 7] = tbl.M22[i]
        end
    end
    for j in 1:n_lock
        ax = obs.lock_cols[j].axis
        for i in obs.lock_members[j]
            A[2i-1, 7+j] = ax == 1 ? tbl.M11[i] : tbl.M12[i]
            A[2i, 7+j] = ax == 1 ? tbl.M21[i] : tbl.M22[i]
        end
    end
    for j in 1:n_drift
        ax = obs.drift_cols[j].axis
        for i in 1:m
            td = obs.tdrift[i]
            A[2i-1, 7+n_lock+j] = (ax == 1 ? tbl.M11[i] : tbl.M12[i]) * td
            A[2i, 7+n_lock+j] = (ax == 1 ? tbl.M21[i] : tbl.M22[i]) * td
        end
    end
    if k > 0
        A[:, n_ast+1:end] .= obs.offset_design
    end

    wμ, wσ = obs.wedge_prior_mean, obs.wedge_prior_sigma
    n_prior = 2 + n_lock + n_drift + k
    prior_rows = zeros(n_prior, n_par)
    prior_res = zeros(n_prior)
    prior_rows[1, 6] = 1 / wσ[1]
    prior_rows[2, 7] = 1 / wσ[2]
    prior_res[1] = wμ[1] / wσ[1]
    prior_res[2] = wμ[2] / wσ[2]
    for j in 1:n_lock
        prior_rows[2+j, 7+j] = 1 / obs.lock_prior_mas
    end
    for j in 1:n_drift
        prior_rows[2+n_lock+j, 7+n_lock+j] = 1 / obs.drift_prior_mas_yr
    end
    if k > 0
        # offsets parameterized as (Δ − Δ̂): zero-mean prior with covariance Σ_Δ
        prior_rows[2+n_lock+n_drift+1:end, n_ast+1:end] .= inv(LowerTriangular(obs.offset_Ld))
    end

    Lt = LowerTriangular(obs.L_noise)
    J = vcat(Lt \ A, prior_rows)
    R = vcat(Lt \ y, prior_res)
    F = qr(J)
    p = F \ R
    Rinv = inv(UpperTriangular(F.R))
    cov = Rinv * Rinv'
    resid = R - J * p
    dof = length(R) - n_par
    chi2red = dof > 0 ? sum(abs2, resid) / dof : NaN

    rs = obs.refsol
    return (;
        ra_deg=rs.ra + p[1] / 3.6e6 / cosd(rs.dec),
        dec_deg=rs.dec + p[2] / 3.6e6,
        parallax_mas=rs.parallax + p[3],
        pmra_mas_yr=rs.pmra + p[4],
        pmdec_mas_yr=rs.pmdec + p[5],
        wedge_mas=p[6:7],
        sigma=sqrt.(diag(cov)[1:7]),
        cov=cov[1:7, 1:7],
        params_offset_grid=p[1:7],
        lock_mas=n_lock > 0 ? p[8:7+n_lock] : nothing,
        lock_sigma_mas=n_lock > 0 ? sqrt.(diag(cov)[8:7+n_lock]) : nothing,
        drift_mas_yr=n_drift > 0 ? p[7+n_lock+1:n_ast] : nothing,
        drift_sigma_mas_yr=n_drift > 0 ? sqrt.(diag(cov)[7+n_lock+1:n_ast]) : nothing,
        offsets_mas=k > 0 ? obs.offset_mean .+ p[n_ast+1:end] : nothing,
        chi2red,
        epoch_jyear=rs.epoch_jyear,
    )
end
export fgs_gls_fit
