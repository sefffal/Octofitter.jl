# ---------------------------------------------------
# Shared sky-path helpers for catalog astrometry
#
# Two things every catalog-level astrometric likelihood needs, and which no
# single likelihood owns:
#
#   1. The five-parameter (and four-parameter) along-scan least-squares
#      refit — the operation a catalog pipeline performed on its own
#      abscissae, and which a model must therefore reproduce before it can
#      compare against published catalog parameters. Pure linear algebra.
#
#   2. A per-epoch offset accumulator: given a target reference and a
#      reference reference, add their sky offsets over a range of an
#      observation's table rows. This is the refs-based replacement for v1's
#      `_simulate_skypath_perturbations!`, whose per-companion
#      `(host + f·planet)/(1+f)` linear coefficient only worked because
#      every companion orbited the same star. A `WeightedPoint` — built by
#      `photocentre(sys, …)` for structurally-fixed membership, or by the
#      likelihood itself for per-draw/per-epoch membership — expresses the
#      same physics exactly and for any topology.
#
# Anything whose blending behaviour is *instrument-specific* (a grating
# response, a resolution taper, a scan-angle-dependent window) does not
# belong here; it belongs to the observation that owns the instrument.
#
# The along-scan convention throughout is the one the Hipparcos/Gaia
# intermediate-data reductions use,
#
#     b_i = Δα*_i · cosϕ_i + Δδ_i · sinϕ_i
#
# with ϕ the scan angle measured so that `cosϕ` multiplies the RA component.
# `GaiaDR4AstromObs` stores ψ = π/2 − ϕ in `scan_pos_angle` and writes the
# same projection as `Δα* sinψ + Δδ cosψ`; pass `cosϕ = sin.(ψ)`,
# `sinϕ = cos.(ψ)` to reuse these routines on such a table.
# ---------------------------------------------------

"""
    prepare_A_5param(cosϕ, sinϕ, epoch, parallax_factor_al,
                     ref_epoch_ra, ref_epoch_dec) -> Matrix{Float64}

Design matrix of the standard five-parameter astrometric solution, one row
per scan, in column order `(Δα*, Δδ, ϖ, μα*, μδ)`.

The parallax column carries the **negative** of the along-scan parallax
factor, matching the sign convention of the Hipparcos/Gaia intermediate
data. Reference epochs are MJD and may differ between RA and Dec (they do in
Hipparcos).

`A` depends only on the scan geometry, which is fixed at construction — so
build it once per observation, not once per sample.
"""
function prepare_A_5param(cosϕ, sinϕ, epoch, parallax_factor_al,
                          ref_epoch_ra, ref_epoch_dec)
    n_obs = length(epoch)
    A = zeros(n_obs, 5)
    @inbounds for i in 1:n_obs
        A[i, 1] = cosϕ[i]                                              # α*
        A[i, 2] = sinϕ[i]                                              # δ
        A[i, 3] = -parallax_factor_al[i]                               # ϖ
        A[i, 4] = cosϕ[i] * (epoch[i] - ref_epoch_ra) / julian_year    # μα*
        A[i, 5] = sinϕ[i] * (epoch[i] - ref_epoch_dec) / julian_year   # μδ
    end
    return A
end

"""
    prepare_A_4param(cosϕ, sinϕ, epoch, ref_epoch_ra, ref_epoch_dec) -> Matrix{Float64}

As [`prepare_A_5param`](@ref) without the parallax column: column order
`(Δα*, Δδ, μα*, μδ)`. For fits in which the parallax is held fixed.
"""
function prepare_A_4param(cosϕ, sinϕ, epoch, ref_epoch_ra, ref_epoch_dec)
    n_obs = length(epoch)
    A = zeros(n_obs, 4)
    @inbounds for i in 1:n_obs
        A[i, 1] = cosϕ[i]
        A[i, 2] = sinϕ[i]
        A[i, 3] = cosϕ[i] * (epoch[i] - ref_epoch_ra) / julian_year
        A[i, 4] = sinϕ[i] * (epoch[i] - ref_epoch_dec) / julian_year
    end
    return A
end

"""
    prepare_pinv_5param(A_prepared, σ) -> Matrix{Float64}

Cached weighted pseudo-inverse `Q = pinv(A ./ σ) ./ σ'` (5×N), so that a
per-epoch-weighted five-parameter fit becomes the single matrix-vector
product `x = Q · b` on the *unweighted* right-hand side. Both `A` and `σ`
are fixed at construction, so this is construction-time work.

Scans a catalog reduction *rejected* carry no information and must be given
`σ = Inf` by the caller — a non-positive formal error would otherwise sneak
the scan back in at weight `1/|σ|`, or produce an Inf-weighted row.
"""
function prepare_pinv_5param(A_prepared, σ)
    all(>(0), σ) || error(
        "prepare_pinv_5param: every σ must be positive. Map rejected scans to " *
        "Inf explicitly (`σ = map(s -> s > 0 ? float(s) : Inf, sres)`) so that " *
        "the weighting is a decision, not an accident.")
    A_scaled = A_prepared ./ σ
    return pinv(A_scaled) ./ permutedims(σ)
end

"""
    fit_5param_prepared(A_prepared, cosϕ, sinϕ, Δα_mas, Δδ_mas,
                        residuals=0.0, σ_formal=0.0;
                        include_chi2=Val(false), buf=Bumper.default_buffer())

Refit the five astrometric parameters to a modelled sky path, exactly as a
catalog pipeline would. Returns `(; parameters)` with `parameters` an
`SVector` in the order `(Δα*, Δδ, μα*, μδ, ϖ)` — note the parallax comes
**last**, while it is the third column of `A_prepared`.

`Δα_mas`/`Δδ_mas` are the modelled offsets per scan; `residuals` is an
optional additive along-scan term (the catalog's own abscissa residuals, for
an IAD channel) and may be a scalar or a vector.

`σ_formal` may be a scalar or a per-epoch vector, and the distinction
matters:

  - **scalar** — the weights cancel out of the least-squares solution
    (`(A/σ) \\ (b/σ) == A \\ b`), so the solve is done unweighted and `1/σ²`
    is folded into the χ² at the end;
  - **vector** — the weights genuinely change the solution, so the solve is
    performed weighted and the residuals come out already whitened.

With `include_chi2=Val(true)` the returned NamedTuple also carries
`chi_squared_astro`, `chi2_reduced` and `dof`; that requires `σ_formal != 0`.

For a fixed `A` and a fixed per-epoch `σ`, prefer
[`prepare_pinv_5param`](@ref) plus [`fit_5param_pinv`](@ref): same answer,
one matvec instead of a factorization per sample.
"""
function fit_5param_prepared(A_prepared, cosϕ, sinϕ, Δα_mas, Δδ_mas,
                             residuals=0.0, σ_formal=0.0;
                             include_chi2=Val(false),
                             buf=Bumper.default_buffer())
    n_obs = size(A_prepared, 1)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    # For a *scalar* weight the explicit σ_formal scaling is redundant: the
    # LSQ x = (A/σ) \ (b/σ) is identical to x = A \ b — the σ cancels — so we
    # solve unweighted and fold 1/σ² into the chi² at the end. For a
    # *per-epoch* σ vector (e.g. Hipparcos `sres`) the weights do NOT cancel
    # and we perform the genuinely weighted solve. Either way we solve
    # against a Bumper-allocated copy of A_prepared (some `\` overloads
    # require an owned matrix rather than a plain view).
    @no_escape buf begin
        b = @alloc(T, n_obs)
        A_dense = @alloc(T, n_obs, size(A_prepared, 2))

        @. b = Δα_mas * cosϕ + Δδ_mas * sinϕ + residuals
        if σ_formal isa Number
            @. A_dense = A_prepared
        else
            @. A_dense = A_prepared / σ_formal
            @. b = b / σ_formal
        end

        x = A_dense \ b

        parameters = @SVector [x[1], x[2], x[4], x[5], x[3]]

        if include_chi2 == Val(true)
            if σ_formal == 0
                error("Asked for `include_chi2=true` but `σ_formal==0`")
            end
            model_predictions = @alloc(T, n_obs)
            residuals_buf = @alloc(T, n_obs)
            mul!(model_predictions, A_dense, x)
            residuals_buf .= b .- model_predictions

            # Scalar σ: residuals are unweighted, apply 1/σ² here.
            # Vector σ: residuals are already whitened by the weighted solve.
            chi_squared_astro = σ_formal isa Number ?
                                dot(residuals_buf, residuals_buf) / (σ_formal * σ_formal) :
                                dot(residuals_buf, residuals_buf)

            n_parameters = 5
            dof = n_obs - n_parameters
            chi2_reduced = chi_squared_astro / dof
        end
    end

    include_chi2 == Val(true) || return (; parameters)
    return (; parameters, chi_squared_astro, chi2_reduced, dof)
end

"""
    fit_5param_pinv(pinv_A, cosϕ, sinϕ, Δα_mas, Δδ_mas, residuals=0.0;
                    buf=Bumper.default_buffer())

As `fit_5param_prepared(…; include_chi2=Val(false))` but using a
pre-computed pseudo-inverse from [`prepare_pinv_5param`](@ref), so the fit
is one matrix-vector product. Returns `(; parameters)` in the same
`(Δα*, Δδ, μα*, μδ, ϖ)` order.
"""
function fit_5param_pinv(pinv_A, cosϕ, sinϕ, Δα_mas, Δδ_mas, residuals=0.0;
                         buf=Bumper.default_buffer())
    n_obs = size(pinv_A, 2)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))
    @no_escape buf begin
        b = @alloc(T, n_obs)
        @. b = Δα_mas * cosϕ + Δδ_mas * sinϕ + residuals
        x_buf = @alloc(T, 5)
        mul!(x_buf, pinv_A, b)
        parameters = @SVector [x_buf[1], x_buf[2], x_buf[4], x_buf[5], x_buf[3]]
    end
    return (; parameters)
end

"""
    fit_4param_prepared(A_factored, cosϕ, sinϕ, Δα_mas, Δδ_mas, σ_formal=0.0;
                        buf=Bumper.default_buffer())

Four-parameter counterpart of [`fit_5param_prepared`](@ref) (parallax held
fixed). `A_factored` may be a factorization of `prepare_A_4param`'s output.
Returns `(; parameters)` as `(Δα*, Δδ, μα*, μδ)`.
"""
function fit_4param_prepared(A_factored, cosϕ, sinϕ, Δα_mas, Δδ_mas, σ_formal=0.0;
                             buf=Bumper.default_buffer())
    n_obs = length(cosϕ)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))
    @no_escape buf begin
        b = @alloc(T, n_obs)
        @. b = Δα_mas * cosϕ + Δδ_mas * sinϕ
        if σ_formal != 0
            @. b *= 1 / σ_formal
        end
        x = A_factored \ b
        parameters = @SVector [x[1], x[2], x[3], x[4]]
    end
    return (; parameters)
end

# ---------------------------------------------------
# Frame offsets for non-Gaia astrometry
#
# Astrometry taken on a *different* frame from the one the model's reference
# point lives on needs its own five-parameter astrometric solution: the
# instrument's own zero point, parallax and proper motion, which are nuisance
# parameters rather than physics. Hipparcos is the case in the tree today
# (`G23HObs`'s `:iad_hip` channel and `HipparcosIADObs`); HST FGS is the next
# one, and it wants this rather than a second copy of it.
#
# What is emphatically *not* wanted is a frame offset per Gaia source. Several
# `G23HObs` in one system share one frame, and that shared frame is what binds
# the system together — the wide pair's relative astrometry constrains the
# wide orbit for free. Inventing per-source offsets would dilute exactly the
# constraint the joint fit exists to exploit. So: offsets for other frames,
# never between sources on the same one.
# ---------------------------------------------------

"""
    FrameOffset(Δra, Δdec, plx, pmra, pmdec)

The five-parameter astrometric solution of an instrument's own frame, in the
tangent plane about that instrument's catalog position: a position offset
[mas], a parallax [mas], and a proper motion [mas/yr].

Built by [`frame_offset`](@ref) from an observation's variables and consumed
by [`frame_offset_alongscan`](@ref). `isbits`, so building one per evaluation
is free.
"""
struct FrameOffset{T<:Number}
    Δra::T
    Δdec::T
    plx::T
    pmra::T
    pmdec::T
end
FrameOffset(Δra, Δdec, plx, pmra, pmdec) =
    FrameOffset(promote(Δra, Δdec, plx, pmra, pmdec)...)

"""
    frame_offset(θ_obs, plx_anchor, ::Type{T}) -> FrameOffset{T}

Read an observation's frame-offset block from its own variables. The names
are fixed, and are the ones `G23HObs` has used since the port:

| variable     | meaning                                                    |
|--------------|------------------------------------------------------------|
| `iad_Δra`    | position offset in α✱ at the instrument's reference epoch  |
| `iad_Δdec`   | position offset in δ                                       |
| `iad_Δplx`   | parallax offset from `plx_anchor`                          |
| `iad_pmra`   | proper motion in α✱ *(absolute, not an offset)*            |
| `iad_pmdec`  | proper motion in δ                                         |

Each is optional and defaults to zero, so an instrument that fits only some
of the five declares only those. The proper motions are absolute because
their natural anchor is the instrument catalog's own solution, which a
`@variables` block expresses directly (`iad_pmra = 4.53 + iad_Δpmra`).

`plx_anchor` is supplied by the caller rather than looked up, because whose
parallax it is, is a modelling decision: `G23HObs` anchors on the Hipparcos
catalog value (its abscissa channel exists for the companion curvature, and
the frame is pure nuisance), while a Hipparcos-only fit anchors on the
system's own `plx` so that the data actually constrain it.
"""
@inline function frame_offset(θ_obs, plx_anchor, ::Type{T}) where {T}
    return FrameOffset{T}(
        _frame_var(θ_obs, :iad_Δra, T),
        _frame_var(θ_obs, :iad_Δdec, T),
        T(plx_anchor) + _frame_var(θ_obs, :iad_Δplx, T),
        _frame_var(θ_obs, :iad_pmra, T),
        _frame_var(θ_obs, :iad_pmdec, T))
end

@inline _frame_var(θ_obs, name::Symbol, ::Type{T}) where {T} =
    hasproperty(θ_obs, name) ? T(getproperty(θ_obs, name)) : zero(T)

"""
    frame_offset_alongscan(off, Δt_yr, cosϕ, sinϕ, parallax_factor_al,
                           Δα=0, Δδ=0)

Project a [`FrameOffset`](@ref) onto one scan: the along-scan coordinate the
instrument would have measured for a source at `Δα`/`Δδ` [mas] off the
frame's own path, `Δt_yr` Julian years after the frame's reference epoch.

    b = (Δra + Δt·pmra + Δα)·cosϕ + (Δdec + Δt·pmdec + Δδ)·sinϕ + ϖ·f_AL

`Δα`/`Δδ` carry whatever the source itself is doing — an orbital reflex, a
photocentre wobble, the Hipparcos grating response. Compare `b` against the
measured abscissa (`proj_meas_alongscan` in a Hipparcos IAD table).
"""
@inline frame_offset_alongscan(off::FrameOffset, Δt_yr, cosϕ, sinϕ, parf,
                               Δα=zero(off.Δra), Δδ=zero(off.Δdec)) =
    (off.Δra + Δt_yr * off.pmra + Δα) * cosϕ +
    (off.Δdec + Δt_yr * off.pmdec + Δδ) * sinϕ +
    off.plx * parf

export FrameOffset, frame_offset, frame_offset_alongscan

# ---------------------------------------------------
# Per-epoch offset accumulation
# ---------------------------------------------------

"""
    accumulate_offsets!(Δα, Δδ, ctx, target, reference)
    accumulate_offsets!(Δα, Δδ, ctx, target, reference, rows)

Add the sky offset of `target` from `reference` — `raoff` into `Δα` and
`decoff` into `Δδ`, both [mas] — for each of this observation's table rows.

`target` and `reference` are resolved reference *values*, not specs: a
`BodyRef`, or a `WeightedPoint` from `ref(ctx, spec)`, `barycentre`,
`photocentre`, or built by the likelihood itself. Resolve them once outside
the loop.

With `rows`, only those table rows are visited, and the results land in
`Δα[k]`/`Δδ[k]` for `k in eachindex(rows)` — so a likelihood whose table
concatenates several instrument channels can give each channel its own
buffer and its own row range. Without it, every row is visited and the
buffers are indexed by table row.

The buffers are **accumulated into**, not overwritten, so a caller can lay
down a reference-point linear motion first (as `GaiaDR4AstromObs` does) and
add the orbital excursion on top.

This is the v2 replacement for `_simulate_skypath_perturbations!`. It is
deliberately one call for the whole source rather than one per companion:
the photocentre of several luminous bodies is the flux-weighted mean of
their *apparent positions*, which is not the superposition of per-companion
photocentres — the point BINARYS makes, and the thing v1's positional
`fluxratio` vector could not represent.
"""
function accumulate_offsets!(Δα, Δδ, ctx::ObsContext, target, reference)
    @inbounds for i in eachindex(ctx.epoch_index)
        sol = solutionat(ctx, i)
        Δα[i] += raoff(sol, target, reference)
        Δδ[i] += decoff(sol, target, reference)
    end
    return (Δα, Δδ)
end

function accumulate_offsets!(Δα, Δδ, ctx::ObsContext, target, reference, rows)
    length(Δα) == length(rows) && length(Δδ) == length(rows) || error(
        "accumulate_offsets!: the output buffers are indexed by position within " *
        "`rows`, so they must have length $(length(rows)); got $(length(Δα)) and " *
        "$(length(Δδ)).")
    @inbounds for (k, i) in enumerate(rows)
        sol = solutionat(ctx, i)
        Δα[k] += raoff(sol, target, reference)
        Δδ[k] += decoff(sol, target, reference)
    end
    return (Δα, Δδ)
end
