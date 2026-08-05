# ---------------------------------------------------
# Kernel phases
#
# Closure phases from a redundant array are linearly dependent: with six
# baselines the four T3 triangles span only three independent directions. The
# kernel-phase basis P₁ is an orthonormal basis of that span, obtained from the
# Cholesky factor of `Tλ Tλᵀ` where `Tλ` is the closure design matrix
# replicated across spectral channels. Projecting the CP residuals onto it
# removes the null direction, which would otherwise be a perfectly-correlated
# row of the covariance and make it singular.
#
# This is the only stage where the merged observation type still branches:
# everything up to and including `cps_model` is shared with the plain
# closure-phase model.
# ---------------------------------------------------

# GRAVITY's four closure triangles over its six baselines. The rows are the
# triangles and the columns the baselines, with the same +1 +1 −1 combination
# `closurephase!` applies.
const GRAVITY_T3_DESIGN = Int8[
    1 -1 0 1 0 0
    1 0 -1 0 1 0
    0 1 -1 0 0 1
    0 0 0 1 -1 1
]

"""
    kernel_phase_basis(design, Λ) -> (Tλ, P₁)

Replicate the closure design matrix `design` across `Λ` spectral channels and
return the wavelength-blocked design `Tλ` together with an orthonormal basis
`P₁` of its row space (the kernel phases), as a `n_kp × (n_T3·Λ)` matrix.

Channels are grouped *within* each baseline block, so `P₁` acts on a residual
vector laid out as `(i_T3 − 1)·Λ + i_wave` — the layout the likelihood fills.
That ordering is what makes the resulting kernel-phase covariance block
diagonal in the kernel-phase index, which the Cholesky path relies on.
"""
function kernel_phase_basis(design::AbstractMatrix, Λ::Integer)
    Tλ = zeros(Int8, Λ * size(design, 1), Λ * size(design, 2))
    for baseline_i in axes(design, 1), baseline_j in axes(design, 2)
        for wavelength_i in 1:Λ
            Tλ[wavelength_i+(baseline_i-1)*Λ, wavelength_i+(baseline_j-1)*Λ] =
                design[baseline_i, baseline_j]
        end
    end

    C, _ = cholesky(Tλ * Tλ')
    P₁ = collect(C) ./ sqrt.(diag(C * C'))
    i_max = findfirst(<=(1e-5), diag(P₁)) - 1
    P₁ = P₁[:, 1:i_max]'
    return Tλ, P₁
end

"""
    _prepare_kernel_phases(row)

Add the per-exposure kernel-phase machinery to a prepared data row: the
blocked design matrix `Tλ`, the basis `P₁`, and the kernel-phase
uncertainties `σ_kp = P₁ · vec(dcps)` (constant, so it is computed once here
rather than every likelihood evaluation).
"""
function _prepare_kernel_phases(row)
    n_T3 = size(row.cps_data, 1)
    n_baselines = size(row.u, 1)
    (n_T3, n_baselines) == size(GRAVITY_T3_DESIGN) || error("""
        The kernel-phase model currently carries only GRAVITY's closure design
        matrix: 6 baselines forming 4 triangles. This exposure has
        $n_baselines baselines and $n_T3 closure phases.

        Use `kernel_phases=false` (plain closure phases with diagonal
        uncertainties), which works for any array.
        """)

    Λ = length(row.eff_wave)
    Tλ, P₁ = kernel_phase_basis(GRAVITY_T3_DESIGN, Λ)
    # GRAVITY's design matrix has rank 3 per channel. The correlation model
    # (`CKP!`) and the block Cholesky both partition the kernel phases into
    # exactly that many blocks of Λ, so a different rank would silently
    # mis-block the covariance rather than fail.
    size(P₁, 1) == 3Λ || error(
        "expected 3 kernel phases per spectral channel, got $(size(P₁, 1)) over " *
        "$Λ channels. The correlation model and the block Cholesky both assume 3.")
    σ_kp = P₁ * vec(row.dcps)
    return (; row..., Tλ, P₁, σ_kp)
end

"""
    _kp_covariance!(Σ, C, σ_kp, jitter)

Build the kernel-phase covariance from the correlation matrix `C`, the
per-kernel-phase uncertainties `σ_kp` and a jitter added along the diagonal:
`Σ = diag(σ) C diag(σ) + jitter² I`.
"""
function _kp_covariance!(Σ, C, σ_kp, jitter)
    @inbounds for j in axes(Σ, 2), i in axes(Σ, 1)
        Σ[i, j] = σ_kp[i] * C[i, j] * σ_kp[j]
    end
    @inbounds for i in axes(Σ, 1)
        Σ[i, i] += jitter^2
    end
    return Σ
end

"""
    _kp_logpdf(Σ, resids, n_blocks)

`logpdf(MvNormal(Σ), resids)` for the block-diagonal kernel-phase covariance,
factorizing each of the `n_blocks` diagonal blocks independently rather than
the whole matrix. Returns `-Inf` (in the residuals' number type) if any block
is not positive definite, which is what a sampler needs from an invalid draw.

The blocks are the `n_blocks` kernel phases, each spanning all `Λ` spectral
channels — the same partition `CKP` populates. v1 sized them by
`length(table.eff_wave)`, which is the number of *exposures* in the table, not
the number of channels in this one; the two coincide only by accident, and
when they do not, the `Cholesky` handed to `PDMat` is not a factorization of
`Σ` at all. `test/runtests.jl` pins this against a plain dense `MvNormal`.
"""
function _kp_logpdf(Σ, resids, n_blocks::Integer)
    T = eltype(resids)
    L = size(Σ, 1)
    Λ = L ÷ n_blocks
    Λ * n_blocks == L || error(
        "kernel-phase covariance is $(L)×$(L), which is not $n_blocks blocks of equal size")
    try
        for k in 1:n_blocks
            rows = (k-1)*Λ+1:k*Λ
            cholesky!(Hermitian(@view Σ[rows, rows]))
        end
        # `Σ` now holds the concatenated upper-triangular factors, and the
        # off-block entries are zero, so this *is* the Cholesky factor of the
        # whole (block-diagonal) matrix. Assembling it by hand skips a second,
        # redundant factorization of the full L×L matrix.
        P = PDMat(Σ, Cholesky(UpperTriangular(Σ)))
        return logpdf(MvNormal(P), resids)
    catch err
        err isa InterruptException && rethrow()
        return convert(T, -Inf)
    end
end
