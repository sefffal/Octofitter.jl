# ---------------------------------------------------
# Banded symmetric positive-definite matrices
#
# The covariance of one RVS epoch spectrum is banded: neighbouring bins share
# detector pixels, and the template-mismatch kernel has compact support. It
# is fixed by the data, so it is factorised once when the observation is
# built. A likelihood evaluation then needs only forward substitutions,
# which are written for any number type so that they run on dual numbers.
# ---------------------------------------------------

"""
    BandedCholesky(A, b)

Cholesky factor `L` of the symmetric positive-definite matrix `A`, whose
entries are zero more than `b` places from the diagonal (`A = L Lᵀ`, `L`
lower triangular with the same bandwidth). Only the lower band of `A` is
read.

`rows[k + 1, i]` holds `L[i, i - k]`, so the entries one substitution step
needs are contiguous.
"""
struct BandedCholesky
    rows::Matrix{Float64}
    logdet::Float64
end

Base.size(F::BandedCholesky) = (size(F.rows, 2), size(F.rows, 2))
bandwidth(F::BandedCholesky) = size(F.rows, 1) - 1

function BandedCholesky(A::AbstractMatrix, b::Integer)
    n = size(A, 1)
    rows = zeros(b + 1, n)
    logdet = 0.0
    for i in 1:n
        for j in max(1, i - b):i
            acc = Float64(A[i, j])
            # L[i, j] = (A[i, j] − Σ_m L[i, m] L[j, m]) / L[j, j]
            for m in max(1, i - b):j-1
                acc -= rows[i-m+1, i] * rows[j-m+1, j]
            end
            if j < i
                rows[i-j+1, i] = acc / rows[1, j]
            else
                acc > 0 || throw(ArgumentError("the matrix is not positive definite (pivot $acc at row $i)"))
                rows[1, i] = sqrt(acc)
                logdet += log(acc)
            end
        end
    end
    return BandedCholesky(rows, logdet)
end

"""
    whiten!(x, F::BandedCholesky) -> x

Overwrite `x` with `L⁻¹ x`. Works for any number type.
"""
function whiten!(x::AbstractVector, F::BandedCholesky)
    R = F.rows
    b = bandwidth(F)
    @inbounds for i in eachindex(x)
        acc = x[i]
        for k in 1:min(b, i - 1)
            acc -= R[k+1, i] * x[i-k]
        end
        x[i] = acc / R[1, i]
    end
    return x
end

"""
    colour!(x, F::BandedCholesky) -> x

Overwrite `x` with `L x`: unit white noise becomes noise with covariance `A`.
"""
function colour!(x::AbstractVector, F::BandedCholesky)
    R = F.rows
    b = bandwidth(F)
    @inbounds for i in reverse(eachindex(x))
        acc = R[1, i] * x[i]
        for k in 1:min(b, i - 1)
            acc += R[k+1, i] * x[i-k]
        end
        x[i] = acc
    end
    return x
end

# In-place Cholesky of a small dense symmetric positive-definite matrix
# (lower triangle), for the mode coefficients. Returns the log determinant,
# or -Inf if the matrix is not positive definite.
function _cholesky_small!(A::AbstractMatrix)
    n = size(A, 1)
    logdet = zero(eltype(A))
    @inbounds for j in 1:n
        d = A[j, j]
        for m in 1:j-1
            d -= A[j, m]^2
        end
        d > 0 || return oftype(logdet, -Inf)
        A[j, j] = sqrt(d)
        logdet += log(d)
        for i in j+1:n
            acc = A[i, j]
            for m in 1:j-1
                acc -= A[i, m] * A[j, m]
            end
            A[i, j] = acc / A[j, j]
        end
    end
    return logdet
end

# x ← L⁻¹ x for the dense lower factor left by `_cholesky_small!`.
function _forward_small!(x::AbstractVector, L::AbstractMatrix)
    @inbounds for i in eachindex(x)
        acc = x[i]
        for m in 1:i-1
            acc -= L[i, m] * x[m]
        end
        x[i] = acc / L[i, i]
    end
    return x
end
