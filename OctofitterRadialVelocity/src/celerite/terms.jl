import Base.+, Base.*, Base.length, Base.product

# Abstract base class and definitions
abstract type Term{T} end
_numtype(::Term) = Float64

function get_terms(term::Term)
    return tuple(term)
end

function get_real_coefficients(term::Term)
    return (Float64[], Float64[])
end

function get_complex_coefficients(term::Term)
    return (Float64[], Float64[], Float64[], Float64[])
end

function get_all_coefficients(term::Term)
    rc = get_real_coefficients(term)
    cc = get_complex_coefficients(term)
    return (rc[1], rc[2], cc[1], cc[2], cc[3], cc[4])
end

function get_value(term::Term, tau)
    coeffs = get_all_coefficients(term)
    t = abs.(tau)
    k = similar(tau)
    fill!(k, 0)
    for i in 1:length(coeffs[1])
        k .+= coeffs[1][i] .* exp.(-coeffs[2][i] .* t)
    end
    for i in 1:length(coeffs[3])
        k .+= (coeffs[3][i] .* cos.(coeffs[6][i] .* t) .+ coeffs[4][i] .* sin.(coeffs[6][i] .* t)) .* exp.(-coeffs[5][i] .* t)
    end
    return k
end
get_psd2(term, omega) =get_psd(term, omega)

function get_psd(term, omega)
    coeffs = get_all_coefficients(term)
    w2 = omega.^2
    w4 = w2.^2
    p = zeros(eltype(w2), size(w2))
    for i in 1:length(coeffs[1])
        a = coeffs[1][i]
        c = coeffs[2][i]
        p = p + a*c ./ (c*c + w2)
    end
    for i in 1:length(coeffs[3])
        a = coeffs[3][i]
        b = coeffs[4][i]
        c = coeffs[5][i]
        d = coeffs[6][i]
        w02 = c*c+d*d
        @. p = p + ((a*c+b*d)*w02+(a*c-b*d).*w2) ./ (w4 + 2.0*(c*c-d*d).*w2+w02*w02)
    end
    return sqrt(2.0 / pi) .* p
end

function length(term::Term)
    return 0
end

function get_parameter_vector(term::Term)
    return Float64[]
end

function set_parameter_vector!(term::Term, vector::Array)
end

# A sum of terms
struct TermSum{T, N, A<:NTuple{N,Term{T}}} <: Term{T}
    terms::A
end

function get_terms(term::TermSum)
    return term.terms
end

function get_all_coefficients(term_sum::TermSum)
    a_real = Float64[]
    c_real = Float64[]
    a_complex = Float64[]
    b_complex = Float64[]
    c_complex = Float64[]
    d_complex = Float64[]
    for term in term_sum.terms
        coeffs = get_all_coefficients(term)
#        a_real = cat(1, a_real, coeffs[1])
#        c_real = cat(1, c_real, coeffs[2])
#        a_complex = cat(1, a_complex, coeffs[3])
#        b_complex = cat(1, b_complex, coeffs[4])
#        c_complex = cat(1, c_complex, coeffs[5])
#        d_complex = cat(1, d_complex, coeffs[6])
        a_real = vcat(a_real, coeffs[1])
        c_real = vcat(c_real, coeffs[2])
        a_complex = vcat(a_complex, coeffs[3])
        b_complex = vcat(b_complex, coeffs[4])
        c_complex = vcat(c_complex, coeffs[5])
        d_complex = vcat(d_complex, coeffs[6])
    end
    return (a_real, c_real, a_complex, b_complex, c_complex, d_complex)
end

function +(t1::Term, t2::Term)
    return TermSum((get_terms(t1)..., get_terms(t2)...))
end

function length(term_sum::TermSum)
    return +(map(length, term_sum.terms)...)
end

function get_parameter_vector(term_sum::TermSum)
    return vcat(map(get_parameter_vector, term_sum.terms)...)
end

function set_parameter_vector!(term_sum::TermSum, vector::Array)
    index = 1
    for term in term_sum.terms
        len = length(term)
        set_parameter_vector!(term, vector[index:index+len-1])
        index = index + len
    end
end

# A product of two terms
struct TermProduct{T,A<:Term{T},B<:Term{T}} <: Term{T}
    term1::A
    term2::B
end

function chain(a...)
    return (el for it in a for el in it)
end

function get_all_coefficients(term_sum::TermProduct)
    c1 = get_all_coefficients(term_sum.term1)
    nr1, nc1 = length(c1[1]), length(c1[3])
    c2 = get_all_coefficients(term_sum.term2)
    nr2, nc2 = length(c2[1]), length(c2[3])

    # First compute real terms
    nr = nr1 * nr2
    ar = Array{Float64}(undef, nr)
    cr = Array{Float64}(undef, nr)
    gen = product(zip(c1[1], c1[2]), zip(c2[1], c2[2]))
    for (i, ((aj, cj), (ak, ck))) in enumerate(gen)
        ar[i] = aj * ak
        cr[i] = cj + ck
    end

    # Then the complex terms
    nc = nr1 * nc2 + nc1 * nr2 + 2 * nc1 * nc2
    ac = Array{Float64}(undef, nc)
    bc = Array{Float64}(undef, nc)
    cc = Array{Float64}(undef, nc)
    dc = Array{Float64}(undef, nc)

    # real * complex
    gen = product(zip(c1[1], c1[2]), zip(c2[3:end]...))
    gen = chain(gen, product(zip(c2[1], c2[2]), zip(c1[3:end]...)))
    for (i, ((aj, cj), (ak, bk, ck, dk))) in enumerate(gen)
        ac[i] = aj * ak
        bc[i] = aj * bk
        cc[i] = cj + ck
        dc[i] = dk
    end

    # complex * complex
    i0 = nr1 * nc2 + nc1 * nr2 + 1
    gen = product(zip(c1[3:end]...), zip(c2[3:end]...))
    for (i, ((aj, bj, cj, dj), (ak, bk, ck, dk))) in enumerate(gen)
        ac[i0 + 2*(i-1)] = 0.5 * (aj * ak + bj * bk)
        bc[i0 + 2*(i-1)] = 0.5 * (bj * ak - aj * bk)
        cc[i0 + 2*(i-1)] = cj + ck
        dc[i0 + 2*(i-1)] = dj - dk

        ac[i0 + 2*(i-1) + 1] = 0.5 * (aj * ak - bj * bk)
        bc[i0 + 2*(i-1) + 1] = 0.5 * (bj * ak + aj * bk)
        cc[i0 + 2*(i-1) + 1] = cj + ck
        dc[i0 + 2*(i-1) + 1] = dj + dk
    end

    return ar, cr, ac, bc, cc, dc
end

function *(t1::Term, t2::Term)
    return TermProduct(t1, t2)
end

function length(term_prod::TermProduct)
    return length(term_prod.term1) + length(term_prod.term2)
end

function get_parameter_vector(term_prod::TermProduct)
    return cat(1, get_parameter_vector(term_prod.term1), get_parameter_vector(term_prod.term2))
end

function set_parameter_vector!(term_prod::TermProduct, vector::Array)
    l = length(term_prod.term1)
    set_parameter_vector!(term_prod.term1, vector[1:l])
    set_parameter_vector!(term_prod.term2, vector[l+1:end])
end

# A real term where b=0 and d=0
struct RealTerm{T} <: Term{T}
    log_a::T
    log_c::T
end

function get_real_coefficients(term::RealTerm)
    return [exp(term.log_a)], [exp(term.log_c)]
end

function length(term::RealTerm)
    return 2
end

function get_parameter_vector(term::RealTerm)
    return [term.log_a, term.log_c]
end

function set_parameter_vector!(term::RealTerm, vector::Array)
    term.log_a = vector[1]
    term.log_c = vector[2]
end

# General Celerite term
struct ComplexTerm{T} <: Term{T}
    log_a::T
    log_b::T
    log_c::T
    log_d::T
end

function get_complex_coefficients(term::ComplexTerm)
    return [exp(term.log_a)], [exp(term.log_b)], [exp(term.log_c)], [exp(term.log_d)]
end

function length(term::ComplexTerm)
    return 4
end

function get_parameter_vector(term::ComplexTerm)
    return [term.log_a, term.log_b, term.log_c, term.log_d]
end

function set_parameter_vector!(term::ComplexTerm, vector::Array)
    term.log_a = vector[1]
    term.log_b = vector[2]
    term.log_c = vector[3]
    term.log_d = vector[4]
end

# A SHO term
struct SHOTerm{T} <: Term{T}
    log_S0::T
    log_Q::T
    log_omega0::T
end

function get_real_coefficients(term::SHOTerm)
    Q = exp(term.log_Q)
    if Q >= 0.5
        return Float64[], Float64[]
    end
    S0 = exp(term.log_S0)
    w0 = exp(term.log_omega0)
    f = sqrt(1.0 - 4.0 * Q^2)
    return (
        0.5*S0*w0*Q*[1.0+1.0/f, 1.0-1.0/f],
        0.5*w0/Q*[1.0-f, 1.0+f]
    )
end

function get_complex_coefficients(term::SHOTerm)
    Q = exp(term.log_Q)
    if Q < 0.5
        return Float64[], Float64[], Float64[], Float64[]
    end
    S0 = exp(term.log_S0)
    w0 = exp(term.log_omega0)
    f = sqrt(4.0 * Q^2 - 1.0)
    return (
        [S0 * w0 * Q],
        [S0 * w0 * Q / f],
        [0.5 * w0 / Q],
        [0.5 * w0 / Q * f],
    )
end

function length(term::SHOTerm)
    return 3
end

function get_parameter_vector(term::SHOTerm)
    return [term.log_S0, term.log_Q, term.log_omega0]
end

function set_parameter_vector!(term::SHOTerm, vector::Array)
    term.log_S0 = vector[1]
    term.log_Q = vector[2]
    term.log_omega0 = vector[3]
end
