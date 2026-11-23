# ---------------------------------------------------
# Image block, part 3: the log-likelihood of a transit with its gradient
#
# Differentiating the likelihood of epoch images with dual numbers costs one
# partial per sampled parameter in every sample of every window, although a
# transit depends on the parameters only through the scene: where each point
# source is in each window, its photon rate and its chromatic moments. So
# the gradient with respect to the scene is carried by hand, in reverse, in
# plain floating point, and the chain rule to the parameters is applied once
# per transit.
#
# With J = −2 ln L of one state, G the Gram matrix of a window (of its mode
# vectors and the residual, weighted by the inverse variance) and H what the
# windows leave for the transit,
#
#     dJ = Σ_windows [ ⟨Ḡ, dG⟩ + Σ_samples d ln v ] ,
#     Ḡ = Z H̄ Zᵀ + (the window's posterior covariance, on its own modes),
#
# where Z maps the transit's modes and the residual onto the window's
# vectors with the window's coefficients eliminated, and H̄ is the transit's
# posterior covariance bordered by its solution. Per sample that is an
# adjoint for each entry of the sample's vector and one for its variance,
# and from there for each source's electrons, position and moments.
# ---------------------------------------------------

import ForwardDiff

# Second derivatives of the Catmull–Rom weights with respect to `t`.
@inline _cubic_ddweights(t) = (2 - 3t, -5 + 9t, 4 - 9t, -1 + 3t)

# One source at one sample with everything the reverse pass needs:
# (f, f_w, f_z, f_ww, f_wz, f_zz), and per shape mode (s, s_w, s_z), the
# subscripts being derivatives with respect to the sample's offset; and
# the stencil (i, j, tx, ty), with i = 0 beyond the table.
@inline function _point2(b::ImageBasis, cvec::NTuple{NT,Float64}, Δw, Δz, ::Val{M}, ::Val{NM}) where {NT,M,NM}
    ux, uy = (Δw - b.x0) / b.hx, (Δz - b.y0) / b.hy
    i, j = floor(Int, ux), floor(Int, uy)
    z3 = ntuple(_ -> 0.0, Val(M))
    (i >= 1 && i <= size(b.terms, 2) - 3 && j >= 1 && j <= size(b.terms, 3) - 3) ||
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), z3, z3, z3, (0, 0, 0.0, 0.0)
    tx, ty = ux - i, uy - j
    wx, wy = _cubic_weights(tx), _cubic_weights(ty)
    dx, dy = _cubic_dweights(tx), _cubic_dweights(ty)
    ex, ey = _cubic_ddweights(tx), _cubic_ddweights(ty)
    C, S = b.terms, b.modes
    f = fw = fz = fww = fwz = fzz = 0.0
    s, sw, sz = z3, z3, z3
    @inbounds for q in 1:4
        jq = j + q - 1
        r0 = r1 = r2 = 0.0
        rs, rsw = z3, z3
        for p in 1:4
            ip = i + p - 1
            c = 0.0
            @simd for n in 1:NT
                c += cvec[n] * C[n, ip, jq]
            end
            r0 += wx[p] * c
            r1 += dx[p] * c
            r2 += ex[p] * c
            node = _mode_nodes(S, cvec, ip, jq, Val(M), Val(NM))
            rs = rs .+ wx[p] .* node
            rsw = rsw .+ dx[p] .* node
        end
        f += wy[q] * r0
        fw += wy[q] * r1
        fz += dy[q] * r0
        fww += wy[q] * r2
        fwz += dy[q] * r1
        fzz += ey[q] * r0
        s = s .+ wy[q] .* rs
        sw = sw .+ wy[q] .* rsw
        sz = sz .+ dy[q] .* rs
    end
    hx, hy = b.hx, b.hy
    return (f, fw / hx, fz / hy, fww / hx^2, fwz / (hx * hy), fzz / hy^2), s, sw ./ hx, sz ./ hy, (i, j, tx, ty)
end

# Add to `cbar` (the adjoint of one source's (1, moments...)) what the
# adjoints `a` of (f, f_w, f_z) and `bm` of the shape values imply.
@inline function _scatter_moments!(cbar::AbstractVector, b::ImageBasis, stencil, a::NTuple{3,Float64}, bm::NTuple{M,Float64},
    ::Val{NM}, ::Val{NT}) where {M,NM,NT}
    i, j, tx, ty = stencil
    i == 0 && return cbar
    wx, wy = _cubic_weights(tx), _cubic_weights(ty)
    dx, dy = _cubic_dweights(tx), _cubic_dweights(ty)
    C, S = b.terms, b.modes
    a2, a3 = a[2] / b.hx, a[3] / b.hy
    @inbounds for q in 1:4, p in 1:4
        ip, jq = i + p - 1, j + q - 1
        wxy = wx[p] * wy[q]
        na = a[1] * wxy + a2 * dx[p] * wy[q] + a3 * wx[p] * dy[q]
        @simd for n in 1:NT
            cbar[n] += na * C[n, ip, jq]
        end
        for m in 1:M
            nb = bm[m] * wxy
            for n in 1:NM
                cbar[n] += nb * S[m, n, ip, jq]
            end
        end
    end
    return cbar
end

# x ← L⁻ᵀ x for the dense lower factor left by `_cholesky_small!`.
function _backward_small!(x::AbstractVector, L::AbstractMatrix)
    @inbounds for i in reverse(eachindex(x))
        acc = x[i]
        for m in i+1:length(x)
            acc -= L[m, i] * x[m]
        end
        x[i] = acc / L[i, i]
    end
    return x
end

# Inverse of the symmetric positive-definite matrix whose lower triangle is
# in `A`, into `Ainv`; `A` is left factorised. Returns log det A, or -Inf.
function _spd_inverse!(Ainv::AbstractMatrix, A::AbstractMatrix)
    n = size(A, 1)
    logdet = _cholesky_small!(A)
    logdet == -Inf && return logdet
    @inbounds for c in 1:n
        col = view(Ainv, :, c)
        fill!(col, 0)
        col[c] = 1
        _forward_small!(col, A)
        _backward_small!(col, A)
    end
    return logdet
end

# As `_eliminate_window!`, keeping what the reverse pass needs: the window's
# posterior covariance `Ainv` and `Y = Ainv G[window modes, transit modes and residual]`.
function _eliminate_window_adjoint!(H::AbstractMatrix, Ainv::AbstractMatrix, Y::AbstractMatrix, G::AbstractMatrix, A::AbstractMatrix,
    d::EpochImages, win::_ImageWindow)
    iw, it = d.window_modes, d.transit_modes
    Kw, Kt, P = length(iw), length(it), size(G, 1)
    inst = d.instrument
    logprec = 0.0
    @inbounds for b in 1:Kw
        for a in b:Kw
            A[a, b] = G[iw[a], iw[b]]
        end
        p = 1 / _window_sigma(inst, win, iw[b])^2
        A[b, b] += p
        logprec += log(p)
    end
    logdet = _spd_inverse!(Ainv, view(A, 1:Kw, 1:Kw))
    logdet == -Inf && return logdet
    @inbounds for c in 1:Kt+1
        jc = c <= Kt ? it[c] : P
        for a in 1:Kw
            acc = 0.0
            for a2 in 1:Kw
                acc += Ainv[a, a2] * G[iw[a2], jc]
            end
            Y[a, c] = acc
        end
    end
    @inbounds for c in 1:Kt+1, c2 in 1:Kt+1
        jc = c <= Kt ? it[c] : P
        j2 = c2 <= Kt ? it[c2] : P
        acc = G[jc, j2]
        for a in 1:Kw
            acc -= G[jc, iw[a]] * Y[a, c2]
        end
        H[c, c2] += acc
    end
    return logdet - logprec
end

# Solve the transit: returns J's part, and fills `Hbar` (dJ = ⟨Hbar, dH⟩) and
# `pbar` (dJ/d precision of each transit mode). `H` is overwritten.
function _transit_solve_adjoint!(Hbar::AbstractMatrix, pbar::AbstractVector, H::AbstractMatrix, A::AbstractMatrix, sigma_t::AbstractVector)
    Kt = size(H, 1) - 1
    logprec = 0.0
    @inbounds for b in 1:Kt
        for a in b:Kt
            A[a, b] = H[a, b]
        end
        p = 1 / sigma_t[b]^2
        A[b, b] += p
        logprec += log(p)
    end
    Ainv = view(Hbar, 1:Kt, 1:Kt)
    logdet = _spd_inverse!(Ainv, view(A, 1:Kt, 1:Kt))
    logdet == -Inf && return Inf
    quad = H[Kt+1, Kt+1]
    @inbounds for a in 1:Kt
        x = 0.0
        for b in 1:Kt
            x += Ainv[a, b] * H[b, Kt+1]
        end
        Hbar[a, Kt+1] = -x
        Hbar[Kt+1, a] = -x
        quad -= H[a, Kt+1] * x
    end
    @inbounds for a in 1:Kt
        for b in 1:Kt
            Hbar[a, b] += Hbar[a, Kt+1] * Hbar[b, Kt+1]
        end
        pbar[a] = Hbar[a, a] - sigma_t[a]^2
    end
    Hbar[Kt+1, Kt+1] = 1
    return quad + logdet - logprec
end

"""
    image_transit_loglike_gradient!(gw, gz, grate, gmoments, gshape, data, t, w, z, rate, moments;
                                    shape = (), flux_scale = 1, jitter = nothing, buf)
        -> (; ll, flux_scale, jitter)

The log-likelihood of transit `t` (the first value of
[`image_transit_loglike`](@ref), same arguments, all `Float64`) and its
gradient: the derivatives with respect to `w`, `z`, `rate`, `moments` and
`shape` are **added** to `gw`, `gz`, `grate`, `gmoments` and `gshape` (as
long as the instrument has shape modes), those with respect to `flux_scale`
and the two `jitter` widths are returned.
"""
function image_transit_loglike_gradient!(gw::AbstractMatrix, gz::AbstractMatrix, grate::AbstractVector, gmoments::AbstractMatrix,
    gshape::AbstractVector, d::EpochImages, t::Integer, w::AbstractMatrix{Float64}, z::AbstractMatrix{Float64},
    rate::AbstractVector{Float64}, moments::AbstractMatrix{Float64}; shape=(), flux_scale::Real=1.0, jitter=nothing,
    buf=Bumper.default_buffer())
    inst = d.instrument
    K, nw = length(rate), length(d.transit_windows[t])
    size(w) == size(z) == size(gw) == size(gz) == (K, nw) || throw(DimensionMismatch("w, z and their gradients must be [source, window] = ($K, $nw)"))
    size(moments) == size(gmoments) == (chromatic_degree(inst), K) || throw(DimensionMismatch("moments must be [degree, source]"))
    length(grate) == K || throw(DimensionMismatch("grate needs one entry per source"))
    length(shape) <= inst.nshapes || throw(DimensionMismatch("$(length(shape)) shape coefficients for $(inst.nshapes) shape modes"))
    length(gshape) == inst.nshapes || throw(DimensionMismatch("gshape needs one entry per shape mode of the instrument"))
    return _image_transit_adjoint!(gw, gz, grate, gmoments, gshape, d, t, w, z, rate, moments, shape, Float64(flux_scale), jitter, buf,
        _static_sizes(d)...)
end

function _image_transit_adjoint!(gw, gz, grate, gmoments, gshape, d::EpochImages, t, w, z, rate, moments, shape, flux_scale, jitter, buf,
    ::Val{M}, ::Val{NM}, ::Val{NT}) where {M,NM,NT}
    inst = d.instrument
    wins = d.transit_windows[t]
    K, nw = length(rate), length(wins)
    iw, it = d.window_modes, d.transit_modes
    Kw, Kt, P = length(iw), length(it), 5 + M
    ron2, bs2, bf2 = inst.read_noise^2, inst.bad_sigma^2, inst.bad_fraction^2
    shp = ntuple(m -> m <= length(shape) ? Float64(shape[m]) : 0.0, Val(M))
    pb = d.p_bad[t]
    nstate = pb == 0 ? 1 : 2
    @no_escape buf begin
        cvecs = @alloc(NTuple{NT,Float64}, K)
        for k in 1:K
            cvecs[k] = _cvec(view(moments, :, k), Val(NT))
        end
        sigma_t = @alloc(Float64, Kt)
        for (c, a) in enumerate(it)
            sigma_t[c] = (jitter !== nothing && (a == 3 || a == 4)) ? jitter[a-2] : _transit_sigma(inst, a)
        end
        electrons, u = @alloc(Float64, K), @alloc(Float64, P)
        G = @alloc(Float64, P, P, 2)
        A = @alloc(Float64, max(Kw, Kt), max(Kw, Kt))
        Ainv, Y = @alloc(Float64, Kw, Kw, nw, 2), @alloc(Float64, Kw, Kt + 1, nw, 2)
        H, Hbar, pbar = @alloc(Float64, Kt + 1, Kt + 1, 2), @alloc(Float64, Kt + 1, Kt + 1, 2), @alloc(Float64, Kt, 2)
        fill!(H, 0)
        fill!(Hbar, 0)
        fill!(pbar, 0)

        # ---- forward: the Gram matrix of every window, its coefficients eliminated
        J1 = J2 = 0.0
        nvalid = 0
        ok = true
        uP = view(u, 1:P-1)
        for (q, iwin) in enumerate(wins)
            win = d.windows[iwin]
            basis = d.bases[win.basis]
            for k in 1:K
                electrons[k] = flux_scale * rate[k] * win.exposure
            end
            ws, zs = view(w, :, q), view(z, :, q)
            fill!(G, 0)
            @inbounds for j in 0:win.n_ac-1, i in 0:win.n_al-1
                y = d.samples[win.offset+1+i+j*win.n_al]
                isnan(y) && continue
                nvalid += 1
                _sample_vector!(uP, basis, win, i, j, ws, zs, electrons, cvecs, Val(M), Val(NM))
                μ, m = _sample_mean(u, win.background, shp)
                var = max(μ, 0.0) + ron2
                varb = var + bs2 + bf2 * m^2
                u[P] = y - μ
                J1 += log(var)
                J2 += log(varb)
                iv, ivb = 1 / var, 1 / varb
                for b in 1:P
                    ub, ubb = u[b] * iv, u[b] * ivb
                    for a in b:P
                        G[a, b, 1] += u[a] * ub
                        G[a, b, 2] += u[a] * ubb
                    end
                end
            end
            @inbounds for s in 1:2, b in 1:P, a in b+1:P
                G[b, a, s] = G[a, b, s]
            end
            for s in 1:nstate
                l = _eliminate_window_adjoint!(view(H, :, :, s), view(Ainv, :, :, q, s), view(Y, :, :, q, s), view(G, :, :, s), A, d, win)
                ok &= l != -Inf
                s == 1 ? (J1 += l) : (J2 += l)
            end
        end
        if ok
            J1 += _transit_solve_adjoint!(view(Hbar, :, :, 1), view(pbar, :, 1), view(H, :, :, 1), A, sigma_t)
            nstate == 2 && (J2 += _transit_solve_adjoint!(view(Hbar, :, :, 2), view(pbar, :, 2), view(H, :, :, 2), A, sigma_t))
        end
        c = nvalid * log(2π)
        ll_good, ll_bad = -(J1 + c) / 2, -(J2 + c) / 2
        result = if !ok || !isfinite(ll_good) || (nstate == 2 && !isfinite(ll_bad))
            (; ll=-Inf, flux_scale=0.0, jitter=(0.0, 0.0))
        else
            ll = nstate == 1 ? ll_good : _logaddexp(log1p(-pb) + ll_good, log(pb) + ll_bad)
            ρ = nstate == 1 ? (1.0, 0.0) : (exp(log1p(-pb) + ll_good - ll), exp(log(pb) + ll_bad - ll))

            # ---- reverse
            Z, Gbar = @alloc(Float64, P, Kt + 1), @alloc(Float64, P, P, 2)
            ZH = @alloc(Float64, P, Kt + 1)
            ubar, Gu = @alloc(Float64, P), @alloc(Float64, P)
            ebar, cbar = @alloc(Float64, K), @alloc(Float64, NT, K)
            pts = @alloc(NTuple{6,Float64}, K)
            shs, shw, shz = @alloc(NTuple{M,Float64}, K), @alloc(NTuple{M,Float64}, K), @alloc(NTuple{M,Float64}, K)
            stencils = @alloc(Tuple{Int,Int,Float64,Float64}, K)
            fill!(cbar, 0)
            gflux = 0.0
            for (q, iwin) in enumerate(wins)
                win = d.windows[iwin]
                basis = d.bases[win.basis]
                for k in 1:K
                    electrons[k] = flux_scale * rate[k] * win.exposure
                end
                fill!(ebar, 0)
                # Ḡ of this window, per state: Z H̄ Zᵀ plus the window's posterior covariance on its own modes
                for s in 1:nstate
                    fill!(Z, 0)
                    @inbounds for c2 in 1:Kt+1
                        Z[c2 <= Kt ? it[c2] : P, c2] += 1
                        for a in 1:Kw
                            Z[iw[a], c2] -= Y[a, c2, q, s]
                        end
                    end
                    @inbounds for c2 in 1:Kt+1, a in 1:P
                        acc = 0.0
                        for c1 in 1:Kt+1
                            acc += Z[a, c1] * Hbar[c1, c2, s]
                        end
                        ZH[a, c2] = acc
                    end
                    @inbounds for b in 1:P, a in 1:P
                        acc = 0.0
                        for c1 in 1:Kt+1
                            acc += ZH[a, c1] * Z[b, c1]
                        end
                        Gbar[a, b, s] = acc
                    end
                    @inbounds for b in 1:Kw, a in 1:Kw
                        Gbar[iw[a], iw[b], s] += Ainv[a, b, q, s]
                    end
                end
                @inbounds for j in 0:win.n_ac-1, i in 0:win.n_al-1
                    y = d.samples[win.offset+1+i+j*win.n_al]
                    isnan(y) && continue
                    sw_ = win.w0 + i * win.dw + j * win.cw
                    sz_ = win.z0 + j * win.dz + i * win.vz
                    f = fw = fz = 0.0
                    shapes = ntuple(_ -> 0.0, Val(M))
                    for k in 1:K
                        pt, s0, s1, s2, st = _point2(basis, cvecs[k], sw_ - w[k, q], sz_ - z[k, q], Val(M), Val(NM))
                        pts[k], shs[k], shw[k], shz[k], stencils[k] = pt, s0, s1, s2, st
                        e = electrons[k]
                        f += e * pt[1]
                        fw += e * pt[2]
                        fz += e * pt[3]
                        shapes = shapes .+ e .* s0
                    end
                    u[1], u[2], u[3], u[4] = 1.0, f, fw, fz
                    for m in 1:M
                        u[4+m] = shapes[m]
                    end
                    μ, m = _sample_mean(u, win.background, shp)
                    var = max(μ, 0.0) + ron2
                    varb = var + bs2 + bf2 * m^2
                    u[P] = y - μ
                    # adjoints of the sample's vector and of its two variances
                    fill!(ubar, 0)
                    vbar = vbarb = 0.0
                    for s in 1:nstate
                        v = s == 1 ? var : varb
                        quad = 0.0
                        for a in 1:P
                            acc = 0.0
                            for b in 1:P
                                acc += Gbar[a, b, s] * u[b]
                            end
                            Gu[a] = acc
                            quad += acc * u[a]
                        end
                        for a in 1:P
                            ubar[a] -= ρ[s] * Gu[a] / v
                        end
                        vb = -ρ[s] / 2 * (1 / v - quad / v^2)
                        s == 1 ? (vbar = vb) : (vbarb = vb)
                    end
                    μbar = (μ > 0 ? vbar + vbarb : 0.0) - ubar[P]
                    mbar = μbar + 2 * bf2 * m * vbarb
                    ubar[2] += mbar
                    for a in 1:M
                        gshape[a] += mbar * u[4+a]
                        ubar[4+a] += shp[a] * mbar
                    end
                    ub = ntuple(a -> ubar[4+a], Val(M))
                    for k in 1:K
                        pt, e = pts[k], electrons[k]
                        ebar[k] += ubar[2] * pt[1] + ubar[3] * pt[2] + ubar[4] * pt[3] + sum(ub .* shs[k]; init=0.0)
                        gw[k, q] -= e * (ubar[2] * pt[2] + ubar[3] * pt[4] + ubar[4] * pt[5] + sum(ub .* shw[k]; init=0.0))
                        gz[k, q] -= e * (ubar[2] * pt[3] + ubar[3] * pt[5] + ubar[4] * pt[6] + sum(ub .* shz[k]; init=0.0))
                        _scatter_moments!(view(cbar, :, k), basis, stencils[k], (e * ubar[2], e * ubar[3], e * ubar[4]), e .* ub, Val(NM), Val(NT))
                    end
                end
                for k in 1:K
                    grate[k] += ebar[k] * flux_scale * win.exposure
                    gflux += ebar[k] * rate[k] * win.exposure
                end
            end
            for k in 1:K, n in 1:NT-1
                gmoments[n, k] += cbar[n+1, k]
            end
            # the transit's widths: dJ/dp = pbar, p = σ⁻², and d ln L = −ρ/2 dJ
            gj = (0.0, 0.0)
            for (c2, a) in enumerate(it)
                (a == 3 || a == 4) || continue
                g = (ρ[1] * pbar[c2, 1] + ρ[2] * pbar[c2, 2]) / sigma_t[c2]^3
                gj = a == 3 ? (g, gj[2]) : (gj[1], g)
            end
            (; ll, flux_scale=gflux, jitter=gj)
        end
        result
    end
end

# ---------------------------------------------------
# Dual numbers in, a dual number out
# ---------------------------------------------------

@inline _value(x::Real) = Float64(x)
@inline _value(x::ForwardDiff.Dual) = ForwardDiff.value(x)
@inline _add_partials(acc, g, x::Real) = acc
@inline _add_partials(acc, g, x::ForwardDiff.Dual) = acc + g * ForwardDiff.partials(x)

"""
    image_transit_loglike_dual(data, t, w, z, rate, moments; shape = (), flux_scale = 1, jitter = nothing, buf)

The first value of [`image_transit_loglike`](@ref) for a scene of
`ForwardDiff` dual numbers, computed with
[`image_transit_loglike_gradient!`](@ref) in plain floating point and the
chain rule applied once. The same number to rounding, at a cost that does
not grow with the number of partials.
"""
function image_transit_loglike_dual(d::EpochImages, t::Integer, w::AbstractMatrix{D}, z::AbstractMatrix{D}, rate::AbstractVector{D},
    moments::AbstractMatrix{D}; shape=(), flux_scale=1.0, jitter=nothing, buf=Bumper.default_buffer()) where {Tag,N,D<:ForwardDiff.Dual{Tag,Float64,N}}
    M = d.instrument.nshapes
    @no_escape buf begin
        wv, zv, ratev, momv = @alloc(Float64, size(w)...), @alloc(Float64, size(z)...), @alloc(Float64, length(rate)), @alloc(Float64, size(moments)...)
        gw, gz, grate, gmom = @alloc(Float64, size(w)...), @alloc(Float64, size(z)...), @alloc(Float64, length(rate)), @alloc(Float64, size(moments)...)
        gshape = @alloc(Float64, M)
        wv .= ForwardDiff.value.(w)
        zv .= ForwardDiff.value.(z)
        ratev .= ForwardDiff.value.(rate)
        momv .= ForwardDiff.value.(moments)
        fill!(gw, 0)
        fill!(gz, 0)
        fill!(grate, 0)
        fill!(gmom, 0)
        fill!(gshape, 0)
        res = image_transit_loglike_gradient!(gw, gz, grate, gmom, gshape, d, t, wv, zv, ratev, momv; shape=map(_value, Tuple(shape)),
            flux_scale=_value(flux_scale), jitter=jitter === nothing ? nothing : map(_value, jitter), buf)
        part = zero(ForwardDiff.Partials{N,Float64})
        if isfinite(res.ll)
            @inbounds for i in eachindex(w)
                part += gw[i] * ForwardDiff.partials(w[i]) + gz[i] * ForwardDiff.partials(z[i])
            end
            @inbounds for k in eachindex(rate)
                part += grate[k] * ForwardDiff.partials(rate[k])
            end
            @inbounds for i in eachindex(moments)
                part += gmom[i] * ForwardDiff.partials(moments[i])
            end
            for (m, s) in enumerate(Tuple(shape))
                part = _add_partials(part, gshape[m], s)
            end
            part = _add_partials(part, res.flux_scale, flux_scale)
            if jitter !== nothing
                part = _add_partials(_add_partials(part, res.jitter[1], jitter[1]), res.jitter[2], jitter[2])
            end
        end
        ForwardDiff.Dual{Tag}(res.ll, part)
    end
end
