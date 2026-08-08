# ---------------------------------------------------
# The forward model: point sources → complex visibility → closure phases
#
# For point sources at sky offsets (Δα✱_j, Δδ_j) with band fluxes f_j,
#
#     V(u,v) = Σ_j f_j exp(−2πi (u Δα✱_j + v Δδ_j)) / Σ_j f_j
#
# with the offsets measured from an arbitrary phase centre. Shifting that
# centre multiplies V by a global phase, which leaves |V|², closure phases and
# kernel phases invariant — so `ref` is free to be whatever reads best.
#
# v1 wrote the j = host term as a literal `cvis_model .+= 1.0`: the host's flux
# was hard-coded to 1 and its position to the phase centre, and the sum was
# normalized by `1/(1 + Σ contrast)`. Everything downstream then had to speak
# in contrast ratios against a body that was not in the sum. Here the host is
# an ordinary entry in `targets` with an ordinary `flux_<band>`, and setting it
# to 1.0 recovers v1 exactly (`test/runtests.jl` asserts that bit-for-bit).
#
# `mas2rad` is spelled the way v1 spelled it — `π / (180 · 3600 · 1000)`,
# applied *after* multiplying by 2π — because the bit-for-bit check runs
# through this line.
# ---------------------------------------------------

"""
    cvis_bin!(cvis; Δdec, Δra, contrast, u, v)

Accumulate the complex visibility of one point source of brightness
`contrast` at sky offset (`Δra`, `Δdec`) [mas] from the phase centre, over the
baselines `u`, `v` [inverse wavelengths].

Accumulates into `cvis` — it is the caller's job to zero it and to divide by
the total brightness afterwards.
"""
function cvis_bin!(cvis; Δdec, Δra, contrast, u, v)
    l2 = contrast
    for I in eachindex(cvis, u, v)
        arg = -2π * (u[I] * Δra + v[I] * Δdec) * π / (180 * 3600 * 1000)
        if !isfinite(arg)
            cvis[I] = NaN
            continue
        end
        s, c = sincos(arg)
        cvis[I] += l2 * (c + s * im)
    end
    return cvis
end

"""
    _cvis_epoch!(cvis, u, v, offsets, weights)

Fill `cvis` with the normalized complex visibility of the whole source list at
one epoch and wavelength. `offsets` is a tuple of `(Δα✱, Δδ)` pairs [mas] and
`weights` the matching tuple of effective brightnesses (band flux × fibre
throughput). Both are tuples, so the loop over sources unrolls and the
per-source types stay concrete.
"""
@inline function _cvis_epoch!(cvis, u, v, offsets::NTuple{N,Any}, weights::NTuple{N,Any}) where {N}
    fill!(cvis, zero(eltype(cvis)))
    norm = zero(eltype(weights))
    for j in 1:N
        Δra, Δdec = offsets[j]
        cvis_bin!(cvis; Δdec, Δra, contrast=weights[j], u, v)
        norm += weights[j]
    end
    # `x * (1/n)`, not `x / n`: v1 spelled the normalization this way and the
    # two differ in the last bit.
    invnorm = one(norm) / norm
    cvis .*= invnorm
    return cvis
end

"""
    closurephase!(cp_out; vis, index_cps1, index_cps2, index_cps3)

Closure phases [degrees] of the complex visibilities `vis`, combining the
three baselines of each triangle as `φ₁ + φ₂ − φ₃`.
"""
function closurephase!(cp_out; vis::AbstractArray, index_cps1::AbstractArray,
                       index_cps2::AbstractArray, index_cps3::AbstractArray)
    visphi(v) = rad2deg.(rem2pi.(atan.(imag.(v), real.(v)), RoundNearest))

    for i_cp in eachindex(index_cps1, index_cps2, index_cps3)
        for i_obs in axes(cp_out, 2)
            visphi1 = visphi(vis[index_cps1[i_cp], i_obs])
            visphi2 = visphi(vis[index_cps2[i_cp], i_obs])
            visphi3 = visphi(vis[index_cps3[i_cp], i_obs])
            cp_out[i_cp, i_obs] = visphi1 + visphi2 - visphi3
        end
    end
    return cp_out
end

# ---------------------------------------------------
# Single-mode fibre injection
#
# v1 computed the pointing as `f·ρ/(1+f)` — the closed-form photocentre of
# *two* bodies — and then applied the resulting throughput to the companion
# while leaving the host at 1.0. Neither half survives: the pointing is now a
# reference (`fiber_pointing`, `Photocentre(:band)` by default, but any body or
# subset works) and every source's throughput is a function of *its own*
# offset from that point, host included. See the note in `InterferometryObs`'s
# docstring — this changes the fibre-coupling numbers relative to v1.
# ---------------------------------------------------

"""
    fiber_coupling_fraction(theta, lambda_w=2.2e-6)

Injection efficiency of a single-mode fibre for a point source `theta` [mas]
off the fibre axis at wavelength `lambda_w` [m], for an 8 m aperture.

Credit: W. Balmer, D. Blakely, and others.
"""
function fiber_coupling_fraction(theta, lambda_w=2.2e-6)
    D = 8
    x = range(-D * 2, D * 2, length=500)
    y = range(-D * 2, D * 2, length=500)
    r = LinearAlgebra.norm.(x, y')
    m = r .< D / 2
    # arcseconds
    phase = reshape(x, :, 1, 1) ./ lambda_w .* reshape(theta, 1, 1, :) * 1e-3 / (180 / pi * 3600) * 2pi
    w_0 = 0.32D
    field_pup = @. m * exp(1im * phase)
    field_fiber = @. exp(-1 * r^2 / (2 * w_0^2))
    Inj = abs2.(sum(field_pup .* field_fiber, dims=(1, 2))) / abs(sum(m .* field_fiber))^2
    return Inj[1, 1, :]
end

"""
    fiber_coupling_interpolator(eff_wave; sep_mas=0:2:100, n_wavelengths=15)

Pre-tabulate [`fiber_coupling_fraction`](@ref) on a (separation, wavelength)
grid spanning `eff_wave` and wrap it in a bilinear interpolator. Sources
outside the grid get zero throughput.
"""
function fiber_coupling_interpolator(eff_wave; sep_mas=0:2:100, n_wavelengths=15)
    λ_lo, λ_hi = extrema(vec(eff_wave))
    # A single spectral channel gives a degenerate axis; widen it by a part in
    # 10⁴ so the grid is still a grid. The throughput varies over ~10% across
    # a real GRAVITY band, so this is far below the interpolation error.
    if λ_lo == λ_hi
        λ_lo *= (1 - 1e-4)
        λ_hi *= (1 + 1e-4)
    end
    λs = range(λ_lo, λ_hi, length=n_wavelengths)
    grid = stack([fiber_coupling_fraction(sep_mas, λ) for λ in λs])
    return Interpolations.linear_interpolation((sep_mas, λs), grid, extrapolation_bc=0.0)
end
