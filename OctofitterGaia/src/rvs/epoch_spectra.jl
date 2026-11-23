# ---------------------------------------------------
# RVS block, part 2: the likelihood of epoch spectra
#
# One published epoch spectrum y (961 bins, normalised, barycentric) is
# modelled as
#
#     y = B(θ) c + ε,      c ~ N(c₀, Λ),      ε ~ N(0, K)
#
# with m(λ) the composite of the emitters' spectra, each at its own velocity
# and weight, normalised to a mean near one, and the columns of B
#
#     T_j(λ) m(λ), j = 0…d     the continuum the pipeline divided by
#     ∂m/∂v                     a velocity zero point for this transit
#     σ² ∂²m/∂v²                a fractional change of the LSF width
#
# The coefficients c are marginalised analytically. K is fixed by the data:
# the published errors, the correlation of neighbouring bins (they are
# interpolated from the same detector pixels), and a compact-support kernel
# for what the template library gets wrong. So K is factorised once, when
# the data are loaded, and a likelihood evaluation costs d + 3 banded
# forward substitutions per transit.
#
# A transit is either like that or bad (anything, with a broad density that
# does not depend on the parameters); the two states are summed.
# ---------------------------------------------------

"Wavelengths [nm, vacuum] of the bins of a DR4 `rvs_epoch_spectrum`."
const DR4_RVS_WAVELENGTH = collect(range(846.0, 870.0; length=961))

"""
    RVSInstrument(; kwargs...)

What the RVS block assumes about the instrument and the pipeline. Nothing
here is sampled; the defaults are placeholders to be calibrated on the
population.

  - `resolving_power = 11500`: of the library spectra, for the width mode.
  - `kappa = RVS_KAPPA` [km s⁻¹ mas⁻¹]: along-scan coupling.
  - `ac_halfwidth = 884` [mas], `ac_profile = GaussianProfile(150)`: the
    window's across-scan half width and the cumulative across-scan
    line-spread function; together they give the fraction of an emitter's
    light the window catches.
  - `continuum_degree = 3`, `continuum_sigma = 0.3`: degree of the
    Chebyshev continuum and the prior width of its coefficients about
    (1, 0, 0, …).
  - `zeropoint_sigma = 0.5` [km/s]: prior width of the per-transit velocity
    zero point.
  - `lsf_width_sigma = 0.1`: prior width of the fractional LSF width change.
  - `interp_corr = 0.25`: correlation of the noise in neighbouring bins.
    Each bin is a linear interpolation between two detector pixels; averaged
    over sub-pixel phase that correlates neighbours at 0.25 and nothing
    beyond.
  - `mismatch_amplitude = 0.025`, `mismatch_length = 110` [km/s]: amplitude
    (in normalised flux) and support of the template-mismatch kernel, a
    Wendland function `(1 - r)⁴ (4r + 1)`, which falls to one half at 0.31
    of its support. The defaults are what two libraries differ by on the
    same star (PHOENIX-ACES against BOSZ 2024: 1–3 % rms for dwarfs of
    3500–12 000 K, correlated over a resolution element). The amplitude
    can be given per bin instead, see [`RVSEpochSpectra`](@ref).
  - `bad_sigma = 1.0`: standard deviation per bin of the bad state, about a
    flat spectrum at one.
  - `limb_darkening = 0.6`: linear coefficient of the rotational profile.
"""
Base.@kwdef struct RVSInstrument{P}
    resolving_power::Float64 = 11500.0
    kappa::Float64 = RVS_KAPPA
    ac_halfwidth::Float64 = 5 * 176.8
    ac_profile::P = GaussianProfile(150.0)
    continuum_degree::Int = 3
    continuum_sigma::Float64 = 0.3
    zeropoint_sigma::Float64 = 0.5
    lsf_width_sigma::Float64 = 0.1
    interp_corr::Float64 = 0.25
    mismatch_amplitude::Float64 = 0.025
    mismatch_length::Float64 = 110.0
    bad_sigma::Float64 = 1.0
    limb_darkening::Float64 = 0.6
end

"Standard deviation of the line-spread function in velocity [km/s]."
lsf_sigma(inst::RVSInstrument) = C_KMS / (inst.resolving_power * 2sqrt(2log(2)))

"""
    window_fraction(inst, Δζ)

Fraction of an emitter's light inside the RVS window, for an emitter `Δζ`
[mas] across scan from the window's centre.
"""
window_fraction(inst::RVSInstrument, Δζ) =
    inst.ac_profile(inst.ac_halfwidth - Δζ) - inst.ac_profile(-inst.ac_halfwidth - Δζ)

"""
    rvs_velocity(inst, v, Δη; Δζ = 0, kappa_ac = 0)

The velocity [km/s] at which an emitter's spectrum appears in an epoch
spectrum: its barycentric line-of-sight velocity `v` (orbit, systemic
velocity, gravitational redshift and template zero point included) plus
the shift from its along-scan offset `Δη` [mas] from the position the
pipeline assumed for the wavelength scale, and optionally a term linear in
the across-scan offset `Δζ`.
"""
rvs_velocity(inst::RVSInstrument, v, Δη; Δζ=0, kappa_ac=0) = v + inst.kappa * Δη + kappa_ac * Δζ

# What is kept per transit: the valid bins, the factor of their covariance
# and the whitened data.
struct _RVSTransit
    bins::Vector{Int}
    chol::BandedCholesky
    white::Vector{Float64}       # L⁻¹ y
    sumsq::Float64               # |L⁻¹ y|²
    ll_bad::Float64              # log density of the bad state
end

"""
    RVSEpochSpectra(flux, flux_error, grid; instrument = RVSInstrument(),
                    wavelength = DR4_RVS_WAVELENGTH, p_bad = 1e-3,
                    mismatch = instrument.mismatch_amplitude)

The epoch spectra of one source, ready for the likelihood: `flux[bin, transit]`
and `flux_error[bin, transit]` as published, with `NaN` in bins that carry
no data. `grid` is the fine grid the emitters' spectra are on (see
[`LogWaveGrid`](@ref)). `p_bad` is the prior probability that a transit is
bad, one number or one per transit.

`mismatch` is the amplitude of the template-mismatch term in normalised
flux: one number, one per bin, or `[bin, transit]`. Template errors sit in
the lines, so a profile (see [`library_mismatch`](@ref)) keeps the continuum
at full weight where one number cannot. It must not depend on the spectra
being fitted.

Spectra must be normalised to a level near one. The pipeline does that
except for double-lined transits (`normalisation_method = 0`); divide those
by their median first.
"""
struct RVSEpochSpectra{I<:RVSInstrument}
    flux::Matrix{Float64}
    flux_error::Matrix{Float64}
    wavelength::Vector{Float64}
    grid::LogWaveGrid
    instrument::I
    p_bad::Vector{Float64}
    mismatch::Matrix{Float64}        # [bin, transit]
    position::Vector{Float64}        # of each bin on the fine grid, in samples from 1, at zero velocity
    cheb::Matrix{Float64}            # [bin, degree + 1]
    transits::Vector{_RVSTransit}
end

ntransits(d::RVSEpochSpectra) = length(d.transits)
nmodes(inst::RVSInstrument) = inst.continuum_degree + 3

function Base.show(io::IO, d::RVSEpochSpectra)
    print(io, "RVSEpochSpectra(", ntransits(d), " transits, ", length(d.wavelength), " bins ",
        d.wavelength[1], "–", d.wavelength[end], " nm, ", sum(t -> length(t.bins), d.transits; init=0), " valid samples)")
end

_wendland(r) = r < 1 ? (1 - r)^4 * (4r + 1) : zero(r)

# Covariance of the valid bins of one spectrum, as a dense lower band:
# `band[k + 1, i]` is the covariance of valid bins i and i - k.
function _rvs_covariance(bins, σ, amplitude, wavelength, inst::RVSInstrument; extra_variance=0.0)
    n = length(bins)
    v = C_KMS .* log.(wavelength[bins])
    b = 1
    for i in 1:n, j in i-1:-1:1
        v[i] - v[j] < inst.mismatch_length || break
        b = max(b, i - j)
    end
    A = zeros(n, n)
    for i in 1:n, j in max(1, i - b):i
        c = amplitude[bins[i]] * amplitude[bins[j]] * _wendland((v[i] - v[j]) / inst.mismatch_length)
        if i == j
            c += σ[bins[i]]^2 + extra_variance
        elseif bins[i] - bins[j] == 1
            c += inst.interp_corr * σ[bins[i]] * σ[bins[j]]
        end
        A[i, j] = c
    end
    return A, b
end

function RVSEpochSpectra(flux::AbstractMatrix, flux_error::AbstractMatrix, grid::LogWaveGrid;
    instrument::RVSInstrument=RVSInstrument(), wavelength=DR4_RVS_WAVELENGTH, p_bad=1e-3,
    mismatch=instrument.mismatch_amplitude)
    size(flux) == size(flux_error) || throw(DimensionMismatch("flux and flux_error differ in size"))
    size(flux, 1) == length(wavelength) || throw(DimensionMismatch("flux has $(size(flux, 1)) bins, wavelength $(length(wavelength))"))
    0 <= instrument.interp_corr < 0.5 || throw(ArgumentError("interp_corr must be in [0, 0.5)"))
    nt = size(flux, 2)
    p_bad = p_bad isa Number ? fill(Float64(p_bad), nt) : collect(Float64, p_bad)
    length(p_bad) == nt || throw(DimensionMismatch("p_bad needs one entry per transit"))
    all(p -> 0 <= p < 1, p_bad) || throw(ArgumentError("p_bad must be in [0, 1)"))
    mismatch = mismatch isa Number ? fill(Float64(mismatch), size(flux)) :
               mismatch isa AbstractVector ? repeat(collect(Float64, mismatch), 1, nt) : Matrix{Float64}(mismatch)
    size(mismatch) == size(flux) || throw(DimensionMismatch("mismatch must be a number, one value per bin, or the size of flux"))
    all(>=(0), mismatch) || throw(ArgumentError("mismatch amplitudes must not be negative"))
    wavelength = collect(Float64, wavelength)
    position = (log.(wavelength) .- grid.ln_wave0) ./ grid.dln_wave .+ 1
    x = (2 .* wavelength .- (wavelength[1] + wavelength[end])) ./ (wavelength[end] - wavelength[1])
    d = instrument.continuum_degree
    cheb = ones(length(wavelength), d + 1)
    d >= 1 && (cheb[:, 2] .= x)
    for j in 2:d
        cheb[:, j+1] .= 2 .* x .* cheb[:, j] .- cheb[:, j-1]
    end
    transits = map(1:nt) do t
        y, σ = view(flux, :, t), view(flux_error, :, t)
        bins = findall(i -> isfinite(y[i]) && isfinite(σ[i]) && σ[i] > 0, eachindex(y))
        length(bins) > nmodes(instrument) || throw(ArgumentError("transit $t has $(length(bins)) valid bins; drop it"))
        amp = view(mismatch, :, t)
        A, b = _rvs_covariance(bins, σ, amp, wavelength, instrument)
        chol = BandedCholesky(A, b)
        white = whiten!(y[bins], chol)
        Abad, _ = _rvs_covariance(bins, σ, amp, wavelength, instrument; extra_variance=instrument.bad_sigma^2)
        cbad = BandedCholesky(Abad, b)
        wbad = whiten!(y[bins] .- 1, cbad)
        ll_bad = -(sum(abs2, wbad) + cbad.logdet + length(bins) * log(2π)) / 2
        _RVSTransit(bins, chol, white, sum(abs2, white), ll_bad)
    end
    return RVSEpochSpectra(Matrix{Float64}(flux), Matrix{Float64}(flux_error), wavelength, grid,
        instrument, p_bad, mismatch, position, cheb, transits)
end

"A copy holding only the transits `inds`."
Base.getindex(d::RVSEpochSpectra, inds::AbstractVector{<:Integer}) =
    RVSEpochSpectra(d.flux[:, inds], d.flux_error[:, inds], d.wavelength, d.grid, d.instrument,
        d.p_bad[inds], d.mismatch[:, inds], d.position, d.cheb, d.transits[inds])

# ---------------------------------------------------
# The mean model
# ---------------------------------------------------

# Add `weight` × one component, at velocity `v`, to the composite `m` and its
# two velocity derivatives at the given bins. False if the shift takes a bin
# off the fine grid.
function _add_component!(m, dm, d2m, c::RVSComponent, v, weight, bins, position, grid::LogWaveGrid)
    shift = log1p(v / C_KMS) / grid.dln_wave
    n = length(c.flux)
    lo = position[first(bins)] - shift
    hi = position[last(bins)] - shift
    (lo >= 2 && hi < n - 1) || return false
    F, D, D2 = c.flux, c.dflux, c.d2flux
    @inbounds for (q, p) in enumerate(bins)
        u = position[p] - shift
        i = floor(Int, u)
        w0, w1, w2, w3 = _cubic_weights(u - i)
        m[q] += weight * (w0 * F[i-1] + w1 * F[i] + w2 * F[i+1] + w3 * F[i+2])
        dm[q] += weight * (w0 * D[i-1] + w1 * D[i] + w2 * D[i+1] + w3 * D[i+2])
        d2m[q] += weight * (w0 * D2[i-1] + w1 * D2[i] + w2 * D2[i+1] + w3 * D2[i+2])
    end
    return true
end

# Fill the design matrix B[valid bin, mode] of transit `t`. False if a
# velocity is out of the range the fine grid covers.
function _rvs_design!(B::AbstractMatrix, d::RVSEpochSpectra, t::Integer,
    components::Tuple, velocities::Tuple, weights::Tuple)
    inst = d.instrument
    bins = d.transits[t].bins
    deg = inst.continuum_degree
    m, dm, d2m = view(B, :, 1), view(B, :, deg + 2), view(B, :, deg + 3)
    fill!(m, 0)
    fill!(dm, 0)
    fill!(d2m, 0)
    level = sum(map((c, w) -> w * c.level, components, weights))
    ok = true
    for k in eachindex(components)
        ok &= _add_component!(m, dm, d2m, components[k], velocities[k], weights[k] / level,
            bins, d.position, d.grid)
    end
    ok || return false
    σ² = lsf_sigma(inst)^2
    @inbounds for q in eachindex(bins)
        d2m[q] *= σ²
        for j in deg:-1:0                       # column 1 is m itself: fill it last
            B[q, j+1] = d.cheb[bins[q], j+1] * m[q]
        end
    end
    return true
end

# Prior precision of mode j; the prior mean is 1 for the first and 0 for the rest.
@inline function _mode_precision(inst::RVSInstrument, j)
    deg = inst.continuum_degree
    σ = j <= deg + 1 ? inst.continuum_sigma : j == deg + 2 ? inst.zeropoint_sigma : inst.lsf_width_sigma
    return 1 / σ^2
end

# ---------------------------------------------------
# Likelihood
# ---------------------------------------------------

"""
    rvs_transit_loglike(data, t, components, velocities, weights; buf) -> (ll, ll_good)

Log-likelihood of epoch spectrum `t` of `data`, given a tuple of
[`RVSComponent`](@ref)s, their velocities [km/s] in this transit (see
[`rvs_velocity`](@ref)) and their weights: the fraction of each emitter's
light that entered the spectrum (window acceptance, times what deblending
left). `ll` is the sum over the good and bad states, `ll_good` the good
state alone.

Works with dual numbers in any argument.
"""
function rvs_transit_loglike(d::RVSEpochSpectra, t::Integer, components::Tuple, velocities::Tuple,
    weights::Tuple; buf=Bumper.default_buffer())
    T = promote_type(map(c -> eltype(c.flux), components)..., map(typeof, velocities)..., map(typeof, weights)..., Float64)
    inst = d.instrument
    tr = d.transits[t]
    n, J = length(tr.bins), nmodes(inst)
    ll_good = @no_escape buf begin
        B = @alloc(T, n, J)
        A = @alloc(T, J, J)
        h = @alloc(T, J)
        if !_rvs_design!(B, d, t, components, velocities, weights)
            T(-Inf)
        else
            @inbounds for j in 1:J
                col = view(B, :, j)
                whiten!(col, tr.chol)
                h[j] = _dot(tr.white, col)
            end
            # A = Λ⁻¹ + ZᵀZ, h = Zᵀ(ỹ − Z c₀), with c₀ the first unit vector
            @inbounds for j in 1:J, i in j:J
                acc = zero(T)
                @simd for q in 1:n
                    acc += B[q, i] * B[q, j]
                end
                A[i, j] = acc
            end
            g11 = A[1, 1]
            rr = tr.sumsq - 2h[1] + g11
            logprec = 0.0
            @inbounds for j in 1:J
                h[j] -= A[j, 1]
                p = _mode_precision(inst, j)
                A[j, j] += p
                logprec += log(p)
            end
            logdetA = _cholesky_small!(A)
            if logdetA == -Inf
                T(-Inf)
            else
                _forward_small!(h, A)
                quad = rr - sum(abs2, h)
                -(quad + tr.chol.logdet + logdetA - logprec + n * log(2π)) / 2
            end
        end
    end
    p = d.p_bad[t]
    ll = p == 0 ? ll_good : _logaddexp(log1p(-p) + ll_good, log(p) + tr.ll_bad)
    return ll, ll_good
end

@inline function _logaddexp(a, b)
    hi, lo = a > b ? (a, b) : (b, a)
    return hi == -Inf ? hi : hi + log1p(exp(lo - hi))
end

"""
    rvs_transit_model(data, t, components, velocities, weights) -> Vector

The expected epoch spectrum of transit `t` at the prior mean of the
marginalised modes (unit continuum, no zero-point or width change), on the
data's bins; `NaN` where the transit has no data.
"""
function rvs_transit_model(d::RVSEpochSpectra, t::Integer, components::Tuple, velocities::Tuple, weights::Tuple)
    bins = d.transits[t].bins
    B = zeros(length(bins), nmodes(d.instrument))
    _rvs_design!(B, d, t, components, velocities, weights) || error("a velocity is outside the range of the spectral grid")
    out = fill(NaN, length(d.wavelength))
    out[bins] .= view(B, :, 1)
    return out
end

"""
    rvs_transit_simulate(rng, data, t, components, velocities, weights; add_noise = true) -> Vector

Draw epoch spectrum `t` from the good state of the model: mode coefficients
from their prior, noise with the covariance of the data. With
`add_noise = false`, the same as [`rvs_transit_model`](@ref).
"""
function rvs_transit_simulate(rng, d::RVSEpochSpectra, t::Integer, components::Tuple, velocities::Tuple,
    weights::Tuple; add_noise::Bool=true)
    add_noise || return rvs_transit_model(d, t, components, velocities, weights)
    inst = d.instrument
    tr = d.transits[t]
    J = nmodes(inst)
    B = zeros(length(tr.bins), J)
    _rvs_design!(B, d, t, components, velocities, weights) || error("a velocity is outside the range of the spectral grid")
    c = [(j == 1) + randn(rng) / sqrt(_mode_precision(inst, j)) for j in 1:J]
    y = B * c .+ colour!(randn(rng, length(tr.bins)), tr.chol)
    out = fill(NaN, length(d.wavelength))
    out[tr.bins] .= y
    return out
end

# ---------------------------------------------------
# Sizing the mismatch term from two libraries
# ---------------------------------------------------

"""
    library_mismatch(lib, other, teff, logg, mh; instrument = RVSInstrument(),
                     wavelength = DR4_RVS_WAVELENGTH, velocity = 0.0, smooth = 60.0, floor = 0.003)
        -> (; amplitude, velocity, residual)

What two spectral libraries disagree on for a star of the given parameters,
as the likelihood would see it: `other`'s spectrum is treated as data and
`lib`'s is fitted to it with the continuum, the line-spread-function width
and a velocity free. `residual` is what is left, per bin, in normalised
flux; `velocity` [km/s] is the fitted shift, the difference between the two
libraries' velocity zero points; `amplitude` is the rms of the residual in
a running window of half width `smooth` [km/s], not less than `floor`, for
use as the `mismatch` of [`RVSEpochSpectra`](@ref).

The profile is in the frame of a star at `velocity` [km/s]; pass the
catalogue velocity of the source and a `smooth` that covers its excursions.
"""
function library_mismatch(lib::RVSLibrary, other::RVSLibrary, teff, logg, mh; instrument::RVSInstrument=RVSInstrument(),
    wavelength=DR4_RVS_WAVELENGTH, velocity=0.0, smooth=60.0, floor=0.003)
    grid = LogWaveGrid(lib)
    grid == LogWaveGrid(other) || throw(ArgumentError("the two libraries are on different wavelength grids"))
    component(l) = begin
        raw = rvs_surface_spectrum!(zeros(length(grid)), l, teff, logg, mh)
        any(isnan, raw) && throw(ArgumentError("($teff, $logg, $mh) is outside $(l.name)"))
        rvs_component!(similar(raw), similar(raw), similar(raw), raw, 0.0, grid)
    end
    n = length(wavelength)
    blank = RVSEpochSpectra(ones(n, 1), fill(1.0, n, 1), grid; instrument, wavelength, mismatch=0.0)
    y = rvs_transit_model(blank, 1, (component(other),), (velocity,), (1.0,))
    c = component(lib)
    J = nmodes(instrument)
    cols = [1:J-2; J]                         # continuum and width; the velocity is fitted, not linearised
    B = zeros(n, J)
    function residual(v)
        _rvs_design!(B, blank, 1, (c,), (velocity + v,), (1.0,)) || return fill(Inf, n)
        X = B[:, cols]
        return y .- X * (X \ y)
    end
    vs = -8.0:0.25:8.0
    χ = [sum(abs2, residual(v)) for v in vs]
    k = clamp(argmin(χ), 2, length(vs) - 1)
    v = vs[k] + step(vs) / 2 * (χ[k-1] - χ[k+1]) / (χ[k-1] - 2χ[k] + χ[k+1])
    r = residual(v)
    lnw = log.(wavelength)
    amplitude = map(1:n) do i
        acc, m = 0.0, 0
        for j in 1:n
            abs(lnw[j] - lnw[i]) * C_KMS <= smooth || continue
            acc += r[j]^2
            m += 1
        end
        max(sqrt(acc / m), floor)
    end
    return (; amplitude, velocity=v, residual=r)
end
