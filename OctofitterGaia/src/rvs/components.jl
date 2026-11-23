# ---------------------------------------------------
# RVS block, part 1: what one emitter contributes to an epoch spectrum
#
# Per likelihood evaluation each emitter's rest-frame spectrum is prepared
# once on the library's fine grid (uniform in ln λ, so uniform in velocity):
# rotationally broadened, with its first and second derivatives in velocity.
# Every transit then only shifts and samples these three arrays.
# ---------------------------------------------------

const C_KMS = C_LIGHT / 1e5

"""
    LogWaveGrid(ln_wave0, dln_wave, n)

A wavelength grid uniform in ln λ: sample `i` (from 1) is at
`exp(ln_wave0 + (i - 1) * dln_wave)` nm, in vacuum. The step in velocity is
`dln_wave * c`.
"""
struct LogWaveGrid
    ln_wave0::Float64
    dln_wave::Float64
    n::Int
end
LogWaveGrid(lib::RVSLibrary) = LogWaveGrid(lib.ln_wave0, lib.dln_wave, nwave(lib))
LogWaveGrid(lib::SeamedLibrary{RVSLibrary}) = LogWaveGrid(lib.cool)

Base.length(g::LogWaveGrid) = g.n
velocity_step(g::LogWaveGrid) = g.dln_wave * C_KMS
wavelengths(g::LogWaveGrid) = exp.(g.ln_wave0 .+ (0:g.n-1) .* g.dln_wave)

# Cumulative distribution of the rotational broadening profile of a
# linearly limb-darkened sphere (Gray 2005, eq. 18.14) in x = Δv / vsini,
# which is non-zero on [-1, 1]. Constant outside, and never evaluated there:
# asin has an infinite slope at ±1 that would turn a zero derivative into NaN.
function _rotation_cdf(x, ε)
    x >= 1 && return one(x)
    x <= -1 && return zero(x)
    return 0.5 + ((1 - ε) * (x * sqrt(1 - x^2) + asin(x)) + (π * ε / 2) * (x - x^3 / 3)) / (π * (1 - ε / 3))
end

"""
    rotational_broaden!(out, spectrum, vsini, grid; limb_darkening = 0.6) -> out

Convolve a spectrum on `grid` with the rotational profile for projected
equatorial velocity `vsini` [km/s]. Each kernel weight is the profile's
integral over one sample, so the result is continuously differentiable in
`vsini` and tends to the input as `vsini → 0` (below half a sample of the
grid it *is* the input). The ends of the grid are extended by their last
value.
"""
function rotational_broaden!(out::AbstractVector, spectrum::AbstractVector, vsini, grid::LogWaveGrid; limb_darkening=0.6)
    h = velocity_step(grid)
    n = length(spectrum)
    if !(vsini > h / 2)                 # the kernel is a single sample
        return copyto!(out, spectrum)
    end
    half = floor(Int, vsini / h + 0.5)
    fill!(out, 0)
    lo = _rotation_cdf((-half - 0.5) * h / vsini, limb_darkening)
    for j in -half:half
        hi = _rotation_cdf((j + 0.5) * h / vsini, limb_darkening)
        w = hi - lo
        lo = hi
        @inbounds for i in 1:n
            out[i] += w * spectrum[clamp(i + j, 1, n)]
        end
    end
    return out
end

"""
    RVSComponent(flux, dflux, d2flux, level)

One emitter's spectrum as the RVS block uses it: `flux` on the fine grid in
the emitter's rest frame, its derivative `dflux` with respect to a velocity
shift [per km/s] (positive velocity moves features to longer wavelengths)
and its second derivative `d2flux` [per km² s⁻²], and `level`, the mean of
`flux`. Only ratios of `flux` between the emitters of a scene matter.
"""
struct RVSComponent{V<:AbstractVector,T}
    flux::V
    dflux::V
    d2flux::V
    level::T
end

"""
    rvs_component!(flux, dflux, d2flux, spectrum, vsini, grid; limb_darkening) -> RVSComponent

Prepare an emitter's rest-frame `spectrum` (on `grid`) for the RVS block:
broaden it by `vsini` [km/s] into `flux` and fill the two derivative arrays.
"""
function rvs_component!(flux::AbstractVector, dflux::AbstractVector, d2flux::AbstractVector,
    spectrum::AbstractVector, vsini, grid::LogWaveGrid; limb_darkening=0.6)
    rotational_broaden!(flux, spectrum, vsini, grid; limb_darkening)
    h = velocity_step(grid)
    n = length(flux)
    total = zero(eltype(flux))
    @inbounds for i in 2:n-1
        dflux[i] = -(flux[i+1] - flux[i-1]) / (2h)
        d2flux[i] = (flux[i+1] - 2flux[i] + flux[i-1]) / h^2
        total += flux[i]
    end
    @inbounds begin
        dflux[1], dflux[n] = dflux[2], dflux[n-1]
        d2flux[1], d2flux[n] = d2flux[2], d2flux[n-1]
        total += flux[1] + flux[n]
    end
    return RVSComponent(flux, dflux, d2flux, total / n)
end

# ---------------------------------------------------
# Velocities and window weights
# ---------------------------------------------------

"""
Along-scan coupling of the RVS wavelength scale [km s⁻¹ mas⁻¹]. RVS is
slitless, so a source displaced along scan from the position the pipeline
assumed is dispersed onto shifted wavelengths: 5 km/s for 34 mas in the
data model's own words. The sign is the one GaiaTwin uses
(Holl et al. 2023, §3.3.1); it is to be confirmed on DR4 data against
`scan_angle_signal_*`.
"""
const RVS_KAPPA = -0.146

"""
    gravitational_redshift(logg, radius)

Gravitational redshift of light leaving a star's surface [km/s], from its
`logg` [cgs] and `radius` [R☉]: `g R / c`, 0.636 km/s for the Sun.
"""
gravitational_redshift(logg, radius) = exp10(logg) * radius * R_SUN_CM / C_LIGHT / 1e5

"""
    GaussianProfile(sigma)

A Gaussian across-scan line-spread function of standard deviation `sigma`
[mas]; calling it gives the cumulative fraction of the light below an
across-scan offset [mas].
"""
struct GaussianProfile
    sigma::Float64
end
(p::GaussianProfile)(ζ) = (1 + erf(ζ / (sqrt(2) * p.sigma))) / 2
