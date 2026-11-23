# ---------------------------------------------------
# Gaia passbands and photon rates
#
# A Gaia passband is a photon response curve S(λ): the fraction of photons
# entering the telescope pupil that end up as photo-electrons. The rate from
# a source with flux density F_λ at the telescope is
#
#     I [e⁻ s⁻¹] = P_A ∫ F_λ S(λ) λ / (h c) dλ
#
# with P_A the pupil area, and its magnitude is -2.5 log10(I) + ZP. Nothing
# else is needed to put a model spectrum on Gaia's scale: no reference
# spectrum, no bolometric correction. The AB zero points published with the
# passbands are reproduced from the curves and P_A to better than a
# millimagnitude, which the tests check.
# ---------------------------------------------------

const H_PLANCK = 6.62607015e-27        # erg s
const C_LIGHT = 2.99792458e10          # cm s⁻¹
const R_SUN_CM = 6.957e10              # IAU 2015 nominal
const PARSEC_CM = 3.0856775814913673e18

"Gaia telescope pupil area, cm² (0.7278 m²; Riello et al. 2021)."
const GAIA_PUPIL_AREA = 0.7278e4

"""
    Passband

A photon response curve on a wavelength grid [nm, vacuum], with the Vega
zero point `zp` [mag] of a rate of 1 e⁻ s⁻¹.
"""
struct Passband
    name::Symbol
    wavelength::Vector{Float64}
    response::Vector{Float64}
    zp::Float64
end
Base.show(io::IO, p::Passband) = print(io, "Passband(:", p.name, ", ",
    first(p.wavelength), "–", last(p.wavelength), " nm, ZP ", p.zp, ")")

_datafile(parts...) = joinpath(dirname(@__DIR__), "data", parts...)

function _read_columns(path)
    rows = [parse.(Float64, split(l)) for l in eachline(path) if !startswith(l, '#') && !isempty(strip(l))]
    return [getindex.(rows, c) for c in 1:length(first(rows))]
end

"""
    gaia_passbands() -> (; G, BP, RP, RVS)

The Gaia (E)DR3 G, BP and RP passbands (Riello et al. 2021) and the DR3
G_RVS passband (Sartoretti et al. 2023), with their VEGAMAG zero points.
"""
function gaia_passbands()
    λ, G, _, BP, _, RP, _ = _read_columns(_datafile("passbands", "GaiaEDR3_passband.dat"))
    zp = parse.(Float64, split(first(eachline(_datafile("passbands", "GaiaEDR3_zeropt.dat"))))[1:6])
    undefined(v) = [x > 90 ? 0.0 : x for x in v]       # 99.99 marks "no response"
    λr, RVS = _read_columns(_datafile("passbands", "GaiaDR3_Grvs_passband.dat"))
    return (;
        G=Passband(:G, λ, undefined(G), zp[1]),
        BP=Passband(:BP, λ, undefined(BP), zp[3]),
        RP=Passband(:RP, λ, undefined(RP), zp[5]),
        RVS=Passband(:RVS, λr ./ 10, RVS, 21.317),
    )
end

# Linear interpolation of the response, zero outside the tabulated range.
function _response_at(p::Passband, λ)
    (λ <= p.wavelength[1] || λ >= p.wavelength[end]) && return 0.0
    i = searchsortedlast(p.wavelength, λ)
    t = (λ - p.wavelength[i]) / (p.wavelength[i+1] - p.wavelength[i])
    return (1 - t) * p.response[i] + t * p.response[i+1]
end

"""
    response_vector(passband, edges_nm; weight = λ -> 1.0) -> Vector

The linear functional that turns a binned spectrum into a photon rate:
`dot(response_vector(p, edges), F)` is the rate in e⁻ s⁻¹ from a spectrum
whose bin-averaged flux density at the telescope is `F`
[erg s⁻¹ cm⁻² cm⁻¹] on bins with the given edges [nm].

`weight(λ)` multiplies the integrand, for moments of the photon
distribution over wavelength.
"""
function response_vector(p::Passband, edges; weight=λ -> 1.0, substep=0.01)
    out = zeros(length(edges) - 1)
    for b in eachindex(out)
        n = max(1, round(Int, (edges[b+1] - edges[b]) / substep))
        h = (edges[b+1] - edges[b]) / n
        acc = 0.0
        for i in 1:n                               # midpoint rule
            λ = edges[b] + (i - 0.5) * h
            acc += _response_at(p, λ) * λ * weight(λ)
        end
        # λ dλ in nm² → cm²
        out[b] = GAIA_PUPIL_AREA / (H_PLANCK * C_LIGHT) * acc * h * 1e-14
    end
    return out
end

"""
    wavenumber_range(passband) -> (ν_lo, ν_hi)

The range of wavenumber [µm⁻¹] over which `passband` carries light: where
its response is above 10⁻⁴ of the peak. The chromatic image model expands
in Legendre polynomials of wavenumber mapped to [-1, 1] over this range, so
the spectra's moments and the instrument's image basis must both use it.
"""
function wavenumber_range(p::Passband)
    support = p.wavelength[p.response.>1e-4*maximum(p.response)]
    return 1e3 / last(support), 1e3 / first(support)
end

"Magnitude of a photon rate [e⁻ s⁻¹] in `passband`'s Vega system."
magnitude(p::Passband, rate) = -2.5 * log10(rate) + p.zp

# ---------------------------------------------------
# Extinction
# ---------------------------------------------------

"""
    extinction_curve(edges_nm) -> Vector

A(λ)/A₀ averaged over each bin, for the Fitzpatrick (1999) curve at
R_V = 3.1. A₀ is the monochromatic extinction at 550 nm, the convention of
the Gaia astrophysical parameters, so a model's `av` is A₀ in magnitudes.
"""
function extinction_curve(edges)
    λ, k = _read_columns(_datafile("extinction", "F99_Rv3.1.dat"))
    (edges[1] >= λ[1] && edges[end] <= λ[end]) ||
        error("extinction curve is tabulated over $(λ[1])–$(λ[end]) nm only")
    at(x) = (i = min(searchsortedlast(λ, x), length(λ) - 1);
    t = (x - λ[i]) / (λ[i+1] - λ[i]);
    (1 - t) * k[i] + t * k[i+1])
    return map(eachindex(edges)[1:end-1]) do b
        n = 20
        h = (edges[b+1] - edges[b]) / n
        sum(at(edges[b] + (i - 0.5) * h) for i in 1:n) / n
    end
end
