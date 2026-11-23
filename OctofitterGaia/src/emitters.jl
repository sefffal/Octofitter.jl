# ---------------------------------------------------
# Emitters
#
# The scene of a system is a list of emitters. An emitter is anything that
# returns a set of weighted sky points and a spectrum: a star returns one
# point, an extended source a quadrature of its surface brightness. The
# image, RVS and BP/RP blocks sum over points and never need a special case.
# ---------------------------------------------------

"""
    AbstractEmitter

A source of light in the scene. The interface is

    scene_points(e, θ_body, θ_system, epoch)   # iterable of (Δα*, Δδ, weight)
    sed!(out, e, θ_body, θ_system)             # observed spectrum on the SED bins

`scene_points` gives offsets [mas] from the body's position and weights that
sum to 1. `sed!` fills `out` with the flux density at the telescope
[erg s⁻¹ cm⁻² cm⁻¹] averaged over each bin, for the whole emitter.
`θ_body` is the body's variable namespace and `θ_system` the system's.

An emitter whose light reaches the RVS epoch spectra also defines

    rvs_grid(e)                                # the LogWaveGrid its spectrum is on
    rvs_spectrum!(out, e, θ_body, θ_system)    # rest-frame spectrum on that grid

with the spectrum at the instrument's resolution, before rotational
broadening, as a flux density at the telescope. Only its ratio to the other
emitters' matters (extinction, common to a system, is left out), and all
emitters of a scene must share one grid.
"""
abstract type AbstractEmitter end

"""
    Star(evol, lib)
    Star(evol, lib, rvs)

A star: one point, with the spectrum of the library `lib` at the body's
`teff`, `logg` and `mh`, diluted by `(radius / distance)²` and reddened by
the system's `av`. `rvs` is the [`RVSLibrary`](@ref) its RVS-window spectrum
comes from, needed only when RVS epoch spectra are modelled; it should be
built from the same model atmospheres as `lib`. Either can be a
[`SeamedLibrary`](@ref).

The body must define `teff` [K], `logg` [cgs], `mh` [dex] and `radius` [R☉],
normally derived from its mass through the evolution model `evol`:

    star = \$evol(mass, system.logage, system.feh)
    teff = star.teff
    logg = star.logg
    mh = star.mh
    radius = star.radius

and the system must define `plx` [mas] and `av`, the monochromatic
extinction at 550 nm in magnitudes (`av = 0.0` for none).
"""
struct Star{E,L<:SEDLibrary,R<:Union{Nothing,RVSWindowLibrary}} <: AbstractEmitter
    evol::E
    lib::L
    rvs::R
    extinction::Vector{Float64}     # A(λ)/A₀ on the library's bins
end
Star(evol, lib::SEDLibrary, rvs::Union{Nothing,RVSWindowLibrary}=nothing) =
    Star(evol, lib, rvs, extinction_curve(bin_edges(lib)))

Base.show(io::IO, s::Star) = print(io, "Star(", s.evol, ", ", libname(s.lib), ")")

nbins(s::Star) = nbins(s.lib)

scene_points(::Star, θ_body, θ_system, epoch) = ((0.0, 0.0, 1.0),)

sed!(out, s::Star, θ_body, θ_system) =
    observed_sed!(out, s, θ_body.teff, θ_body.logg, θ_body.mh, θ_body.radius, θ_system.plx, θ_system.av)

rvs_grid(s::Star{<:Any,<:Any,<:RVSWindowLibrary}) = LogWaveGrid(s.rvs)

function rvs_spectrum!(out, s::Star{<:Any,<:Any,<:RVSWindowLibrary}, θ_body, θ_system)
    rvs_surface_spectrum!(out, s.rvs, θ_body.teff, θ_body.logg, θ_body.mh)
    dilution = (θ_body.radius * R_SUN_CM * θ_system.plx / (1000 * PARSEC_CM))^2
    out .*= dilution
    return out
end

"""
    observed_sed!(out, star, teff, logg, mh, radius, plx, av) -> out

The spectrum of a star at the telescope: the library spectrum times
`(radius / distance)²` times `10^(-0.4 · av · A(λ)/A₀)`. `radius` in R☉,
`plx` in mas.
"""
function observed_sed!(out, s::Star, teff, logg, mh, radius, plx, av)
    surface_sed!(out, s.lib, teff, logg, mh)
    dilution = (radius * R_SUN_CM * plx / (1000 * PARSEC_CM))^2
    k = s.extinction
    @inbounds @simd for n in eachindex(out)
        out[n] *= dilution * exp10(-0.4 * av * k[n])
    end
    return out
end

# ---------------------------------------------------
# Linear functionals of the SED
# ---------------------------------------------------

_legendre(n, x) = n == 0 ? one(x) : n == 1 ? x :
                  ((2n - 1) * x * _legendre(n - 1, x) - (n - 1) * _legendre(n - 2, x)) / n

"""
    SEDFunctionals(edges_nm; passbands = gaia_passbands(), degree = 6)

The fixed vectors that turn an observed SED on bins with the given edges
into what the instrument model needs. Each is a linear functional, so each
is additive over the components of a blend:

  - `G`, `BP`, `RP`, `RVS`: photon rates [e⁻ s⁻¹], `dot(f.G, sed)`;
  - `chroma[n + 1, :]`, n = 0…`degree`: the G-band photon rate weighted by
    the Legendre polynomial `P_n` of wavenumber, mapped to [-1, 1] over
    `wavenumber_range` [µm⁻¹]. The ratio to the n = 0 row is the
    photon-weighted moment ⟨P_n(ν)⟩ that the chromatic image model uses.

See [`photon_rate`](@ref) and [`chromatic_moments`](@ref).
"""
struct SEDFunctionals
    edges::Vector{Float64}
    G::Vector{Float64}
    BP::Vector{Float64}
    RP::Vector{Float64}
    RVS::Vector{Float64}
    chroma::Matrix{Float64}
    wavenumber_range::Tuple{Float64,Float64}
    passbands::NamedTuple{(:G, :BP, :RP, :RVS),NTuple{4,Passband}}
end

function SEDFunctionals(edges::AbstractVector; passbands=gaia_passbands(), degree::Integer=6)
    edges = collect(Float64, edges)
    G = passbands.G
    ν_lo, ν_hi = wavenumber_range(G)
    x(λ) = 2 * (1e3 / λ - ν_lo) / (ν_hi - ν_lo) - 1
    chroma = zeros(degree + 1, length(edges) - 1)
    for n in 0:degree
        chroma[n+1, :] .= response_vector(G, edges; weight=λ -> _legendre(n, x(λ)))
    end
    return SEDFunctionals(edges, chroma[1, :], response_vector(passbands.BP, edges),
        response_vector(passbands.RP, edges), response_vector(passbands.RVS, edges),
        chroma, (ν_lo, ν_hi), passbands)
end
SEDFunctionals(lib::SEDLibrary; kwargs...) = SEDFunctionals(bin_edges(lib); kwargs...)

function Base.show(io::IO, f::SEDFunctionals)
    print(io, "SEDFunctionals(", length(f.G), " bins ", f.edges[1], "–", f.edges[end],
        " nm; G, BP, RP, RVS rates; chromatic moments to degree ", size(f.chroma, 1) - 1, ")")
end

"""
    photon_rate(f::SEDFunctionals, band, sed)

Photon rate [e⁻ s⁻¹] of the observed spectrum `sed` in `band`, one of
`:G`, `:BP`, `:RP`, `:RVS`.
"""
photon_rate(f::SEDFunctionals, band::Symbol, sed) = _dot(getfield(f, band), sed)

function _dot(v::AbstractVector{Float64}, sed)
    acc = zero(eltype(sed))
    @inbounds @simd for n in eachindex(v, sed)
        acc += v[n] * sed[n]
    end
    return acc
end

"""
    chromatic_moments!(c, f::SEDFunctionals, sed) -> rate

Fill `c[n]`, n = 1…degree, with the G-band photon-weighted Legendre moments
of wavenumber ⟨P_n(ν)⟩ of `sed`, and return the G photon rate (the n = 0
moment is 1 by construction and is not stored).

The moments of a blend are the rate-weighted means of its components'
moments, which a single colour is not.
"""
function chromatic_moments!(c::AbstractVector, f::SEDFunctionals, sed)
    rate = _dot(f.G, sed)
    @inbounds for n in eachindex(c)
        c[n] = _dot(view(f.chroma, n + 1, :), sed) / rate
    end
    return rate
end

"""
    BandRate(star, band; functionals = SEDFunctionals(star.lib))

A function `(teff, logg, mh, radius, plx, av) -> photon rate [e⁻ s⁻¹]` of a
[`Star`](@ref) in `band` (`:G`, `:BP`, `:RP` or `:RVS`), for a body's
`flux_<band>` variable:

    rate_G = BandRate(star, :G)
    …
    flux_G = \$rate_G(teff, logg, mh, radius, system.plx, system.av)

This is the same number the raw-data likelihood computes from the same
spectrum, so observations that use `flux_G` (a `Photocentre`, for instance)
and the pixel model agree by construction.
"""
struct BandRate{S<:Star}
    star::S
    band::Symbol
    response::Vector{Float64}
end
BandRate(star::Star, band::Symbol; functionals=SEDFunctionals(star.lib)) =
    BandRate(star, band, getfield(functionals, band))

function (r::BandRate)(teff, logg, mh, radius, plx, av)
    T = float(promote_type(typeof(teff), typeof(logg), typeof(mh), typeof(radius), typeof(plx), typeof(av)))
    buf = Bumper.default_buffer()
    @no_escape buf begin
        sed = @alloc(T, nbins(r.star))
        observed_sed!(sed, r.star, teff, logg, mh, radius, plx, av)
        _dot(r.response, sed)
    end
end
