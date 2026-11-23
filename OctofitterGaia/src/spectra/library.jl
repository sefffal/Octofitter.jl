# ---------------------------------------------------
# Spectral library: low-resolution SEDs on (Teff, log g, [M/H]) nodes
#
# The node product is built once per library by `products/build_sed_nodes.jl`
# and holds log10 of the bin-averaged surface flux density. Everything the
# image and BP/RP blocks need from a spectrum is a linear functional of this
# vector (see `SEDFunctionals`), so a single interpolation per emitter per
# likelihood evaluation serves all of them.
# ---------------------------------------------------

"""
    SpectralLibrary(path)

Low-resolution model spectra on a `Teff × log g × [M/H]` node grid, read
from a node-product file. `surface_sed!` interpolates them.

Fields: `teff` [K], `logg` [cgs], `mh` [dex] node axes; `edges` bin edges
[nm, vacuum]; `log_flux[bin, teff, logg, mh]`, log₁₀ of the surface flux
density [erg s⁻¹ cm⁻² cm⁻¹] averaged over the bin; `present[teff, logg, mh]`,
false where the library has no model and the node was filled by
extrapolation.
"""
struct SpectralLibrary
    name::String
    teff::Vector{Float64}
    logg::Vector{Float64}
    mh::Vector{Float64}
    edges::Vector{Float64}
    log_flux::Array{Float32,4}
    present::Array{Bool,3}
    path::String
end

function SpectralLibrary(path::AbstractString)
    h5open(path, "r") do f
        SpectralLibrary(read(attributes(f)["library"]), read(f["teff"]), read(f["logg"]),
            read(f["mh"]), read(f["edges_nm"]), read(f["log_flux"]),
            read(f["present"]) .!= 0, String(path))
    end
end

function Base.show(io::IO, lib::SpectralLibrary)
    print(io, "SpectralLibrary(", lib.name, "; ", count(lib.present), " models, ", nbins(lib),
        " bins ", lib.edges[1], "–", lib.edges[end], " nm; Teff ", extrema(lib.teff),
        ", log g ", extrema(lib.logg), ", [M/H] ", extrema(lib.mh), ")")
end

"Number of wavelength bins of the library's spectra."
nbins(lib::SpectralLibrary) = size(lib.log_flux, 1)

"""
    surface_sed!(out, lib, teff, logg, mh) -> out

Fill `out` with the surface flux density [erg s⁻¹ cm⁻² cm⁻¹] in each of the
library's bins at (`teff`, `logg`, `mh`), by Catmull–Rom interpolation of
log flux in the three parameters. The result is continuous with continuous
first derivatives, and works with `ForwardDiff` duals.

The point must lie in a cell whose eight corners are models the library
has; otherwise `out` is filled with `NaN`.
"""
function surface_sed!(out::AbstractVector, lib::SpectralLibrary, teff, logg, mh)
    length(out) == nbins(lib) || throw(DimensionMismatch("out has $(length(out)) bins, the library $(nbins(lib))"))
    sT = _stencil(lib.teff, teff)
    sg = _stencil(lib.logg, logg)
    sm = _stencil(lib.mh, mh)
    if sT === nothing || sg === nothing || sm === nothing
        return fill!(out, NaN)
    end
    (iT, wT, _), (ig, wg, _), (im, wm, _) = sT, sg, sm
    # indices 2 and 3 of a stencil are the cell the point is in
    @inbounds for c in (2, 3), b in (2, 3), a in (2, 3)
        lib.present[iT[a], ig[b], im[c]] || return fill!(out, NaN)
    end
    fill!(out, 0)
    lf = lib.log_flux
    @inbounds for c in 1:4, b in 1:4, a in 1:4
        w = wT[a] * wg[b] * wm[c]
        iszero(w) && continue
        i, j, k = iT[a], ig[b], im[c]
        @simd for n in eachindex(out)
            out[n] += w * lf[n, i, j, k]
        end
    end
    @inbounds @simd for n in eachindex(out)
        out[n] = exp10(out[n])
    end
    return out
end
