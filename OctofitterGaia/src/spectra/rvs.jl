# ---------------------------------------------------
# Spectral library: the RVS window at the instrument's resolution
#
# The node product is built once per library by `products/build_rvs_nodes.jl`:
# every model, convolved with the RVS line-spread function, on a grid uniform
# in ln λ, compressed to a few principal components of log flux. Between
# models only the component weights are interpolated, so a spectrum at any
# (Teff, log g, [M/H]) costs one small stencil and one matrix–vector product.
#
# The spectra are absolute (surface flux density), so the relative line depth
# of two blended stars follows from their radii and temperatures and is not
# a free parameter.
# ---------------------------------------------------

"""
    RVSLibrary(path)

Model spectra of the RVS window (845–871.5 nm) at the instrument's
resolution, as a principal-component emulator on a `Teff × log g × [M/H]`
node grid. [`rvs_surface_spectrum!`](@ref) evaluates it.

The wavelength grid is uniform in ln λ: sample `i` (from 1) is at
`exp(ln_wave0 + (i - 1) * dln_wave)` nm in vacuum, and `dln_wave * c` is the
step in velocity. `rms_wave` is the rms over the library's models of what
the truncation to `ncomponents` leaves, in log₁₀ flux at each sample: the
emulator's own error, to be carried into the likelihood's variance.
"""
struct RVSLibrary
    name::String
    teff::Vector{Float64}
    logg::Vector{Float64}
    mh::Vector{Float64}
    ln_wave0::Float64
    dln_wave::Float64
    mean::Vector{Float64}            # [wave] log10 flux
    basis::Matrix{Float64}           # [component, wave]
    weights::Array{Float32,4}        # [component, teff, logg, mh]
    present::Array{Bool,3}
    rms_wave::Vector{Float64}
    resolving_power::Float64
    path::String
end

function RVSLibrary(path::AbstractString)
    h5open(path, "r") do f
        a = attributes(f)
        RVSLibrary(read(a["library"]), read(f["teff"]), read(f["logg"]), read(f["mh"]),
            read(f["ln_wave0"]), read(f["dln_wave"]), read(f["mean"]),
            permutedims(Float64.(read(f["basis"]))),
            read(f["weights"]), read(f["present"]) .!= 0, Float64.(read(f["rms_wave"])),
            read(a["resolving_power"]), String(path))
    end
end

nwave(lib::RVSLibrary) = length(lib.mean)
ncomponents(lib::RVSLibrary) = size(lib.basis, 1)

"Vacuum wavelengths [nm] of the library's samples."
wavelengths(lib::RVSLibrary) = exp.(lib.ln_wave0 .+ (0:nwave(lib)-1) .* lib.dln_wave)

function Base.show(io::IO, lib::RVSLibrary)
    λ = wavelengths(lib)
    print(io, "RVSLibrary(", lib.name, "; ", count(lib.present), " models, ", nwave(lib),
        " samples ", round(λ[1]; digits=2), "–", round(λ[end]; digits=2), " nm at ",
        round(lib.dln_wave * C_LIGHT / 1e5; digits=2), " km/s, R = ", lib.resolving_power, ", ",
        ncomponents(lib), " components)")
end

"""
    rvs_weights!(w, lib, teff, logg, mh) -> w

The emulator's component weights at (`teff`, `logg`, `mh`), by Catmull–Rom
interpolation between the library's models: continuous with continuous
first derivatives. Filled with `NaN` unless the point lies in a cell whose
eight corners are models the library has.
"""
function rvs_weights!(w::AbstractVector, lib::RVSLibrary, teff, logg, mh)
    length(w) == ncomponents(lib) || throw(DimensionMismatch("w has $(length(w)) entries, the library $(ncomponents(lib)) components"))
    sT = _stencil(lib.teff, teff)
    sg = _stencil(lib.logg, logg)
    sm = _stencil(lib.mh, mh)
    if sT === nothing || sg === nothing || sm === nothing
        return fill!(w, NaN)
    end
    (iT, wT, _), (ig, wg, _), (im, wm, _) = sT, sg, sm
    @inbounds for c in (2, 3), b in (2, 3), a in (2, 3)
        lib.present[iT[a], ig[b], im[c]] || return fill!(w, NaN)
    end
    fill!(w, 0)
    @inbounds for c in 1:4, b in 1:4, a in 1:4
        x = wT[a] * wg[b] * wm[c]
        iszero(x) && continue
        i, j, k = iT[a], ig[b], im[c]
        @simd for n in eachindex(w)
            w[n] += x * lib.weights[n, i, j, k]
        end
    end
    return w
end

"""
    rvs_surface_spectrum!(out, lib, w) -> out
    rvs_surface_spectrum!(out, lib, teff, logg, mh) -> out

Fill `out` with the surface flux density [erg s⁻¹ cm⁻² cm⁻¹] on the
library's wavelength grid, in the star's rest frame and at the instrument's
resolution, from component weights `w` or directly from the stellar
parameters.
"""
function rvs_surface_spectrum!(out::AbstractVector, lib::RVSLibrary, w::AbstractVector)
    length(out) == nwave(lib) || throw(DimensionMismatch("out has $(length(out)) samples, the library $(nwave(lib))"))
    B = lib.basis
    @inbounds for i in eachindex(out)
        acc = zero(eltype(out)) + lib.mean[i]
        @simd for n in eachindex(w)
            acc += B[n, i] * w[n]
        end
        out[i] = exp10(acc)
    end
    return out
end

function rvs_surface_spectrum!(out::AbstractVector, lib::RVSLibrary, teff, logg, mh)
    @no_escape Bumper.default_buffer() begin
        w = @alloc(eltype(out), ncomponents(lib))
        rvs_weights!(w, lib, teff, logg, mh)
        rvs_surface_spectrum!(out, lib, w)
        nothing             # `out` may itself live in the arena, and must not be the block's value
    end
    return out
end
