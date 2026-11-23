# ---------------------------------------------------
# Image block, part 1: the chromatic image basis
#
# The image of a point source in an AF window is linear in its spectrum, so
# it is a sum over a small basis,
#
#     image(Δw, Δz) = Σ_n c_n C_n(Δw, Δz),      c_n = ⟨P_n(ν)⟩
#
# with c_n the photon-weighted Legendre moments of wavenumber over the
# source's spectrum in the G band (`chromatic_moments!`) and C_n fixed
# functions of the offset from the source: the instrument. C_n already
# includes the integration over the sample and the smearing of TDI, so the
# model evaluates it at sample centres. A blend is the sum of its
# components' images, each with its own moments, which a single colour per
# source cannot express.
#
# The C_n are a calibration product, one set per device (field of view, CCD,
# gate, window class). Until DR4 they come from GaiaTwin's wavefront optics
# (`products/build_image_basis_gaiatwin.jl`).
# ---------------------------------------------------

"""
    ImageBasis(path, device)

The chromatic image basis of one device, read from an instrument product.
`terms[n, i, j]` is `C_{n-1}` at along-scan offset `x0 + (i - 1) hx` and
across-scan offset `y0 + (j - 1) hy` [mas] from the source: the fraction of
the light in a sample centred there, for the n-th Legendre moment of the
spectrum. Between nodes it is interpolated bicubically (Catmull–Rom),
continuous with continuous first derivatives; beyond the table it is zero.

`wavenumber_range` [µm⁻¹] is the range the Legendre polynomials are mapped
over; it must be that of the [`SEDFunctionals`](@ref) the moments come from.
"""
struct ImageBasis
    name::String
    x0::Float64
    hx::Float64
    y0::Float64
    hy::Float64
    terms::Array{Float32,3}
    wavenumber_range::Tuple{Float64,Float64}
end

function ImageBasis(path::AbstractString, device::AbstractString)
    h5open(path, "r") do f
        haskey(f, device) || throw(ArgumentError("no device \"$device\" in $path; it has $(join(keys(f), ", "))"))
        g = f[device]
        ν = read(attributes(f)["wavenumber_range"])
        ImageBasis(string(read(attributes(f)["instrument"]), ", ", device), read(g["x0"]), read(g["hx"]),
            read(g["y0"]), read(g["hy"]), read(g["terms"]), (ν[1], ν[2]))
    end
end

"Number of basis terms: the degree of the chromatic expansion plus one."
nterms(b::ImageBasis) = size(b.terms, 1)

function Base.show(io::IO, b::ImageBasis)
    print(io, "ImageBasis(", b.name, "; ", nterms(b), " terms, ±", round(Int, -b.x0), " × ±", round(Int, -b.y0),
        " mas at ", b.hx, " × ", b.hy, " mas)")
end

# Catmull–Rom weights along one uniform axis and the array index of the first
# of the four nodes they apply to; `nothing` when the four are not all inside
# the table.
@inline function _axis_stencil(x, x0, h, n)
    u = (x - x0) / h
    i = floor(Int, u)
    (i >= 1 && i <= n - 3) || return nothing
    return i, _cubic_weights(u - i)
end

"""
    image_terms!(out, basis, Δw, Δz) -> out

Fill `out[n]` with `C_{n-1}` at the offset (`Δw`, `Δz`) [mas] of a sample
centre from the source, along and across scan. Zeros beyond the table.
"""
function image_terms!(out::AbstractVector, b::ImageBasis, Δw, Δz)
    fill!(out, 0)
    sx = _axis_stencil(Δw, b.x0, b.hx, size(b.terms, 2))
    sy = _axis_stencil(Δz, b.y0, b.hy, size(b.terms, 3))
    (sx === nothing || sy === nothing) && return out
    (i, wx), (j, wy) = sx, sy
    C = b.terms
    @inbounds for q in 1:4, p in 1:4
        w = wx[p] * wy[q]
        for n in eachindex(out)
            out[n] += w * C[n, i+p-1, j+q-1]
        end
    end
    return out
end

"""
    image_fraction(basis, moments, Δw, Δz)

Fraction of a source's G-band photons in a sample centred (`Δw`, `Δz`) [mas]
from it, for a source whose chromatic moments ⟨P_n(ν)⟩, n = 1…degree, are
`moments` (the zeroth is one; see [`chromatic_moments!`](@ref)).
"""
function image_fraction(b::ImageBasis, moments::AbstractVector, Δw, Δz)
    length(moments) == nterms(b) - 1 ||
        throw(DimensionMismatch("$(length(moments)) moments for a basis of degree $(nterms(b) - 1)"))
    T = promote_type(eltype(moments), typeof(Δw), typeof(Δz), Float64)
    sx = _axis_stencil(Δw, b.x0, b.hx, size(b.terms, 2))
    sy = _axis_stencil(Δz, b.y0, b.hy, size(b.terms, 3))
    (sx === nothing || sy === nothing) && return zero(T)
    (i, wx), (j, wy) = sx, sy
    C = b.terms
    acc = zero(T)
    @inbounds for q in 1:4, p in 1:4
        c = T(C[1, i+p-1, j+q-1])
        for n in eachindex(moments)
            c += moments[n] * C[n+1, i+p-1, j+q-1]
        end
        acc += wx[p] * wy[q] * c
    end
    return acc
end

"""
    window_mean!(μ, basis, w, z, background, sources) -> μ

Expected electrons in the samples of one window: `μ[k]` for the sample whose
centre is at (`w[k]`, `z[k]`) [mas], the window's along- and across-scan
coordinates. `sources` is a tuple with one `(; w, z, electrons, moments)`
per point of the scene: its position in the same coordinates, the electrons
it delivers in the exposure, and its chromatic moments. `background` is in
electrons per sample.
"""
function window_mean!(μ::AbstractArray, b::ImageBasis, w::AbstractArray, z::AbstractArray, background, sources::Tuple)
    size(μ) == size(w) == size(z) || throw(DimensionMismatch("μ, w and z must have the same size"))
    fill!(μ, background)
    for s in sources
        @inbounds for k in eachindex(μ, w, z)
            μ[k] += s.electrons * image_fraction(b, s.moments, w[k] - s.w, z[k] - s.z)
        end
    end
    return μ
end
