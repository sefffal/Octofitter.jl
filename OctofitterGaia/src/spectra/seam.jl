# ---------------------------------------------------
# Joining two libraries in effective temperature
#
# No single grid of model atmospheres covers every star. Where two meet they
# disagree (by a few hundredths of a magnitude and a few per cent in the
# lines, for PHOENIX and BOSZ at 10 000–12 000 K), and a step in flux at a
# fixed temperature is a step in the posterior. So the join is a range, not
# a point: below it one library, above it the other, and inside it their
# geometric mean with a weight that rises smoothly from 0 to 1.
# ---------------------------------------------------

"""
    SeamedLibrary(cool, hot; teff = (10_000, 12_000))

Two libraries of the same kind (two [`SpectralLibrary`](@ref)s or two
[`RVSLibrary`](@ref)s on the same wavelength bins) used as one: `cool`
below `teff[1]`, `hot` above `teff[2]`, and in between

    log F = (1 - w) log F_cool + w log F_hot,    w = t² (3 - 2t),

with `t` the fraction of the way across the range. The result is
continuous with a continuous first derivative in `teff`. Inside the range a
point must be covered by both libraries.
"""
struct SeamedLibrary{L}
    cool::L
    hot::L
    teff_lo::Float64
    teff_hi::Float64
    function SeamedLibrary(cool::L, hot::L; teff=(10_000.0, 12_000.0)) where {L<:Union{SpectralLibrary,RVSLibrary}}
        teff[1] < teff[2] || throw(ArgumentError("the seam needs teff[1] < teff[2]"))
        _same_sampling(cool, hot) || throw(ArgumentError("the two libraries are on different wavelength bins"))
        (cool.teff[end] >= teff[2] && hot.teff[1] <= teff[1]) ||
            throw(ArgumentError("the seam $(teff) K is not inside both libraries ($(cool.teff[end]) K, $(hot.teff[1]) K)"))
        new{L}(cool, hot, teff[1], teff[2])
    end
end

_same_sampling(a::SpectralLibrary, b::SpectralLibrary) = a.edges == b.edges
_same_sampling(a::RVSLibrary, b::RVSLibrary) = (a.ln_wave0, a.dln_wave, nwave(a)) == (b.ln_wave0, b.dln_wave, nwave(b))

"Either kind of library, single or seamed."
const SEDLibrary = Union{SpectralLibrary,SeamedLibrary{SpectralLibrary}}
const RVSWindowLibrary = Union{RVSLibrary,SeamedLibrary{RVSLibrary}}

libname(lib::Union{SpectralLibrary,RVSLibrary}) = lib.name
libname(lib::SeamedLibrary) = string(lib.cool.name, " below ", lib.teff_lo, " K, ", lib.hot.name, " above ", lib.teff_hi, " K")
bin_edges(lib::SpectralLibrary) = lib.edges
bin_edges(lib::SeamedLibrary{SpectralLibrary}) = lib.cool.edges
nbins(lib::SeamedLibrary{SpectralLibrary}) = nbins(lib.cool)
nwave(lib::SeamedLibrary{RVSLibrary}) = nwave(lib.cool)

Base.show(io::IO, lib::SeamedLibrary) = print(io, "SeamedLibrary(", libname(lib), ")")

@inline function _seam_weight(lib::SeamedLibrary, teff)
    t = (teff - lib.teff_lo) / (lib.teff_hi - lib.teff_lo)
    return t^2 * (3 - 2t)
end

# out ← out^(1 - w) · other^w, the geometric mean
function _geometric_blend!(out, other, w)
    @inbounds for n in eachindex(out, other)
        out[n] = exp((1 - w) * log(out[n]) + w * log(other[n]))
    end
    return out
end

function _seamed!(f!::F, out::AbstractVector, lib::SeamedLibrary, teff, logg, mh) where {F}
    teff <= lib.teff_lo && return f!(out, lib.cool, teff, logg, mh)
    teff >= lib.teff_hi && return f!(out, lib.hot, teff, logg, mh)
    f!(out, lib.cool, teff, logg, mh)
    @no_escape Bumper.default_buffer() begin
        other = @alloc(eltype(out), length(out))
        f!(other, lib.hot, teff, logg, mh)
        _geometric_blend!(out, other, _seam_weight(lib, teff))
        nothing
    end
    return out
end

surface_sed!(out::AbstractVector, lib::SeamedLibrary{SpectralLibrary}, teff, logg, mh) =
    _seamed!(surface_sed!, out, lib, teff, logg, mh)
rvs_surface_spectrum!(out::AbstractVector, lib::SeamedLibrary{RVSLibrary}, teff, logg, mh) =
    _seamed!(rvs_surface_spectrum!, out, lib, teff, logg, mh)
