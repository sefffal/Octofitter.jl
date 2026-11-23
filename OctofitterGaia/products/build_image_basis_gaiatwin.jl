# Build the chromatic image basis of the AF windows from GaiaTwin's wavefront
# optics: the stand-in instrument until the basis can be calibrated on DR4
# epoch images.
#
# For each device the pixel-integrated monochromatic image PSF_ν(Δw, Δz) is
# computed at wavenumbers across the G band and expanded, sample offset by
# sample offset, in Legendre polynomials of wavenumber:
#
#     PSF_ν(Δw, Δz) ≈ Σ_n C_n(Δw, Δz) P_n(x(ν)),   x ∈ [-1, 1] over `wavenumber_range`.
#
# The image of a source with photon-weighted moments c_n = ⟨P_n⟩ is then
# Σ_n c_n C_n. The fit is weighted by the G passband, since that is where the
# photons are.
#
# Output: one HDF5 file that `OctofitterGaia.ImageBasis` reads,
#
#   attributes   instrument, wavenumber_range [µm⁻¹], degree, the optics' parameters
#   <device>/x0, hx, y0, hy     the table's first node and steps, mas (along, across scan)
#   <device>/terms[n, i, j]     C_{n-1}
#
# with devices named "fov<f>_ccd<k>" (2-D AF windows: one pixel per sample,
# no across-scan motion smear).
#
# Usage:
#   julia build_image_basis_gaiatwin.jl <GaiaTwin checkout> <out.h5> [half_al=900] [half_ac=1800] [hx=4] [hy=12] [fov1_ccd5 ...]

using Pkg
Pkg.activate(; temp=true)
Pkg.develop(; path=ARGS[1])
Pkg.add("HDF5")
using GaiaTwin, HDF5, LinearAlgebra
include(joinpath(@__DIR__, "..", "src", "passbands.jl"))       # the package's own passband and wavenumber range

const DEGREE = 6
const NNU = 33

legendre(n, x) = n == 0 ? one(x) : n == 1 ? x : ((2n - 1) * x * legendre(n - 1, x) - (n - 1) * legendre(n - 2, x)) / n

# The rung-2 optics used in GaiaTwin's own scripts (scripts/rung2/unequal_pairs.jl).
standin_optics() = GaiaTwin.WFEOptics(wfe=GaiaTwin.WFE(coma_al=12.0), sigma_al=22.4, sigma_ac=290.0)

function main(out, half_al, half_ac, hx, hy, devices)
    G = gaia_passbands().G
    ν_lo, ν_hi = wavenumber_range(G)
    ν = collect(range(ν_lo, ν_hi; length=NNU))
    x = 2 .* (ν .- ν_lo) ./ (ν_hi - ν_lo) .- 1
    weight = [sqrt(max(_response_at(G, 1e3 / v), 0.0)) for v in ν]
    V = [legendre(n, xi) for xi in x, n in 0:DEGREE]             # [wavenumber, term]
    fit = pinv(weight .* V) .* weight'                          # [term, wavenumber]: coefficients from the monochromatic values
    optics = standin_optics()
    h5open(out, "w") do f
        a = attributes(f)
        a["instrument"] = "GaiaTwin WFEOptics (coma_al = 12 nm, sigma_al = 22.4 mas, sigma_ac = 290 mas)"
        a["wavenumber_range"] = [ν_lo, ν_hi]
        a["degree"] = DEGREE
        a["n_wavenumbers"] = NNU
        for dev in devices
            m = match(r"^fov(\d)_ccd(\d)$", dev)
            m === nothing && error("device \"$dev\" is not of the form fov<f>_ccd<k>")
            fov, ccd = parse(Int, m.captures[1]), parse(Int, m.captures[2])
            mono = map(ν) do v
                GaiaTwin.psf2d(optics; fov, ccd, weights=([1 / v], [1.0]), half_al, half_ac, hx, hy)
            end
            nx, ny = size(mono[1].v)
            terms = zeros(Float32, DEGREE + 1, nx, ny)
            for j in 1:ny, i in 1:nx, k in 1:NNU, n in 1:DEGREE+1
                terms[n, i, j] += fit[n, k] * mono[k].v[i, j]
            end
            g = create_group(f, dev)
            g["x0"], g["hx"], g["y0"], g["hy"] = mono[1].x0, mono[1].hx, mono[1].y0, mono[1].hy
            g["terms"] = terms
            println(dev, ": ", nx, " × ", ny, " nodes; C_0 sums to ", round(sum(terms[1, :, :]) * hx * hy / (GaiaTwin.PIX_AL * GaiaTwin.PIX_AC); digits=4)); flush(stdout)
        end
    end
    println("wrote ", out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    arg(i, d) = length(ARGS) >= i ? parse(Float64, ARGS[i]) : d
    devices = length(ARGS) >= 7 ? ARGS[7:end] : ["fov$(f)_ccd$(k)" for f in 1:2 for k in 1:9]
    main(ARGS[2], arg(3, 900.0), arg(4, 1800.0), arg(5, 4.0), arg(6, 12.0), devices)
end
