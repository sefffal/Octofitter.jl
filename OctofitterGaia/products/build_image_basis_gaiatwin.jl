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
# Shape modes. A calibration is uncertain in a few directions, and the model
# marginalises over them (`ImageInstrument`). Until a template library exists
# they are the principal components of an ensemble of perturbed optics about
# the stand-in: wavefront defocus along and across scan, coma, and the width
# of the detector's smearing in both directions (`MODE_SCATTER`), drawn in
# antithetic pairs. Each member's change of the chromatic terms up to
# `mode_degree` is one row; the components are scaled so that a coefficient
# of one is one standard deviation of the ensemble.
#
# Output: one HDF5 file that `OctofitterGaia.ImageBasis` reads,
#
#   attributes   instrument, wavenumber_range [µm⁻¹], degree, the optics' parameters
#   <device>/x0, hx, y0, hy     the table's first node and steps, mas (along, across scan)
#   <device>/terms[n, i, j]     C_{n-1}
#   <device>/modes[m, n, i, j]  shape mode m: its change of C_{n-1}, n − 1 ≤ mode_degree
#   <device>/mode_variance[m]   share of the ensemble's variance in mode m
#
# with devices named "fov<f>_ccd<k>" (2-D AF windows, no across-scan motion
# smear; one pixel per sample, two across scan for AF1, which is `ccd1`).
#
# Usage:
#   julia build_image_basis_gaiatwin.jl <GaiaTwin checkout> <out.h5> [half_al=900] [half_ac=1800] [hx=4] [hy=12] [fov1_ccd5 ...]
#
# The number of shape modes is `IMAGE_BASIS_NMODES` in the environment (3).

using Pkg
Pkg.activate(; temp=true)
Pkg.develop(; path=ARGS[1])
Pkg.add("HDF5")
using GaiaTwin, HDF5, LinearAlgebra, Random
include(joinpath(@__DIR__, "..", "src", "passbands.jl"))       # the package's own passband and wavenumber range

const DEGREE = 6
const NNU = 33
const NNU_MODES = 9
const MODE_SCATTER = (defocus_al=10.0, defocus_ac=10.0, coma_al=2.0, sigma_al=0.03, sigma_ac=0.03)   # nm rms, nm rms, nm rms, fractions

legendre(n, x) = n == 0 ? one(x) : n == 1 ? x : ((2n - 1) * x * legendre(n - 1, x) - (n - 1) * legendre(n - 2, x)) / n

# The rung-2 optics used in GaiaTwin's own scripts (scripts/rung2/unequal_pairs.jl).
standin_optics() = GaiaTwin.WFEOptics(wfe=GaiaTwin.WFE(coma_al=12.0), sigma_al=22.4, sigma_ac=290.0)

"The stand-in optics displaced by `δ` standard deviations of `MODE_SCATTER`: the optics and the extra wavefront to image them with."
function perturbed_optics(δ)
    o, s = standin_optics(), MODE_SCATTER
    optics = GaiaTwin.WFEOptics(wfe=o.wfe, sigma_al=o.sigma_al * (1 + s.sigma_al * δ[4]), sigma_ac=o.sigma_ac * (1 + s.sigma_ac * δ[5]))
    return optics, GaiaTwin.WFE(defocus_al=s.defocus_al * δ[1], defocus_ac=s.defocus_ac * δ[2], coma_al=s.coma_al * δ[3])
end

# Wavenumbers across the G band and the matrix that turns monochromatic values at them into Legendre coefficients.
function chromatic_fit(G, n, degree)
    ν_lo, ν_hi = wavenumber_range(G)
    ν = collect(range(ν_lo, ν_hi; length=n))
    x = 2 .* (ν .- ν_lo) ./ (ν_hi - ν_lo) .- 1
    weight = [sqrt(max(_response_at(G, 1e3 / v), 0.0)) for v in ν]
    V = [legendre(k, xi) for xi in x, k in 0:degree]            # [wavenumber, term]
    return ν, pinv(weight .* V) .* weight'                      # [term, wavenumber]
end

# The chromatic terms of one device for the given optics: terms[n, i, j].
function chromatic_terms(optics, extra, fov, ccd, ν, fit, half_al, half_ac, hx, hy)
    ac_box = ccd == 1 ? 2 * GaiaTwin.PIX_AC : GaiaTwin.PIX_AC   # AF1 samples are two pixels across scan
    mono = map(ν) do v
        GaiaTwin.psf2d(optics; fov, ccd, extra, weights=([1 / v], [1.0]), ac_box, half_al, half_ac, hx, hy)
    end
    nx, ny = size(mono[1].v)
    terms = zeros(Float64, size(fit, 1), nx, ny)
    for j in 1:ny, i in 1:nx, k in eachindex(ν), n in axes(fit, 1)
        terms[n, i, j] += fit[n, k] * mono[k].v[i, j]
    end
    return terms, mono[1]
end

function main(out, half_al, half_ac, hx, hy, devices; nmodes=3, mode_degree=2, nensemble=16, seed=1)
    G = gaia_passbands().G
    ν_lo, ν_hi = wavenumber_range(G)
    ν, fit = chromatic_fit(G, NNU, DEGREE)
    νm, fitm = chromatic_fit(G, NNU_MODES, mode_degree)
    optics = standin_optics()
    rng = Random.Xoshiro(seed)
    half = [randn(rng, 5) for _ in 1:nensemble÷2]
    members = [half; .-half]
    h5open(out, "w") do f
        a = attributes(f)
        a["instrument"] = "GaiaTwin WFEOptics (coma_al = 12 nm, sigma_al = 22.4 mas, sigma_ac = 290 mas)"
        a["wavenumber_range"] = [ν_lo, ν_hi]
        a["degree"] = DEGREE
        a["n_wavenumbers"] = NNU
        a["mode_degree"] = mode_degree
        a["mode_ensemble"] = "$(length(members)) members, antithetic; rms " * join(("$k = $v" for (k, v) in pairs(MODE_SCATTER)), ", ")
        for dev in devices
            m = match(r"^fov(\d)_ccd(\d)$", dev)
            m === nothing && error("device \"$dev\" is not of the form fov<f>_ccd<k>")
            fov, ccd = parse(Int, m.captures[1]), parse(Int, m.captures[2])
            terms, tab = chromatic_terms(optics, GaiaTwin.WFE(), fov, ccd, ν, fit, half_al, half_ac, hx, hy)
            nx, ny = size(terms, 2), size(terms, 3)
            g = create_group(f, dev)
            g["x0"], g["hx"], g["y0"], g["hy"] = tab.x0, tab.hx, tab.y0, tab.hy
            g["terms"] = Float32.(terms)
            line = string(dev, ": ", nx, " × ", ny, " nodes; C_0 sums to ",
                round(sum(terms[1, :, :]) * hx * hy / (GaiaTwin.PIX_AL * (ccd == 1 ? 2 : 1) * GaiaTwin.PIX_AC); digits=4))
            if nmodes > 0
                nominal, _ = chromatic_terms(optics, GaiaTwin.WFE(), fov, ccd, νm, fitm, half_al, half_ac, hx, hy)
                D = zeros(length(members), length(nominal))
                for (p, δ) in enumerate(members)
                    o, extra = perturbed_optics(δ)
                    D[p, :] .= vec(first(chromatic_terms(o, extra, fov, ccd, νm, fitm, half_al, half_ac, hx, hy))) .- vec(nominal)
                end
                E = eigen(Symmetric(D * D'))                     # ascending
                order = sortperm(E.values; rev=true)[1:nmodes]
                modes = zeros(Float32, nmodes, size(nominal)...)
                for (k, e) in enumerate(order)
                    v = D' * E.vectors[:, e] ./ sqrt(length(members))
                    v .*= sign(v[argmax(abs.(v))])               # a fixed sign, for reproducible products
                    modes[k, :, :, :] .= reshape(v, size(nominal))
                end
                g["modes"] = modes
                g["mode_variance"] = E.values[order] ./ sum(E.values)
                line *= string("; modes hold ", join(round.(E.values[order] ./ sum(E.values); digits=3), ", "), " of the ensemble variance")
            end
            println(line); flush(stdout)
        end
    end
    println("wrote ", out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    arg(i, d) = length(ARGS) >= i ? parse(Float64, ARGS[i]) : d
    devices = length(ARGS) >= 7 ? ARGS[7:end] : ["fov$(f)_ccd$(k)" for f in 1:2 for k in 1:9]
    main(ARGS[2], arg(3, 900.0), arg(4, 1800.0), arg(5, 4.0), arg(6, 12.0), devices; nmodes=parse(Int, get(ENV, "IMAGE_BASIS_NMODES", "3")))
end
