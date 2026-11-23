# Make `gaiatwin_image_fixture.h5`: a small chromatic image basis of one
# device with its shape modes, in the layout of the instrument product, and
# what GaiaTwin itself computes for two stars in AF windows, for the tests
# of the image model. This script is how the committed fixture was produced.
#
#   julia make_gaiatwin_image_fixture.jl <GaiaTwin checkout> <out.h5>
#
# Besides the basis (group `fov1_ccd5`), the group `reference` holds, for an
# 18 × 12 window of one-pixel samples:
#   w, z [sample]                   sample centres, mas (along, across scan)
#   teff, position[2, source]       the two sources and where they are
#   lam [µm], weights[lam, source]  GaiaTwin's G-band photon weights of each
#   image[sample…, source]          GaiaTwin's polychromatic PSF (`psf2d`) at the sample
#                                   centres, from a 1 × 3 mas table: GaiaTwin interpolates
#                                   bilinearly, which on its usual 4 × 12 mas table is
#                                   itself in error by 10⁻³ of the peak
#
# and the group `transits` the windows of the same pair in a few transits,
# as an instrument that is not the calibrated one records them: in each
# transit the optics are displaced from the stand-in by one draw of the
# ensemble the shape modes were built from (`delta[:, transit]`, in standard
# deviations), every window has its own sub-pixel phase, and the samples
# carry photon and read noise:
#   samples[al, ac, ccd, transit]   electrons, background included
#   w0, z0[ccd, transit]            centre of sample (0, 0), mas
#   position[2, source, transit]    where the sources are, mas (along, across scan)
#   rate[source] [e⁻ s⁻¹], exposure [s], background [e⁻ per sample], read_noise [e⁻]

include(joinpath(@__DIR__, "..", "..", "products", "build_image_basis_gaiatwin.jl"))

out = ARGS[2]
main(out, 660.0, 1320.0, 12.0, 30.0, ["fov1_ccd5"]; mode_degree=1)

optics = standin_optics()
teff = [5800.0, 3600.0]
position = [0.0 131.7; 0.0 -260.4]                              # the second star 131.7 mas along and −260.4 mas across scan
w = [(i - 9.5 + 0.31) * GaiaTwin.PIX_AL for i in 1:18, j in 1:12]
z = [(j - 6.5 - 0.17) * GaiaTwin.PIX_AC for i in 1:18, j in 1:12]
lam = first(GaiaTwin.sed_weights(teff[1]))
weights = reduce(hcat, [last(GaiaTwin.sed_weights(T)) for T in teff])
image = zeros(18, 12, 2)
for k in 1:2
    p = GaiaTwin.psf2d(optics; fov=1, ccd=5, teff=teff[k], half_al=800.0, half_ac=1500.0, hx=1.0, hy=3.0)
    for j in 1:12, i in 1:18
        image[i, j, k] = p(w[i, j] - position[1, k], z[i, j] - position[2, k])
    end
end
h5open(out, "r+") do f
    g = create_group(f, "reference")
    attributes(g)["source"] = "GaiaTwin.psf2d on a 1 x 3 mas table, PHOENIX photon weights (" * GaiaTwin.sed_source() * ")"
    g["w"], g["z"], g["teff"], g["position"], g["lam"], g["weights"], g["image"] = w, z, teff, position, collect(lam), weights, image
end

# The pair as six transits of nine windows record it: G = 11.5 and 14.5, gate 9 (1024 lines).
rng = Random.Xoshiro(20261123)
nt, nccd = 6, 9
exposure = 1024 * GaiaTwin.LINE_S
rate = [10^(0.4 * (GaiaTwin.ZP_G - G)) for G in (11.5, 14.5)]
background, read_noise = 3.0, 4.48
delta = randn(rng, 5, nt)
samples = zeros(Float32, 18, 12, nccd, nt)
w0, z0 = zeros(nccd, nt), zeros(nccd, nt)
positions = zeros(2, 2, nt)
for t in 1:nt
    θ = 2π * rand(rng)                                           # the pair's orientation to the scan
    primary = (40 * randn(rng), 120 * randn(rng))
    positions[:, 1, t] .= primary
    positions[:, 2, t] .= primary .+ 180 .* (cos(θ), sin(θ))
    o, extra = perturbed_optics(delta[:, t])
    tables = [GaiaTwin.psf2d(o; fov=1, ccd=5, teff=T, extra, half_al=800.0, half_ac=1500.0, hx=1.0, hy=3.0) for T in teff]
    for k in 1:nccd
        w0[k, t] = primary[1] + (-8.5 + rand(rng) - 0.5) * GaiaTwin.PIX_AL
        z0[k, t] = primary[2] + (-5.5 + rand(rng) - 0.5) * GaiaTwin.PIX_AC
        for j in 0:11, i in 0:17
            x, y = w0[k, t] + i * GaiaTwin.PIX_AL, z0[k, t] + j * GaiaTwin.PIX_AC
            E = background + exposure * sum(rate[s] * tables[s](x - positions[1, s, t], y - positions[2, s, t]) for s in 1:2)
            samples[i+1, j+1, k, t] = E + sqrt(E + read_noise^2) * randn(rng)
        end
    end
end
h5open(out, "r+") do f
    g = create_group(f, "transits")
    attributes(g)["source"] = "GaiaTwin.psf2d on a 1 x 3 mas table; optics displaced per transit by `delta` standard deviations of the mode ensemble"
    g["samples"], g["w0"], g["z0"], g["position"], g["delta"] = samples, w0, z0, positions, delta
    g["rate"], g["exposure"], g["background"], g["read_noise"] = rate, exposure, background, read_noise
    g["sample_size"] = [GaiaTwin.PIX_AL, GaiaTwin.PIX_AC]
end
println("wrote ", out, " (", round(filesize(out) / 1e3; digits=1), " kB); spectra from ", GaiaTwin.sed_source())
