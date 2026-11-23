# Make `gaiatwin_image_fixture.h5`: a small chromatic image basis of one
# device, in the layout of the instrument product, and what GaiaTwin itself
# computes for two stars in an AF window, for the test of the image model's
# mean. This script is how the committed fixture was produced.
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

include(joinpath(@__DIR__, "..", "..", "products", "build_image_basis_gaiatwin.jl"))

out = ARGS[2]
main(out, 660.0, 1320.0, 12.0, 30.0, ["fov1_ccd5"])

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
println("wrote ", out, " (", round(filesize(out) / 1e3; digits=1), " kB); spectra from ", GaiaTwin.sed_source())
