# Make `gaiatwin_rvs_fixture.h5`: RVS epoch spectra of a double-lined pair
# simulated by GaiaTwin, with the truth that made them. GaiaTwin is not a
# dependency of the package or of its tests; this script is how the
# committed fixture was produced.
#
#   julia make_gaiatwin_rvs_fixture.jl <GaiaTwin checkout> <out.h5>
#
# Contents:
#   ln_wave0 [ln nm], dln_wave, spectra[sample, component]
#                                   GaiaTwin's rest-frame library spectra (PHOENIX-ACES at
#                                   R = 11 500, pixel-integrated, unit mean) on its fine grid
#   velocity[component, transit]    km/s
#   electrons[component, transit]   per CCD, inside the window
#   pixel_wave [nm], composite[pixel, transit]
#                                   `composite_pixels`: the noise-free blend on detector pixels
#   wave [nm], flux[bin, transit], flux_error[bin, transit]
#                                   `dr4_epoch_spectrum`: three CCDs with noise, interpolated to
#                                   the DR4 bins, averaged and normalised

using Pkg
Pkg.activate(; temp=true)
Pkg.develop(; path=ARGS[1])
Pkg.add("HDF5")
using GaiaTwin, HDF5, Random

lib = GaiaTwin.RVSLibrary()
inst = GaiaTwin.RVSInstrument()
teff, logg = (5800.0, 4000.0), (4.5, 4.5)
comps = [GaiaTwin.library_spectrum(lib, teff[k], logg[k]) for k in 1:2]

nt = 16
rng = Xoshiro(20261123)
phase = 2π .* rand(rng, nt)
γ = 23.0
velocity = vcat((γ .+ 31.0 .* sin.(phase))', (γ .- 52.0 .* sin.(phase))')
grvs = (8.5, 10.3)
total(g) = 4.4167 * 10^(-0.4 * (g - 21.317))               # electrons per CCD
electrons = vcat(fill(total(grvs[1]), nt)', (total(grvs[2]) .* (0.85 .+ 0.15 .* rand(rng, nt)))')

pixel_wave = exp.(GaiaTwin.pixel_grid(inst)) ./ 10
composite = reduce(hcat, [GaiaTwin.composite_pixels(lib, inst, comps, velocity[:, t], electrons[:, t]) for t in 1:nt])
epochs = [GaiaTwin.dr4_epoch_spectrum(lib, comps, velocity[:, t], electrons[:, t]; rng) for t in 1:nt]

h5open(ARGS[2], "w") do f
    attributes(f)["source"] = "GaiaTwin.jl, src/rvs (composite_pixels, dr4_epoch_spectrum)"
    attributes(f)["teff"] = collect(teff)
    attributes(f)["logg"] = collect(logg)
    f["ln_wave0"] = lib.lnλ0 - log(10)
    f["dln_wave"] = lib.dlnλ
    f["spectra"] = Float32.(reduce(hcat, comps))
    f["velocity"] = velocity
    f["electrons"] = electrons
    f["pixel_wave"] = pixel_wave
    f["composite"] = composite
    f["wave"] = collect(epochs[1].λ_nm)
    f["flux"] = Float32.(reduce(hcat, [e.flux for e in epochs]))
    f["flux_error"] = Float32.(reduce(hcat, [e.sigma for e in epochs]))
end
println("wrote ", ARGS[2], " (", round(filesize(ARGS[2]) / 1e3; digits=1), " kB)")
