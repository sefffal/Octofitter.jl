# Building the node products

The package loads two kinds of file per spectral library, both on the same
`Teff × log g × [M/H]` nodes. Neither is in the repository; these scripts
make them. They have their own environment (`Project.toml` here).

| Step | Script | Output |
|---|---|---|
| 1 | `fetch_phoenix_aces.jl <out>` | `<out>/staging/aces/`: one small file per model with the spectrum in 0.1 nm bins over 300–1100 nm and the native-resolution cut-out 842–875 nm. Downloads 21 GB, keeps 550 MB; resumable. |
| 1 | `fetch_bosz2024.jl <out>` | `<out>/staging/bosz/`, the same for BOSZ 2024 (solar-scaled, the R = 50 000 files, seven metallicities from −1 to +0.5). Downloads 20 GB, keeps 310 MB. BOSZ is on air wavelengths and gives the Eddington flux; the script converts both. Where a node exists more than once the product is MARCS rather than ATLAS9 (7500–8000 K) and microturbulence 1 km/s for log g ≥ 3.5, 2 km/s below; the other choice of each goes to `<out>/staging/bosz-alt/`, for comparison only. |
| 2 | `build_sed_nodes.jl <staging> <out.h5>` | Low-resolution SEDs (5 nm bins, 320–1100 nm), read by `SpectralLibrary`. |
| 3 | `build_rvs_nodes.jl <staging> <out.h5>` | The RVS window at R = 11 500 on a 2 km/s grid as a principal-component emulator, read by `RVSLibrary`. Prints the truncation error for several numbers of components. |

A new library needs only its own fetch script writing the same staging
layout: `grid.h5` with `sed_edges_angstrom` and `rvs_wave_angstrom` and the
attributes `library`, `reference`, `source` and `resolving_power` (`Inf` for
spectra that are not broadened; the RVS builder convolves with the
difference in quadrature), and per model an HDF5 file with `sed`, `rvs` and
the attributes `PHXTEFF`, `PHXLOGG`, `PHXM_H`, `PHXALPHA`. Fluxes are surface
flux densities per unit wavelength in erg s⁻¹ cm⁻² cm⁻¹ on vacuum
wavelengths. All RVS products must be built on the same wavelength grid
(the defaults) to be compared or joined (`library_mismatch`,
`SeamedLibrary`).

Libraries are not rectangular in the node axes. Both builders fill the
missing nodes of the tensor grid by linear extrapolation along log g, then
Teff, and record which nodes are real in `present`. The package interpolates
only inside cells whose eight corners are real; the filled values give the
cubic stencil its end slopes.

## The image basis

`build_image_basis_gaiatwin.jl <GaiaTwin checkout> <out.h5>` writes the
chromatic image basis that `ImageBasis` reads, for the 2-D AF windows of
both fields of view and the nine AF strips, from GaiaTwin's wavefront
optics. It is the stand-in instrument until the basis is calibrated on DR4
epoch images, and it needs a GaiaTwin checkout (the script sets up its own
temporary environment). 68 MB at the default 4 × 12 mas sampling over
±900 × ±1800 mas.

## Fixtures

The test fixtures are cut from the finished products by
`../test/fixtures/make_fixtures.jl`, and those that hold GaiaTwin's own
output by the `make_gaiatwin_*_fixture.jl` scripts beside it.

The MIST evolution grid is not built here: it is a DataDep of the package,
reduced from the distributed isochrones on first use
(`OctofitterGaia.build_mist_grid`).
