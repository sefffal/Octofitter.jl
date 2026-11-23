# Building the node products

The package loads two kinds of file per spectral library, both on the same
`Teff × log g × [M/H]` nodes. Neither is in the repository; these scripts
make them. They have their own environment (`Project.toml` here).

| Step | Script | Output |
|---|---|---|
| 1 | `fetch_phoenix_aces.jl <out>` | `<out>/staging/aces/`: one small file per model with the spectrum in 0.1 nm bins over 300–1100 nm and the native-resolution cut-out 842–875 nm. Downloads 21 GB, keeps 550 MB; resumable. |
| 2 | `build_sed_nodes.jl <staging> <out.h5>` | Low-resolution SEDs (5 nm bins, 320–1100 nm), read by `SpectralLibrary`. |
| 3 | `build_rvs_nodes.jl <staging> <out.h5>` | The RVS window at R = 11 500 on a 2 km/s grid as a principal-component emulator, read by `RVSLibrary`. Prints the truncation error for several numbers of components. |

A new library needs only its own fetch script writing the same staging
layout: `grid.h5` with `sed_edges_angstrom` and `rvs_wave_angstrom`, and per
model an HDF5 file with `sed`, `rvs` and the attributes `PHXTEFF`, `PHXLOGG`,
`PHXM_H`, `PHXALPHA`. Fluxes are surface flux densities per unit wavelength in
erg s⁻¹ cm⁻² cm⁻¹ on vacuum wavelengths.

Libraries are not rectangular in the node axes. Both builders fill the
missing nodes of the tensor grid by linear extrapolation along log g, then
Teff, and record which nodes are real in `present`. The package interpolates
only inside cells whose eight corners are real; the filled values give the
cubic stencil its end slopes.

The test fixtures are cut from the finished products by
`../test/fixtures/make_fixtures.jl`.

The MIST evolution grid is not built here: it is a DataDep of the package,
reduced from the distributed isochrones on first use
(`OctofitterGaia.build_mist_grid`).
