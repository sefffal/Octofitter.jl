# Data shipped with OctofitterGaia

Small instrument tables. Everything larger (evolution grids, spectral
library node products) is a DataDep.

| File | Content | Source |
|---|---|---|
| `passbands/GaiaEDR3_passband.dat` | G, BP and RP photon response curves with uncertainties, 320–1100 nm in 1 nm steps; `99.99` where a band is undefined | Gaia (E)DR3 passbands, version 2: Riello et al. 2021, A&A 649, A3; https://www.cosmos.esa.int/web/gaia/edr3-passbands |
| `passbands/GaiaEDR3_zeropt.dat` | Zero points of G, BP, RP with uncertainties, VEGAMAG and AB | same |
| `passbands/GaiaDR3_Grvs_passband.dat` | G_RVS response, wavelength in Å; zero point 21.317 ± 0.002 mag (Vega) | Sartoretti et al. 2023, A&A 674, A6; copy served by the SVO Filter Profile Service (`GAIA/GAIA3.Grvs`) |
| `extinction/F99_Rv3.1.dat` | A(λ)/A(550 nm) for the Fitzpatrick (1999) curve at R_V = 3.1, 300–1100 nm in 1 nm steps | tabulated with `dust_extinction` 1.5, `parameter_averages.F99` |
