"""
    OctofitterGaia

A forward model of the rawest Gaia DR4 products for Octofitter: epoch images
(SM and AF windows), RVS epoch spectra and BP/RP epoch spectra, predicted
from one scene of emitters whose spectra are tied to stellar evolution
models.
"""
module OctofitterGaia

using Octofitter
using PlanetOrbits
using DataDeps
using HDF5
using LinearAlgebra
import Bumper
using Bumper: @no_escape, @alloc
import p7zip_jll

include("interp.jl")
include("passbands.jl")
include("evolution/mist.jl")
include("spectra/library.jl")
include("emitters.jl")

export MISTGrid, SpectralLibrary, Star, SEDFunctionals, BandRate
export gaia_passbands, photon_rate, chromatic_moments!, surface_sed!, sed!, scene_points

function __init__()
    register(DataDep("MIST_v1.2_vvcrit0.4",
        """
        Dataset: MESA Isochrones and Stellar Tracks (MIST) v1.2 isochrones, v/vcrit = 0.4,
                 with UBV(RI)c, 2MASS, Kepler, Hipparcos, Tycho and Gaia photometry
        Author: Choi et al.; Dotter
        Website: https://mist.science
        Publications: Choi et al. 2016, ApJ 823, 102; Dotter 2016, ApJS 222, 8;
                      Paxton et al. 2011, 2013, 2015 (MESA)

        Isochrones for [Fe/H] from -4 to +0.5 and ages from 10^5 to 10^10.3 yr.
        After the download they are reduced to a 40 MB grid file and the
        text tables are deleted.

        File size: 152MiB (1.1GiB while it is unpacked)
        """,
        "https://mist.science/data/tarballs_v1.2/MIST_v1.2_vvcrit0.4_UBVRIplus.txz",
        "c8425286fc0320f52fdbb0c5f7bbed775dd2be91e048dded8ee97e299236fd6f",
        post_fetch_method=_mist_post_fetch
    ))
    return
end

end
