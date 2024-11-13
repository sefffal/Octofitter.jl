module OctofitterRadialVelocity

using Octofitter
using PlanetOrbits
using Tables, TypedTables
using Distributions
using DataDeps
using AbstractGPs
# using TemporalGPs
using FITSIO
using StaticArrays
using LinearAlgebra
using Bumper

# mjd2jd = Octofitter.mjd2jd
# jd2mjd = Octofitter.jd2mjd


include("rv-absolute.jl")
include("rv-absolute-margin.jl")
include("rv-relative.jl")

include("data-sources/harps_rvbank.jl")
include("data-sources/harps_dr1.jl")
include("data-sources/hires.jl")
include("data-sources/lick.jl")
include("data-sources/ces.jl")
include("compat/radvel.jl")


rvpostplot = Octofitter.rvpostplot



function __init__()

    register(DataDep("ESOHARPS_DR1_rvs",
        """
        Dataset:     ESO/HARPS Radial Velocities Catalog
        Author:      Barbieri, M.
        License:     CC0-1.0
        Publication: https://arxiv.org/abs/2312.06586
        Website:     https://archive.eso.org/dataset/ADP.2023-12-04T15:16:53.464
        
        The first public data release of the HARPS radial velocities catalog. This data release aims to provide the astronomical community with a catalog of radial velocities obtained with spectroscopic observations acquired from 2003 to 2023 with the High Accuracy Radial Velocity Planet Searcher (HARPS) spectrograph installed at the ESO 3.6m telescope in La Silla Observatory (Chile), and spanning wavelengths from 3800 to 6900 Angstrom. The catalog comprises 289843 observations of 6488 unique astronomical objects.
        Radial velocities reported in this catalog are obtained using the HARPS pipeline, with a typical precision of 0.5 m/s, which is essential for the search and validation of exoplanets. Additionally, independent radial velocities measured on the Hα spectral line are included, with a typical error of around 300 m/s suitable for various astrophysical applications where high precision is not critical. This catalog includes 282294 radial velocities obtained through the HARPS pipeline and 288972 derived from the Hα line, collectively forming a time-series dataset that provides a historical record of measurements for each object.
        Further, each object has been cross-referenced with the SIMBAD astronomical database to ensure accurate identification, enabling users to locate and verify objects with existing records in astronomical literature. Information provided for each object includes: astrometric parameters (coordinates, parallaxes, proper motions, radial velocities), photometric parameters (apparent magnitudes in the visible and near-infrared), spectral types and object classifications.
        
        File size: 160MiB
        """,
        "https://dataportal.eso.org/dataPortal/file/ADP.2023-12-04T15:16:53.464",
        "9cff9058cb126e76eb9841d2e3fe3f385c1ebe386662633f21e7db78d2ba6b14"
    ))

    register(DataDep("HARPS_RVBank",
        """
        Dataset:     A public HARPS radial velocity database corrected for systematic errors
        Author:      Trifonov et al.
        License:     CC0-1.0
        Publication: https://www.aanda.org/articles/aa/full_html/2020/04/aa36686-19/aa36686-19.html
        Website:     https://www2.mpia-hd.mpg.de/homes/trifonov/HARPS_RVBank.html

        A public HARPS radial velocity database corrected for systematic errors. (2020)

        Updated in 2023 to coincide with the ESO/HARPS Radial Velocities Catalog release.
        Publication: https://arxiv.org/abs/2312.06586
        Website:     https://archive.eso.org/dataset/ADP.2023-12-04T15:16:53.464
        
        File size: 38MiB
        """,
        "https://github.com/3fon3fonov/HARPS_RVBank/raw/master/HARPS_RVBank_ver02.csv.zip",
        "9218ebd833f8971dcf304c7a6d9835de1c988dc5faae131f3eb939c7b9682586",
        post_fetch_method=unpack

    ))

    register(DataDep("HIRES_rvs",
        """
        Dataset:     A public HIRES radial velocity database corrected for systematic errors
        Author:      Butler et al.
        License:     
        Publication: https://ui.adsabs.harvard.edu/abs/2017yCat..51530208B/abstract
        Website:     https://ebps.carnegiescience.edu/data/hireskeck-data

        File size: 3.7MiB
        """,
        "https://drive.google.com/uc?id=10xCy8UIH8wUAnNJ8zCN8kfFzazWnw-f_&export=download",
        "ad68c2edb69150318e8d47e34189fe104f2a5194a4fcd363c78c741755893251",
        post_fetch_method=unpack
    ))


    register(DataDep("Lick_rvs",
        """
        Dataset:     The Twenty-Five Year Lick Planet Search
        Author:      Fischer et al.
        License:     
        Publication: https://iopscience.iop.org/article/10.1088/0067-0049/210/1/5#apjs488421t2

        A public Lick radial velocity database.
        
        File size: 780k
        """,
        "https://content.cld.iop.org/journals/0067-0049/210/1/5/revision1/apjs488421t2_mrt.txt?Expires=1698868925&Signature=YyKJ4p64PeQg2sh~VAYj6aXxH8b-lH0F0lS6GF0YP07V7oaZWzM4sthpMRldUE7cHQZbMkwoW0R-Jq2FymIYqIlAnT1-qs-y~JifD1A1WThaBOEP2gl5JGgDOGXXMCLK4VuKM3ZucSUu9TWIb3vbNqrG7l~V9LIs-K2bW~KcM-syfRzJ1YC6TSiej1PHJVhoxN-SUQRAw2lkLVQ-eea30IFOw9RSmYFqrqUQGnwx7fdkbTd5ZSvQ~BmB0HZjsav890rZpEVWlCs8ITLpKab3aEysIptlezpS90boNDi3CR-p7We2M9WfibcsemIa72HH7cZS~S1Ri8QTQra5nTY8eQ__&Key-Pair-Id=KL1D8TIY3N7T8",
    ))


    register(DataDep("CES_rvs",
        """
        Dataset:     The planet search programme at the ESO Coude Echelle spectrometer and HARPS.IV. The search for Jupiter analogues around solar-like stars
        Author:      
        License:     
        Publication: 10.1051/0004-6361/201116551 (https://www.aanda.org/10.1051/0004-6361/201116551)
        The planet search programme at the ESO Coude Echelle spectrometer and HARPS.IV. The search for Jupiter analogues around solar-like stars. (2013)
    
        Extracted file size: 2.3M
        """,
        "http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J/A+A/552/A78",
        "24d1ce0356fb6b62ec6a131c3d9f55ef3aced37660c635d87f33d697e074cdfb",
        post_fetch_method=unpack
    ))

    register(DataDep("UVES_rvs",
        """
        Dataset: A Reanalysis of the UVES M Dwarf Planet Search Program
        Author:  R. P. Butler, H. R. A. Jones, F. Feng, M. Tuomi, G. Anglada-Escudé, and Sandy Keiser
        License:     
        Publication: 10.3847/1538-3881/ab4905 (https://iopscience.iop.org/article/10.3847/1538-3881/ab4905)
        
        The UVES (Ultraviolet and Visible Spectrometer) M Dwarf Planet Search program surveyed 40 M dwarfs and 1 M giant from 2000 through 2007 March. Two of the M dwarfs were double-lined spectroscopic binaries. The 38 single-lined M dwarfs in this survey are among the nearest and brightest M dwarfs. Starting with the reduced 1D spectra provided by the UVES team, we reanalyzed the UVES velocities of Proxima Cen as part of the “Pale Red Dot” program. The velocity rms decreased from 3.6 to 2.3 m s−1. Motivated by this result, we have harvested all of the raw data from the UVES M Dwarf Planet Search from the European Southern Observatory (ESO) archives and have written custom packages to generate 1D spectra from the raw data, and velocities from the 1D spectra. The median improvement in the velocity rms from the new analysis is 1.8 m s−1. Six of the 38 M dwarfs from the original study had a velocity rms &lt; 4 m s−1. In the reanalysis presented here, 22 of these stars have a velocity rms &lt; 4 m s−1. We improve the upper limits on possible planets orbiting these stars by a factor of typically two to three. For many of these M dwarfs, these observations represent the first epoch of high-precision velocity measurements.
    
        """,
        "https://content.cld.iop.org/journals/1538-3881/158/6/251/revision1/ajab4905t2_mrt.txt",
    ))

    # 
    
    return
end

end
