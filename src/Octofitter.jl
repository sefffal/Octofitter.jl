module Octofitter


using Printf
using Tables, TypedTables
using Distributions, DistributionsAD
using Bijectors
using AbstractMCMC
using AdvancedHMC
using NamedTupleTools
using ForwardDiff
using Logging
using Statistics
using StatsBase
using NamedTupleTools
using OrderedCollections
using KernelDensity

# Many users are unfamiliar with Julia, and they want to load their data from CSV.
# We export the CSV package to help them on their journey.
using CSV
export CSV

using Reexport
@reexport using PlanetOrbits

export KernelDensity

# Re-export from TypedTables
export Table, FlexTable

using Base.Threads: @threads
using StaticArrays 
using MCMCChains: MCMCChains, Chains
using Random
using DataDeps
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

const mjup2msol = PlanetOrbits.mjup2msol_IAU

# Re-export the Chains constructor.
export Chains 

include("units.jl")
include("orbit-models.jl")
include("distributions.jl")
include("variables.jl")
include("parameterizations.jl")

include("likelihoods/system.jl")
include("likelihoods/relative-astrometry.jl")
include("likelihoods/photometry.jl")
include("likelihoods/hgca.jl")
include("likelihoods/gaia.jl")
include("likelihoods/gaia-new3.jl")
include("likelihoods/hipparcos.jl")
include("likelihoods/observable.jl")


include("logdensitymodel.jl")
include("initialization.jl")
include("optimization.jl")
include("sampling.jl")

include("analysis.jl")
include("macros.jl")
include("sonora.jl")
include("BHAC.jl")

include("io.jl")
include("io-orbitize.jl")

include("sbc.jl")
include("predictive-distributions.jl")
include("cross-validation.jl")


function octofit_pigeons end

export octofit_pigeons

function __init__()

    # Some octofitter likelihoods can cause errors when nansafe_mode isn't set to true.
    set_preferences!(ForwardDiff, "nansafe_mode" => true)

    # List DataDeps here.
    # These are optional, automatic downloads required for certain
    # functionality.
    register(DataDep("HGCA_eDR3",
        """
        Dataset: Hipparcos Gaia Catalog of Accelerations
        Author: Brandt et al
        License: arXiv.org Non-exclusive license to distribute
        Website: https://arxiv.org/abs/2105.11662

        A catalog by Brandt et al containing cross calibrated proper motions
        between Hipparcos and GAIA eDR3.

        File size: 19MiB
        """,
        # "http://physics.ucsb.edu/~tbrandt/HGCA_vEDR3.fits", # file is now 404
        "https://raw.githubusercontent.com/t-brandt/orvara/master/HGCA_vEDR3.fits",
        "23684d583baaa236775108b360c650e79770a695e16914b1201f290c1826065c"
    ))

    register(DataDep("SonoraBobcatEvoPhot",
        """
        Dataset: Sonora Bobcat: cloud-free, substellar atmosphere models, spectra, photometry, evolution, and chemistry
        Author: Marley et al
        License: Creative Commons Attribution 4.0 International
        Website: https://zenodo.org/record/5063476

        ``Presented here are models for non-irradiated, substellar mass objects belonging to the Sonora model series, described in Marley et al. (2021).''
        
        This download contains just the thermal evolution and photometry tables.

        File size: 1MiB
        """,
        "https://zenodo.org/record/5063476/files/evolution_and_photometery.tar.gz?download=1",
        "2198426d1ca0e410fda7b63c3b7f45f3890a8d9f2fcf0a3a1e36e14185283ca5",
        post_fetch_method=unpack
    ))

    register(DataDep("Whereistheplanet",
        """
        Dataset:     Planet astrometry and orbit fits from whereistheplanet.com
        Author:      Wang et al.
        License:     BSD-3 Clause
        Website:     https://github.com/semaphoreP/whereistheplanet

        File size: 10MiB
        """,
        "https://github.com/semaphoreP/whereistheplanet/archive/refs/heads/master.zip",
        # "c02e7c601dc94d7acd0c58398b518038b036d1507f790f3419b574b39d515197",
        post_fetch_method=unpack
    ))

    register(DataDep("Hipparcos_IAD",
        """
        Dataset: The Hipparcos 2 Catalog -- intermediate astrometry data
        Author: Van Leeuwen and Michalik

        Website: https://www.cosmos.esa.int/web/hipparcos/hipparcos-2

        > Van Leeuwen and Michalik (2021) have provided a human readable version of the IAD of the Java tool in a 
        > zip file [warning: ~350 MB]. The data is provided as one ASCII file per star. An eleven line header gives
        > auxiliary information together with the astrometric reference parameters. The rest of each file is mostly
        > identical to the IAD format of the DVD: it provides the individual observations of this star expressed as
        > abscissa residuals against the astrometric reference solution. This data is currently not available from
        > the ESA legacy archive.

        File size: 332MiB
        """,
        "https://www.cosmos.esa.int/documents/532822/6470227/ResRec_JavaTool_2014.zip/a58ad12e-cffb-f959-0ed5-2ae26899f61a?t=1631109433177&download=true",
        "db850403b396ebfa493a5f457530edfac2c2fab33ad2c8795eb70c0e5a828b59",
        post_fetch_method=unpack
    ))

    return
end

include("precompile.jl")
end
