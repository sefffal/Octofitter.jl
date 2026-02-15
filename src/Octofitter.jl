module Octofitter


using Printf
using Tables, TypedTables
using Distributions, DistributionsAD
using Bijectors
using AbstractMCMC
using AdvancedHMC
using NamedTupleTools
using ForwardDiff
using DifferentiationInterface
using Logging
using Statistics
using StatsBase
using NamedTupleTools
using OrderedCollections
using KernelDensity
using LinearAlgebra

# Many users are unfamiliar with Julia, and they want to load their data from CSV.
# We export the CSV package to help them on their journey.
using CSV
export CSV

using Reexport
@reexport using PlanetOrbits

export KernelDensity

# Re-export from TypedTables
export Table

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
include("macros.jl")

# Helper for checking tables are well-formed
equal_length_cols(tab) = allequal(length(getproperty(tab, col)) for col in Tables.columnnames(tab))

include("likelihoods/system.jl")
include("likelihoods/relative-astrometry.jl")
include("likelihoods/photometry.jl")
include("likelihoods/hgca.jl")
include("likelihoods/hipparcos.jl")
include("likelihoods/hgca-linfit.jl")
include("likelihoods/g23h.jl")
include("likelihoods/gaia-dr4.jl")

include("likelihoods/prior-observable.jl")
include("likelihoods/prior-planet-order.jl")
include("likelihoods/prior-non-crossing.jl")


include("logdensitymodel.jl")
include("initialization.jl")
include("sampling.jl")

include("analysis.jl")
include("sonora.jl")
include("BHAC.jl")

include("io.jl")
include("io-orbitize.jl")

include("sbc.jl")
include("cross-validation.jl")

"""
    using Pigeons
    octofit_pigeons(model; nrounds, n_chains=16, n_chains_variational=16)

Use Pigeons.jl to sample from intractable posterior distributions.
`Pigeons` must be loaded by the user.

```julia
using Pigeons
model = Octofitter.LogDensityModel(System, autodiff=:ForwardDiff, verbosity=4)
chain, pt = octofit_pigeons(model)
```
"""
function octofit_pigeons end

export octofit_pigeons


include_dependency("../Project.toml")
const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const OCTO_VERSION = VersionNumber(TOML.parsefile(joinpath(ROOT_DIR, "Project.toml"))["version"])
const OCTO_VERSION_STR = "v$(string(OCTO_VERSION))"

function __init__()

    if isinteractive() && get(ENV, "CI", "") != "true" 
        @info """\
Welcome to Octofitter $(OCTO_VERSION_STR) 🐙
Check for new releases: https://github.com/sefffal/Octofitter.jl/releases/
Read the documentation: https://sefffal.github.io/Octofitter.jl/$(OCTO_VERSION_STR)
"""
    end

    if isinteractive() && get(ENV, "CI", "") != "true" && Threads.nthreads() == 1
        @info """\
Note: Julia was started with only one thread. Some models may run faster if you supply multiple threads.
To enable multithreading, run:

        julia --threads=auto

or set the environment variable `JULIA_NUM_THREADS=auto` and restart julia.
"""
    end


    # if running on a Mac, prompt users to load AppleAccelerate
    if Sys.isapple()
        apple_accelerate_loaded = 
            any(contains("Accelerate.framework"), getfield.(BLAS.get_config().loaded_libs, :libname))
        if !apple_accelerate_loaded && isinteractive() && get(ENV, "CI", "") != "true" 
            @info """\
Note for macOS users of Octofitter: Some models may run faster if you use Apple's optimized math libraries.
To enable these optimizations, run:

    julia> using AppleAccelerate

Place `using AppleAccelerate` at the start of your script to suppress this message.
            """
        end
    end

    # Some octofitter likelihoods can cause errors when nansafe_mode isn't set to true.
    nansafe_mode = load_preference(ForwardDiff, "nansafe_mode")
    if isnothing(nansafe_mode) || !nansafe_mode
        set_preferences!(ForwardDiff, "nansafe_mode" => true)
    end

    # JPL Horizons seems to allow their SSL certifcates to expire all the time.
    # We have to disable certificate checking in order to be able to query it
    # consistently
    if !haskey(ENV, "JULIA_NO_VERIFY_HOSTS")
        ENV["JULIA_NO_VERIFY_HOSTS"] = "**.nasa.gov"
    end

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

    register(DataDep("HGCA_DR2",
        """
        Dataset: Hipparcos Gaia Catalog of Accelerations
        Author: Brandt et al
        License: arXiv.org Non-exclusive license to distribute
        Website: https://arxiv.org/abs/2105.11662

        A catalog by Brandt et al containing cross calibrated proper motions
        between Hipparcos and GAIA eDR3.

        File size: 19MiB
        """,
        "https://content.cld.iop.org/journals/0067-0049/241/2/39/revision1/apjsab13b2fits.tar.gz",
        "a42b2046ba572c5cd74121a2d9c09b2baba1a63c2ab406b9c673c2544b8ee47a",
        post_fetch_method=unpack
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

    register(DataDep("DE440_Ephemeris",
        """
        Dataset: JPL Development Ephemeris DE440
        Author: NASA Jet Propulsion Laboratory
        License: Public Domain
        Website: https://naif.jpl.nasa.gov/

        High-precision planetary ephemeris from NASA JPL covering years 1550-2650 CE.
        Used for calculating Earth's barycentric position relative to the Solar System Barycenter.
        Replaces online NASA Horizons queries with offline SPICE calculations.

        File size: 128MiB
        """,
        [
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
        ],
        "c340a944068f6ffdb3b2ce755cf736895917b0251a175c4578a36abb9ffdc72e"
    ))

    register(DataDep("BHAC15_GAIA",
        """
        Dataset: BHA15 Isochrone models
        Author: Isabelle Baraffe, Derek Homeier, France Allard, and Gilles Chabrier
        Publication: "New evolutionary models for pre-main sequence and main sequence low-mass stars down to the hydrogen-burning limit"

        """,
        "https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_iso.GAIA",
        "43ba70b5ae87d32fdc2cd8b1346ad705b24c97d17b223510aa8a7bc8d753ab76"
    ))

    register(DataDep("G23H_Catalog",
        """
        Dataset: G23H - Calibrated Gaia DR2, DR3, and Hipparcos Catalog
        Author: Thompson et al.
        License: Creative Commons
        Website: http://dx.doi.org/10.11570/26.0002

        A composite catalog combining calibrated proper motions from Hipparcos, Gaia DR2, and Gaia DR3,
        along with astrometric excess noise calibration and RV variability data from the 'paired' catalog.

        This catalog enables joint modeling of:
        - Hipparcos proper motions and intermediate astrometric data (IAD)
        - Hipparcos-Gaia proper motion anomaly (from HGCA)
        - Calibrated Gaia DR2 proper motions
        - DR3-DR2 scaled position differences
        - Gaia DR3 proper motions
        - Gaia astrometric excess noise (RUWE/UEVA)
        - Gaia RV variability constraints

        For use with GaiaHipparcosUEVAJointObs likelihood objects.

        WARNING: This is a large file (~14 GB). Download may take considerable time.
        For local testing, you can pass a local file path directly to the catalog parameter.

        File format: Apache Arrow (feather)
        File size: ~14 GB
        """,
        "https://www.canfar.net/storage/vault/file/AstroDataCitationDOI/CISTI.CANFAR/26.0002/data/G23H-v1.0.feather",
        # No hash verification for this large file - users should verify integrity manually
    ))

    return
end
end
