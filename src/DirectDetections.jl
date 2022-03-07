module DirectDetections


const __SOURCE_DIR__ = dirname(Base.source_path())
is_developed() = occursin(joinpath(".julia", "dev"), __SOURCE_DIR__)
if is_developed()
    @info "Packaged is dev'd, skipping precompilation and include(\"precompile.jl\")" maxlog=1
    __precompile__(false)
end

using Printf
using ComponentArrays # TODO: remove last uses of component arrays
using Tables, TypedTables
using Distributions
using Bijectors
using AbstractMCMC
using AdvancedHMC
using NamedTupleTools
using ForwardDiff
using Logging
using Statistics
using StatsBase
using NamedTupleTools
using DirectImages
using DirectOrbits

# Re-export these from DirectOrbits
export mjd, KeplerianElements, KeplerianElementsDeg
# Re-export from TypedTables
export Table, FlexTable

using Base.Threads: @threads
using StaticArrays 
using MCMCChains: MCMCChains, Chains
using Random
using DataDeps
using RecipesBase
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

const mjup2msol = 0.0009543

# Re-export the Chains constructor.
export Chains 

include("types.jl")
include("distributions.jl")
include("models.jl")
include("sampling.jl")
include("analysis.jl")
include("macros.jl")
include("hgca.jl")
include("sonora.jl")

function __init__()
    init_plots()

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
        "http://physics.ucsb.edu/~tbrandt/HGCA_vEDR3.fits",
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

    return
end

if !is_developed()
    include("precompile.jl")
end

end
