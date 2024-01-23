module OctofitterRadialVelocity

using Octofitter
using PlanetOrbits
using Tables, TypedTables
using Distributions
using DataDeps
using RecipesBase
using AbstractGPs
# using TemporalGPs
using Distances: Distances
using FITSIO
using StaticArrays
using LinearAlgebra

export jd2mjd, mjd2jd

# Radial Velocity data type
const rv_cols = (:epoch, :rv, :σ_rv)

function rvpostplot end
function rvpostplot! end

"""
    StarAbsoluteRVLikelihood(
        (;inst_idx=1, epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;inst_idx=1, epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;inst_idx=1, epoch=5100.2,  rv=7.90,  σ_rv=.11),
    )

Represents a likelihood function of relative astometry between a host star and a secondary body.
`:epoch` (mjd), `:rv` (km/s), and `:σ_rv` (km/s), and `:inst_idx` are all required.

`:inst_idx` is used to distinguish RV time series between insturments so that they may optionally
be fit with different zero points and jitters.
In addition to the example above, any Tables.jl compatible source can be provided.
"""
struct StarAbsoluteRVLikelihood{TTable<:Table,GP} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_names::Vector{String}
    gaussian_process::GP
    function StarAbsoluteRVLikelihood(observations...; instrument_names=nothing, gaussian_process=nothing)
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
        # Check instrument indexes are contiguous
        if length(unique(table.inst_idx)) != maximum(table.inst_idx)
            error("instrument indexes must run from 1:N without gaps")
        end
        # We sort the data first by instrument index then by epoch to make some later code faster
        m = maximum(table.epoch)
        ii = sortperm(table.inst_idx .* (10m) .+ table.epoch)
        table = table[ii,:]
        if isnothing(instrument_names)
            instrument_names = string.(1:maximum(table.inst_idx))
        end
        return new{typeof(table),typeof(gaussian_process)}(table, instrument_names, gaussian_process)
    end
end
StarAbsoluteRVLikelihood(observations::NamedTuple...;kwargs...) = StarAbsoluteRVLikelihood(observations; kwargs...)
export StarAbsoluteRVLikelihood

"""
Radial velocity likelihood.
"""
function Octofitter.ln_like(
    rvlike::StarAbsoluteRVLikelihood,
    θ_system,
    elements,
    num_epochs::Val{L}=Val(length(rvlike.table))
) where L

    T = typeof(θ_system.M)
    ll = zero(T)

    # Each RV instrument index can have it's own barycentric RV offset and jitter.
    # We support up to 10 insturments. Grab the offsets (if defined) and store in 
    # a tuple. Then we can index by inst_idxs to look up the right offset and jitter.
    barycentric_rv_inst = (
        hasproperty(θ_system, :rv0_1) ? θ_system.rv0_1 : convert(T, NaN),
        hasproperty(θ_system, :rv0_2) ? θ_system.rv0_2 : convert(T, NaN),
        hasproperty(θ_system, :rv0_3) ? θ_system.rv0_3 : convert(T, NaN),
        hasproperty(θ_system, :rv0_4) ? θ_system.rv0_4 : convert(T, NaN),
        hasproperty(θ_system, :rv0_5) ? θ_system.rv0_5 : convert(T, NaN),
        hasproperty(θ_system, :rv0_6) ? θ_system.rv0_6 : convert(T, NaN),
        hasproperty(θ_system, :rv0_7) ? θ_system.rv0_7 : convert(T, NaN),
        hasproperty(θ_system, :rv0_8) ? θ_system.rv0_8 : convert(T, NaN),
        hasproperty(θ_system, :rv0_9) ? θ_system.rv0_9 : convert(T, NaN),
        hasproperty(θ_system, :rv0_10) ? θ_system.rv0_10 : convert(T, NaN),
    )
    jitter_inst = (
        hasproperty(θ_system, :jitter_1) ? θ_system.jitter_1 : convert(T, NaN),
        hasproperty(θ_system, :jitter_2) ? θ_system.jitter_2 : convert(T, NaN),
        hasproperty(θ_system, :jitter_3) ? θ_system.jitter_3 : convert(T, NaN),
        hasproperty(θ_system, :jitter_4) ? θ_system.jitter_4 : convert(T, NaN),
        hasproperty(θ_system, :jitter_5) ? θ_system.jitter_5 : convert(T, NaN),
        hasproperty(θ_system, :jitter_6) ? θ_system.jitter_6 : convert(T, NaN),
        hasproperty(θ_system, :jitter_7) ? θ_system.jitter_7 : convert(T, NaN),
        hasproperty(θ_system, :jitter_8) ? θ_system.jitter_8 : convert(T, NaN),
        hasproperty(θ_system, :jitter_9) ? θ_system.jitter_9 : convert(T, NaN),
        hasproperty(θ_system, :jitter_10) ? θ_system.jitter_10 : convert(T, NaN),
    )

    # Vector of radial velocity of the star at each epoch. Go through and sum up the influence of
    # each planet and put it into here. 
    rv_star_buf_all_inst = zeros(T, (L))
    noise_var_all_inst =  zeros(T, (L))

    # Work through RV data and model one instrument at a time
    for inst_idx in 1:maximum(rvlike.table.inst_idx)
        if isnan(jitter_inst[inst_idx]) || isnan(barycentric_rv_inst[inst_idx])
            @warn("`jitter_$inst_idx` and `rv0_$inst_idx` must be provided (were NaN)", maxlog=1)
        end
        istart = findfirst(==(inst_idx), rvlike.table.inst_idx)
        iend = findlast(==(inst_idx), rvlike.table.inst_idx)
        epochs = vec(@view rvlike.table.epoch[istart:iend])
        σ_rvs = vec(@view rvlike.table.σ_rv[istart:iend])
        rvs = vec(@view rvlike.table.rv[istart:iend])
        istart = istart[1]
        iend = iend[1]
        noise_var = @view noise_var_all_inst[istart:iend]
        rv_star_buf = @view rv_star_buf_all_inst[istart:iend] 
        rv_star_buf .+= rvs
        rv_star_buf .+= barycentric_rv_inst[inst_idx] 
        
        # # Zygote version:
        # rv_star_buf = @view(rv_star_buf_all_inst[istart:iend]) .+ barycentric_rv_inst[inst_idx]

        # rv_star_buf .+= barycentric_rv_inst[inst_idx] 
        # rv_star_buf += barycentric_rv_inst[inst_idx] 

        # Go through all "influences" on the RV signal, and subtract their models.
        # You could consider `rv_star` as the residuals after subtracting these.
        for planet_i in eachindex(elements)
            orbit = elements[planet_i]
            # Need to structarrays orbit if we want this to SIMD
            planet_mass = θ_system.planets[planet_i].mass
            for epoch_i in eachindex(epochs)
                rv_star_buf[epoch_i] -= radvel(orbit, epochs[epoch_i], planet_mass*Octofitter.mjup2msol)
            end
            # Zygote version:
            # rv_star_buf -= radvel.(orbit, epochs, planet_mass*Octofitter.mjup2msol)
        end
        

        noise_var .= σ_rvs.^2 .+ jitter_inst[inst_idx]^2

        if !isnothing(rvlike.gaussian_process)
            # Fit a GP

            local gp
            try
                gp = @inline rvlike.gaussian_process(θ_system)
                
            catch err
                if err isa PosDefException
                    @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    return -Inf
                elseif err isa ArgumentError
                    @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    return -Inf
                else
                    rethrow(err)
                end
            end                

            fx = gp(epochs, noise_var)

        else
            # Don't fit a GP
            fx = MvNormal(Diagonal((noise_var)))
        end

        try
            ll += logpdf(fx, rv_star_buf)
        catch err
            if err isa PosDefException || err isa DomainError
                return -Inf
            else
                rethrow(err)
            end
        end
    end
    
    return ll
end


        # We have inst: Normal(inst jitter)
        # Or we have  : GP(inst parameters)
        # Can we fold these into one interface?

        # The GP has to be fit to all the residuals and then evaluated at each epoch.
        # The Normal(jitter) could be fit to all data as an MVNormal()
        # (epochs, residuals) -> logpdf(Normal(0, (σ_rvs[i]^2 + jitter_inst^2)))


        # Now calculate the likelihoods (Gaussian Proccess over residuals)
    #     Per((len_scale, per_scale)) = 
    #     # with_lengthscale(SqExponentialKernel(), len_scale) ∘ PeriodicTransform(1 / per_scale)
    #     # with_lengthscale(SqExponentialKernel(), len_scale) ∘ PeriodicTransform(1 / per_scale)
    #     with_lengthscale(SqExponentialKernel(), len_scale) ∘ PeriodicTransform(1 / per_scale)
    # # f = GP(PeriodicKernel() * SqExponentialKernel())
    # f = GP(Per([3,1.2]))

        # Way to unify this interface:
        # RVLike(data, noise_models=(
        #   (system, epochs) ->  Normal()
        # ))


# Generate new radial velocity observations for a planet
function Octofitter.generate_from_params(like::StarAbsoluteRVLikelihood, θ_planet, elem::PlanetOrbits.AbstractOrbit)

    # Get epochs and uncertainties from observations
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return StarAbsoluteRVLikelihood(radvel_table)
end




# Generate new radial velocity observations for a star
function Octofitter.generate_from_params(like::StarAbsoluteRVLikelihood, θ_system, orbits::Vector{<:Visual{KepOrbit}})

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = radvel.(reshape(orbits, :, 1), epochs, transpose(planet_masses))
    # TODO: Question: is adding jitter like this appropriate in a generative model? I think so.
    noise = randn(length(epochs)) .* θ_system.jitter
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return StarAbsoluteRVLikelihood(radvel_table)
end

mjd2jd(jd) = jd + 2400000.5
jd2mjd(jd) = jd - 2400000.5


include("harps_rvbank.jl")
include("harps_dr1.jl")
include("hires.jl")
include("lick.jl")
include("radvel.jl")


# Plot recipe for astrometry data
@recipe function f(rv::StarAbsoluteRVLikelihood)
   
    xguide --> "time (mjd)"
    yguide --> "radvel (m/s)"

    multiple_instruments = hasproperty(rv.table,:inst_idx) && 
                           length(unique(rv.table.inst_idx)) > 1
    if !multiple_instruments
        @series begin
            color --> :black
            label := nothing
            seriestype := :scatter
            markersize--> 0
            yerr := rv.table.σ_rv
            rv.table.epoch, rv.table.rv
        end
    else
        for inst_idx in sort(unique(rv.table.inst_idx))
            @series begin
                label := nothing
                seriestype := :scatter
                markersize--> 0
                color-->inst_idx
                markerstrokecolor-->inst_idx
                yerr := rv.table.σ_rv[rv.table.inst_idx.==inst_idx]
                rv.table.epoch[rv.table.inst_idx.==inst_idx], rv.table.rv[rv.table.inst_idx.==inst_idx]
            end
        end
    end


end




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
        "bc32ee889fba6cb0871efe8b4278b0c994cf302c7538aeb6aca630e27bdbb6d8",
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
    
    return
end

end
