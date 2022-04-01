# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

# General
using StatsBase
using StaticArrays

# Data
using CSV
using JSON
using TypedTables
using DataFrames 

# Plots
using Plots
using PairPlots

# Orbit fitting
using DirectOrbits
using DirectDetections

using ImageTransformations
using CoordinateTransformations

# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------
# Load data in .json file to dictionary 
function loadJSON(filepath::String)
    dict = open(filepath, "r") do f
        JSON.parse(f)
    end
    return dict
end

# Sample system parameters from prior distributions
function drawfrompriors(system::System)
    θ = DirectDetections.sample_priors(system)
    arr2nt = DirectDetections.make_arr2nt(system)
    θnt = arr2nt(θ)
    return θnt
end

# Generate new astrometry observations
function newobs(obs::Astrometry, elem::KeplerianElements, θ_planet)

    # Get epochs and uncertainties from observations
    epochs = obs.table.epoch
    σ_ras = obs.table.σ_ra 
    σ_decs = obs.table.σ_dec

    # Generate now astrometry data
    ras = DirectOrbits.raoff.(elem, epochs)
    decs = DirectOrbits.decoff.(elem, epochs)
    astrometry_table = Table(epoch=epochs, ra=ras, dec=decs, σ_ra=σ_ras, σ_dec=σ_decs)

    return Astrometry(astrometry_table)
end

# Generate new radial velocity observations for a planet
function newobs(obs::RadialVelocity, elem::KeplerianElements, θ_planet)

    # Get epochs and uncertainties from observations
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocity(radvel_table)
end

# Generate new radial velocity observations for a star
function newobs(obs::RadialVelocity, elems::Vector{<:KeplerianElements}, θ_system)

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = DirectOrbits.radvel.(reshape(elems, :, 1), epochs, transpose(planet_masses))
    noise = randn(length(epochs)) .* θ_system.jitter
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocity(radvel_table)
end


# Generate new images
function newobs(obs::Images, elem::KeplerianElements, θ_planet)


    newrows = map(obs.table) do row
        (;band, image, platescale, epoch, psf) = row

        # Generate new astrometry point
        os = DirectOrbits.orbitsolve(elem, epoch)
        ra = DirectOrbits.raoff(os)
        dec = DirectOrbits.decoff(os)

        phot = θ_planet[band]

        dx = ra/platescale
        dy = -dec/platescale
        translation_tform = Translation(dx + mean(axes(psf,1)), dy + mean(axes(psf,2)))
        # TBD if we want to support rotations for handling negative sidelobes.

        psf_positioned = warp(psf, translation_tform, axes(image), fillvalue=0)

        psf_scaled = psf_positioned .* phot ./ maximum(psf)
        
        injected = image .+ psf_scaled

        return merge(row, (;image=injected))
    end

    return Images(newrows)
end

# Need images function

# Generate calibration data
function calibrationsystem(system::System)

    # Get parameter values sampled from priors
    θ_newsystem = drawfrompriors(system)

    # Generate new orbits for each planet in the system
    elements = map(eachindex(system.planets)) do i
        planet = system.planets[i]
        θ_newplanet = θ_newsystem.planets[i]

        if (hasproperty(θ_newplanet, :a) && θ_newplanet.a <= 0) ||
            (hasproperty(θ_newplanet, :e) && !(0 <= θ_newplanet.e < 1))
            out_of_bounds[] = true
        end

        neworbit = DirectDetections.construct_elements(DirectDetections.orbittype(planet), θ_newsystem, θ_newplanet)

        return neworbit
    end

    # Generate new observations for each planet in the system
    newplanets = map(1:length(system.planets)) do i
        planet = system.planets[i]
        elem = elements[i]

        newplanet_obs = map(planet.observations) do obs
            return newobs(obs, elem, planet)
        end
        newplanet = Planet{DirectDetections.orbittype(planet)}(planet.priors, planet.derived, newplanet_obs..., name=planet.name)
    end

    # Generate new observations for the star
    newstar_obs = map(system.observations) do obs
        return newobs(obs, collect(elements), θ_newsystem)
    end

    # Generate new system
    newsystem = System(system.priors, system.derived, newstar_obs..., newplanets..., name=system.name)

    return θ_newsystem, newsystem
end

# Run chains on new model
function calibrationhmc(system::System, target_accept, num_chains, adaptation, iterations, thinning, tree_depth, verbosity)

    # Get parameter values sampled from priors and generate a new system
    θ_newsystem, newsystem = calibrationsystem(system)
    θ_array = DirectDetections.result2mcmcchain(newsystem, [θ_newsystem])

    # Run chains
    @time chains = DirectDetections.hmc(
        newsystem, target_accept,
        num_chains = num_chains,
        adaptation = adaptation,
        iterations = iterations,
        thinning = thinning,
        tree_depth = tree_depth,
        verbosity = verbosity
    )

    # Compute rank statistic of each parameter in the chains
    chainkeys = string.(keys(chains))
    priorsampledict = Dict()
    rdict = Dict()

    for key in chainkeys
        posteriorsamples = vec(chains[key])
        priorsample = θ_array[key][1]
        r = sum(posteriorsamples .< priorsample)
        priorsampledict[key] = priorsample
        rdict[key] = r
    end

    return priorsampledict, rdict, chains
end

# Calibrate model by running chains and save results 
function calibrate(system::System, chainparams::Dict, saveas::String)

    # Get chain parameters 
    target_accept = chainparams["target_accept"]
    num_chains = chainparams["num_chains"]
    adaptation = chainparams["adaptation"]
    iterations = chainparams["iterations"]
    thinning = chainparams["thinning"]
    tree_depth = chainparams["tree_depth"]
    verbosity = chainparams["verbosity"]

    # Run chains
    calib = calibrationhmc(system, target_accept, num_chains, adaptation, iterations, thinning, tree_depth, verbosity)
    priorsampledict, rdict, chains = calib
    
    # Save chain parameters
    chainparams_data = JSON.json(chainparams)
    open("$(saveas)_chain_params.json", "w") do f
        write(f, chainparams_data)
    end

    # Save sampled parameter values
    priorsampledict_data = JSON.json(priorsampledict)
    open("$(saveas)_prior_samples.json", "w") do f
        write(f, priorsampledict_data)
    end

    # Save cumulative probabilities of sampled values w.r.t. the chains
    rdict_data = JSON.json(rdict)
    open("$(saveas)_rank_stats.json", "w") do f
        write(f, rdict_data)
    end

    # Save chains
    CSV.write("$(saveas)_chains.csv", chains)

    return nothing
end
export calibrate

# Generates histograms and corner plot of calibration results
function calibrationplots(datadir::String, plotsdir::String; histdpi::Int=300,
                          histcolour::Union{Symbol,String}="#1E90FF", filetype::String="png")

    # Get all files with r values
    files = readdir(datadir)
    rankfiles = datadir .* [i for i ∈ files if occursin("rank_stats", i)]

    # Get named tuples of CDF values from files
    rankarr = map(rankfiles) do i
        dict = loadJSON(i)
        return (; (Symbol(k) => v for (k,v) in dict)...)
    end

    # Convert to data frame
    rankdata = DataFrame(rankarr)
    colnames = names(rankdata)

    # Plot histograms
    for name in colnames
        data = rankdata[!,name]
        nbins = floor(Int, sqrt(length(data)))
        hist = histogram(data, label="", dpi=histdpi, color=histcolour)
        xlabel!(name)
        savefig(hist, plotsdir * "$name.$filetype")
    end

    # Plot corner plot 
    c = corner(rankdata, hist_kwargs=(;nbins=5), hist2d_kwargs=(;nbins=5))
    savefig(c, plotsdir * "corner.$filetype")

    return nothing 
end
export calibrationplots

# ----------------------------------------------------------------------------------------------------------------------