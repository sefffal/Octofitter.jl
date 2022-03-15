# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

# General
using StatsBase
using StaticArrays

# Data
using JSON
using JLD2
using TypedTables
using DataFrames 

# Plots
using Plots
using PairPlots

# Orbit fitting
using DirectOrbits
using DirectDetections

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
function newobs(obs::Astrometry, elem::KeplerianElements)

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
function newobs(obs::RadialVelocity, elem::KeplerianElements)

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
            return newobs(obs, elem)
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
function calibrationhmc(system::System, acceptance, num_chains, adaptation, iterations, tree_depth, verbosity)

    # Get parameter values sampled from priors and generate a new system
    θ_newsystem, newsystem = calibrationsystem(system)
    θ_array = DirectDetections.result2mcmcchain(newsystem, [θ_newsystem])

    # Run chains
    @time chains = DirectDetections.hmc(
        newsystem, acceptance,
        num_chains = num_chains,
        adaptation = adaptation,
        iterations = iterations,
        tree_depth = tree_depth,
        verbosity = verbosity
    )

    # Calculate cumulative probabilities of sampled values w.r.t. the chains
    chainkeys = string.(keys(chains))
    sampleddict = Dict()
    cdfdict = Dict()

    for key in chainkeys

        # Generate empirical CDF using chains
        paramdist = vec(chains[key])
        paramcdf = ecdf(paramdist)

        # Get cumulative probability of sampled value
        sampledval = θ_array[key][1]
        cdfval = paramcdf(sampledval)

        # Set dictionary elements
        sampleddict[key] = sampledval
        cdfdict[key] = cdfval
    end

    return sampleddict, cdfdict, chains
end

# Calibrate model by running chains and save results 
function calibrate(system::System, chainparams::Dict, saveas::String)

    # Get chain parameters 
    acceptance = chainparams["acceptance"]
    num_chains = chainparams["num_chains"]
    adaptation = chainparams["adaptation"]
    iterations = chainparams["iterations"]
    tree_depth = chainparams["tree_depth"]
    verbosity = chainparams["verbosity"]

    # Run chains
    sampleddict, cdfdict, chains = calibrationhmc(system, acceptance, num_chains, adaptation, iterations, tree_depth, verbosity)
    
    # Save chain parameters
    chainparams_data = JSON.json(chainparams)
    open("$(saveas)_chain_params.json", "w") do f
        write(f, chainparams_data)
    end

    # Save sampled parameter values
    sampleddict_data = JSON.json(sampleddict)
    open("$(saveas)_sampled_vals.json", "w") do f
        write(f, sampleddict_data)
    end

    # Save cumulative probabilities of sampled values w.r.t. the chains
    cdfdict_data = JSON.json(cdfdict)
    open("$(saveas)_cdf_vals.json", "w") do f
        write(f, cdfdict_data)
    end

    # Save chains
    jldsave("$(saveas)_chains.jld2"; chains)

    return nothing
end
export calibrate

# Generates histograms and corner plot of calibration results
function calibrationplots(datadir::String, plotsdir::String; histdpi::Int=300,
                          histcolour::Union{Symbol,String}="#1E90FF", filetype::String="png")

    # Get all files with CDF values
    files = readdir(datadir)
    cdffiles = datadir .* [i for i ∈ files if occursin("cdf_vals", i)]

    # Get named tuples of CDF values from files
    cdfarr = map(cdffiles) do i
        dict = loadJSON(i)
        return (; (Symbol(k) => v for (k,v) in dict)...)
    end

    # Convert to data frame
    cdfdata = DataFrame(cdfarr)
    colnames = names(cdfdata)

    # Plot histograms
    for name in colnames
        data = cdfdata[!,name]
        hist = histogram(data, label="", dpi=histdpi, color=histcolour)
        xlabel!(name)
        savefig(hist, saveas * "$name.$filetype")
    end

    # Plot corner plot 
    c = corner(df)
    savefig(c, saveas * "corner.$filetype")

    return nothing 
end
export calibrationplots

# ----------------------------------------------------------------------------------------------------------------------