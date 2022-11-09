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
using MCMCChains
using PlanetOrbits
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
export drawfrompriors

# Generate new astrometry observations
function newobs(obs::Astrometry, elem::VisualOrbit, θ_planet)

    # Get epochs and uncertainties from observations
    epochs = obs.table.epoch
    σ_ras = obs.table.σ_ra 
    σ_decs = obs.table.σ_dec

    # Generate now astrometry data
    ras = raoff.(elem, epochs)
    decs = decoff.(elem, epochs)
    astrometry_table = Table(epoch=epochs, ra=ras, dec=decs, σ_ra=σ_ras, σ_dec=σ_decs)

    return Astrometry(astrometry_table)
end

# Generate new radial velocity observations for a planet
function newobs(obs::RadialVelocity, elem::VisualOrbit, θ_planet)

    # Get epochs and uncertainties from observations
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocity(radvel_table)
end

# Generate new radial velocity observations for a star
function newobs(obs::RadialVelocity, elems::Vector{<:VisualOrbit}, θ_system)

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = radvel.(reshape(elems, :, 1), epochs, transpose(planet_masses))
    noise = randn(length(epochs)) .* θ_system.jitter
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocity(radvel_table)
end

# Generate new images
function newobs(obs::Images, elements::Vector{<:VisualOrbit}, θ_system)

    newrows = map(obs.table) do row
        (;band, image, platescale, epoch, psf) = row

        injected = copy(image)
        
        for i in eachindex(elements)
            θ_planet = θ_system.planets[i]
            elem = elements[i]
            # Generate new astrometry point
            os = orbitsolve(elem, epoch)

            # TODO: this does not consider the shift to the images due to the motion of the star
            ra = raoff(os)
            dec = decoff(os)

            phot = θ_planet[band]

            # TODO: verify signs
            dx = ra/platescale
            dy = -dec/platescale
            translation_tform = Translation(
                mean(axes(psf,1))-mean(axes(image,1))+mean(dims(image,1))+dx,
                mean(axes(psf,2))-mean(axes(image,2))+mean(dims(image,2))+dy
            )
            # TBD if we want to support rotations for handling negative sidelobes.

            psf_positioned = warp(arraydata(psf), translation_tform, axes(image), fillvalue=0)
            psf_positioned[.! isfinite.(psf_positioned)] .= 0
            psf_scaled = psf_positioned .* phot ./ maximum(filter(isfinite, psf_positioned))
            injected .+= psf_scaled
        end

        return merge(row, (;image=injected))
    end

    return Images(newrows)
end

"""
Specific HGCA proper motion modelling. Model the GAIA-Hipparcos/Δt proper motion
using 5 position measurements averaged at each of their epochs.
"""
function newobs(obs::ProperMotionAnomHGCA, elements, θ_system)
    ll = 0.0

    # This observation type just wraps one row from the HGCA (see hgca.jl)
    hgca = obs.table
    # Roughly over what time period were the observations made?
    dt_gaia = 1038 # EDR3: days between  Date("2017-05-28") - Date("2014-07-25")
    dt_hip = 4*365
    # How many points over Δt should we average the proper motion and stellar position
    # at each epoch? This is because the PM is not an instantaneous measurement.
    N_ave = 25 #5

    # Look at the position of the star around both epochs to calculate 
    # our modelled delta-position proper motion

    # First epoch: Hipparcos
    ra_hip_model = 0.0
    dec_hip_model = 0.0
    pmra_hip_model = 0.0
    pmdec_hip_model = 0.0
    # The model can support multiple planets
    for i in eachindex(elements)
        θ_planet = θ_system.planets[i]
        orbit = elements[i]
        if θ_planet.mass < 0
            return -Inf
        end
        # Average multiple observations over a timescale +- dt/2
        # to approximate what HIPPARCOS would have measured.
        for δt = range(-dt_hip/2, dt_hip/2, N_ave)
            # RA and dec epochs are usually slightly different
            # Note the unit conversion here from jupiter masses to solar masses to 
            # make it the same unit as the stellar mass (element.mu)
            # TODO: we can't yet use the orbitsolve interface here for the pmra calls,
            # meaning we calculate the orbit 2x as much as we need.
            o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_hip[1])+δt)
            o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_hip[1])+δt)
            ra_hip_model += -raoff(o_ra) * θ_planet.mass*mjup2msol/orbit.M
            dec_hip_model += -decoff(o_dec) * θ_planet.mass*mjup2msol/orbit.M
            pmra_hip_model += pmra(o_ra, θ_planet.mass*mjup2msol)
            pmdec_hip_model += pmdec(o_dec, θ_planet.mass*mjup2msol)
        end
    end
    ra_hip_model/=N_ave
    dec_hip_model/=N_ave
    pmra_hip_model/=N_ave
    pmdec_hip_model/=N_ave

    # Last epoch: GAIA
    ra_gaia_model = 0.0
    dec_gaia_model = 0.0
    pmra_gaia_model = 0.0
    pmdec_gaia_model = 0.0
    # The model can support multiple planets
    for i in eachindex(elements)
        θ_planet = θ_system.planets[i]
        orbit = elements[i]
        if θ_planet.mass < 0
            return -Inf
        end
        # Average multiple observations over a timescale +- dt
        # to approximate what HIPPARCOS and GAIA would have measured.
        for δt = range(-dt_gaia/2, dt_gaia/2, N_ave)
            # RA and dec epochs are usually slightly different
            # Note the unit conversion here from jupiter masses to solar masses to 
            # make it the same unit as the stellar mass (element.M)
            o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_gaia[1])+δt)
            o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_gaia[1])+δt)
            ra_gaia_model += -raoff(o_ra) * θ_planet.mass*mjup2msol/orbit.M
            dec_gaia_model += -decoff(o_dec) * θ_planet.mass*mjup2msol/orbit.M
            pmra_gaia_model += pmra(o_ra, θ_planet.mass*mjup2msol)
            pmdec_gaia_model += pmdec(o_dec, θ_planet.mass*mjup2msol)
        end
    end
    ra_gaia_model/=N_ave
    dec_gaia_model/=N_ave
    pmra_gaia_model/=N_ave
    pmdec_gaia_model/=N_ave

    # Model the GAIA-Hipparcos delta-position velocity
    pmra_hg_model = (ra_gaia_model - ra_hip_model)/(years2mjd(hgca.epoch_ra_gaia[1]) - years2mjd(hgca.epoch_ra_hip[1]))
    pmdec_hg_model = (dec_gaia_model - dec_hip_model)/(years2mjd(hgca.epoch_dec_gaia[1]) - years2mjd(hgca.epoch_dec_hip[1]))

    return ProperMotionAnomHGCA(merge(hgca[1], (;
        pmra_hip=pmra_hip_model,
        pmdec_hip=pmdec_hip_model,
        pmra_gaia=pmra_gaia_model,
        pmdec_gaia=pmdec_gaia_model,
        pmra_hg=pmra_hg_model,
        pmdec_hg=pmdec_hg_model,
    )))

end

# Generate calibration data
function generate(system::System, θ_newsystem = drawfrompriors(system))

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

    return newsystem
end
export generate

# Run chains on new model
function calibrationhmc(system::System, target_accept, num_chains, adaptation, iterations, thinning, tree_depth, verbosity)

    # Get parameter values sampled from priors and generate a new system
    θ_newsystem = drawfrompriors(system)
    newsystem = generate(system, θ_newsystem)
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
    priorsampledict = OrderedDict()
    rdict = OrderedDict()

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
function calibrate(system::System, chainparams::AbstractDict, saveas::AbstractString)

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

# Get data required for rank statistic plots
function getplotdata(datadir::String, plotsdir::String, grmax)

    # Load files and create file arrays
    files = readdir(datadir)
    rankfiles = datadir .* [i for i ∈ files if occursin("rank_stats.json", i)]
    usablerankfiles = []
    chainfiles = datadir .* [i for i ∈ files if occursin("chains.csv", i)]
    pairedfiles = [(rankfiles[i], chainfiles[i]) for i ∈ 1:length(rankfiles)]

    # Determine if each chain is usable based on Gelman-Rubin diagnostic
    for pair ∈ pairedfiles 

        # Load chains data as matrix
        rankfile, chainfile = pair 
        df = CSV.read(chainfile, DataFrame)
        m = Matrix(df)

        # Get number of subchains and number of iterations per chain
        nsubchains = Int(m[:,2][end])
        maxval = Int(m[:,1][end])

        # Convert subchains from matrix to Chains object
        subchains = []
        for i ∈ 0:(nsubchains - 1)
            startidx = i*maxval + 1
            stopidx = (i + 1)*maxval 
            push!(subchains, m[startidx:stopidx,3:end])
        end
        subchains = Tuple(subchains)
        c = Chains(cat(subchains..., dims=nsubchains), propertynames(df)[3:end])

        # Execute Gelman-Rubin diagnostic
        grdiagnostic = gelmandiag(c)
        psrf = grdiagnostic.nt.psrf 

        # Save rank statistic file if Gelman-Rubin condition is satisfied
        if count(i -> i > grmax, psrf) == 0
            push!(usablerankfiles, rankfile)
        end
    end

    # Save numbers of usable chains
    open(plotsdir * "usable_chains.txt", "w") do f
        write(f, "Usable chains (Gelman-Rubin ≤ $grmax):\n")
        for file ∈ usablerankfiles 
            filenumber = split(split(file, "\\")[end], "_")[1]
            write(f, "$filenumber\n")
        end
        write(f, "Total: $(length(usablerankfiles))/$(length(rankfiles))")
    end

    # Get maximum histogram values
    m = Matrix(CSV.read(chainfiles[1], DataFrame))
    maxval = Int(m[:,1][end] * Int(m[:,2][end]))

    return maxval, usablerankfiles
end

# Generate histograms and corner plot of calibration results
function calibrationplots(datadir::String, plotsdir::String; grmax::Number=1.2, histdpi::Int=300,
    histcolour::Union{Symbol,String}="#1E90FF", filetype::String="png")

    # Get all rank statistic files and max histogram value
    maxval, rankfiles = getplotdata(datadir, plotsdir, grmax)

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
        @info "Plotting rank statistic histogram for $name"
        data = rankdata[!,name]
        nbins = floor(Int, sqrt(length(data)))
        histvals = fit(StatsBase.Histogram, vec(data), range(0, stop=maxval, length=nbins+1))
        hist = plot(histvals, label="", dpi=histdpi, color=histcolour, bins=())
        xlabel!(name)
        savefig(hist, plotsdir * "$name.$filetype")
    end

    # Plot corner plot 
    @info "Plotting rank statistic corner plot"
    c = PairPlots.corner(rankdata, hist_kwargs=(;nbins=5), hist2d_kwargs=(;nbins=5), plotcontours=false)
    savefig(c, plotsdir * "corner.$filetype")

    return nothing 
end
export calibrationplots

# ----------------------------------------------------------------------------------------------------------------------