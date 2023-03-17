# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

# General
using StatsBase
# using StaticArrays

# Data
using TOML
# using CSV
# using JSON
# using TypedTables
# using DataFrames 

# Plots
# using Plots
# using PairPlots


# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------
# Load data in .json file to dictionary 
# function loadJSON(filepath::String)
#     dict = open(filepath, "r") do f
#         TOML.parse(f)
#     end
#     return dict
# end


# Run chains on new model
function calibrationhmc(system::System, target_accept, num_chains, adaptation, iterations, thinning, tree_depth, verbosity)

    # Get parameter values sampled from priors and generate a new system
    θ_newsystem = drawfrompriors(system)

    msg = stringify_nested_named_tuple(θ_newsystem)
    println("="^80)
    println("Drew parameters $msg")

    θ_array = Octofitter.result2mcmcchain([θ_newsystem])
    newsystem = generate_from_params(system, θ_newsystem)

    model = Octofitter.LogDensityModel(newsystem; autodiff=:ForwardDiff, verbosity=4)

    # Run chains
    chains = Octofitter.advancedhmc(
        model, target_accept,
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
    target_accept = chainparams[:target_accept]
    num_chains = chainparams[:num_chains]
    adaptation = chainparams[:adaptation]
    iterations = chainparams[:iterations]
    thinning = chainparams[:thinning]
    tree_depth = chainparams[:tree_depth]
    verbosity = chainparams[:verbosity]

    # Run chains
    calib = calibrationhmc(system, target_accept, num_chains, adaptation, iterations, thinning, tree_depth, verbosity)
    priorsampledict, rdict, chains = calib

    
    
    # Save chain parameters
    open("$(saveas)_sampler_parameters.toml", "w") do f
        TOML.print(f, chainparams)
    end

    # Save sampled parameter values
    open("$(saveas)_parameters.toml", "w") do f
        TOML.print(f, priorsampledict)
    end

    # Save cumulative probabilities of sampled values w.r.t. the chains
    open("$(saveas)_rank_stats.toml", "w") do f
        TOML.print(f, rdict)
    end

    return nothing
end
export calibrate

# Get data required for rank statistic plots
function getplotdata(datadir::String, plotsdir::String, grmax)

    # Load files and create file arrays
    files = readdir(datadir)
    rankfiles = datadir .* [i for i ∈ files if occursin("rank_stats.toml", i)]
    usablerankfiles = []
    chainfiles = datadir .* [i for i ∈ files if occursin("chains.fits", i)]
    pairedfiles = [(rankfiles[i], chainfiles[i]) for i ∈ 1:length(rankfiles)]

    # Determine if each chain is usable based on Gelman-Rubin diagnostic
    for (rankfile, chainfile) ∈ pairedfiles 

        # Load chains data as matrix
        df = loadchain("$(saveas)_chains.fits")
        m = Array(df)

        # Get number of subchains and number of iterations per chain
        # TODO: check this, I think this assumes the data was flattend
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
    # TODO: rather than discarding, we need an error message. To be representative, we 
    # need to increase convergence not discard.
    open(plotsdir * "usable_chains.txt", "w") do f
        write(f, "Usable chains (Gelman-Rubin ≤ $grmax):\n")
        for file ∈ usablerankfiles 
            filenumber = split(split(file, "\\")[end], "_")[1]
            write(f, "$filenumber\n")
        end
        write(f, "Total: $(length(usablerankfiles))/$(length(rankfiles))")
    end

    # Get maximum histogram values
    # TODO: I don't quite understand what this code is doing yet
    m = loadchain(chainfiles[1])
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
        dict = TOML.parsefile(i)
        return (; (Symbol(k) => v for (k,v) in dict)...)
    end

    # Convert to data frame
    rankdata = Table(rankarr)
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
