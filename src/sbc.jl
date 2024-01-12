# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

# General
using StatsBase
# using StaticArrays

# Data
using TOML
using Random
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
function calibrationhmc(
    system::System;
    rng=Random.default_rng(),
    verbosity=0,
    target_accept,
    θ = sample_priors(rng, system),
    kwargs...
)

    θ_newsystem_flat = θ
    # Get parameter values sampled from priors and generate a new system
    # θ_newsystem = drawfrompriors(system)
    arr2nt = make_arr2nt(system)
    θ_newsystem = arr2nt(θ_newsystem_flat)

    verbosity > 0 && println("= SBC ="*"="^70)
    verbosity > 0 && println("Drew parameters: ", stringify_nested_named_tuple(θ_newsystem))

    newsystem = generate_from_params(system, θ_newsystem)

    model = Octofitter.LogDensityModel(newsystem; autodiff=:ForwardDiff, verbosity=4)

    # Run chains
    chains = octofit(
        model, target_accept;
        verbosity,
        kwargs...
    )

    # Also calculate the log posterior and log likelihood of these parameters (Modrak et al 2022)

    ln_like_function = make_ln_like(newsystem, θ_newsystem)
    loglike = ln_like_function(model.system, θ_newsystem)

    ln_prior_function = make_ln_prior(newsystem)
    logprior = ln_prior_function(θ_newsystem_flat)

    logpost = logprior + loglike

    
    θ_array = Octofitter.result2mcmcchain([(;loglike, logpost, θ_newsystem...)])

    # θ_t = sample.z.θ
    # θ = model.invlink(θ_t)
    # resolved_namedtuple = model.arr2nt(θ)
    # # Add log posterior, tree depth, and numerical error reported by
    # # the sampler.
    # # Also recompute the log-likelihood and add that too.
    # loglike = ln_like(model.system, resolved_namedtuple)
    # @show (;
    #     loglike = loglike,
    #     logpost = sample.z.ℓπ.value,
    

    # Compute a good guess for how much thinning we should do based on 
    # the effective sample size
    # ess = median(filter(isfinite, MCMCChains.ess_rhat(chains)[:,:ess]))
    # keep_ii = range(
    #     start=1,
    #     stop=size(chains,1),
    #     step=round(Int, size(chains,1)/ess, RoundUp)
    # )
    # verbosity > 1 && @info "Thinning to keep $(length(keep_ii)) independent samples."
    # chains_thinned = chains[keep_ii,:,:]

    # Compute rank statistic of each parameter in the chains
    # that is also in θ_array (otherwise we can't compare)
    chainkeys = string.(intersect(keys(θ_array), keys(chains)))
    priorsampledict = OrderedDict()
    rdict = OrderedDict()

    for key in chainkeys
        posteriorsamples = vec(chains[key])
        priorsample = only(θ_array[key])
        # Find rank
        r = count(<(priorsample), posteriorsamples)
        priorsampledict[key] = priorsample
        rdict[key] = r/length(chains[key])*100 # Should be on scale of 1-100
    end

    return priorsampledict, rdict, chains
end

# Calibrate model by running chains and save results 
function sbctrial(system::System, chainparams, saveas::AbstractString)
    chainparams = Dict{Symbol,Any}(pairs(chainparams))
    
    # Get verbosity from chain parameters to decide how much we should log
    # Rest can just be forwarded
    verbosity = chainparams[:verbosity] = get(chainparams, :verbosity, 2)
    target_accept = get(chainparams, :target_accept, 0.85)
    delete!(chainparams, :target_accept)

    # Run chains
    calib = calibrationhmc(system; target_accept, verbosity, chainparams...)
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

    Octofitter.savechain("$(saveas)_chains.fits", chains)

    return
end

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

    return
end

# ----------------------------------------------------------------------------------------------------------------------
