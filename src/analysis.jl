
# This file contains funtions for analysing results from chains.
using Dates

"""
    projectpositions(chains.planets[1], mjd("2020-02-02"))

Given the posterior for a particular planet in the model and a modified julian date(s),
return `ra` and `dec` offsets in mas for each sampling in the posterior.
"""
function projectpositions(chains, planet_key, times)

    ras = zeros(size(chains,1),size(chains,3),length(times))
    decs = zeros(size(chains,1),size(chains,3),length(times))
    
    for is in collect(Iterators.partition(1:size(chains,1), 5000))
        for j in 1:size(chains,3)
            els = Octofitter.construct_elements(chains, planet_key, is .* j)
            for (i,el) in zip(is,els)
                for (k, t) in enumerate(times)
                    o = orbitsolve(el, t)
                    ras[i,j,k] = raoff(o)
                    decs[i,j,k] = decoff(o)
                end
            end
        end
    end
    return ras, decs
end
export projectpositions



function bayesfactor(chain, planet, property)
    prior = model.system.planets[planet].priors.priors[property]
    post = chain["$planet.$property"]

    # Fit a histogram to the posterior.
    # TODO: in future this could be a KDE
    nbins = floor(Int,sqrt(length(post))/2)
    h = fit(Histogram, vec(post), range(0,maximum(post),length=nbins))
    hn = normalize(h, mode=:pdf)
    
    i_first_nonzero_bin = findfirst(>(0), hn.weights)
    bin_centre = mean(hn.edges[1][i_first_nonzero_bin:i_first_nonzero_bin+1])
    
    lnbf = logpdf(prior, bin_centre) - log(hn.weights[i_first_nonzero_bin])

    # Now fit a Gaussian to the posterior in order to estimate bayes factors for
    # very strong detections but assuming the posterior is normally distributed
    post = fit(Normal, vec(post))
    lnbf_exrap = logpdf(prior, 1e-6) - logpdf(post, 1e-6)

    return (;lnbf, islowerlimit=i_first_nonzero_bin>1, lnbf_exrap)
end





function plotchains end
function plotchains! end
function timeplotgrid end
export plotchains, plotchains!, timeplotgrid
