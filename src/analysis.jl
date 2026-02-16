
# This file contains functions for analysing results from chains.
using Dates

"""
    octoplot(
        system::Octofitter.Model,
        chain::MCMCChains.Chains;
        fname="\$(system.name)-plot-grid",
        color = "\$(first(keys(system.planets)))_e",
        colorbartitle=string(color),
        clims = quantile(vec(chain[color]),(0.01, 0.99)),
        cmap = :plasma,
        dpi=200,
        kwargs...
    )

Given a  `LogDensityModel` and an MCMC Chain (sampled either from 
the posterior or the prior), produce a panel of 9 plots visualizing the orbits.
The output is saved to a file based on the name of the system (`fname`).
"""
function octoplot(::Octofitter.LogDensityModel, ::Any)
    error("You must load the Makie.jl package (`using CairoMakie`, or `using GLMakie`) before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function octoplot!(::Octofitter.LogDensityModel, ::Any)
    error("You must load the Makie.jl package (`using CairoMakie`, or `using GLMakie`) before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end

function octocorner(::Octofitter.LogDensityModel, ::Any)
    error("You must load the Makie package (eg `using CairoMakie`) and PairPlots package (eg `using PairPlots`) before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function octocorner(model::LogDensityModel, args...; kwargs...)
    return octocorner(model.system, args...; kwargs...)
end


function astromplot(args...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function astromplot!(args...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end



function hgcaplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function hgcaplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


function masspostplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function masspostplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


function rvpostplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function rvpostplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end

function rvpostplot_animated( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


function gaiatimeplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function gaiatimeplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


function gaiastarplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function gaiastarplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


function skytrackplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function skytrackplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


function dotplot( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end
function dotplot!( args...; kwargs...)
    error("You must load the Makie package (eg `using CairoMakie`)  before calling this function. Then, pass your model and chain result as arguments.  If you're seeing this message despite loading those packages, check that you are passing the correct argument types.")
end


export octoplot, octocorner
