#=
Need plots for:

orbits.#
orbits shaded by property.
orbits shaded by z position.
multiple planets^
overtop images^
phot vs a plots.
histograms of all properties.

astrometry#
pma plots scattered by property.
pma time plots colored by property.

Spectra/photometry

function to generate all.
=#

import Contour as Contourlib
# using CairoMakie
using GLMakie
using DirectOrbits
using StatsBase

function plot_orbits!(gp::GridPosition, args...; kwargs...)
    ax_sky = Axis(
        gp,
        xreversed=true,
        autolimitaspect=1,
        xlabel="ΔRA (mas)",
        ylabel="ΔDEC (mas)"
    )
    plot_orbits!(ax_sky, args...; kwargs...)
end
function plot_orbits!(ax::Axis, chains, planet_key;
    color,
    N=500,
    ii = rand(1:size(chains,1) * size(chains,3), N),
    colorrange=quantile(filter(isfinite, color),(0.05,0.95))
)
    orbits = DirectDetections.construct_elements(chains, planet_key, ii)
    posns = DirectOrbits.kep2cart_ν.(orbits, range(0, 2π, length=90)')
    ras = getproperty.(posns, :x)
    decs = getproperty.(posns, :y)
    color = color[ii]
    for (ras, decs, c) in zip(eachrow(ras), eachrow(decs), color)
        lines!(ax, ras, decs, linewidth=0.5; color=fill(c, size(ras)), colorrange, colormap=:plasma)
    end
    lines!(ax, [0], [0],  color=:black, label=string(planet_key))
    # Makie.scatter!(ax, [0], [0], marker="⋆", color=:white, markersize=30)
    # Makie.scatter!(ax, [0], [0], marker="⋆", color=:black, markersize=20)
end

function plot_image!(ax::Axis, image, platescale;
    colorrange=quantile(filter(isfinite, image),(0.05,0.95)),
    colormap=:magma,
)
    Makie.heatmap!(
        ax,
        (collect(axes(img,1)) .* platescale),
        (collect(axes(img,2)) .* platescale),
        collect(img[end:-1:begin,:]),
        colormap,
        colorrange,
    )
end

function plot_astrom!(ax::Axis, system::System)
    for planet in system.planets
        if !isnothing(planet.astrometry)
            errorbars!(ax,
                collect(planet.astrometry.ra),
                collect(planet.astrometry.dec),
                collect(planet.astrometry.σ_dec),
                direction=:y,
                color=:white,
                linewidth=4,
            )
            errorbars!(ax,
                collect(planet.astrometry.ra),
                collect(planet.astrometry.dec),
                collect(planet.astrometry.σ_ra),
                direction=:x,
                color=:white,
                linewidth=4,
            )
            errorbars!(ax,
                collect(planet.astrometry.ra),
                collect(planet.astrometry.dec),
                collect(planet.astrometry.σ_dec),
                direction=:y,
                color=:red,
            )
            errorbars!(ax,
                collect(planet.astrometry.ra),
                collect(planet.astrometry.dec),
                collect(planet.astrometry.σ_ra),
                direction=:x,
                color=:red,
            )
        end
    end
end

function plot_pma_epoch!(ax::Axis, system::System, ind::Integer)
    system_pma = system.propermotionanom
    Makie.errorbars!(
        ax,
        [system_pma.pm_ra[ind]],
        [system_pma.pm_dec[ind]],
        [system_pma.σ_pm_ra[ind]],
        direction=:x,
        color=:red,
    )
    Makie.errorbars!(
        ax,
        [system_pma.pm_ra[ind]],
        [system_pma.pm_dec[ind]],
        [system_pma.σ_pm_dec[ind]],
        direction=:y,
        color=:red,
    )
    Makie.hlines!(ax, 0; color=:black, label="")
    Makie.vlines!(ax, 0; color=:black, label="")
end

function plot_pma_epoch!(ax::Axis, chains, ind::Integer;
    color=:black,
    N=500,
    ii = rand(1:size(chains,1) * size(chains,3), N),
    colorrange=quantile(filter(isfinite, color),(0.05,0.95))
)
    system_pma = chains.info.model.propermotionanom
    vx = zeros(size(ii))
    vy = zeros(size(ii))
    for j in keys(chains.info.model.planets)
        elements = DirectDetections.construct_elements(chains, j, ii)
        mass = vec(chains["$j[mass]"][ii])
        vx .+= getindex.(propmotionanom.(elements, system_pma.ra_epoch[ind], mass.*DirectDetections.mjup2msol),1)
        vy .+= getindex.(propmotionanom.(elements, system_pma.dec_epoch[ind], mass.*DirectDetections.mjup2msol),2)
    end
    if typeof(color) <: AbstractArray
        color = color[ii]
    end
    Makie.scatter!(ax, vx, vy; color, markersize=5, colorrange, colormap=:plasma)
end

function plot_pma_time!(ax::Axis, system::System, direction; orientation=:horizontal)
    system_pma = system.propermotionanom
    if direction == :ra
        y = collect(system_pma.pm_ra)
        σ = collect(system_pma.σ_pm_ra)
        epoch = collect(system_pma.ra_epoch)
    elseif direction == :dec
        y = collect(system_pma.pm_dec)
        σ = collect(system_pma.σ_pm_dec)
        epoch = collect(system_pma.dec_epoch)
    else
        error("Unsupported direction. Must be one of :ra or :dec")
    end
    if orientation == :horizontal
        Makie.errorbars!(ax, epoch, y, σ, direction=:y, color=:white, linewidth=4)
        Makie.errorbars!(ax, epoch, y, σ, direction=:y, color=:red)
    else
        Makie.errorbars!(ax, y, epoch, σ, direction=:x, color=:white, linewidth=4)
        Makie.errorbars!(ax, y, epoch, σ, direction=:x, color=:red)
    end
end
function plot_pma_time!(ax::Axis, chains, direction;
    orientation=:horizontal,
    color,
    N=100,
    ii = rand(1:size(chains,1)*size(chains,3), N),
    colormap=:plasma,
    colorrange=quantile(filter(isfinite, color),(0.05,0.95)),
    linewidth=0.5
)
    system_pma = chains.info.model.propermotionanom
    if direction == :ra
        dir_ind = 1
    elseif direction == :dec
        dir_ind = 2
    else
        error("Unsupported direction. Must be one of :ra or :dec")
    end

    epoch_min = minimum([system_pma.ra_epoch; system_pma.dec_epoch;])
    epoch_max = maximum([system_pma.ra_epoch; system_pma.dec_epoch;])
    epoch_dt = epoch_max - epoch_min
    epoch = range(epoch_min - 0.2epoch_dt, epoch_max + 0.2epoch_dt, length=200)

    vy = zeros(length(ii), length(epoch))
    for j in keys(chains.info.model.planets)
        elements = DirectDetections.construct_elements(chains, j, ii)
        mass = vec(chains["$j[mass]"][ii])
        vy .+= getindex.(propmotionanom.(elements, epoch', mass.*DirectDetections.mjup2msol), dir_ind)
    end
    # if typeof(color) <: AbstractArray
        color = color[ii]
    # end
    if orientation == :horizontal
        for (row,c) in zip(eachrow(vy), color)
            Makie.lines!(ax, epoch, row; color=fill(c, size(epoch)), colormap, colorrange,  linewidth)
        end
    else
        for (row,c) in zip(eachrow(vy), color)
            Makie.lines!(ax, row, epoch; color=fill(c, size(epoch)), colormap, colorrange,  linewidth)
        end
    end
end

function hist_prop!(ax::Axis, prop; color=:blue, label="")
    h = StatsBase.normalize(fit(Histogram, vec(prop)), mode=:pdf)
    Makie.stairs!(ax, h; color, label)
    Makie.ylims!(ax, low=0)
end

function plot_phot!(gp::GridPosition, planet::Planet; axis=(;), kwargs...)
    bands = collect(string.(planet.photometry.band))
    ax = Axis(gp;
        xlabel="band",
        xticks=(collect(eachindex(bands)), bands),
        axis...
    )
    plot_phot!(ax, planet)
    return ax
end
function plot_phot!(ax::Axis, planet::Planet)
    bands = string.(planet.photometry.band)
    Makie.errorbars!(ax, 
        eachindex(bands),
        collect(planet.photometry.phot),
        collect(planet.photometry.σ_phot),
        color=:red,
    )
end

# TODO
function plot_phot!(ax::Axis, chains, planet_key;
    N,
    ii = rand(1:size(chains,1)*size(chains,3), N),
    color=:black
)
    planet = chains.info.model.planets[planet_key]
    bands = string.(planet.photometry.band)
    Makie.errorbars!(ax, 
        eachindex(bands),
        collect(planet.photometry.phot),
        collect(planet.photometry.σ_phot),
        color=:red,
    )
end

function compute_levels(H)
    # Calculate levels for contours
   levels = 1 .- exp.(-0.5 .* (0.5:0.5:2.1).^2)
   ii = sortperm(reshape(H,:))
   h2flat = H[ii]
   sm = cumsum(h2flat)
   sm /= sm[end]
   if all(isnan, sm) || length(H) <= 1
       @warn "Could not compute valid contours"
       plotcontours = false
       V = [0]
   else
       V = sort(map(v0 -> h2flat[sm .≤ v0][end], levels))
       if any(==(0), diff(V))
           @warn "Too few points to create valid contours"
       end
   end
   return V
end
function margin_hist(
    chains,
    key1,
    key2,
    range1=range(extrema(chains[key1])...,length=30),
    range2=range(extrema(chains[key2])...,length=30);
    color=:black
)

    # TODO: this needs to be generalized somehow
    priors = NamedTuple(chains.info.model.planets[1].priors.priors)    
    
    ha = fit(Histogram, vec(chains[key1]), range1)
    hl = fit(Histogram, vec(chains[key2]), range2)

    ## Main 2D histogram
    hla = fit(
        Histogram,
        (vec(chains[key1]), vec(chains[key2])),
        (range1, range2)
    )

    fig = Figure(
        resolution = (650*2, 400*2),
        figure_padding=1,
        font="Arial",
        fontsize=20,
    )
    
    ax = Axis(
        fig[2,1],
        xlabel=key1,
        ylabel=key2,

    )
    heatmap!(ax, hla, colormap=cgrad([:white, color]))

    l = compute_levels(hla.weights)[1:1]
    ax1 = hla.edges[1][1:end-1] .+ step(hla.edges[1])/2
    ax2 = hla.edges[2][1:end-1] .+ step(hla.edges[2])/2
    for cl in Contourlib.levels(Contourlib.contours(ax2, ax1, hla.weights', l))
        lvl = Contourlib.level(cl) # the z-value of this contour level
        thefirst=true
        for line in Contourlib.lines(cl)
            xs, ys = coordinates(line) # coordinates of this line segment
#             p = Makie.GeometryBasics.Polygon(
#                 Makie.GeometryBasics.Point2.(zip(ys, xs)),
#             )
#             poly!(ax, p, color=RGBAf0(model_colours[cnum],0.4))
            lines!(ax, ys, xs; color)

            thefirst=false
        end
    end


    # Maginal histograms
    ax_l = Axis(
        fig[2,2],
    )
    hidexdecorations!(ax_l)
    hideydecorations!(ax_l)
    
    # Also plot the prior pdf scaled to the same height
    # x = range(extrema(range2)..., length=250)
    # prior_pdf = pdf.(priors.mass,x)
    # Makie.lines!(ax_l,  prior_pdf ./ maximum(prior_pdf) .* maximum(hl.weights),x;color=:black, label="prior",linewidth=1)

    stairs!(ax_l, hl.weights, hl.edges[1][1:end-1].+step(hl.edges[1]); color)


    ax_a = Axis(
        fig[1,1],
    )
    hideydecorations!(ax_a)
    hidexdecorations!(ax_a)

    # Also plot the prior pdf scaled to the same height
    # x = range(extrema(range1)..., length=250)
    # prior_pdf = pdf.(priors.a, x)
    # Makie.lines!(ax_a, x, prior_pdf ./ maximum(prior_pdf) .* maximum(ha.weights), color=:black, label="prior",linewidth=1)
    
    stairs!(ax_a, ha; color)

    linkxaxes!(ax, ax_a)
    linkyaxes!(ax, ax_l)
    xlims!(ax, extrema(range1))
    ylims!(ax, extrema(range2))
    xlims!(ax_l, low=0, high=maximum(hl.weights)*1.1)
    ylims!(ax_a, low=0, high=maximum(ha.weights)*1.1)
    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)

    rowgap!(fig.layout,0)
    colgap!(fig.layout,0)

    # colsize!(fig.layout, 2, )
    # rowsize!(fig.layout, 1, 60)
    rowsize!(fig.layout, 2, Auto(3))
    colsize!(fig.layout, 2, Aspect(1,1.0))
    
    
    # Legend(fig[1,2], ax_a, tellheight=false, tellwidth=false, valign=:top, halign=:right, patchsize=(10,10), framevisible=false)
    return fig
end



# N= 1500
# ii = rand(1:size(chain,1)*size(chain,3),N)
# Zs = chain["b[Z]"][ii]
# Js = chain["b[J]"][ii]
# Ls = chain["b[L]"][ii]
# plot!(vcat(Zs', Js', Ls'), color=:black, lw=1,  alpha=0.05, label="")
# # scatter!(Ls', color=:black, lw=1,  alpha=0.05, label="")
# scatter!(vcat(Zs', Js', Ls'), color=:black, lw=1, ms=4, alpha=0.01, label="")
# s

##


N = 500
ii = rand(1:size(chains,1)*size(chains,3), N)

fig =  Figure(resolution=(1600,1400))
ax_sky = Axis(fig[1,1:2], xreversed=true,autolimitaspect=1, xlabel="ΔRA (mas)", ylabel="ΔDEC (mas)")

if !isnothing(chains.info.model.images)
    plot_image!(ax_sky, chains.info.model.images.image[end], chains.info.model.images.platescale[end])
end

plot_orbits!(ax_sky, chains, :b; color=chains["b[mass]"], ii)
plot_astrom!(ax_sky, chains.info.model)

# pma_grid = GridLayout()
# fig.layout[2,1:3] = pma_grid

# ax_pma1 = Axis(pma_grid[1:2,1], xlabel="Δv_ra (mas/yr)", ylabel="Δv_dec (mas/yr)", title="HIPPARCOS", autolimitaspect=1,)
# plot_pma_epoch!(ax_pma1, chains, 1, color=chains["b[mass]"],  N=N)
# plot_pma_epoch!(ax_pma1, chains.info.model, 1)

# ax_pma_time_ra = Axis(pma_grid[1,2], xlabel="epoch (mjd)", ylabel="Δv_ra (mas/yr)")
# plot_pma_time!(ax_pma_time_ra, chains, :ra, color=chains["b[mass]"], N=N)
# plot_pma_time!(ax_pma_time_ra, chains.info.model, :ra)

# ax_pma_time_dec = Axis(pma_grid[2,2], xlabel="epoch (mjd)", ylabel="Δv_dec (mas/yr)")
# plot_pma_time!(ax_pma_time_dec, chains, :dec, color=chains["b[mass]"], N=N)
# plot_pma_time!(ax_pma_time_dec, chains.info.model, :dec)

# ax_pma2 = Axis(pma_grid[1:2,3],xlabel="Δv_ra (mas/yr)", ylabel="Δv_dec (mas/yr)", title="GAIA", autolimitaspect=1,)
# plot_pma_epoch!(ax_pma2, chains, 2, color=chains["b[mass]"], N=N)
# plot_pma_epoch!(ax_pma2, chains.info.model, 2)

ax_pma1 = Axis(fig[2,1], xlabel="Δv_ra (mas/yr)", ylabel="Δv_dec (mas/yr)", title="HIPPARCOS", autolimitaspect=1,)
plot_pma_epoch!(ax_pma1, chains, 1; color=chains["b[mass]"], ii)
plot_pma_epoch!(ax_pma1, chains.info.model, 1)

ax_pma2 = Axis(fig[2,2],xlabel="Δv_ra (mas/yr)", ylabel="Δv_dec (mas/yr)", title="GAIA", autolimitaspect=1,)
plot_pma_epoch!(ax_pma2, chains, 2; color=chains["b[mass]"], ii)
plot_pma_epoch!(ax_pma2, chains.info.model, 2)

ax_pma_time_ra = Axis(fig[3,1:2], xlabel="epoch (mjd)", ylabel="Δv_ra (mas/yr)")
plot_pma_time!(ax_pma_time_ra, chains, :ra; color=chains["b[mass]"], ii)
plot_pma_time!(ax_pma_time_ra, chains.info.model, :ra)

ax_pma_time_dec = Axis(fig[4,1:2], xlabel="epoch (mjd)", ylabel="Δv_dec (mas/yr)")
plot_pma_time!(ax_pma_time_dec, chains, :dec; color=chains["b[mass]"], ii)
plot_pma_time!(ax_pma_time_dec, chains.info.model, :dec)

# Histograms

ax_hist = Axis(
    fig[1,3],
    ylabel="posterior density",
    xlabel="mass (Mjup)",
)
hideydecorations!(ax_hist)
hist_prop!(ax_hist, chains["b[mass]"], color=:black)

ax_hist = Axis(
    fig[2,3],
    ylabel="posterior density",
    xlabel="a (au)",
)
hideydecorations!(ax_hist)
hist_prop!(ax_hist, chains["b[a]"], color=:black)

ax_hist = Axis(
    fig[3,3],
    ylabel="posterior density",
    xlabel="e",
)
hideydecorations!(ax_hist)
hist_prop!(ax_hist, chains["b[e]"], color=:black)

plot_phot!(fig[4,3], chains.info.model.planets[1])

linkaxes!(ax_pma1, ax_pma2)
linkaxes!(ax_pma_time_ra, ax_pma_time_dec)

rowsize!(fig.layout, 1, Auto(2))
# rowsize!(fig.layout, 2, Auto(2))
# rowsize!(fig.layout, 3, Auto(1))
# colsize!(fig.layout, 1, Auto(2/3))
# colsize!(fig.layout, 2, Auto(1/3))

fig
##
