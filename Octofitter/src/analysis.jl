# This file contains funtions for analysing results from chains.
using Dates
using Measurements

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

# """
#     sampleorbits(chains, planet_key, 100)

# Given the posterior for a particular planet in the model, and the number of orbits to sample,
# return a random subset of length N.
# """
# function sampleorbits(chains, planetnum, N)
#     return map(rand(eachindex(chains.planets[planetnum].a), N)) do i
#         return VisualOrbit(;
#             a = chains.planets[planetnum].a[i],
#             i = chains.planets[planetnum].i[i],
#             e = chains.planets[planetnum].e[i],
#             ω = chains.planets[planetnum].ω[i],
#             Ω = chains.planets[planetnum].Ω[i],
#             μ = chains.μ[i],
#             plx = chains.plx[i],
#             τ = chains.planets[planetnum].τ[i],
#         )
#     end
# end
# export sampleorbits



function bayesfactor(chain, planet, property)
    prior = chain.info.model.planets[planet].priors.priors[property]
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


"""
    plotchains(
        chain, planet_key;
        N=1500,
        ii = rand(1:size(chain,1)*size(chain,3), N),
        color=length(chain.info.model.planets) == 0 || !haskey(chain, string(planet_key)*"_a") ? nothing : string(planet_key)*"_a",
        colorbartitle=color,
        clims=nothing,
        cmap=:plasma,
        alpha=30/length(ii),
        attime=nothing,
        kwargs...,
    )

Draw samples from a posterior chain for a given planet given by name `planet_key` and visualize them in some way.
Use `kind` to control what plot is made. A few options: :astrometry, :radvel, :trueanom, :meananom, :eccanom, :x, :y, :z, (:x, :y), :raoff, :decoff, :pmra, :pmdec, :accra, :accdec, :radvel, :posangle, :projectedseparation.
See PlanetOrbits documentation for more details.

Inputs:
* chain                   The chain to draw from
* planet_key              Planet name in the model (symbol)
* N=1500                  Number of samples to draw for the plot 
* kind=nothing            Specify what kind of plot to make. 
* ii=...                  Specific row numbers to use, if you want to e.g. plot the same 100 samples in a few different plots
* color="\$planet_key_a"   Column name to to map colors to. Semi-major axis by default but can be any column or an arbitrary array.
* colorbartitle=color     Name for colourbar
* clims=nothing           Tuple of colour limits (min and max)
* cmap=:plasma            Colormap
* alpha=...               Transparency of the lines
"""
function plotchains end
"""
    plotchains!(
        plot, ...;
        kwargs...,
    )

See plotchains
"""
function plotchains! end
export plotchains, plotchains!

function plotposterior end
export plotposterior
function plotposterior! end
export plotposterior!
function plotmodel end
export plotmodel
function plotmodel! end
export plotmodel!
export timeplot
export timeplotgrid
function plotpos3d end
function plotpos3d! end
export plotpos3d
export plotpos3d!
function plotposviews end
export plotposviews


# Optionally depend on Plots. We can't get the control we need using Plots
# recipes alone, unfortunately.
using Requires
function init_plots()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin


        function plotchains(args...;kwargs...)
            p = Plots.plot()
            plotchains!(p, args...;kwargs...)
        end
        function plotchains!(
            p, chain, planet_key;
            N=1500,
            ii = rand(1:size(chain,1)*size(chain,3), N),
            color=length(chain.info.model.planets) == 0 || !haskey(chain, string(planet_key)*"_a") ? nothing : string(planet_key)*"_a",
            colorbartitle=color,
            clims=nothing,
            cmap=:plasma,
            alpha=30/length(ii),
            attime=nothing,
            mass=nothing,
            kwargs...,
        )
            # Construct orbits from chain
            orbits = Octofitter.construct_elements(chain, planet_key, ii)
            if !isnothing(attime)
                orbits = orbitsolve.(orbits, attime)
            end

            # Handle the colours of each element
            if isnothing(color)
                colours = nothing
            elseif typeof(color) <: Union{Symbol,AbstractString}
                if haskey(chain, Symbol(color))
                    colours = chain[color][ii]
                else
                    colours = nothing
                    kwargs = (;kwargs..., color)
                end
            elseif typeof(color) <: AbstractArray # They passed in the values directly
                colours = color[ii]
            else
                error("unkown color option")
            end

            if !isnothing(colours)
                if isnothing(clims)
                    clims=quantile(Iterators.flatten(
                        extrema(colours)
                    ), (0.01,0.99))
                end
                minc = first(clims)
                dc = last(clims) - minc
                clims = (minc, dc+minc)
            end


            cmap= Plots.cgrad(cmap)
            for i in eachindex(orbits)
                if !isnothing(colours)
                    c = colours[i]
                    kwargs_loc = (;color=cmap[(c-minc)/dc], kwargs...)
                else
                    kwargs_loc = kwargs
                end
                m = nothing
                if !isnothing(mass)
                    if length(mass) > 1
                        m = mass[ii][i]*mjup2msol
                    else
                        m = mass*mjup2msol
                    end
                end
                Plots.plot!(p, orbits[i]; label="", alpha, mass=m, kwargs_loc...)
            end
            # Colorbar
            # Plots.scatter!(p, [NaN],[NaN]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)
            Plots.scatter!(p, [0],[0]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)
            Plots.scatter!(p, [0],[0], marker=(:star, :black, 5), markerstrokewidth=1, markerstrokecolor=:white,  label="")
            return p
            # # Star at centre

        end




        # function plotposterior(args...;kwargs...)
        #     Plots.plot()
        #     plotposterior!(args...;kwargs...)
        # end
        # function plotposterior!(
        #     chain,
        #     planet_key,
        #     keys::Nothing;
        #     N=1500,
        #     ii = rand(1:size(chain,1)*size(chain,3), N),
        #     alpha=0.02,
        #     lw=0.3,
        #     color=1,
        #     kwargs...
        # )
        #     elements = construct_elements(chain, planet_key, ii)

        #     Plots.plot!(;
        #         size=(700,700),
        #         dpi=200,
        #         fmt=:png,
        #         framestyle=:box,
        #         grid=:none,
        #         minorticks=true,
        #         aspectratio=1,
        #         fontfamily="Arial",
        #         margin=8Plots.mm,
        #         kwargs...
        #     )        
        #     for i in eachindex(elements)
        #         Plots.plot!(elements[i], label="",color=color, lw=lw, alpha=alpha,)
        #     end
        #     Plots.scatter!([0],[0], marker=(:star,:black, 5), markerstrokewidth=1, markerstrokecolor=:white, label="")
        # end
        # function plotposterior!(
        #     chain,
        #     planet_key,
        #     property;
        #     N=1500,
        #     ii = rand(1:size(chain,1)*size(chain,3), N),
        #     alpha=0.02,
        #     cmap=:turbo,
        #     rev=true,
        #     colorbartitle= property isa Symbol ? string(property) : "",#"semi-major axis (au)",
        #     clims=nothing,
        #     lw=0.3,
        #     kwargs...
        # )
        #     elements = construct_elements(chain, planet_key, ii)
        #     N = length(ii)

        #     Plots.plot!(;
        #         size=(550,550),
        #         dpi=200,
        #         fmt=:png,
        #         framestyle=:box,
        #         grid=:none,
        #         minorticks=true,
        #         aspectratio=1,
        #         fontfamily="Arial",
        #         margin=8Plots.mm,
        #         kwargs...
        #     )
        #     Plots.xlabel!("Δ right ascension (mas)")
        #     Plots.ylabel!("Δ declination (mas)")

        #     if typeof(property) <: Union{Symbol,AbstractString}
        #         # k = string(planet_key)*"[$property]"
        #         k = property
        #         colours = chain[k][ii]
        #     else # Assume they passed in the values directly
        #         colours = property[ii]
        #     end

        #     if isnothing(clims)
        #         clims = extrema(colours)
        #     end
            
        #     minc = first(clims)
        #     dc = last(clims) - minc
        #     clims = (minc, dc+minc)

        #     # Lines to periastron
        #     if !haskey(kwargs, :kind)  || kwargs[:kind] == :astrometry
        #         xs = hcat(zeros(N), raoff.(elements, periastron.(elements)), fill(NaN, N))'[:]
        #         ys = hcat(zeros(N), decoff.(elements, periastron.(elements)), fill(NaN, N))'[:]
        #         Plots.plot!(xs, ys; lw, alpha, color="#777", label="")
        #     end
        #     cmap= Plots.cgrad(cmap; rev)
        #     for i in eachindex(colours)#sortperm(colours, rev=false)
        #         c = colours[i]
        #         Plots.plot!(elements[i]; label="",color=cmap[(c-minc)/dc], lw, alpha, kwargs...)
        #     end
        #     # Star at centre
        #     Plots.scatter!([0],[0], marker=(:star, :black, 5), markerstrokewidth=1, markerstrokecolor=:white,  label="")
        #     # Colorbar
        #     Plots.scatter!([0],[0]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)
        # end

        # function plotpos3d(args...;kwargs...)
        #     Plots.plot()
        #     plotpos3d!(args...;kwargs...)
        # end
        # function plotpos3d!(
        #     chain,
        #     planet_key,
        #     property;
        #     N=1500,
        #     ii = rand(1:size(chain,1)*size(chain,3), N),
        #     alpha=0.05,
        #     cmap=:turbo,
        #     rev=true,
        #     colorbartitle= property isa Symbol ? string(property) : "",#"semi-major axis (au)",
        #     clims=nothing,
        #     lw=0.3,
        #     kwargs...
        # )
        #     elements = construct_elements(chain, planet_key, ii)
        #     N = length(ii)

        #     Plots.plot!(;
        #         size=(550,550),
        #         dpi=200,
        #         fmt=:png,
        #         framestyle=:box,
        #         grid=:none,
        #         minorticks=true,
        #         aspectratio=1,
        #         fontfamily="Arial",
        #         margin=8Plots.mm,
        #         kwargs...
        #     )

        #     if typeof(property) <: Union{Symbol,AbstractString}
        #         # k = string(planet_key)*"[$property]"
        #         k = property
        #         colours = chain[k][ii]
        #     else # Assume they passed in the values directly
        #         colours = property[ii]
        #     end

        #     if isnothing(clims)
        #         clims = extrema(colours)
        #     end
            
        #     minc = first(clims)
        #     dc = last(clims) - minc
        #     clims = (minc, dc+minc)

        #     # Lines to periastron
        #     xs = hcat(zeros(N), PlanetOrbits.posx.(elements, periastron.(elements)), fill(NaN, N))'[:]
        #     ys = hcat(zeros(N), PlanetOrbits.posy.(elements, periastron.(elements)), fill(NaN, N))'[:]
        #     zs = hcat(zeros(N), PlanetOrbits.posz.(elements, periastron.(elements)), fill(NaN, N))'[:]
        #     Plots.plot!(xs, ys, zs; lw, alpha, color="#777", label="")
        #     cmap= Plots.cgrad(cmap; rev)
        #     for i in eachindex(colours)
        #         c = colours[i]
        #         Plots.plot!(elements[i]; label="",color=cmap[(c-minc)/dc], lw, alpha, kind=(:x,:z,:y))
        #     end
        #     # Star at centre
        #     Plots.scatter!([0],[0],[0], marker=(:star, :black, 5), markerstrokewidth=1, markerstrokecolor=:white,  label="")
        #     # Colorbar
        #     Plots.scatter!([0],[0],[0]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)
        # end

        # function plotposviews(
        #     chain, planet_key;
        #     N=1500,
        #     alpha=(N <= 15 ? 1 : 30/N),
        #     ii = rand(1:size(chain,1)*size(chain,3), N),
        #     color=length(chain.info.model.planets) == 0 ? nothing : string(planet_key)*"_a",
        #     colorbartitle=color,
        #     cmap=:plasma,
        #     clims=isnothing(color) ? nothing : quantile(Iterators.flatten(
        #         extrema(chain[color])
        #         ), (0.01,0.99)),
        #     kwargs...,
        # )
        #     p1 = plotposterior(
        #         chain, planet_key, color;
        #         ii, lw=1, alpha,
        #         cmap, rev=false,
        #         clims,
        #         colorbartitle,
        #         kind=(:x,:y),
        #         size=(500,500),
        #         margin=0Plots.mm,
        #         kwargs...
        #     )
        #     if length(chain.info.model.planets[planet_key].observations) > 0 && chain.info.model.planets[planet_key].observations[1] isa Astrometry
        #         el = Octofitter.construct_elements(chain, :b, 1)
        #         scalefacty = PlanetOrbits.posy(el,0) / decoff(el, 0)
        #         scalefactx = PlanetOrbits.posx(el,0) / raoff(el, 0)
        #         Plots.scatter!(
        #             p1,
        #             scalefactx .* chain.info.model.planets[planet_key].observations[1].table.ra,
        #             scalefacty .* chain.info.model.planets[planet_key].observations[1].table.dec,
        #             color=:black,label=""
        #         )
        #     end
        #     p2 = plotposterior(
        #         chain, planet_key, color;
        #         ii, lw=1, alpha,
        #         cmap, rev=false,
        #         clims,
        #         colorbartitle,
        #         kind=(:z,:y),
        #         size=(500,500),
        #         margin=0Plots.mm,
        #         kwargs...
        #     )
        #     if length(chain.info.model.planets[planet_key].observations) > 0 && chain.info.model.planets[planet_key].observations[1] isa Astrometry
        #         el = Octofitter.construct_elements(chain, :b, 1)
        #         Plots.hline!(p2, scalefacty .* chain.info.model.planets[planet_key].observations[1].table.dec, color=:black,label="")
        #     end
        #     p3 = plotposterior(
        #         chain, planet_key, color;
        #         ii, lw=1, alpha,
        #         cmap, rev=false,
        #         clims,
        #         colorbartitle,
        #         kind=(:x,:z),
        #         size=(500,500),
        #         kwargs...
        #     )
        #     if length(chain.info.model.planets[planet_key].observations) > 0 && chain.info.model.planets[planet_key].observations[1] isa Astrometry
        #         el = Octofitter.construct_elements(chain, :b, 1)
        #         Plots.vline!(p3, scalefactx .* chain.info.model.planets[planet_key].observations[1].table.ra, color=:black,label="")
        #     end
        #     # Set limits correctly
        #     Plots.ylims!(p2, Plots.ylims(p1))
        #     Plots.xlims!(p3, Plots.xlims(p1))
        #     Plots.ylims!(p3, Plots.xlims(p2))


        #     Plots.plot(p1,p2,p3; layout=(2,2))#, size=(650,650))
        # end
        
        # function plotmodel(args...; kwargs...)
        #     Plots.plot()
        #     plotmodel!(args...;kwargs...)
        # end
        # function plotmodel!(images::Images,i; kwargs...)
        #     img = DirectImages.DirectImage(images.image[i])
        #     img.PLATESCALE = images.platescale[i]
        #     imshow!(img;skyconvention=true,kwargs...)
        #     # plot!(sampled,color=:white,alpha=5/N,label="")
        #     Plots.xlims!(img.PLATESCALE.*extrema(axes(img,1)))
        #     Plots.ylims!(img.PLATESCALE.*extrema(axes(img,2)))

        # end

        # function plotmodel!(
        #     chain,
        #     N=1500,
        #     system=chain.info.model;
        #     alpha=(N <= 15 ? 1 : 30/N),
        #     ii = rand(1:size(chain,1)*size(chain,3), N),
        #     color=length(system.planets) == 0 ? nothing : string(first(keys(system.planets)))*"_a",
        #     colorbartitle=color,
        #     plotpma=!isnothing(propermotionanom(system)),
        #     # TODO: ideally this is based on if there is a mass variable
        #     plotmass=!isnothing(propermotionanom(system)),
        #     plotimages=!isnothing(images(system)),
        #     cmap=:plasma,
        #     imagecmap=:Greys,
        #     clims=isnothing(color) ? nothing : quantile(Iterators.flatten(
        #         # extrema(planet.a)
        #         # extrema([sample.planets[key].a for sample in chain])
        #         # for key in keys(chain[1].planets)
        #         extrema(chain[color])
        #     ), (0.01,0.99)),
        #     lims=nothing,
        #     kwargs...
        # )

        #     # Plot orbits
        #     p_orbits = Plots.plot(
        #         xlims=:symmetric,
        #         ylims=:symmetric,
        #     )

        #     # plot images?
        #     if plotimages
        #         if length(unique(images(system).table.platescale)) != 1
        #             @warn "Plotmodel does not yet support images with multiple platescales. Taking first image only."
        #             img = first(images(system).image)
        #         else
        #             img = mean(cat(images(system).table.image...,dims=3),dims=3)[:,:,1]
        #         end
        #         # # To get the colour scales to work out
        #         # img ./= only(quantile(filter(isfinite, arraydata(img)), [0.98]))#[0.9995]))
        #         # img .*= maximum(clims) - minimum(clims)
        #         # img .+= minimum(clims)
        #         # implot!(img; color=imagecmap, skyconvention=true, lims)
        #     end
            
        #     for planet_key in eachindex(system.planets)
        #         plotposterior!(
        #             chain, planet_key, color; ii, lw=1, alpha=alpha,
        #             cmap=cmap, rev=false,
        #             # cmap=:turbo,
        #             clims=clims,
        #             colorbartitle=colorbartitle,
        #             kwargs...
        #         )

        #         # Planet astrometry?
        #         astrom = Octofitter.astrometry(system.planets[planet_key])
        #         if !isnothing(astrom)
        #             # Put black error bars over thicker white error bars to make them easier to see
        #             # when plotted over e.g. images, orbits.
        #             Plots.scatter!(astrom,marker=(:white,:circle,0),label="",linewidth=0,color=:white,markerstrokewidth=3, markerstrokecolor=:white)
        #             Plots.scatter!(astrom,marker=(:black,:circle,0),label="",linewidth=0,color=:black,markerstrokewidth=1, markerstrokecolor=:black)
        #         end
        #     end
        #     # We will override this if we have more information available further down
        #     final_plot = p_orbits

        #     # astrometric acceleration?
        #     if plotpma && !isnothing(propermotionanom(system))
        #         pma_plots = pmaplot(chain; ii, color)
                
        #         # # Crete new final plot
        #         l = eval(:(Plots.@layout [
        #             A{0.65h}
        #             B
        #         ]))
        #         final_plot = Plots.plot(p_orbits, pma_plots, layout=l, size=(550,750),
        #             margin=5Plots.mm, guidefontsize=9, titlefontsize=9, top_margin=0Plots.mm)
        #     end


        #     if plotmass
        #         p = Plots.plot(;xlabel="", linewidth=3, label="")
        #         for p in keys(system.planets)
        #             m = vec(chain["$p_mass"])
        #             h = fit(Histogram, m, nbins=round(Int,sqrt(length(m))))#, nbins=100)
        #             Plots.plot!(h.edges[1][begin:end-1],seriestype=:step, h.weights, label="$p_mass",color=1,lw=2)
        #         end

        #         layout = eval(:(Plots.@layout [
        #             A{0.8h}
        #             B
        #         ]))
        #         final_plot = Plots.plot(final_plot, p, layout=layout, size=(550,850))
        #     end

        #     return final_plot
        # end


        """
        Plot a parameter against time.

        Example:
        ```
        plot(
            timeplot(chains_pma, :b, :mass, :ra),
            timeplot(chains_pma, :b, :mass, :dec),
            timeplot(chains_pma, :b, :mass, :rv),
            timeplot(chains_pma, :b, :mass, :pmra),
            timeplot(chains_pma, :b, :mass, :pmdec),
            layout = @layout([
                A
                B
                C
                D E
            ]),
            framestyle=:box,
            size=(500,700)
        )
        ```
        """
        function timeplot(args...; kwargs...)
            Plots.plot()
            timeplot!(args...; kwargs...)
        end
        function timeplot!(
            chain,
            planet_key,
            color,
            prop;
            N = 1500,
            ii = rand(1:size(chain,1)*size(chain,3), N),
            alpha = 0.05,
            cmap=:plasma,
            clims = quantile(vec(chain[color]),(0.01, 0.99)),
            kwargs...
        )
            # Handle the colours of each element
            if typeof(color) <: Union{Symbol,AbstractString}
                color = chain[color][ii]
            else # They passed in the values directly
                color = color[ii]
            end
            planet = chain.info.model.planets[planet_key]
            if prop == :ra
                ylabel = "RA (mas)"
            elseif prop == :dec
                ylabel = "DEC (mas)"
            elseif prop == :sep
                ylabel = "SEP (mas)"
            elseif prop == :pa
                ylabel = "PA (°)"
            elseif prop == :rv
                ylabel = "RV m/s"
            elseif prop == :pmra
                ylabel = "∂RA mas/yr"
            elseif prop == :pmdec
                ylabel = "∂DEC mas/yr"
            else
                error("Unsupported property. Choose :ra, :dec, :sep, :pa, :rv, :pmra, or :pmdec")
            end
            p1 = Plots.plot!(;
                ylabel,
                legend=:none
            )
            xerr = nothing
            if planet_key isa String || planet_key isa Symbol
                planet_keys = tuple(planet_key)
            else
                planet_keys = planet_key
            end
            all_epochs_planet = mapreduce(vcat, keys(chain.info.model.planets)) do planet_key
                planet = chain.info.model.planets[planet_key]
                reduce(vcat, Any[hasproperty(t.table, :epoch) ? t.table.epoch : t.table.ra_epoch for t in planet.observations if hasproperty(t.table, :epoch) || hasproperty(t.table, :ra_epoch)], init=[])
            end
            all_epochs_star = mapreduce(vcat, chain.info.model.observations, init=Float64[]) do obs
                if hasproperty(obs.table, :epoch)
                    return obs.table.epoch
                else
                    return Float64[]
                end
            end
            all_epochs = [all_epochs_planet; all_epochs_star]
            if isempty(all_epochs) 
                all_epochs = mjd() .+ [-365*2, +365*2]
            end
            t = range((extrema(all_epochs) .+ [-365, 365])..., length=100)
            # t = range((extrema(all_epochs) .+ [-365, 365])..., length=500)
            y = nothing
            ribbon = nothing
            zcolor = nothing
            if prop == :ra
                if !isnothing(astrometry(planet))
                    y = astrometry(planet).table.ra
                    yerr = astrometry(planet).table.σ_ra
                    x = astrometry(planet).table.epoch
                end
                elements = Octofitter.construct_elements(chain, planet_key, ii)
                fit = raoff.(elements, t')
            elseif prop == :dec
                if !isnothing(astrometry(planet))
                    y = astrometry(planet).table.dec
                    yerr = astrometry(planet).table.σ_dec
                    x = astrometry(planet).table.epoch
                end
                elements = Octofitter.construct_elements(chain, planet_key, ii)
                fit = decoff.(elements, t')
            elseif prop == :sep
                if !isnothing(astrometry(planet))
                    if hasproperty(astrometry(planet).table, :sep)
                        y = astrometry(planet).table.sep
                        yerr = astrometry(planet).table.σ_sep
                    else
                        xx = astrometry(planet).table.ra
                        yy = astrometry(planet).table.dec
                        xxerr = astrometry(planet).table.σ_ra
                        yyerr = astrometry(planet).table.σ_dec
                        y_unc = sqrt.((xx .± xxerr).^2 .+ (yy .± yyerr).^2)
                        y = Measurements.value.(y_unc)
                        yerr = Measurements.uncertainty.(y_unc)
                    end
                    x = astrometry(planet).table.epoch
                end
                elements = Octofitter.construct_elements(chain, planet_key, ii)
                fit = projectedseparation.(elements, t')
                if Symbol("$(planet_key)_mass") in keys(chain)
                    fit  = fit# .- projectedseparation.(elements, t', chain["$(planet_key)_mass"][ii].*Octofitter.mjup2msol)
                end
            elseif prop == :pa
                if !isnothing(astrometry(planet))
                    if hasproperty(astrometry(planet).table, :pa)
                        y = rad2deg.(rem2pi.(astrometry(planet).table.pa, RoundNearest))
                        yerr = rad2deg.(astrometry(planet).table.σ_pa)
                    else
                        xx = astrometry(planet).table.ra
                        yy = astrometry(planet).table.dec
                        xxerr = astrometry(planet).table.σ_ra
                        yyerr = astrometry(planet).table.σ_dec
                        y_unc = atand.((xx .± xxerr), (yy .± yyerr))
                        y = Measurements.value.(y_unc)
                        yerr = Measurements.uncertainty.(y_unc)
                    end
                    x = astrometry(planet).table.epoch
                end
                elements = Octofitter.construct_elements(chain, planet_key, ii)
                fit = rad2deg.(posangle.(elements, t'))
                # Avoid ugly wrapping
                diffs = diff(fit, dims=1)
                fit[findall(abs.(diffs) .> 100)] .= NaN
            elseif prop == :pmra
                if !isnothing(propermotionanom(chain.info.model))
                    hgca = propermotionanom(chain.info.model).table[1]
                    t = range(years2mjd(hgca.epoch_ra_hip)-365*5, years2mjd(hgca.epoch_ra_gaia)+365*5, length=100)
                    y = [hgca.pmra_hip, hgca.pmra_hg, hgca.pmra_gaia]
                    yerr = [hgca.pmra_hip_error, hgca.pmra_hg_error, hgca.pmra_gaia_error]
                    x = years2mjd.([hgca.epoch_ra_hip, (hgca.epoch_ra_hip + hgca.epoch_ra_gaia)/2, hgca.epoch_ra_gaia])
                    xerr = [4*365/2, 25*365/2, 3*365/2]
                end
                fit = 0
                for planet_key in planet_keys
                    elements = Octofitter.construct_elements(chain, planet_key, ii)
                    fit = fit .+ pmra.(elements, t', collect(chain["$(planet_key)_mass"][ii]).*Octofitter.mjup2msol)
                end
                if :pmra in keys(chain)
                    fit = fit .+ chain["pmra"][ii]
                end
            elseif prop == :pmdec
                if !isnothing(propermotionanom(chain.info.model))
                    hgca = propermotionanom(chain.info.model).table[1]
                    t = range(years2mjd(hgca.epoch_dec_hip)-365*5, years2mjd(hgca.epoch_dec_gaia)+365*5, length=100)
                    y = [hgca.pmdec_hip, hgca.pmdec_hg, hgca.pmdec_gaia]
                    yerr = [hgca.pmdec_hip_error, hgca.pmdec_hg_error, hgca.pmdec_gaia_error]
                    x = years2mjd.([hgca.epoch_dec_hip, (hgca.epoch_dec_hip + hgca.epoch_dec_gaia)/2, hgca.epoch_dec_gaia])
                    xerr = [4*365/2, 25*365/2, 3*365/2]
                end
                fit = 0
                for planet_key in planet_keys
                    elements = Octofitter.construct_elements(chain, planet_key, ii)
                    fit = fit .+ pmdec.(elements, t', collect(chain["$(planet_key)_mass"][ii]).*Octofitter.mjup2msol)
                end
                if :pmdec in keys(chain)
                    fit = fit .+ chain["pmdec"][ii]
                end
            elseif prop == :rv
                # ribbon = chain[jitter_inst][ii]
                fit = 0
                for planet_key in planet_keys
                    elements = Octofitter.construct_elements(chain, planet_key, ii)
                    fit = fit .+ radvel.(elements, t', collect(chain["$(planet_key)_mass"][ii]).*mjup2msol)
                end
                for obs in chain.info.model.observations
                    # TODO: make this pluggable instead of this hacky workaround
                    if startswith(string(typeof(obs)), "RadialVelocity")
                        if haskey(chain,:rv0_1)
                            barycentric_rv_inst_1 = median(vec(chain["rv0_1"]))
                            jitter = barycentric_rv_inst_1 = median(vec(chain["jitter_1"]))
                        end
                        if haskey(chain,:rv0_2)
                            barycentric_rv_inst_2 = median(vec(chain["rv0_2"]))
                            jitter = barycentric_rv_inst_2 = median(vec(chain["jitter_2"]))
                        end
                        if haskey(chain,:rv0_3)
                            barycentric_rv_inst_3 = median(vec(chain["rv0_3"]))
                            jitter = barycentric_rv_inst_3 = median(vec(chain["jitter_3"]))
                        end
                        if haskey(chain,:rv0_4)
                            barycentric_rv_inst_4 = median(vec(chain["rv0_4"]))
                            jitter = barycentric_rv_inst_4 = median(vec(chain["jitter_4"]))
                        end

                        idxes = length(unique(obs.table.inst_idx))
                        y = [[] for _ in 1:idxes]
                        x = [[] for _ in 1:idxes]
                        yerr = [[] for _ in 1:idxes]
                        for row in obs.table
                            if row.inst_idx == 1
                                barycentric_rv_inst = barycentric_rv_inst_1
                            elseif row.inst_idx == 2
                                barycentric_rv_inst = barycentric_rv_inst_2
                            elseif row.inst_idx == 3
                                barycentric_rv_inst = barycentric_rv_inst_3
                            elseif row.inst_idx == 4
                                barycentric_rv_inst = barycentric_rv_inst_4
                            end
                            push!(x[row.inst_idx], row.epoch)
                            push!(y[row.inst_idx], row.rv - barycentric_rv_inst)
                            push!(yerr[row.inst_idx], row.σ_rv + jitter)
                        end
                    end
                end
            end
            for i in axes(fit, 1)
                c = (color[i] - first(clims))/(last(clims)-first(clims))
                Plots.plot!(
                    mjd2date.(t), view(fit, i, :);
                    # line_z=chain[color][ii],
                    color=Plots.cgrad(cmap)[c],
                    alpha,
                    fillalpha=alpha/2,
                    ribbon=isnothing(ribbon) ? nothing : ribbon[i],
                    kwargs...,
                )
            end
            if !isnothing(y)
                if eltype(y) <: AbstractArray
                    for j in eachindex(y)
                        Plots.scatter!(
                            p1,
                            mjd2date.(x[j]),
                            y[j];
                            yerr=!isnothing(yerr) ? yerr[j] : nothing,
                            xerr=!isnothing(xerr) ? xerr[j] : nothing,
                            markerstrokewidth=1.5,
                            markerstrokecolor=:auto,
                            markersize=2.5,
                            color=j
                        )
                    end
                else
                    Plots.scatter!(
                        p1,
                        mjd2date.(x),
                        y; yerr, xerr,
                        color=1,
                        markerstrokewidth=1.5,
                        markerstrokecolor=:auto,
                        markersize=2.5,
                    )
                end
            end
            p1
        end

        function timeplotgrid(
            chains;
            color = "$(first(keys(chains.info.model.planets)))_e",
            clims = quantile(vec(chains[color]),(0.01, 0.99)),
            cmap = :plasma,
            alpha=0.1,
            N=1500,
            ii = rand(1:size(chains,1)*size(chains,3), N),
            kwargs...
        )
            kwargs = (; 
                ii,
                clims,
                cmap,
                alpha,
                kwargs...
            )
            planet_keys = keys(chains.info.model.planets)

            ppost = Plots.plot()
            for (i,planet_key) in enumerate(planet_keys)
                # plotchains!(ppost, chains, planet_key; color, body=:primary, mass=chains["$(planet_key)_mass"], rev=false, colorbar=nothing, kwargs...)
                plotchains!(ppost, chains, planet_key; color, rev=false, colorbar=nothing, kwargs...)
                Plots.scatter!([0],[0],marker=(:star, :white, :black, 5),label="")
                astrom = astrometry(chains.info.model.planets[planet_key])
                if !isnothing(astrom) && hasproperty(astrom.table, :ra)
                    els = construct_elements(chains, planet_key, :)
                    sols = orbitsolve.(els, astrom.table.epoch')
                    xs = median(raoff.(sols,chains["$(planet_key)_mass"].*mjup2msol), dims=1)'
                    ys = median(decoff.(sols,chains["$(planet_key)_mass"].*mjup2msol), dims=1)'
                    Plots.scatter!(
                        ppost,
                        astrom.table.ra  .- xs,
                        astrom.table.dec .- ys,
                        label="",
                        color=i,
                        markerstrokewidth=1.5,
                        markerstrokecolor=:auto,
                        markersize=3.5,
                    )
                end
            end
            psep = Plots.plot()
            for planet_key in planet_keys
                timeplot!(chains, planet_key, color, :sep; clims, kwargs...)
            end
            ppa = Plots.plot()
            for planet_key in planet_keys
                timeplot!(chains, planet_key, color, :pa; clims, kwargs...)
            end
            
            ppmra = Plots.plot()
            timeplot!(chains, planet_keys, color, :pmra;clims,  kwargs...)
            
            ppmdec = Plots.plot()
            timeplot!(chains, planet_keys, color, :pmdec;clims,  kwargs...)
            
            prv = Plots.plot()
            timeplot!(chains, planet_keys, color, :rv; clims, kwargs...)
            

            # Vertical
            # Plots.plot(
            #     psep,
            #     ppa,
            #     ppmra,
            #     ppmdec,
            #     ppost,
            #     prv,
            #     layout = (3,2),
            #     framestyle=:box,
            #     grid=false,
            #     size=(1000,1200),
            #     margin=4Plots.mm
            # )
            # Horizontal
            Plots.plot(
                ppost,
                psep,
                ppa,
                prv,
                ppmra,
                ppmdec,
                layout = (2,3),
                framestyle=:box,
                grid=false,
                size=(1400,900),
                margin=4Plots.mm
            )
        end

        function pmaplot(chains; kwargs...)
            pmaplot(chains, propermotionanom(chains.info.model); kwargs...)
        end
        function pmaplot(chains, pma::ProperMotionAnomHGCA; color, N=1500, ii =rand(1:size(chains,1)*size(chains,3), N), kwargs...)
            
            hgca = pma.table

            # Roughly over what time period were the observations made?
            dt_gaia = 1038 # EDR3: days between  Date("2017-05-28") - Date("2014-07-25")
            dt_hip = 4*365
            # How many points over Δt should we average the proper motion and stellar position
            # at each epoch? This is because the PM is not an instantaneous measurement.
            N_ave = 25
        
            ra_gaia_model = zeros(size(ii))
            dec_gaia_model = zeros(size(ii))
            pmra_gaia_model = zeros(size(ii))
            pmdec_gaia_model = zeros(size(ii))
            ra_hip_model = zeros(size(ii))
            dec_hip_model = zeros(size(ii))
            pmra_hip_model = zeros(size(ii))
            pmdec_hip_model = zeros(size(ii))
            pmra_hg_model = zeros(size(ii))
            pmdec_hg_model = zeros(size(ii))

            planets = keys(chains.info.model.planets)
        
            for (j,jj) in enumerate(ii)
                elements = [
                    Octofitter.construct_elements(chains, p, jj)
                    for p in planets
                ]

                # First epoch: Hipparcos
                for i in eachindex(elements)
                    orbit = elements[i]
                    p = planets[i]
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
                        ra_hip_model[j] += -raoff(o_ra) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
                        dec_hip_model[j] += -decoff(o_dec) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
                        pmra_hip_model[j] += pmra(o_ra, chains["$(p)_mass"][jj]*mjup2msol)
                        pmdec_hip_model[j] += pmdec(o_dec, chains["$(p)_mass"][jj]*mjup2msol)
                    end
                end
                ra_hip_model[j]/=N_ave
                dec_hip_model[j]/=N_ave
                pmra_hip_model[j]/=N_ave
                pmdec_hip_model[j]/=N_ave
            
                # First epoch: Hipparcos
                for i in eachindex(elements)
                    orbit = elements[i]
                    p = planets[i]
                    # Average multiple observations over a timescale +- dt/2
                    # to approximate what HIPPARCOS would have measured.
                    for δt = range(-dt_gaia/2, dt_gaia/2, N_ave)
                        # RA and dec epochs are usually slightly different
                        # Note the unit conversion here from jupiter masses to solar masses to 
                        # make it the same unit as the stellar mass (element.mu)
                        # TODO: we can't yet use the orbitsolve interface here for the pmra calls,
                        # meaning we calculate the orbit 2x as much as we need.
                        o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_gaia[1])+δt)
                        o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_gaia[1])+δt)
                        ra_gaia_model[j] += -raoff(o_ra) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
                        dec_gaia_model[j] += -decoff(o_dec) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
                        pmra_gaia_model[j] += pmra(o_ra, chains["$(p)_mass"][jj]*mjup2msol)
                        pmdec_gaia_model[j] += pmdec(o_dec, chains["$(p)_mass"][jj]*mjup2msol)
                    end
                end
                ra_gaia_model[j]/=N_ave
                dec_gaia_model[j]/=N_ave
                pmra_gaia_model[j]/=N_ave
                pmdec_gaia_model[j]/=N_ave
            
            
                # Model the GAIA-Hipparcos delta-position velocity
                pmra_hg_model[j] = (ra_gaia_model[j] - ra_hip_model[j])/(years2mjd(hgca.epoch_ra_gaia[1]) - years2mjd(hgca.epoch_ra_hip[1]))
                pmdec_hg_model[j] = (dec_gaia_model[j] - dec_hip_model[j])/(years2mjd(hgca.epoch_dec_gaia[1]) - years2mjd(hgca.epoch_dec_hip[1]))
            end
            # # Compute the likelihood at all three epochs (Hipparcos, GAIA-Hip, GAIA)
            # pmra_model = (pmra_hip_model, pmra_hg_model, pmra_gaia_model)
            # pmdec_model = (pmdec_hip_model, pmdec_hg_model, pmdec_gaia_model)
            # pmra_meas = (hgca.pmra_hip[1], hgca.pmra_hg[1], hgca.pmra_gaia[1])
            # pmdec_meas = (hgca.pmdec_hip[1], hgca.pmdec_hg[1], hgca.pmdec_gaia[1])
            # σ_pmra = (hgca.pmra_hip_error[1], hgca.pmra_hg_error[1], hgca.pmra_gaia_error[1])
            # σ_pmdec = (hgca.pmdec_hip_error[1], hgca.pmdec_hg_error[1], hgca.pmdec_gaia_error[1])
            
            # return Plots.plot(
            #     Plots.histogram2d(pmra_hip_model, pmdec_hip_model),
            #     Plots.histogram2d(pmra_hg_model, pmdec_hg_model),
            #     Plots.histogram2d(pmra_gaia_model, pmdec_gaia_model),
            # )

            pma_extrema_x = extrema(vcat(pmra_hip_model, pmra_hg_model, pmra_gaia_model, ))
            pma_extrema_y = extrema(vcat(pmdec_hip_model, pmdec_hg_model, pmdec_gaia_model, ))

            if !isnothing(color)
                subplots = [
                    Plots.scatter(
                        pmra_hip_model .+ chains["pmra"][ii],
                        pmdec_hip_model .+ chains["pmdec"][ii],
                        marker_z=chains[color][ii],
                        title="Hipparcos",
                    ),
                    Plots.scatter(
                        pmra_hg_model .+ chains["pmra"][ii],
                        pmdec_hg_model .+ chains["pmdec"][ii],
                        marker_z=chains[color][ii],
                        title="Gaia-Hipparcos",
                    ),
                    Plots.scatter(
                        pmra_gaia_model .+ chains["pmra"][ii],
                        pmdec_gaia_model .+ chains["pmdec"][ii],
                        marker_z=chains[color][ii],
                        title="Gaia",
                    ),
                ]
            else
                error("Not implemented")
            end

            for p in subplots
                # Plots.scatter!(p, [0], [0], marker=(5, :circle, :red),label="")
                # Plots.hline!(p, [0], color=:black, label="", linewidth=1)
                # Plots.vline!(p, [0], color=:black, label="", linewidth=1)
                Plots.xlabel!(p, "Δμ ra - mas/yr")
            end
            Plots.ylabel!(subplots[1], "Δμ dec - mas/yr")


            pmra_meas = (hgca.pmra_hip[1], hgca.pmra_hg[1], hgca.pmra_gaia[1])
            pmdec_meas = (hgca.pmdec_hip[1], hgca.pmdec_hg[1], hgca.pmdec_gaia[1])
            σ_pmra = (hgca.pmra_hip_error[1], hgca.pmra_hg_error[1], hgca.pmra_gaia_error[1])
            σ_pmdec = (hgca.pmdec_hip_error[1], hgca.pmdec_hg_error[1], hgca.pmdec_gaia_error[1])
            for i in 1:3
                Plots.scatter!(
                    subplots[i],
                    [pmra_meas[i]], [pmdec_meas[i]],
                    xerr=[σ_pmra[i]], yerr=[σ_pmdec[i]],
                    markersize=6,
                    markerstrokewidth=3,
                    color=:white
                )
                Plots.scatter!(
                    subplots[i],
                    [pmra_meas[i]], [pmdec_meas[i]],
                    xerr=[σ_pmra[i]], yerr=[σ_pmdec[i]],
                    markersize=6,
                    markerstrokewidth=1,
                    color=:black
                )
            end

            pma_plot =  Plots.plot(
                subplots...,
                alpha=1,
                legend=false,
                colorbar=false,
                label="",
                markerstrokewidth=0,
                ms=1,
                aspectratio=1,
                layout=(1,3),
                framestyle=:box,
                minorticks=true,
                grid=false,
                # xlims=pma_extrema_x,
                # ylims=pma_extrema_y,
                # clims=clims_pma
            )

            return pma_plot


        end
    end
end
