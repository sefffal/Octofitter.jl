# This file contains funtions for analysing results from chains.

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
            els = DirectDetections.construct_elements(chains, planet_key, is .* j)
            for (i,el) in zip(is,els)
                for (k, t) in enumerate(times)
                    ra, dec, _ = kep2cart(el, t)
                    ras[i,j,k] = ra
                    decs[i,j,k] = dec
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
#         return KeplerianElements(;
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

@recipe function f(astrom::Astrometry) where T
    xerror := astrom.σ_ra
    yerror := astrom.σ_dec
    xflip --> true
    xguide --> "Δ right ascension (mas)"
    yguide --> "Δ declination (mas)"

    return astrom.ra, astrom.dec
end

function plotposterior end
export plotposterior
function plotposterior! end
export plotposterior!
function plotmodel end
export plotmodel
function plotmodel! end
export plotmodel!

# Optionally depend on Plots. If the user imports it, this code will be run to set up
# our `imshow` function.
using Requires
function init_plots()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin

        function plotposterior(args...;kwargs...)
            Plots.plot()
            plotposterior!(args...;kwargs...)
        end
        function plotposterior!(
            chain,
            planet_key,
            keys::Nothing,
            N=1500;
            alpha=0.02,
            lw=0.3,
            color=1,
            kwargs...
        )
            ii = rand(1:size(chain,1)*size(chain,3), N)
            elements = construct_elements(chain, planet_key, ii)

            Plots.plot!(;
                size=(700,700),
                dpi=200,
                fmt=:png,
                framestyle=:box,
                grid=:none,
                minorticks=true,
                aspectratio=1,
                fontfamily="Arial",
                margin=8Plots.mm,
                kwargs...
            )
            Plots.xlabel!("Δ right ascension (mas)")
            Plots.ylabel!("Δ declination (mas)")

        
            for i in eachindex(elements)
                Plots.plot!(elements[i], label="",color=color, lw=lw, alpha=alpha,)
            end
            Plots.scatter!([0],[0], marker=(:star,:black, 5), markerstrokewidth=1, markerstrokecolor=:white, label="")
        end
        function plotposterior!(
            chain,
            planet_key,
            property,
            N=1500;
            alpha=0.02,
            cmap=:turbo,
            rev=true,
            colorbartitle= property isa Symbol ? "$planet_key[$property]" : "",#"semi-major axis (au)",
            clims=nothing,
            lw=0.3,
            kwargs...
        )
            ii = rand(1:size(chain,1)*size(chain,3), N)
            elements = construct_elements(chain, planet_key, ii)

            Plots.plot!(;
                size=(550,550),
                dpi=200,
                fmt=:png,
                framestyle=:box,
                grid=:none,
                minorticks=true,
                aspectratio=1,
                fontfamily="Arial",
                margin=8Plots.mm,
                kwargs...
            )
            Plots.xlabel!("Δ right ascension (mas)")
            Plots.ylabel!("Δ declination (mas)")

            if property isa Symbol
                k = string(planet_key)*"[$property]"
                colours = chain[k][ii]
            else # Assume they passed in the values directly
                colours = property[ii]
            end

            if isnothing(clims)
                clims = extrema(colours)
            end
            
            minc = first(clims)
            dc = last(clims) - minc
            clims = (minc, dc+minc)

            cmap= Plots.cgrad(cmap; rev)
            for i in eachindex(colours)#sortperm(colours, rev=false)
                c = colours[i]
                Plots.plot!(elements[i]; label="",color=cmap[(c-minc)/dc], lw, alpha)
            end
            Plots.scatter!([0],[0], marker=(:star,:black, 5), markerstrokewidth=1, markerstrokecolor=:white,  label="")
            Plots.scatter!([0],[0]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)

        end
        
        function plotmodel(args...; kwargs...)
            Plots.plot()
            plotmodel!(args...;kwargs...)
        end
        function plotmodel!(images::Images,i; kwargs...)
            img = DirectImages.DirectImage(images.image[i])
            img.PLATESCALE = images.platescale[i]
            imshow!(img;skyconvention=true,kwargs...)
            # plot!(sampled,color=:white,alpha=5/N,label="")
            Plots.xlims!(img.PLATESCALE.*extrema(axes(img,1)))
            Plots.ylims!(img.PLATESCALE.*extrema(axes(img,2)))

        end

        function plotmodel!(
            chain,
            N=1500,
            system=chain.info.model;
            alpha=(N <= 15 ? 1 : 30/N),
            color= isnothing(system.propermotionanom) ? :a : :mass,
            plotpma=!isnothing(system.propermotionanom),
            # TODO: ideally this is based on if there is a mass variable
            plotmass=!isnothing(system.propermotionanom),
            cmap=:plasma,
            imagecmap=:Greys,
            pma_scatter=nothing,
            clims=quantile(Iterators.flatten(
                # extrema(planet.a)
                # extrema([sample.planets[key].a for sample in chain])
                # for key in keys(chain[1].planets)
                extrema(chain["$pk[$color]"] for pk in keys(system.planets))
            ), (0.01,0.99)),
            lims=nothing,
            kwargs...
        )

            # Plot orbits
            p_orbits = Plots.plot(
                xlims=:symmetric,
                ylims=:symmetric,
            )

            # plot images?
            if !isnothing(system.images)
                if length(unique(system.images.platescale)) != 1
                    @warn "Plotmodel does not yet support images with multiple platescales. Taking first image only."
                    img = DirectImages.DirectImage(first(system.images.image))
                else
                    img = DirectImages.DirectImage(
                        DirectImages.stack(maximum,system.images.image)
                    )
                end
                img.PLATESCALE = system.images.platescale[1]
                # To get the colour scales to work out
                img ./= only(quantile(filter(isfinite, arraydata(img)), [0.98]))#[0.9995]))
                img .*= maximum(clims) - minimum(clims)
                img .+= minimum(clims)
                imshow!(img; color=imagecmap, skyconvention=true, lims)
            end
            
            for planet_key in eachindex(system.planets)
                plotposterior!(
                    chain, planet_key, color, N; lw=1, alpha=alpha,
                    cmap=cmap, rev=false,
                    # cmap=:turbo,
                    clims=clims,
                    kwargs...
                )

                # Planet astrometry?
                astrom = DirectDetections.astrometry(system.planets[planet_key])
                if !isnothing(astrom)
                    # Put black error bars over thicker white error bars to make them easier to see
                    # when plotted over e.g. images, orbits.
                    Plots.scatter!(astrom,marker=(:white,:circle,0),label="",linewidth=0,color=:white,markerstrokewidth=3, markerstrokecolor=:white)
                    Plots.scatter!(astrom,marker=(:black,:circle,0),label="",linewidth=0,color=:black,markerstrokewidth=1, markerstrokecolor=:black)
                end
            end
            # We will override this if we have more information available further down
            final_plot = p_orbits

            # astrometric acceleration?
            if plotpma


                lpma =length(system.propermotionanom.ra_epoch)
                pma_ii = sortperm(system.propermotionanom.ra_epoch)
                if lpma == 2 && 
                    isapprox(minimum(system.propermotionanom.ra_epoch), mjd("1991"), atol=365) &&
                    isapprox(maximum(system.propermotionanom.ra_epoch), mjd("2016"), atol=365)

                    titles=["Hipparcos", "GAIA EDR3"]
                else
                    titles = string.(system.propermotionanom.ra_epoch)
                end
                system_pma = system.propermotionanom
                
                pma_extrema_x = [-0.0,0.0]
                pma_extrema_y = [-0.0,0.0]
                pma_plots = map(sortperm(system_pma.ra_epoch)) do i
                    vx = zeros(prod(size(chain,[1,3])))
                    vy = zeros(prod(size(chain,[1,3])))
                    for j in keys(system.planets)
                        elements = construct_elements(chain, j, 1:prod(size(chain,[1,3])))
                        # mass = [sample.planets[j].mass for sample in chain]
                        mass = reshape(chain["$j[mass]"],:)
                        vx .+= getindex.(propmotionanom.(elements, system_pma.ra_epoch[i], mass.*mjup2msol),1)
                        vy .+= getindex.(propmotionanom.(elements, system_pma.dec_epoch[i], mass.*mjup2msol),2)
                    end
                    Plots.plot(
                        framestyle=:box,
                        minorticks=true,
                        aspectratio=1,
                        grid=false,
                        xlims=:symmetric,
                        ylims=:symmetric,
                    )
                    if !isnothing(pma_scatter)
                        if length(system.planets) > 1
                            @warn "Multiplanet PMA scatter plots not yet implemented"
                        end
                        ii = rand(1:size(chain,1), N)
                        planet_key = first(system.planets).name
                        prop = chain["$planet_key[$pma_scatter]"][ii]
                        # if pma_scatter == :mass
                        #     prop ./=mjup2msol
                        # end
                        if color == pma_scatter
                            clims_pma = clims
                        else
                            clims_pma = quantile(prop, (0.01, 0.99))
                        end
                        Plots.scatter!(vx[ii], vy[ii], marker_z=prop, alpha=1, color=:plasma, legend=false, colorbar=false, label="", markerstrokewidth=0, ms=1, clims=clims_pma)
                    else
                        h = fit(Histogram, (vx, vy))#, (-1:0.05:1, -1:0.05:1))
                        Plots.plot!(h, color=Plots.cgrad([Plots.RGBA(0,0,0,0), Plots.RGBA(0,0,0,1)]), colorbar=false)
                    end
                    pma_plot = Plots.scatter!(
                        [system_pma.pm_ra[i]],
                        [system_pma.pm_dec[i]],
                        xerror=[system_pma.σ_pm_ra[i]],
                        yerror=[system_pma.σ_pm_dec[i]],
                        label="",
                        color=:red,
                        markercolor=:red,
                        markerstrokecolor=:red,
                    )
                    # Plots.scatter!([0], [0], marker=(5, :circle, :red),label="")
                    Plots.hline!(pma_plot, [0], color=:black, label="")
                    Plots.vline!(pma_plot, [0], color=:black, label="")
                    Plots.title!(pma_plot, titles[i])
                    # Plots.xlims!(-1,1)
                    # Plots.ylims!(-1,1)
                    Plots.xlabel!(pma_plot, "Δμ ra - mas/yr")
                    Plots.ylabel!(pma_plot, "Δμ dec - mas/yr")
                    pma_extrema_x[1] = min(pma_extrema_x[1], minimum(vx), system_pma.pm_ra[i] - system_pma.σ_pm_ra[i], system_pma.pm_ra[i] + system_pma.σ_pm_ra[i])
                    pma_extrema_x[2] = max(pma_extrema_x[2], maximum(vx), system_pma.pm_ra[i] - system_pma.σ_pm_ra[i], system_pma.pm_ra[i] + system_pma.σ_pm_ra[i])
                    pma_extrema_y[1] = min(pma_extrema_y[1], minimum(vy), system_pma.pm_dec[i] - system_pma.σ_pm_dec[i], system_pma.pm_dec[i] + system_pma.σ_pm_dec[i])
                    pma_extrema_y[2] = max(pma_extrema_y[2], maximum(vy), system_pma.pm_dec[i] - system_pma.σ_pm_dec[i], system_pma.pm_dec[i] + system_pma.σ_pm_dec[i])
                    return pma_plot
                end
                for plot in pma_plots
                    # Plots.xlims!(plot, (pma_extrema_x...,))
                    # Plots.ylims!(plot, (pma_extrema_y...,))
                    lim = maximum(abs, ([pma_extrema_x; pma_extrema_y]))
                    Plots.xlims!(plot, -lim, lim)
                    Plots.ylims!(plot, -lim, lim)
                end
                # pma_plot = plot(pma_plots..., size=(700,300),margin=5Plots.mm, guidefontsize=9, titlefontsize=9, top_margin=0Plots.mm)
                # pma_plot = Plots.plot(pma_plots, link=:both)

                # Crete new final plot
                l = eval(:(Plots.@layout [
                    A{0.65h}
                    Plots.grid(1, $lpma)
                ]))
                final_plot = Plots.plot(p_orbits, pma_plots..., layout=l, size=(550,750),
                    margin=5Plots.mm, guidefontsize=9, titlefontsize=9, top_margin=0Plots.mm)
            end

            if plotmass
                p = Plots.plot(;xlabel="", linewidth=3, label="")
                for p in keys(system.planets)
                    m = vec(chain["$p[mass]"])
                    h = fit(Histogram, m, nbins=round(Int,sqrt(length(m))))#, nbins=100)
                    Plots.plot!(h.edges[1][begin:end-1],seriestype=:step, h.weights, label="$p[mass]")
                end

                layout = eval(:(Plots.@layout [
                    A{0.8h}
                    B
                ]))
                final_plot = Plots.plot(final_plot, p, layout=layout, size=(550,850))
            end

            return final_plot
        end
    end
end