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
            els = DirectDetections.construct_elements(chains, planet_key, is .* j)
            for (i,el) in zip(is,els)
                for (k, t) in enumerate(times)
                    o = kep2cart(el, t)
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
    xerror := astrom.table.σ_ra
    yerror := astrom.table.σ_dec
    xflip --> true
    xguide --> "Δ right ascension (mas)"
    yguide --> "Δ declination (mas)"

    return astrom.table.ra, astrom.table.dec
end

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

            # Lines to periastron
            xs = hcat(zeros(1500), raoff.(elements, periastron.(elements)), fill(NaN, 1500))'[:]
            ys = hcat(zeros(1500), decoff.(elements, periastron.(elements)), fill(NaN, 1500))'[:]
            Plots.plot!(xs, ys; lw, alpha, color="#777", label="")
            cmap= Plots.cgrad(cmap; rev)
            for i in eachindex(colours)#sortperm(colours, rev=false)
                c = colours[i]
                Plots.plot!(elements[i]; label="",color=cmap[(c-minc)/dc], lw, alpha)
            end
            # Star at centre
            Plots.scatter!([0],[0], marker=(:star, :black, 5), markerstrokewidth=1, markerstrokecolor=:white,  label="")
            # Colorbar
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
            color= isnothing(propermotionanom(system)) ? :a : :mass,
            plotpma=false,#!isnothing(propermotionanom(system)),
            # TODO: ideally this is based on if there is a mass variable
            plotmass=!isnothing(propermotionanom(system)),
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
            if !isnothing(images(system))
                if length(unique(images(system).table.platescale)) != 1
                    @warn "Plotmodel does not yet support images with multiple platescales. Taking first image only."
                    img = DirectImages.DirectImage(first(images(system).image))
                else
                    img = DirectImages.DirectImage(
                        DirectImages.stack(maximum,images(system).table.image)
                    )
                end
                img.PLATESCALE = images(system).table.platescale[1]
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


                lpma =length(propermotionanom(system).table.ra_epoch)
                # pma_ii = sortperm(propermotionanom(system).table.ra_epoch)
                if lpma == 2 && 
                    isapprox(minimum(propermotionanom(system).table.ra_epoch), mjd("1991"), atol=365) &&
                    isapprox(maximum(propermotionanom(system).table.ra_epoch), mjd("2016"), atol=365)

                    titles=["Hipparcos", "GAIA"]
                elseif lpma == 3 && 
                    isapprox(minimum(propermotionanom(system).table.ra_epoch), mjd("1991"), atol=365) &&
                    isapprox(minimum(propermotionanom(system).table.ra_epoch), (mjd("1991")+mjd("2016"))/2, atol=365) &&
                    isapprox(maximum(propermotionanom(system).table.ra_epoch), mjd("2016"), atol=365)

                    titles=["Hipparcos", "GAIA-Hipparcos", "GAIA"]
                else
                    titles = string.(round.(Int, propermotionanom(system).table.ra_epoch))
                end
                system_pma = propermotionanom(system)
                
                pma_extrema_x = [-0.0,0.0]
                pma_extrema_y = [-0.0,0.0]
                pma_plots = map(sortperm(system_pma.table.ra_epoch)) do i
                    vx = zeros(prod(size(chain,[1,3])))
                    vy = zeros(prod(size(chain,[1,3])))
                    for j in keys(system.planets)
                        elements = construct_elements(chain, j, 1:prod(size(chain,[1,3])))
                        # mass = [sample.planets[j].mass for sample in chain]
                        mass = reshape(chain["$j[mass]"],:)
                        vx .+= pmra.(elements, system_pma.table.ra_epoch[i], mass.*mjup2msol) .+ vec(chain["pmra"])
                        vy .+= pmdec.(elements, system_pma.table.dec_epoch[i], mass.*mjup2msol) .+ vec(chain["pmdec"])
                    end
                    Plots.plot(
                        framestyle=:box,
                        minorticks=true,
                        aspectratio=1,
                        grid=false,
                        # xlims=:symmetric,
                        # ylims=:symmetric,
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
                        [system_pma.table.pm_ra[i]],
                        [system_pma.table.pm_dec[i]],
                        xerror=[system_pma.table.σ_pm_ra[i]],
                        yerror=[system_pma.table.σ_pm_dec[i]],
                        label="",
                        color=:red,
                        markercolor=:red,
                        markerstrokecolor=:red,
                    )
                    # Plots.scatter!([0], [0], marker=(5, :circle, :red),label="")
                    # Plots.hline!(pma_plot, [0], color=:black, label="")
                    # Plots.vline!(pma_plot, [0], color=:black, label="")
                    Plots.title!(pma_plot, titles[i])
                    # Plots.xlims!(-1,1)
                    # Plots.ylims!(-1,1)
                    Plots.xlabel!(pma_plot, "Δμ ra - mas/yr")
                    Plots.ylabel!(pma_plot, "Δμ dec - mas/yr")
                    pma_extrema_x[1] = min(minimum(vx), system_pma.table.pm_ra[i] - system_pma.table.σ_pm_ra[i], system_pma.table.pm_ra[i] + system_pma.table.σ_pm_ra[i])
                    pma_extrema_x[2] = max(maximum(vx), system_pma.table.pm_ra[i] - system_pma.table.σ_pm_ra[i], system_pma.table.pm_ra[i] + system_pma.table.σ_pm_ra[i])
                    pma_extrema_y[1] = min(minimum(vy), system_pma.table.pm_dec[i] - system_pma.table.σ_pm_dec[i], system_pma.table.pm_dec[i] + system_pma.table.σ_pm_dec[i])
                    pma_extrema_y[2] = max(maximum(vy), system_pma.table.pm_dec[i] - system_pma.table.σ_pm_dec[i], system_pma.table.pm_dec[i] + system_pma.table.σ_pm_dec[i])
                    return pma_plot
                end
                # for plot in pma_plots
                #     Plots.xlims!(plot, pma_extrema_x...)
                #     Plots.ylims!(plot, pma_extrema_y...)
                # end
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
        function timeplot(
            chain,
            planet_key,
            color,
            prop, 
            N = 500;
            ii = rand(1:size(chain,1)*size(chain,3), N),
            kwargs...
        )
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
            p1 = Plots.plot(;
                ylabel,
                legend=:none
            )
            xerr = nothing
            elements = DirectDetections.construct_elements(chain, planet_key, ii)
            all_epochs = reduce(vcat, Any[t.table.epoch for t in planet.observations if haskey(t.table, :epoch)], init=[mjd()])
            t = range((extrema(all_epochs) .+ [-365, 365])..., length=100)
            y = nothing
            if prop == :ra
                y = astrometry(planet).table.ra
                yerr = astrometry(planet).table.σ_ra
                fit = raoff.(elements, t')'
                x = astrometry(planet).table.epoch
            elseif prop == :dec
                y = astrometry(planet).table.dec
                yerr = astrometry(planet).table.σ_dec
                fit = decoff.(elements, t')'
                x = astrometry(planet).table.epoch
            elseif prop == :sep
                xx = astrometry(planet).table.ra
                yy = astrometry(planet).table.dec
                xxerr = astrometry(planet).table.σ_ra
                yyerr = astrometry(planet).table.σ_dec
                y_unc = sqrt.((xx .± xxerr).^2 .+ (yy .± yyerr).^2)
                y = Measurements.value.(y_unc)
                yerr = Measurements.uncertainty.(y_unc)
                fit = projectedseparation.(elements, t')'
                x = astrometry(planet).table.epoch
            elseif prop == :pa
                xx = astrometry(planet).table.ra
                yy = astrometry(planet).table.dec
                xxerr = astrometry(planet).table.σ_ra
                yyerr = astrometry(planet).table.σ_dec
                y_unc = atand.((xx .± xxerr), (yy .± yyerr))
                y = Measurements.value.(y_unc)
                yerr = Measurements.uncertainty.(y_unc)
                fit = rad2deg.(posangle.(elements, t')')
                x = astrometry(planet).table.epoch
            elseif prop == :pmra
                hgca = propermotionanom(chain.info.model).table[1]
                t = range(years2mjd(hgca.epoch_ra_hip)-365*5, years2mjd(hgca.epoch_ra_gaia)+365*5, length=100)
                y = [hgca.pmra_hip, hgca.pmra_hg, hgca.pmra_gaia]
                yerr = [hgca.pmra_hip_error, hgca.pmra_hg_error, hgca.pmra_gaia_error]
                fit = chain["pmra"][ii]' .+ pmra.(elements, t', collect(chain["$planet_key[mass]"][ii]).*DirectDetections.mjup2msol)'
                x = years2mjd.([hgca.epoch_ra_hip, (hgca.epoch_ra_hip + hgca.epoch_ra_gaia)/2, hgca.epoch_ra_gaia])
                xerr = [4*365/2, 25*365/2, 3*365/2]
            elseif prop == :pmdec
                hgca = propermotionanom(chain.info.model).table[1]
                t = range(years2mjd(hgca.epoch_dec_hip)-365*5, years2mjd(hgca.epoch_dec_gaia)+365*5, length=100)
                y = [hgca.pmdec_hip, hgca.pmdec_hg, hgca.pmdec_gaia]
                yerr = [hgca.pmdec_hip_error, hgca.pmdec_hg_error, hgca.pmdec_gaia_error]
                fit = chain["pmdec"][ii]' .+ pmdec.(elements, t', collect(chain["$planet_key[mass]"][ii]).*DirectDetections.mjup2msol)'
                x = years2mjd.([hgca.epoch_dec_hip, (hgca.epoch_dec_hip + hgca.epoch_dec_gaia)/2, hgca.epoch_dec_gaia])
                xerr = [4*365/2, 25*365/2, 3*365/2]
            elseif prop == :rv
                # model lines
                fit = radvel.(elements, t', collect(chain["$planet_key[mass]"][ii]))'
                for obs in planet.observations
                    if obs isa RadialVelocity
                        # plot rv data over model
                        y = obs.table.rv
                        yerr = obs.table.σ_rv
                        x = obs.table.epoch
                    end
                end
            end
            Plots.plot!(
                mjd2date.(t), fit,
                line_z=repeat(
                    collect(chain["$planet_key[$color]"][ii]),
                    1, length(t)
                )',
                alpha=0.05;
                kwargs...
            )
            if !isnothing(y)
                Plots.scatter!(
                    p1,
                    mjd2date.(x),
                    y; yerr, xerr,
                    # markersize=1.5,
                    # color=:black,
                    color=1,
                    # markerstrokewidth=3,
                )
            end
            p1
        end

        function timeplotgrid(
            chains;
            color="b[e]",
            clims = quantile(vec(chains[color]),(0.01, 0.99)),
            N=1500,
            ii = rand(1:size(chains,1)*size(chains,3), N)
        )
            kwargs = (; 
                ii,
                clims,
                cmap=:plasma,
                alpha=0.1
            )
            ppost = plotposterior(chains, :b, :e; rev=false, colorbar=nothing, kwargs...)
            for (i,planet_key) in enumerate(keys(chains.info.model.planets))
                astrom = astrometry(chains.info.model.planets[planet_key])
                if !isnothing(astrom)
                    Plots.scatter!(ppost, astrom, label="", color=i)
                end
            end
            Plots.plot(
                timeplot(chains, :b, :e, :sep; kwargs...),
                timeplot(chains, :b, :e, :pa; kwargs...),
                timeplot(chains, :b, :e, :pmra; kwargs...),
                timeplot(chains, :b, :e, :pmdec; kwargs...),
                ppost,
                timeplot(chains, :b, :e, :rv; kwargs...),
                layout = (3,2),
                framestyle=:box,
                grid=false,
                size=(1000,1200),
            )
        end
    end
end