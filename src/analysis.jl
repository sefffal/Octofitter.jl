"""
    projectpositions(chains.planets[1], mjd("2020-02-02"))

Given the posterior for a particular planet in the model and a modified julian date(s),
return `ra` and `dec` offsets in mas for each sampling in the posterior.
"""
function projectpositions(chain, planet_key, times)

    ras = zeros(size(chain,1) * length(times))
    decs = zeros(size(chain,1) * length(times))

    els = construct_elements(chain, planet_key, 1:size(chain,1))
    
    for (j,el) in enumerate(els)

        for (k, t) in enumerate(times)
            i = j * length(times) + k - 1
            ra, dec, _ = kep2cart(el, t)
            ras[i] = ra
            decs[i] = dec
        end
    end
    return ras, decs
end
export projectpositions

"""
    sampleorbits(chains, planet_key, 100)

Given the posterior for a particular planet in the model, and the number of orbits to sample,
return a random subset of length N.
"""
function sampleorbits(chains, planetnum, N)
    return map(rand(eachindex(chains.planets[planetnum].a), N)) do i
        return KeplerianElements(;
            a = chains.planets[planetnum].a[i],
            i = chains.planets[planetnum].i[i],
            e = chains.planets[planetnum].e[i],
            ω = chains.planets[planetnum].ω[i],
            Ω = chains.planets[planetnum].Ω[i],
            μ = chains.μ[i],
            plx = chains.plx[i],
            τ = chains.planets[planetnum].τ[i],
        )
    end
end
export sampleorbits


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
            ii = rand(1:size(chain,1), N)
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
            Plots.scatter!([0],[0], marker=(:star,:black, 5), label="")
        end
        function plotposterior!(
            chain,
            planet_key,
            property,
            N=1500;
            alpha=0.02,
            cmap=:turbo,
            rev=true,
            colorbartitle="",#"semi-major axis (au)",
            clims=nothing,
            lw=0.3,
            kwargs...
        )
            ii = rand(1:size(chain,1), N)
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
            for i in sortperm(colours, rev=false)
                c = colours[i]
                Plots.plot!(elements[i]; label="",color=cmap[(c-minc)/dc], lw, alpha)
            end
            Plots.scatter!([0],[0], marker=(:star,:black, 5), label="")
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
            system,
            N=1500,
            alpha=0.02;
            plotpma=true,
            cmap=:plasma,
            pma_scatter=nothing,
            clims=extrema(Iterators.flatten(
                # extrema(planet.a)
                # extrema([sample.planets[key].a for sample in chain])
                # for key in keys(chain[1].planets)
                extrema(chain["$pk[a]"] for pk in keys(system.planets))
            )),
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
                    @warn "Plotmodel does not yet support images with multiple platescales"
                    img = DirectImages.DirectImage(first(system.images.image))
                else
                    img = DirectImages.DirectImage(
                        DirectImages.stack(maximum,system.images.image)
                    )
                end
                img.PLATESCALE = system.images.platescale[1]
                # To get the colour scales to work out
                img ./= only(quantile(filter(isfinite, arraydata(img)), [0.9995]))
                img .*= maximum(clims) - minimum(clims)
                img .+= minimum(clims)
                imshow!(img; color=:Greys, skyconvention=true, lims)
            end
            
            for planet_key in eachindex(system.planets)
                plotposterior!(
                    chain, planet_key, :a, N; lw=1, alpha=alpha, colorbartitle="semi-major axis (au)",
                    cmap=cmap, rev=false,
                    # cmap=:turbo,
                    clims=clims,
                    kwargs...
                )

                # Planet astrometry?
                astrom = DirectDetections.astrometry(system.planets[planet_key])
                if !isnothing(astrom)
                    Plots.scatter!(astrom,marker=(:black,:circle,3),label="")
                end
            end
            # We will override this if we have more information available further down
            final_plot = p_orbits

            # astrometric acceleration?
            if !isnothing(system.propermotionanom) && plotpma


                titles=["GAIA EDR3", "Hipparcos",]
                system_pma = system.propermotionanom
                pma_plots = map(sortperm(system_pma.ra_epoch)) do i
                    vx = zeros(length(chain))
                    vy = zeros(length(chain))
                    for j in keys(system.planets)
                        ii = rand(eachindex(chain), N)
                        elements = construct_elements(chain, j, 1:size(chain,1))
                        # mass = [sample.planets[j].mass for sample in chain]
                        mass = chain["$j[mass]"]
                        vx .+= getindex.(propmotionanom.(elements, system_pma.ra_epoch[i], getproperty.(elements, :μ), mass),1)
                        vy .+= getindex.(propmotionanom.(elements, system_pma.dec_epoch[i], getproperty.(elements, :μ), mass),2)
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
                        prop = getproperty(chains.planets[1], pma_scatter)[ii]
                        if pma_scatter == :mass
                            prop ./=mjup2msol
                        end
                        ii_sub = rand(ii, 4000)
                        Plots.scatter!(vx[ii_sub], vy[ii_sub], marker_z=prop[ii_sub], alpha=1, color=:plasma, colorbar=true, label="", markerstrokewidth=0, ms=1)
                    else
                        h = fit(Histogram, (vx, vy))#, (-1:0.05:1, -1:0.05:1))
                        Plots.plot!(h, color=Plots.cgrad([Plots.RGBA(0,0,0,0), Plots.RGBA(0,0,0,1)]), colorbar=false)
                    end
                    Plots.scatter!(
                        [system_pma.pm_ra[i]],
                        [system_pma.pm_dec[i]],
                        xerror=[system_pma.σ_pm_ra[i]],
                        yerror=[system_pma.σ_pm_dec[i]],
                        label="",
                        color=1,
                        markerstrokecolor=1
                    )
                    # Plots.scatter!([0], [0], marker=(5, :circle, :red),label="")
                    Plots.hline!([0], color=:black, label="")
                    Plots.vline!([0], color=:black, label="")
                    Plots.title!(titles[i])
                    # Plots.xlims!(-1,1)
                    # Plots.ylims!(-1,1)
                    Plots.xlabel!("Δμ ra - mas/yr")
                    Plots.ylabel!("Δμ dec - mas/yr")
                end     
                # pma_plot = plot(pma_plots..., size=(700,300),margin=5Plots.mm, guidefontsize=9, titlefontsize=9, top_margin=0Plots.mm)

                # Crete new final plot
                l = eval(:(Plots.@layout [
                    A{0.65h}
                    B C
                ]))
                final_plot = Plots.plot(p_orbits, pma_plots..., layout=l, size=(550,750),
                    margin=5Plots.mm, guidefontsize=9, titlefontsize=9, top_margin=0Plots.mm)
            end

            return final_plot
        end
    end
end