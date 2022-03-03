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
            keys::Nothing;
            N=1500,
            ii = rand(1:size(chain,1)*size(chain,3), N),
            alpha=0.02,
            lw=0.3,
            color=1,
            kwargs...
        )
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
            property;
            N=1500,
            ii = rand(1:size(chain,1)*size(chain,3), N),
            alpha=0.02,
            cmap=:turbo,
            rev=true,
            colorbartitle= property isa Symbol ? string(property) : "",#"semi-major axis (au)",
            clims=nothing,
            lw=0.3,
            kwargs...
        )
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

            if typeof(property) <: Union{Symbol,AbstractString}
                # k = string(planet_key)*"[$property]"
                k = property
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
            ii = rand(1:size(chain,1)*size(chain,3), N),
            color=length(system.planets) == 0 ? nothing : string(first(keys(system.planets)))*"[a]",
            colorbartitle=color,
            plotpma=!isnothing(propermotionanom(system)),
            # TODO: ideally this is based on if there is a mass variable
            plotmass=!isnothing(propermotionanom(system)),
            cmap=:plasma,
            imagecmap=:Greys,
            clims=isnothing(color) ? nothing : quantile(Iterators.flatten(
                # extrema(planet.a)
                # extrema([sample.planets[key].a for sample in chain])
                # for key in keys(chain[1].planets)
                extrema(chain[color])
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
                    chain, planet_key, color; ii, lw=1, alpha=alpha,
                    cmap=cmap, rev=false,
                    # cmap=:turbo,
                    clims=clims,
                    colorbartitle=colorbartitle,
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
            if plotpma && !isnothing(propermotionanom(system))
                pma_plots = pmaplot(chain; ii, color)
                
                # # Crete new final plot
                l = eval(:(Plots.@layout [
                    A{0.65h}
                    B
                ]))
                final_plot = Plots.plot(p_orbits, pma_plots, layout=l, size=(550,750),
                    margin=5Plots.mm, guidefontsize=9, titlefontsize=9, top_margin=0Plots.mm)
            end


            if plotmass
                p = Plots.plot(;xlabel="", linewidth=3, label="")
                for p in keys(system.planets)
                    m = vec(chain["$p[mass]"])
                    h = fit(Histogram, m, nbins=round(Int,sqrt(length(m))))#, nbins=100)
                    Plots.plot!(h.edges[1][begin:end-1],seriestype=:step, h.weights, label="$p[mass]",color=1,lw=2)
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
            elements = DirectDetections.construct_elements(chain, planet_key, ii)
            all_epochs = mapreduce(vcat, keys(chain.info.model.planets)) do planet_key
                planet = chain.info.model.planets[planet_key]
                reduce(vcat, Any[hasproperty(t.table, :epoch) ? t.table.epoch : t.table.ra_epoch for t in planet.observations if hasproperty(t.table, :epoch) || hasproperty(t.table, :ra_epoch)], init=[])
            end
            if isempty(all_epochs) 
                all_epochs = mjd() .+ [-365*2, +365*2]
            end
            t = range((extrema(all_epochs) .+ [-365, 365])..., length=100)
            y = nothing
            ribbon = nothing
            if prop == :ra
                if !isnothing(astrometry(planet))
                    y = astrometry(planet).table.ra
                    yerr = astrometry(planet).table.σ_ra
                    x = astrometry(planet).table.epoch
                end
                fit = raoff.(elements, t')
            elseif prop == :dec
                if !isnothing(astrometry(planet))
                    y = astrometry(planet).table.dec
                    yerr = astrometry(planet).table.σ_dec
                    x = astrometry(planet).table.epoch
                end
                fit = decoff.(elements, t')
            elseif prop == :sep
                if !isnothing(astrometry(planet))
                    xx = astrometry(planet).table.ra
                    yy = astrometry(planet).table.dec
                    xxerr = astrometry(planet).table.σ_ra
                    yyerr = astrometry(planet).table.σ_dec
                    y_unc = sqrt.((xx .± xxerr).^2 .+ (yy .± yyerr).^2)
                    y = Measurements.value.(y_unc)
                    yerr = Measurements.uncertainty.(y_unc)
                    x = astrometry(planet).table.epoch
                end
                fit = projectedseparation.(elements, t')
            elseif prop == :pa
                if !isnothing(astrometry(planet))
                    xx = astrometry(planet).table.ra
                    yy = astrometry(planet).table.dec
                    xxerr = astrometry(planet).table.σ_ra
                    yyerr = astrometry(planet).table.σ_dec
                    y_unc = atand.((xx .± xxerr), (yy .± yyerr))
                    y = Measurements.value.(y_unc)
                    yerr = Measurements.uncertainty.(y_unc)
                    x = astrometry(planet).table.epoch
                end
                fit = rad2deg.(posangle.(elements, t'))
            elseif prop == :pmra
                hgca = propermotionanom(chain.info.model).table[1]
                t = range(years2mjd(hgca.epoch_ra_hip)-365*5, years2mjd(hgca.epoch_ra_gaia)+365*5, length=100)
                y = [hgca.pmra_hip, hgca.pmra_hg, hgca.pmra_gaia]
                yerr = [hgca.pmra_hip_error, hgca.pmra_hg_error, hgca.pmra_gaia_error]
                fit = chain["pmra"][ii] .+ pmra.(elements, t', collect(chain["$planet_key[mass]"][ii]).*DirectDetections.mjup2msol)
                x = years2mjd.([hgca.epoch_ra_hip, (hgca.epoch_ra_hip + hgca.epoch_ra_gaia)/2, hgca.epoch_ra_gaia])
                xerr = [4*365/2, 25*365/2, 3*365/2]
            elseif prop == :pmdec
                hgca = propermotionanom(chain.info.model).table[1]
                t = range(years2mjd(hgca.epoch_dec_hip)-365*5, years2mjd(hgca.epoch_dec_gaia)+365*5, length=100)
                y = [hgca.pmdec_hip, hgca.pmdec_hg, hgca.pmdec_gaia]
                yerr = [hgca.pmdec_hip_error, hgca.pmdec_hg_error, hgca.pmdec_gaia_error]
                fit = chain["pmdec"][ii] .+ pmdec.(elements, t', collect(chain["$planet_key[mass]"][ii]).*DirectDetections.mjup2msol)
                x = years2mjd.([hgca.epoch_dec_hip, (hgca.epoch_dec_hip + hgca.epoch_dec_gaia)/2, hgca.epoch_dec_gaia])
                xerr = [4*365/2, 25*365/2, 3*365/2]
            elseif prop == :rv
                # TODO: how do we indicate jitter?
                if :rv in keys(chain)
                    sys_rv = chain["rv"][ii]
                else
                    sys_rv = 0
                end
                fit = sys_rv .+ radvel.(elements, t', collect(chain["$planet_key[mass]"][ii]).*mjup2msol)
                if :jitter in keys(chain)
                    ribbon = chain["jitter"][ii]
                end
                for obs in chain.info.model.observations
                    if obs isa RadialVelocity
                        # plot rv data over model
                        y = obs.table.rv
                        yerr = obs.table.σ_rv
                        x = obs.table.epoch
                    end
                end
            end
            # if haskey(kwargs, :cmap)
                # kwargs = (;kwargs..., color=kwargs[:cmap])
            # end

            # Plots.plot!(
            #     mjd2date.(t), fit;
            #     line_z=repeat(
            #         collect(chain[color][ii]),
            #         1, length(t)
            #     )',
            #     alpha,
            #     fillalpha=alpha,
            #     ribbon=ribbon,
            #     kwargs...,
            # )
            for i in axes(fit, 1)
                c = (chain[color][ii[i]] - first(clims))/(last(clims)-first(clims))
                Plots.plot!(
                    mjd2date.(t), view(fit, i, :);
                    # line_z=chain[color][ii],
                    color=Plots.cgrad(cmap)[c],
                    alpha,
                    fillalpha=alpha,
                    ribbon=isnothing(ribbon) ? nothing : ribbon[i],
                    kwargs...,
                )
            end
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

            ppost = Plots.plot()
            for (i,planet_key) in enumerate(keys(chains.info.model.planets))
                plotposterior!(chains, planet_key, color; rev=false, colorbar=nothing, kwargs...)
                astrom = astrometry(chains.info.model.planets[planet_key])
                if !isnothing(astrom)
                    Plots.scatter!(ppost, astrom, label="", color=i)
                end
            end
            psep = Plots.plot()
            for planet_key in keys(chains.info.model.planets)
                timeplot!(chains, planet_key, color, :sep; kwargs...)
            end
            ppa = Plots.plot()
            for planet_key in keys(chains.info.model.planets)
                timeplot!(chains, planet_key, color, :pa; kwargs...)
            end
            ppmra = Plots.plot()
            for planet_key in keys(chains.info.model.planets)
                timeplot!(chains, planet_key, color, :pmra; kwargs...)
            end
            ppmdec = Plots.plot()
            for planet_key in keys(chains.info.model.planets)
                timeplot!(chains, planet_key, color, :pmdec; kwargs...)
            end
            prv = Plots.plot()
            for planet_key in keys(chains.info.model.planets)
                timeplot!(chains, planet_key, color, :rv; kwargs...)
            end


            Plots.plot(
                psep,
                ppa,
                ppmra,
                ppmdec,
                ppost,
                prv,
                layout = (3,2),
                framestyle=:box,
                grid=false,
                size=(1000,1200),
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
                    DirectDetections.construct_elements(chains, p, j)
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
                        ra_hip_model[j] += -raoff(o_ra) * chains["$p[mass]"][jj]*mjup2msol/orbit.M
                        dec_hip_model[j] += -decoff(o_dec) * chains["$p[mass]"][jj]*mjup2msol/orbit.M
                        pmra_hip_model[j] += pmra(o_ra, chains["$p[mass]"][jj]*mjup2msol)
                        pmdec_hip_model[j] += pmdec(o_dec, chains["$p[mass]"][jj]*mjup2msol)
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
                        ra_gaia_model[j] += -raoff(o_ra) * chains["$p[mass]"][jj]*mjup2msol/orbit.M
                        dec_gaia_model[j] += -decoff(o_dec) * chains["$p[mass]"][jj]*mjup2msol/orbit.M
                        pmra_gaia_model[j] += pmra(o_ra, chains["$p[mass]"][jj]*mjup2msol)
                        pmdec_gaia_model[j] += pmdec(o_dec, chains["$p[mass]"][jj]*mjup2msol)
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
                        marker_z=chains[color],
                        title="Hipparcos",
                    ),
                    Plots.scatter(
                        pmra_hg_model .+ chains["pmra"][ii],
                        pmdec_hg_model .+ chains["pmdec"][ii],
                        marker_z=chains[color],
                        title="Gaia-Hipparcos",
                    ),
                    Plots.scatter(
                        pmra_gaia_model .+ chains["pmra"][ii],
                        pmdec_gaia_model .+ chains["pmdec"][ii],
                        marker_z=chains[color],
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
                    lw=2,
                )
            end

            pma_plot =  Plots.plot(
                subplots...,
                alpha=1,
                color=:plasma,
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