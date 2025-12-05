module OctofitterRadialVelocityMakieExt
using AbstractGPs
using AstroImages
using Dates
using Makie
using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using Printf
using Statistics
using StatsBase
using LinearAlgebra: normalize

## Need three panels
# 1) Mean model (orbit + GP) and data (- mean instrument offset)
# 2) Residuals of above
# 3) Phase folded curve
function Octofitter.rvpostplot(
    model,
    results,
    args...;
    fname="$(model.system.name)-rvpostplot.png",
    kwargs...,
)
    fig = Figure(
        size=(600, 305 + 135*length(model.system.planets))
    )
    Octofitter.rvpostplot!(fig.layout, model, results,args...;kwargs...)

    # Show warning about pathfinder approximation
    if hasproperty(results.info, :sampler) && results.info.sampler == "pathfinder"
        Label(fig[0,1], Makie.rich("pathfinder approximation",color=:darkred, font=:bold), tellwidth=false)
    end


    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function Octofitter.rvpostplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains,
    sample_idx = argmax(results["logpost"][:]);
    show_perspective=true,
    show_summary=false,
    show_legend=true,
    show_phase=true,
    show_hist=true,
    ts=nothing
)
    if gridspec_or_fig isa Figure
        gs = GridLayout(gridspec_or_fig[1,1])
    else
        gs = gridspec_or_fig
    end
    n_planet_count = length(model.system.planets)


    rv_likes_all = filter(model.system.observations) do obs
        obs isa StarAbsoluteRVLikelihood || 
        obs isa OctofitterRadialVelocity.MarginalizedStarAbsoluteRVLikelihood# ||
        # obs isa OctofitterRadialVelocity.StarAbsoluteRVLikelihood_Celerite_Shared
    end

    #  filter RV objects to exclude any that are outside the plotted region
    if isnothing(ts)
        rv_likes = rv_likes_all
    else
        rv_likes = filter(rv_likes_all) do obs
            f,l = extrema(obs.table.epoch)
            range_overlaps = f <= last(ts) && l >= first(ts)
            if range_overlaps
                return any(epoch -> first(ts) <= epoch <= last(ts), obs.table.epoch)
            end
            return false
        end
    end

    # if length(rv_likes) > 1
    #     error("`rvpostplot` requires a system with only one StarAbsoluteRVLikelihood. Combine the data together into a single likelihood object.")
    # end
    # if length(rv_likes) != 1
    #     error("`rvpostplot` requires a system with a StarAbsoluteRVLikelihood.")
    # end
    # Start filling the RV plot

    els_by_planet = map(keys(model.system.planets)) do planet_key
        Octofitter.construct_elements(model, results, planet_key, sample_idx)
    end
    M_by_planet = [
        results[string(planet_key)*"_mass"][sample_idx] .* Octofitter.mjup2msol
        for planet_key in keys(model.system.planets)
    ]

    # Sometimes we want the true RV, including eg perspecive acceleration, and
    # sometimes we want just the perturbation orbit. 
    # This function returns an orbit-element's parent orbit-element if its
    # an AbsoluteVisualOrbit, otherwise, just returns the passed in element.
    nonabsvis_parent(element::AbsoluteVisual) = element.parent 
    nonabsvis_parent(element::AbstractOrbit) = element


    # Model plot vs raw data
    all_epochs = vec(vcat((rvs.table.epoch for rvs in rv_likes)...))
    tmin, tmax = extrema(all_epochs)
    delta = tmax-tmin
    if delta == 0
        delta = 1
    end
    if !isnothing(ts)
        ts_grid = ts
    else
        ts_grid = range(tmin-0.015delta, tmax+0.015delta,length=10000)
    end
    # Ensure the curve has points at exactly our data points. Otherwise for fine structure
    # we might miss them unless we add a very very fine grid.
    ts = sort(vcat(ts_grid, all_epochs))
    # RV = radvel.(els[ii], ts', M[ii])
    # RV_map = radvel.(els, ts, M)

    # For secondary date axis on top
    date_start = mjd2date(ts[begin])
    date_end = mjd2date(ts[end])
    date_start = Date(year(date_start), month(date_start))
    date_end = Date(year(date_end), month(date_end))
    dates = range(date_start, date_end, step=Year(1))
    dates_str = string.(year.(dates))
    if length(dates) == 1
        dates = range(date_start, date_end, step=Month(1))
        dates_str = map(d->string(year(d),"-",lpad(month(d),2,'0')),dates)
    else
        year_step = 1
        while length(dates) > 8
            year_step += 1
            dates = range(date_start, date_end, step=Year(year_step))
        end
        dates_str = string.(year.(dates))
    end
    ax_fit = Axis(
        gs[1,1],
        ylabel="RV [m/s]",
        xaxisposition=:top,
        xticks=(
            mjd.(string.(dates)),
            dates_str
        ),
        xgridvisible=false,
        ygridvisible=false,
        alignmode=Makie.Mixed(;
            left=nothing,
            right=nothing,
            top=0,
            bottom=nothing,
        )
    )
    ax_resid = Axis(
        gs[2,1],
        xlabel="time [MJD]",
        ylabel="Residuals",
        xgridvisible=false,
        ygridvisible=false,
    )
    linkxaxes!(ax_fit, ax_resid)
    Makie.xlims!(ax_resid, extrema(ts))

    if show_hist
        ax_resid_hist = Axis(
            gs[2,2],
        )
        hidedecorations!(ax_resid_hist)
        linkyaxes!(ax_resid, ax_resid_hist)
    end
    rowgap!(gs, 1, 0)


    # Horizontal zero line
    hlines!(ax_resid, 0, color=:black, linewidth=2)
    show_hist && hlines!(ax_resid_hist, 0, color=:black, linewidth=2)

    # Perspective acceleration line
    if first(els_by_planet) isa AbsoluteVisual && show_perspective
        lines!(ax_fit, ts_grid, radvel.(first(els_by_planet), ts_grid, 0.0), color=:orange)
    end
        
    nt_format = Octofitter.mcmcchain2result(model, results, sample_idx:sample_idx)[1]


    any_models_have_a_gp = false
    for rvs in rv_likes
        any_models_have_a_gp |= hasproperty(rvs, :gaussian_process) && !isnothing(rvs.gaussian_process)
    end


    # Main blue orbit line in top panel
    RV = sum(map(els_by_planet,M_by_planet) do els, M
        radvel.(els, ts, M)
    end)
    # Use a narrow line if we're overplotting a complicated GP
    # if !any_models_have_a_gp
        lines!(ax_fit, ts, RV, color=:blue, linewidth=0.5)
    # end



    # Create an axis for each planet, with orbit model vs. phase
    ax_phase_by_planet= Axis[]
    i_planet = 0
    if show_phase
        for (planet, els, M) in zip(model.system.planets, els_by_planet, M_by_planet)
            i_planet += 1
            ax_phase = Axis(
                gs[2+i_planet,1],
                xlabel="Phase",
                ylabel="RV [m/s]",
                xticks=-0.5:0.1:0.5,
                xgridvisible=false,
                ygridvisible=false,
            )
            Makie.xlims!(ax_phase, -0.5,0.5)
            if i_planet != n_planet_count
                hidexdecorations!(ax_phase,  grid=false)
            end
            push!(ax_phase_by_planet, ax_phase)

            text = Makie.rich(string(planet.name), font=:bold, fontsize=16)
            Label(
                gs[2+i_planet,2,],
                text,
                padding = (5, 0, 0, 0),
                valign=:top,
                halign=:left,
                justification=:left,
                tellwidth=false,
                tellheight=false,
            )
            if show_summary
                el = Octofitter.construct_elements(model, results,planet.name,:)
                credP = cred_internal_format(quantile(period.(el), (0.14, 0.5,0.84))...)
                text = Makie.rich( 
                    "P", Makie.subscript(string(planet.name)),
                    " = $credP d",
                    font="Monaco",
                    fontsize=10
                )
                credE = cred_internal_format(quantile(eccentricity.(el), (0.14, 0.5,0.84))...)
                text = Makie.rich(text, "\n", 
                    "e", Makie.subscript(string(planet.name)),
                    " = $credE",
                    font="Monaco",
                    fontsize=10
                )
                credMass = cred_internal_format(quantile(results["$(planet.name)_mass"][:], (0.14, 0.5,0.84))...)
                text = Makie.rich(text, "\n", 
                    "m", Makie.subscript(string(planet.name)),
                    " = $credMass M", Makie.subscript("jup"),
                    font="Monaco",
                    fontsize=10
                )
                try
                    inclination(el[1])
                catch
                    text = Makie.rich(
                        text, " (min)",
                        font="Monaco",
                        fontsize=10
                    )
                end

                Label(
                    gs[2+i_planet,2,],
                    text,
                    padding = (5, 0, 0, 0),
                    valign=:bottom,
                    halign=:left,
                    justification=:left,
                    tellwidth=false,
                    tellheight=false,
                )
            end


            # Plot orbit models
            T = period(els)
            phases = -0.5:0.005:0.5
            ts_phase_folded = phases .* T
            # Shift to RV zero crossing
            rv_zero_cross = 1
            for i in LinearIndices(ts_phase_folded)[2:end]
                a = radvel(els, ts_phase_folded[i-1])
                b = radvel(els, ts_phase_folded[i])
                if a <= 0 <= b
                    rv_zero_cross = i
                    break
                end
            end
            ts_phase_folded = ts_phase_folded .+ ts_phase_folded[rv_zero_cross]
            # Don't include any perspective acceleration
            RV = radvel.(nonabsvis_parent(els), ts_phase_folded, M)
            Makie.lines!(
                ax_phase,
                phases,
                RV,
                color=:blue,
                linewidth=5
            )
        end
    end


    marker_symbols = [
        :circle,
        :rect,
        :diamond,
        :hexagon,
        :cross,
        :xcross,
        :utriangle,
        :dtriangle,
        :ltriangle,
        :rtriangle,
        :pentagon,
        :star4,
        :star5,
        :star6,
        :star8,
        :vline,
        :hline,
    ]

    

    # Calculate RVs minus the median instrument-specific offsets.
    
    # For the phase-folded binned data, 
    # we also collect all data minus perspective acceleration and any GP,
    # as well as the data uncertainties, and data + jitter + GP quadrature uncertainties
    N = 0
    masks = []
    for rvs in rv_likes
        push!(masks, N .+ (1:length(rvs.table.rv)))
        N += length(rvs.table.rv)
    end
    epochs_all = zeros(N)
    resids_all = zeros(N)
    uncer_all = zeros(N)
    uncer_jitter_gp_all = zeros(N)

    rv_like_idx = 0
    for rvs in rv_likes
        rv_like_idx += 1 
        mask = masks[rv_like_idx]
        rvs_off_sub = collect(rvs.table.rv)
        jitters = zeros(length(rvs_off_sub))

        color = Makie.wong_colors()[mod1(rv_like_idx, length(Makie.wong_colors()))]
        marker_symbol = marker_symbols[mod1(rv_like_idx, length(marker_symbols))]

        like_obj_name = Octofitter.normalizename(likelihoodname(rvs))

        if rvs isa StarAbsoluteRVLikelihood
            jitter = nt_format.observations[like_obj_name].jitter
            barycentric_rv_inst = nt_format.observations[like_obj_name].offset
        elseif rvs isa MarginalizedStarAbsoluteRVLikelihood
            # TODO: marginalized RV likelihood

            barycentric_rv_inst = _find_rv_zero_point_maxlike(rvs, nt_format, els_by_planet)
            jitter = nt_format.observations[like_obj_name].jitter
        else
            error("plotting not yet implemented for this type of data")
        end
        if ismissing(barycentric_rv_inst )
            barycentric_rv_inst = 0.
        end
        if ismissing(jitter)
            jitter = 0.
        end

        # Apply barycentric rv offset correction for this instrument
        # using the MAP parameters

        # if rvs isa OctofitterRadialVelocity.StarAbsoluteRVLikelihood_Celerite_Shared
        #     for i in eachindex(rvs_off_sub)
        #         inst_idx = rvs.table.inst_idx[i]
        #         rvs_off_sub[i] -= barycentric_rv_inst[inst_idx]
        #         rvs_off_sub[i] -= rvs.trend_function(nt_format, rvs.table.epoch[i], inst_idx)
        #         jitters[i] = jitter[inst_idx]
        #     end
        # else
            rvs_off_sub .-= barycentric_rv_inst
            θ_obs = nt_format.observations[like_obj_name]
            trend = rvs.trend_function.((θ_obs,), rvs.table.epoch)
            trend = map(trend) do t
                if ismissing(t) 
                    0.
                else
                    t
                end
            end
            rvs_off_sub .-= trend
            jitters .= jitter
        # end

        # Calculate the residuals minus the orbit model and any perspecive acceleration
        model_at_data_by_planet = map(els_by_planet,M_by_planet) do els, M
            radvel.(els, rvs.table.epoch, M)
        end
        model_at_data = sum(model_at_data_by_planet)
        resids = rvs_off_sub .- model_at_data 
        data_minus_off_and_gp  = zeros(length(resids))
        perspective_accel_to_remove = radvel.(first(els_by_planet), rvs.table.epoch, 0.0)

        data = rvs.table

        ts_inst = sort(vcat(
            vec(data.epoch),
            range((extrema(data.epoch) )...,step=step(ts_grid)
        )))


        # Plot a gaussian process per-instrument
        # If not using a GP, we fit a GP with a "ZeroKernel"
        map_gp = nothing
        if hasproperty(rvs, :gaussian_process) && !isnothing(rvs.gaussian_process)
            θ_obs = nt_format.observations[Octofitter.normalizename(likelihoodname(rvs))]
            map_gp = rvs.gaussian_process(θ_obs)
        end
        if isnothing(map_gp)
            map_gp = GP(ZeroKernel())
        elseif isdefined(Main, :TemporalGPs) && map_gp isa Main.TemporalGPs.LTISDE
            # Unwrap the temporal GPs wrapper so that we can calculate mean_and_var
            # We don't need the speed up provided by LTISDE for plotting once.
            map_gp = map_gp.f
        end

        if map_gp isa OctofitterRadialVelocity.Celerite.CeleriteGP

            # Compute GP model
            OctofitterRadialVelocity.Celerite.compute!(map_gp, vec(rvs.table.epoch), vec(
                sqrt.(rvs.table.σ_rv.^2 .+ jitters.^2)
            ))

            y_inst, var = OctofitterRadialVelocity.Celerite.predict(map_gp, vec(resids), collect(ts_inst); return_var=true)
            y_at_dat, var_at_dat = OctofitterRadialVelocity.Celerite.predict(map_gp, vec(resids), collect(vec(data.epoch)); return_var=true)
            var = max.(0, var)
            var_at_dat = max.(0, var_at_dat)

            resids = resids .-= y_at_dat

            # errs_data_jitter = sqrt.(
            #     data.σ_rv.^2 .+
            #     jitter.^2
            # )
            errs_data_jitter_gp = sqrt.(
                data.σ_rv.^2 .+
                jitters.^2 .+
                var_at_dat
            )
        else

            fx = map_gp(
                # x
                vec(rvs.table.epoch),
                # y-err
                vec(
                    sqrt.(rvs.table.σ_rv.^2 .+ jitters.^2)
                )
            )
            # condition GP on residuals (data - orbit - inst offsets)
            map_gp_posterior = posterior(fx, vec(resids))
            y_inst, var = mean_and_var(map_gp_posterior, ts_inst)

            # Subtract MAP GP from residuals
            resids = resids .-= mean(map_gp_posterior, vec(data.epoch))
            data_minus_off_and_gp .= rvs_off_sub .- mean(map_gp_posterior, vec(data.epoch))

            # errs_data_jitter = sqrt.(
            #     data.σ_rv.^2 .+
            #     jitters.^2
            # )
            errs_data_jitter_gp = sqrt.(
                data.σ_rv.^2 .+
                jitters.^2 .+
                mean_and_var(map_gp_posterior, vec(data.epoch))[2]
            )
        end


        epochs_all[mask] .= vec(rvs.table.epoch)
        resids_all[mask] .= resids .- perspective_accel_to_remove
        uncer_all[mask] .= data.σ_rv
        uncer_jitter_gp_all[mask] .= errs_data_jitter_gp
        # rvs_all_minus_accel_minus_perspective[mask] .= rvs_off_sub .- mean(map_gp_posterior, vec(data.epoch)) .- perspective_accel_to_remove
        # errs_all_data_jitter_gp[mask] .= errs_data_jitter_gp

        # RV_sample_idxnst =  radvel.(els, ts_inst, M)
        RV_sample = sum(map(els_by_planet,M_by_planet) do els, M
            radvel.(els, ts_inst, M)
        end)
        obj = band!(ax_fit, ts_inst,
            vec(y_inst .+ RV_sample .- sqrt.(var)),# .+ jitter^2)),
            vec(y_inst .+ RV_sample .+ sqrt.(var)),# .+ jitter^2)),
            color=(color, 0.35)
        )
        # Try to put bands behind everything else
        translate!(obj, 0, 0, -10)

        # Draw the full model ie. RV + perspective + GP
        # We darken the colour by plotting a faint black line under it
        if any_models_have_a_gp
            lines!(
                ax_fit,
                ts_inst,
                RV_sample .+ y_inst,
                color=(:black,1),
                linewidth=0.1
            )
            lines!(
                ax_fit,
                ts_inst,
                RV_sample .+ y_inst,
                color=(color,0.8),
                # color=:blue,
                linewidth=0.1
            )
        end
        


        # Model plot vs raw data
        errorbars!(
            ax_fit,
            data.epoch,
            rvs_off_sub,
            errs_data_jitter_gp,
            linewidth=0.5,
            color="#CCC",
        )
        errorbars!(
            ax_fit,
            data.epoch,
            rvs_off_sub,
            data.σ_rv,
            linewidth=0.5,
            color=color
        )

        errorbars!(
            ax_resid,
            data.epoch,
            resids,
            errs_data_jitter_gp,
            linewidth=0.5,
            color="#CCC",
        )
        errorbars!(
            ax_resid,
            data.epoch,
            resids,
            data.σ_rv,
            linewidth=0.5,
            color=color
        )

        

        Makie.scatter!(
            ax_fit,
            data.epoch,
            rvs_off_sub,
            color=color,
            markersize=4,
            strokecolor=:black,
            strokewidth=0.1,
            marker=marker_symbol
        )

        Makie.scatter!(
            ax_resid,
            data.epoch,
            resids,
            color=color,
            markersize=4,
            strokecolor=:black,
            strokewidth=0.1,
            marker=marker_symbol
        )

        
        if show_hist
            h = fit(Histogram, resids)
            Makie.stairs!(
                ax_resid_hist,
                h.weights,
                h.edges[1][1:end-1] .+ step(h.edges[1])/2,
                step=:center,
                color=color,
                linewidth=0.5
            )
            xlims!(ax_resid_hist, low=0)
        end

        for (ax_phase, els,rv_model_this_planet) in zip(ax_phase_by_planet, els_by_planet, model_at_data_by_planet)
            # Plot orbit models


            # Phase-folded plot
            T = period(els)
            # Shift to RV zero crossing
            # Find orbit fraction at zero-crossing
            ts_phase_folded = ( -0.5:0.005:0.5) .* T
            rv_zero_cross = 1
            for i in LinearIndices(ts_phase_folded)[2:end]
                a = radvel(els, ts_phase_folded[i-1])
                b = radvel(els, ts_phase_folded[i])
                if a <= 0 <= b
                    rv_zero_cross = i
                    break
                end
            end
            t_offset = ts_phase_folded[rv_zero_cross]            
            phase_folded = mod.(data.epoch .- t_offset .- T/2, T)./T .- 0.5
             

            errorbars!(
                ax_phase,
                phase_folded,
                resids .+ rv_model_this_planet .- perspective_accel_to_remove,
                errs_data_jitter_gp,
                linewidth=0.5,
                color="#CCC",
            )
            errorbars!(
                ax_phase,
                phase_folded,
                resids .+ rv_model_this_planet .- perspective_accel_to_remove,
                data.σ_rv,
                linewidth=0.5,
                color=Makie.wong_colors()[mod1(rv_like_idx, length(Makie.wong_colors()))]
            )
            
            Makie.scatter!( 
                ax_phase,
                phase_folded,
                resids .+ rv_model_this_planet .- perspective_accel_to_remove,
                color=Makie.wong_colors()[mod1(rv_like_idx, length(Makie.wong_colors()))],
                markersize=4,
                strokecolor=:black,
                strokewidth=0.1,
                marker=marker_symbol
            )
        end


    end


    # println(epochs_all)
    # println(resids_all)
    # println(uncer_all)

    # Now that we have gone through all datasets and have the residuals, go through 
    # each planet and plot binned residuals
    for (ax_phase, planet, els, M) in zip(ax_phase_by_planet, model.system.planets, els_by_planet, M_by_planet)

        T = period(els)
        # Shift to RV zero crossing
        # Find orbit fraction at zero-crossing
        ts_phase_folded = ( -0.5:0.005:0.5) .* T
        rv_zero_cross = 1
        for i in LinearIndices(ts_phase_folded)[2:end]
            a = radvel(els, ts_phase_folded[i-1])
            b = radvel(els, ts_phase_folded[i])
            if a <= 0 <= b
                rv_zero_cross = i
                break
            end
        end
        t_offset = ts_phase_folded[rv_zero_cross]            

        # Binned values on phase folded plot
        # Noise weighted (including jitter and GP)
        # bins = -0.45:0.1:0.45
        bins = -0.495:0.05:0.495
        binned = zeros(length(bins))
        binned_unc = zeros(length(bins))
        phase_folded = mod.(epochs_all .- t_offset .- T/2, T)./T .- 0.5
        
        # the rv component we plot is the residual RV, with the signal of this *particular*
        # planet added back in.
        rv = collect(resids_all)
        rv .+= radvel.(els, epochs_all, M)

        for (i,bin_cent) in enumerate(bins)
            mask = bin_cent - step(bins)/2 .<= phase_folded .<= bin_cent + step(bins/2)
            if count(mask) == 0
                binned[i] = NaN
                continue
            end
            binned[i] = mean(
                rv[mask],
                ProbabilityWeights(1 ./ uncer_jitter_gp_all[mask].^2)
            )
            binned_unc[i] = std(
                rv[mask],
                ProbabilityWeights(1 ./ uncer_jitter_gp_all[mask].^2)
            )
        end
        errorbars!(
            ax_phase,
            bins,
            binned,
            binned_unc,
            color=:black,
            linewidth=2,
        )
        scatter!(
            ax_phase,
            bins,
            binned,
            color=:red,
            markersize=10,
            strokecolor=:black,
            strokewidth=2,
        )
    end


 
    if show_legend
        markers =  [
            # Group 1
            [
                MarkerElement(color = Makie.wong_colors()[mod1(i,length(Makie.wong_colors()))], marker=marker_symbols[mod1(i,length(marker_symbols))], markersize = 15)
                for i in 1:length(rv_likes_all)
            ],
            # Group 2
            [
                [
                    LineElement(color = Makie.wong_colors()[mod1(i,length(Makie.wong_colors()))], linestyle = :solid,
                    points = Point2f[((-1+i)/length(rv_likes_all), 0), ((-1+i)/length(rv_likes_all), 1)])
                    for i in 1:length(rv_likes_all)
                ],
                LineElement(color = "#CCC", linestyle = :solid,
                    points = Point2f[(0.0, 0), (0.0, 1)]),
                LineElement(color = :blue,linewidth=4,),
                MarkerElement(color = :red, strokecolor=:black, strokewidth=2, marker=:circle, markersize = 15),   
            ]
        ]
        labels = [
            [likelihoodname(rv) for rv in rv_likes_all],
            [
                "data uncertainty",
                any_models_have_a_gp ? "uncertainty,\njitter, and GP" : "uncertainty and\njitter",
                "orbit model",
                "binned",
            ]
        ]
        if first(els_by_planet) isa AbsoluteVisual && show_perspective
            push!(markers[end], LineElement(color = :orange,linewidth=4,))
            push!(labels[end], "perspective")
        end

        xs = show_phase ? (1:3) : (1:2)
        ys = show_hist ? 3 : 2
        Legend(
            gs[xs,ys],
            markers,
            labels,
            ["instrument", "legend"],
            valign=:top,
            framevisible=false,
        )
    end
    Makie.rowsize!(gs, 1, Auto(2.4))
    if show_phase
        Makie.rowsize!(gs, 2, Auto(1))
        for i in 1:n_planet_count
            Makie.rowsize!(gs, 2+i, Auto(2))
        end
    end
    if show_hist
        Makie.colgap!(gs, 1, 0)
        Makie.colsize!(gs, 2, Aspect(2,1.0))
    end

end

function Octofitter.rvpostplot_animated(model, chain; framerate=4,compression=1, fname="rv-posterior.mp4", N=min(size(chain,1),50), kwargs...)
    return Octofitter.rvpostplot_animated(identity, model, chain; framerate,compression, fname, N, kwargs...)
end
function Octofitter.rvpostplot_animated(callback::Base.Callable, model, chain; framerate=4,compression=1, fname="rv-posterior.mp4", N=min(size(chain,1),50), kwargs...)
    imgs = []
    print("generating plots:")
    for i in rand(1:size(chain,1), N)
        print(".")
        f = tempname()*".png"
        fig = Octofitter.rvpostplot(model,chain,i; fname=f, kwargs...)
        callback(fig)
        save(f,fig)
        push!(imgs, load(f))
        rm(f)
    end
    println()
    fig = Figure(
        size=reverse(size(first(imgs))),
        figure_padding=(0,0,0,0)
    )
    ax = Axis(fig[1,1],
        autolimitaspect=1,
        leftspinevisible=false,
        rightspinevisible=false,
        topspinevisible=false,
        bottomspinevisible=false,
    )
    hidedecorations!(ax)
    i = Observable(imgs[1])
    image!(ax, @lift(rotr90($i)))
    print("animating       :")
    Makie.record(fig, fname, imgs; framerate, compression) do img
        print(".")
        i[] = img
    end
    println()
    return fname
end


function _find_rv_zero_point_maxlike(
    rvlike,#::MarginalizedStarAbsoluteRVLikelihood,
    θ_system,
    planet_orbits::Tuple,
)
    T = Octofitter._system_number_type(θ_system)

    # Data for this instrument:
    epochs = rvlike.table.epoch
    σ_rvs = rvlike.table.σ_rv
    rvs = rvlike.table.rv

    like_obj_name = Octofitter.normalizename(likelihoodname(rvlike))
    jitter = θ_system.observations[like_obj_name].jitter

    if hasproperty(θ_system, :observations) && hasproperty(θ_system.observations, like_obj_name)
        θ_obs = θ_system.observations[like_obj_name]
        jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : 0
    end

    # RV residual calculation: measured RV - model
    resid = zeros(T, length(rvs))
    resid .+= rvs
    # Start with model ^

    # Go through all planets and subtract their modelled influence on the RV signal:
    for planet_i in eachindex(planet_orbits)
        orbit = planet_orbits[planet_i]
        planet_mass = θ_system.planets[planet_i].mass
        for i_epoch in eachindex(epochs)
            sol = orbitsolve(orbit, epochs[i_epoch])
            resid[i_epoch] -= radvel(sol, planet_mass*Octofitter.mjup2msol)
        end
    end
    
    # Marginalize out the instrument zero point using math from the Orvara paper
    A = zero(T)
    B = zero(T)
    for i_epoch in eachindex(epochs)
        # The noise variance per observation is the measurement noise and the jitter added
        # in quadrature
        var = σ_rvs[i_epoch]^2 + jitter^2
        A += 1/var
        B -= 2resid[i_epoch]/var
    end

    rv0 = B/2A

    return -rv0
end




# From PairPlots.jl
function cred_internal_format(low,mid,high)
    largest_error = max(abs(high-mid), abs(mid-low))
    # Fallback for series with no variance
    if largest_error == 0
        return "$mid"
    end

    digits_after_dot = max(0, 1 - round(Int, log10(largest_error)))
    return @sprintf(
        "%.*f ± %.*f",
        digits_after_dot, mid,
        digits_after_dot, largest_error,
    )

    return title 
end

end