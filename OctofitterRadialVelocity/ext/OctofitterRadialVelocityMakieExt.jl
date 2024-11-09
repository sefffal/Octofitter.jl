module OctofitterRadialVelocityMakieExt
using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using Makie
using Statistics
using StatsBase
using AbstractGPs
# using TemporalGPs
using Dates
using AstroImages

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

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function Octofitter.rvpostplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains,
    sample_idx = argmax(results["logpost"][:]);
    show_perspective=true
)
    gs = gridspec_or_fig
    n_planet_count = length(model.system.planets)


    rv_likes = filter(model.system.observations) do obs
        obs isa StarAbsoluteRVLikelihood || obs isa OctofitterRadialVelocityMakieExt.MarginalizedStarAbsoluteRVLikelihood
    end
    # if length(rv_likes) > 1
    #     error("`rvpostplot` requires a system with only one StarAbsoluteRVLikelihood. Combine the data together into a single likelihood object.")
    # end
    # if length(rv_likes) != 1
    #     error("`rvpostplot` requires a system with a StarAbsoluteRVLikelihood.")
    # end
    # Start filling the RV plot

    els_by_planet = map(keys(model.system.planets)) do planet_key
        Octofitter.construct_elements(results,planet_key, :)
    end
    M_by_planet = [
        results[string(planet_key)*"_mass"] .* Octofitter.mjup2msol
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
    ts_grid = range(tmin-0.015delta, tmax+0.015delta,length=10000)
    # Ensure the curve has points at exactly our data points. Otherwise for fine structure
    # we might miss them unless we add a very very fine grid.
    ts = sort(vcat(ts_grid, all_epochs))
    # RV = radvel.(els[ii], ts', M[ii])
    # RV_map = radvel.(els[sample_idx], ts, M[sample_idx])

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
        )
    )
    ax_resid = Axis(
        gs[2,1],
        xlabel="time [MJD]",
        ylabel="Residuals",
    )
    linkxaxes!(ax_fit, ax_resid)
    Makie.xlims!(ax_resid, extrema(ts))

    ax_resid_hist = Axis(
        gs[2,2],
    )
    hidedecorations!(ax_resid_hist)
    linkyaxes!(ax_resid, ax_resid_hist)
    rowgap!(gs, 1, 0)


    # Horizontal zero line
    hlines!(ax_resid, 0, color=:black, linewidth=2)
    hlines!(ax_resid_hist, 0, color=:black, linewidth=2)

    # Perspective acceleration line
    if first(els_by_planet)[sample_idx] isa AbsoluteVisual && show_perspective
        lines!(ax_fit, ts_grid, radvel.(first(els_by_planet)[sample_idx], ts_grid, 0.0), color=:orange)
    end
        
    nt_format = Octofitter.mcmcchain2result(model, results)


    any_models_have_a_gp = false
    for rvs in rv_likes
        any_models_have_a_gp |= hasproperty(rvs, :gaussian_process) && !isnothing(rvs.gaussian_process)
    end


    # Main blue orbit line in top panel
    RV = sum(map(els_by_planet,M_by_planet) do els, M
        radvel.(els[sample_idx], ts, M[sample_idx])
    end)
    # Use a narrow line if we're overplotting a complicated GP
    if !any_models_have_a_gp
        lines!(ax_fit, ts, RV, color=:blue, linewidth=3)
    end



    # Create an axis for each planet, with orbit model vs. phase
    ax_phase_by_planet= Axis[]
    i_planet = 0
    for (planet, els, M) in zip(model.system.planets, els_by_planet, M_by_planet)
        i_planet += 1
        ax_phase = Axis(
            gs[2+i_planet,1],
            xlabel="Phase",
            ylabel="RV [m/s]",
            xticks=-0.5:0.1:0.5,
        )
        Makie.xlims!(ax_phase, -0.5,0.5)
        if i_planet != n_planet_count
            hidexdecorations!(ax_phase,  grid=false)
        end
        push!(ax_phase_by_planet, ax_phase)

        Label(
            gs[2+i_planet,2,],
            string(planet.name),
            fontsize = 16,
            font = :bold,
            padding = (5, 0, 0, 0),
            valign=:top,
            halign=:left,
            tellwidth=false,
            tellheight=false,
        )


        # Plot orbit models
        t_peri = periastron(els[sample_idx])
        T = period(els[sample_idx])
        phases = -0.5:0.005:0.5
        ts_phase_folded = ((phases .+ 0.5) .* T) .+ t_peri .+ T/2
        # Don't include any perspective acceleration
        RV = radvel.(nonabsvis_parent(els[sample_idx]), ts_phase_folded, M[sample_idx])
        Makie.lines!(
            ax_phase,
            phases,
            RV,
            color=:blue,
            linewidth=5
        )
    end


    # for (planet_key, els, ax_phase) in zip(keys(model.system.planets), els_by_planet, ax_phase_by_planet)

       

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


        if hasproperty(rvs,:offset_symbol)
            barycentric_rv_inst = nt_format[sample_idx][rvs.offset_symbol]
            jitter = nt_format[sample_idx][rvs.jitter_symbol]
        else
            barycentric_rv_inst = _find_rv_zero_point_maxlike(rvs, nt_format[sample_idx], getindex.(els_by_planet, sample_idx))
            jitter = nt_format[sample_idx][rvs.jitter_symbol]
        end

        # Apply barycentric rv offset correction for this instrument
        # using the MAP parameters
        rvs_off_sub .-= barycentric_rv_inst
        jitters .= jitter

        # Calculate the residuals minus the orbit model and any perspecive acceleration
        model_at_data_by_planet = map(els_by_planet,M_by_planet) do els, M
            radvel.(els[sample_idx], rvs.table.epoch, M[sample_idx])
        end
        model_at_data = sum(model_at_data_by_planet)
        resids = rvs_off_sub .- model_at_data 
        data_minus_off_and_gp  = zeros(length(resids))
        perspective_accel_to_remove = radvel.(first(els_by_planet)[sample_idx], rvs.table.epoch, 0.0)

        data = rvs.table

        ts_inst = sort(vcat(
            vec(data.epoch),
            range((extrema(data.epoch) )...,step=step(ts_grid)
        )))


        # Plot a gaussian process per-instrument
        # If not using a GP, we fit a GP with a "ZeroKernel"
        map_gp = nothing
        if hasproperty(rvs, :gaussian_process) && !isnothing(rvs.gaussian_process)
            row = results[sample_idx,:,:];
            nt = (Table((row)))[1]
            map_gp = rvs.gaussian_process(nt)
        end
        if isnothing(map_gp)
            map_gp = GP(ZeroKernel())
        # Drop TemporalGPs for now due to compilation failures
        # elseif map_gp isa TemporalGPs.LTISDE
        #     # Unwrap the temporal GPs wrapper so that we can calculate mean_and_var
        #     # We don't need the speed up provided by LTISDE for plotting once.
        #     map_gp = map_gp.f
        end

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
        y, var = mean_and_var(map_gp_posterior, ts_inst)

        # Subtract MAP GP from residuals
        resids = resids .-= mean(map_gp_posterior, vec(data.epoch))
        data_minus_off_and_gp .= rvs_off_sub .- mean(map_gp_posterior, vec(data.epoch))
        y_inst, var = mean_and_var(map_gp_posterior, ts_inst)

        # errs_data_jitter = sqrt.(
        #     data.σ_rv.^2 .+
        #     jitter.^2
        # )
        errs_data_jitter_gp = sqrt.(
            data.σ_rv.^2 .+
            jitter.^2 .+
            mean_and_var(map_gp_posterior, vec(data.epoch))[2]
        )


        epochs_all[mask] .= vec(rvs.table.epoch)
        resids_all[mask] .= resids
        uncer_all[mask] .= data.σ_rv
        uncer_jitter_gp_all[mask] .= errs_data_jitter_gp
        # rvs_all_minus_accel_minus_perspective[mask] .= rvs_off_sub .- mean(map_gp_posterior, vec(data.epoch)) .- perspective_accel_to_remove
        # errs_all_data_jitter_gp[mask] .= errs_data_jitter_gp

        # RV_sample_idxnst =  radvel.(els[sample_idx], ts_inst, M[sample_idx])
        RV_sample = sum(map(els_by_planet,M_by_planet) do els, M
            radvel.(els[sample_idx], ts_inst, M[sample_idx])
        end)
        obj = band!(ax_fit, ts_inst,
            vec(y_inst .+ RV_sample .- sqrt.(var)),# .+ jitter^2)),
            vec(y_inst .+ RV_sample .+ sqrt.(var)),# .+ jitter^2)),
            color=(Makie.wong_colors()[rv_like_idx], 0.35)
        )
        # Try to put bands behind everything else
        translate!(obj, 0, 0, -10)

        # Draw the full model ie. RV + perspective + GP
        # We darken the colour by plotting a faint black line under it
        if any_models_have_a_gp
            lines!(
                ax_fit,
                ts_inst,
                RV_sample .+ y,
                color=(:black,1),
                linewidth=0.3
            )
            lines!(
                ax_fit,
                ts_inst,
                RV_sample .+ y,
                color=(Makie.wong_colors()[rv_like_idx],0.8),
                # color=:blue,
                linewidth=0.3
            )
        end
        


        # Model plot vs raw data
        errorbars!(
            ax_fit,
            data.epoch,
            rvs_off_sub,
            errs_data_jitter_gp,
            linewidth=1,
            color="#CCC",
        )
        errorbars!(
            ax_fit,
            data.epoch,
            rvs_off_sub,
            data.σ_rv,
            # linewidth=1,
            color=Makie.wong_colors()[rv_like_idx]
        )

        errorbars!(
            ax_resid,
            data.epoch,
            resids,
            errs_data_jitter_gp,
            linewidth=1,
            color="#CCC",
        )
        errorbars!(
            ax_resid,
            data.epoch,
            resids,
            data.σ_rv,
            color=Makie.wong_colors()[rv_like_idx]
        )

        

        Makie.scatter!(
            ax_fit,
            data.epoch,
            rvs_off_sub,
            color=Makie.wong_colors()[rv_like_idx],
            markersize=4,
            strokecolor=:black,strokewidth=0.1,
        )

        Makie.scatter!(
            ax_resid,
            data.epoch,
            resids,
            color=Makie.wong_colors()[rv_like_idx],
            markersize=4,
            strokecolor=:black,strokewidth=0.1,
        )

        h = fit(Histogram, resids)
        Makie.stairs!(
            ax_resid_hist,
            h.weights,
            h.edges[1][1:end-1] .+ step(h.edges[1])/2,
            step=:center,
            color=Makie.wong_colors()[rv_like_idx],
        )
        xlims!(ax_resid_hist, low=0)

        for (ax_phase, els,rv_model_this_planet) in zip(ax_phase_by_planet, els_by_planet, model_at_data_by_planet)
            # Plot orbit models
            t_peri = periastron(els[sample_idx])
            T = period(els[sample_idx])

            # Don't include any perspective acceleration
            # Phase-folded plot
            phase_folded = mod.(data.epoch .- t_peri .- T/2, T)./T .- 0.5
            errorbars!(
                ax_phase,
                phase_folded,
                resids .+ rv_model_this_planet,
                errs_data_jitter_gp,
                linewidth=1,
                color="#CCC",
            )
            errorbars!(
                ax_phase,
                phase_folded,
                resids .+ rv_model_this_planet,
                data.σ_rv,
                # linewidth=1,
                color=Makie.wong_colors()[rv_like_idx]
            )
            
            Makie.scatter!( 
                ax_phase,
                phase_folded,
                resids .+ rv_model_this_planet,
                color=Makie.wong_colors()[rv_like_idx],
                markersize=4,
                strokecolor=:black,strokewidth=0.1,
            )
        end


    end


    # Now that we have gone through all datasets and have the residuals, go through 
    # each planet and plot binned residuals
    for (ax_phase, planet, els, M) in zip(ax_phase_by_planet, model.system.planets, els_by_planet, M_by_planet)

    # epochs_all
    # resids_all
    # uncer_all
    # uncer_jitter_gp_all
    # RV_sample = sum(map(els_by_planet,M_by_planet) do els, M
    #             radvel.(els[sample_idx], ts_inst, M[sample_idx])
    #         end)

        t_peri = periastron(els[sample_idx])
        T = period(els[sample_idx])

        # Binned values on phase folded plot
        # Noise weighted (including jitter and GP)
        # bins = -0.45:0.1:0.45
        bins = -0.495:0.05:0.495
        binned = zeros(length(bins))
        binned_unc = zeros(length(bins))
        phase_folded = mod.(epochs_all .- t_peri .- T/2, T)./T .- 0.5
        
        # the rv component we plot is the residual RV, with the signal of this *particular*
        # planet added back in.
        rv = collect(resids_all)
        rv .+= radvel.(els[sample_idx], epochs_all, M[sample_idx])

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


 
    markers =  [
        # Group 1
        [
            MarkerElement(color = Makie.wong_colors()[i], marker=:circle, markersize = 15)
            for i in 1:length(rv_likes)
        ],
        # Group 2
        [
            [
                LineElement(color = Makie.wong_colors()[i], linestyle = :solid,
                points = Point2f[((-1+i)/length(rv_likes), 0), ((-1+i)/length(rv_likes), 1)])
                for i in 1:length(rv_likes)
            ],
            LineElement(color = "#CCC", linestyle = :solid,
                points = Point2f[(0.0, 0), (0.0, 1)]),
            LineElement(color = :blue,linewidth=4,),
            MarkerElement(color = :red, strokecolor=:black, strokewidth=2, marker=:circle, markersize = 15),   
        ]
    ]
    labels = [
        [rv.instrument_name for rv in rv_likes],
        [
            "data uncertainty",
            any_models_have_a_gp ? "uncertainty,\njitter, and GP" : "uncertainty and\njitter",
            "orbit model",
            "binned",
        ]
    ]
    if show_perspective
        push!(markers[end], LineElement(color = :orange,linewidth=4,))
        push!(labels[end], "perspective")
    end

    Legend(
        gs[1:3,3],
        markers,
        labels,
        ["instrument", ""],
        valign=:top,
        framevisible=false,
    )
    Makie.rowsize!(gs, 1, Auto(2))
    Makie.rowsize!(gs, 2, Auto(1))
    for i in 1:n_planet_count
        Makie.rowsize!(gs, 2+i, Auto(2))
    end
    Makie.colgap!(gs, 1, 0)
    Makie.colsize!(gs, 2, Aspect(2,1.0))


end

function Octofitter.rvpostplot_animated(model, chain; framerate=4,compression=1, fname="rv-posterior.mp4", N=min(size(chain,1),50))
    imgs = []
    print("generating plots")
    for i in rand(1:size(chain,1), N)
        print(".")
        f = tempname()*".png"
        Octofitter.rvpostplot(model,chain,i,fname=f)
        push!(imgs, load(f))
        rm(f)
    end
    println()
    fig = Figure()
    ax = Axis(fig[1,1],autolimitaspect=1)
    hidedecorations!(ax)
    i = Observable(imgs[1])
    image!(ax, @lift(rotr90($i)))
    print("animating")
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

    jitter = getproperty(θ_system, rvlike.jitter_symbol)

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


end