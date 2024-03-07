module OctofitterRadialVelocityMakieExt
using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using Makie
using Statistics
using StatsBase
using AbstractGPs
using TemporalGPs
using Dates

## Need three panels
# 1) Mean model (orbit + GP) and data (- mean instrument offset)
# 2) Residuals of above
# 3) Phase folded curve
function OctofitterRadialVelocity.rvpostplot(
    model,
    results,
    fname="$(model.system.name)-rvpostplot.png",
    args...
)
    fig = Figure()
    OctofitterRadialVelocity.rvpostplot!(fig.layout, model, results,args...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function OctofitterRadialVelocity.rvpostplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains,
    sample_idx = argmax(results["logpost"][:]),
    planet_key = first(keys(model.system.planets))
)
    gs = gridspec_or_fig


    rvs = only(filter(obs->obs isa StarAbsoluteRVLikelihood, model.system.observations))
    ii = rand(1:size(results,1),100)


    # Start filling the RV plot
    els = Octofitter.construct_elements(results,planet_key, :)
    M = (results[string(planet_key)*"_mass"] .* PlanetOrbits.mjup2msol)

    # For phase-folded plot
    t_peri = periastron(els[sample_idx])
    T = period(els[sample_idx])

    # Model plot vs raw data
    ts_grid = range((extrema(rvs.table.epoch) )...,length=10000)
    # Ensure the curve has points at exactly our data points. Otherwise for fine structure
    # we might miss them unless we add a very very fine grid.
    ts = sort(vcat(ts_grid, vec(rvs.table.epoch)))
    RV = radvel.(els[ii], ts', M[ii])
    RV_map = radvel.(els[sample_idx], ts, M[sample_idx])

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
    # hidexdecorations!(ax_fit,grid=false)

    ax_phase = Axis(
        gs[3,1],
        xlabel="Phase",
        ylabel="RV [m/s]",
        xticks=-0.5:0.1:0.5,
    )
    # xlims!(ax_phase, -0.5,0.5)
    rowgap!(gs, 1, 0)

    # Horizontal zero line
    hlines!(ax_resid, 0, color=:blue, linewidth=5)

    # for rv_post in eachrow(RV)
    #     lines!(
    #         ax_fit,
    #         ts,
    #         rv_post,
    #         color=(:blue, 0.05)
    #     )
    # end

    # Calculate RVs minus the median instrument-specific offsets.
    # Use the MAP parameter values
    rvs_off_sub = collect(rvs.table.rv)
    jitters_all = zeros(length(rvs_off_sub))
    for inst_idx in 1:maximum(rvs.table.inst_idx)
        barycentric_rv_inst = results["rv0_$inst_idx"][sample_idx]
        jitter = results["jitter_$inst_idx"][sample_idx]
        thisinst_mask = vec(rvs.table.inst_idx.==inst_idx)
        # Apply barycentric rv offset correction for this instrument
        # using the MAP parameters
        rvs_off_sub[thisinst_mask] .+= barycentric_rv_inst
        jitters_all[thisinst_mask] .= jitter
    end
    # Calculate the residuals minus the orbit model
    # Use the MAP values
    model_at_data = radvel.(els[sample_idx], rvs.table.epoch, M[sample_idx]) 
    resids_all = rvs_off_sub .- model_at_data 
    rvs_off_gp_sub = zeros(length(resids_all))
    errs_all = zeros(length(resids_all))
    data_minus_off_and_gp  = zeros(length(resids_all))

    # Model plot vs phase-folded data
    phases = -0.5:0.02:0.5
    ts_phase_folded = ((phases .+ 0.5) .* T) .+ t_peri .+ T/4
    RV = radvel.(els[sample_idx], ts_phase_folded, M[sample_idx])
    Makie.lines!(
        ax_phase,
        phases,
        RV,
        color=:blue,
        linewidth=5
    )
    Makie.xlims!(ax_phase, -0.5,0.5)
    
    # Main blue orbit line in top panel
    # lines!(ax_fit, ts, RV_map, color=:darkblue, linewidth=0.2)


    # plot!(ax_fit, ts, posterior_gp; bandscale=1, color=(:black,0.3))
    for inst_idx in 1:maximum(rvs.table.inst_idx)
        # barycentric_rv_inst = results["rv0_$inst_idx"][sample_idx]
        thisinst_mask = vec(rvs.table.inst_idx.==inst_idx)
        jitters = jitters_all[thisinst_mask]
        jitter = only(unique(jitters))
        data = rvs.table[thisinst_mask]
        resid = resids_all[thisinst_mask]
        rvs_off_sub_this = rvs_off_sub[thisinst_mask]

        ts_inst = sort(vcat(
            vec(data.epoch),
            range((extrema(data.epoch) )...,step=step(ts_grid)
        )))


        ## Gaussian Process. We have one per instrument
        # I think this is implementable as 
        # η₁ = h # the amplitude of the covariance function
        # η₂  # equivalent to the evolution timescale of features in the stellar surface which produce activity-induced RV variations
        # η₃  # is equivalent to the stellar rotation period
        # η₄   # gives a measure of the level of high-frequency variability structure in the GP model
        # η₃? = θ  # the period of the correlated noise 
        # η₂? = λ  # characteristic decay timescale of the correlation (spot lifetime)
        # η₄? = ω # coherence scale (sometimes called the structure parameter)


        # η₁ = results["gp_η₁"][sample_idx]
        # η₂ = results["gp_η₂"][sample_idx]
        # η₃ = results["gp_η₃"][sample_idx]
        # η₄ = results["gp_η₄"][sample_idx]
        # kernel = η₁^2 * 
        #             # (Matern52Kernel() ∘ ScaleTransform(2/η₂)) *  
        #             # (ApproxPeriodicKernel{7}(r=η₄) ∘ ScaleTransform(1/η₃)) #
        #             # This is a closer match to what other packages use
        #             (SqExponentialKernel() ∘ ScaleTransform(1/(η₂))) *
        #             (PeriodicKernel(r=[η₄]) ∘ ScaleTransform(1/(η₃)))
        # gp_naive = GP(kernel)
        # map_gp = gp_naive
        # # map_gp = to_sde(gp_naive, SArrayStorage(Float64))

        map_gp = nothing
        if !isnothing(rvs.gaussian_process)
            row = results[sample_idx,:,:];
            nt = (Table((row)))[1]
            map_gp = rvs.gaussian_process(nt)
        end
        if isnothing(map_gp)
            map_gp = GP(ZeroKernel())
        elseif map_gp isa TemporalGPs.LTISDE
            # Unwrap the temporal GPs wrapper so that we can calculate mean_and_var
            # We don't need the speed up provided by LTISDE for plotting once.
            map_gp = map_gp.f
        end

        fx = map_gp(
            # x
            vec(rvs.table.epoch[thisinst_mask]),
            # y-err
            vec(
                sqrt.(rvs.table.σ_rv[thisinst_mask].^2 .+ jitters_all[thisinst_mask].^2)
            )
        )
        # condition GP on residuals (data - orbit - inst offsets)
        map_gp_posterior = posterior(fx, vec(resids_all[thisinst_mask]))
        y, var = mean_and_var(map_gp_posterior, ts_inst)

        # We have a single GP fit. We plot the mean directly but not the std.
        # We want to show separate uncertainty bands per instrument by adding
        # in the jitter in quadrature

        # Subtract MAP GP from residuals
        resid = resids_all[thisinst_mask] .-= mean(map_gp_posterior, vec(data.epoch))
        # rvs_off_gp_sub[thisinst_mask] .= rvs_off_sub[thisinst_mask] .- mean(map_gp_posterior, vec(data.epoch))
        data_minus_off_and_gp[thisinst_mask] .= rvs_off_sub_this .- mean(map_gp_posterior, vec(data.epoch))
        y_inst, var = mean_and_var(map_gp_posterior, ts_inst)

        errs = sqrt.(
            data.σ_rv.^2 .+
            jitter.^2 .+
            mean_and_var(map_gp_posterior, vec(data.epoch))[2]
        )
        errs_all[thisinst_mask] .= errs

        RV_sample_idxnst =  radvel.(els[sample_idx], ts_inst, M[sample_idx])
        band!(ax_fit, ts_inst,
            vec(y_inst .+ RV_sample_idxnst .- sqrt.(var .+ jitter^2)),
            vec(y_inst .+ RV_sample_idxnst .+ sqrt.(var .+ jitter^2)),
            color=(Makie.wong_colors()[inst_idx], 0.35)
        )
        lines!(
            ax_fit,
            ts_inst,
            radvel.(els[sample_idx], ts_inst, M[sample_idx]) .+ y,
            color=(:black,1),
            linewidth=0.3
        )
        lines!(
            ax_fit,
            ts_inst,
            radvel.(els[sample_idx], ts_inst, M[sample_idx]) .+ y,
            color=(Makie.wong_colors()[inst_idx],0.8),
            # color=:blue,
            linewidth=0.3
        )
        


        # Model plot vs raw data
        errorbars!(
            ax_fit,
            data.epoch,
            rvs_off_sub_this,
            data.σ_rv,
            # linewidth=1,
            color=Makie.wong_colors()[inst_idx]
        )
        # scatter!(
        #     ax_fit,
        #     data.epoch,
        #     rvs_off_sub_this,
        #     color=Makie.wong_colors()[inst_idx],
        #     markersize=4,
        # )

        errorbars!(
            ax_resid,
            data.epoch,
            resid,
            errs,
            linewidth=1,
            color="#CCC",
        )
        errorbars!(
            ax_resid,
            data.epoch,
            resid,
            data.σ_rv,
            # linewidth=1,
            color=Makie.wong_colors()[inst_idx]
        )
        # scatter!(
        #     ax_resid,
        #     data.epoch,
        #     resid,
        #     color=Makie.wong_colors()[inst_idx],
        #     markersize=4
        # )

        # Phase-folded plot
        phase_folded = mod.(data.epoch .- t_peri .- T/4, T)./T .- 0.5
        errorbars!(
            ax_phase,
            phase_folded,
            data_minus_off_and_gp[thisinst_mask],
            errs,
            linewidth=1,
            color="#CCC",
        )
        errorbars!(
            ax_phase,
            phase_folded,
            data_minus_off_and_gp[thisinst_mask],
            data.σ_rv,
            # linewidth=1,
            color=Makie.wong_colors()[inst_idx]
        )
        # scatter!(
        #      ax_phase,
        #     phase_folded,
        # #     rvs_off_sub_this .- mean(map_gp_posterior, vec(data.epoch)),
        #     data_minus_off_and_gp[thisinst_mask],
        #     color=Makie.wong_colors()[inst_idx],
        #     markersize=4
        # )
    end
    for inst_idx in 1:maximum(rvs.table.inst_idx)
        # barycentric_rv_inst = results["rv0_$inst_idx"][sample_idx]
        thisinst_mask = vec(rvs.table.inst_idx.==inst_idx)
        jitter = jitters_all[thisinst_mask]
        data = rvs.table[thisinst_mask]
        resid = resids_all[thisinst_mask]
        rvs_off_sub_this = rvs_off_sub[thisinst_mask]

        Makie.scatter!(
            ax_fit,
            data.epoch,
            rvs_off_sub_this,
            color=Makie.wong_colors()[inst_idx],
            markersize=4,
            strokecolor=:black,strokewidth=0.1,
        )

        Makie.scatter!(
            ax_resid,
            data.epoch,
            resid,
            color=Makie.wong_colors()[inst_idx],
            markersize=4,
            strokecolor=:black,strokewidth=0.1,
        )
        phase_folded = mod.(data.epoch .- t_peri .- T/4, T)./T .- 0.5
        Makie.scatter!(
            ax_phase,
            phase_folded,
            data_minus_off_and_gp[thisinst_mask],
            color=Makie.wong_colors()[inst_idx],
            markersize=4,
            strokecolor=:black,strokewidth=0.1,
        )
    end


    Makie.xlims!(ax_resid, extrema(ts))


    # Binned values on phase folded plot
    # Noise weighted (including jitter and GP)
    bins = -0.45:0.1:0.45
    binned = zeros(length(bins))
    binned_unc = zeros(length(bins))
    phase_folded = mod.(rvs.table.epoch[:] .- t_peri .- T/4, T)./T .- 0.5
    jitters = map(eachrow(rvs.table)) do row
        inst_idx = row[].inst_idx
        results["jitter_$inst_idx"][sample_idx]
    end
    for (i,bin_cent) in enumerate(bins)
        mask = bin_cent - step(bins)/2 .<= phase_folded .<= bin_cent + step(bins/2)
        if count(mask) == 0
            binned[i] = NaN
            continue
        end
        binned[i] = mean(
            data_minus_off_and_gp[mask],
            ProbabilityWeights(1 ./ errs_all[mask].^2)
        )
        binned_unc[i] = sem(
            data_minus_off_and_gp[mask],
            ProbabilityWeights(1 ./ errs_all[mask].^2)
        )
    end
    errorbars!(
        ax_phase,
        bins,
        binned,
        binned_unc,
        color=:black,
        linewidth=1,
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
    Legend(
        gs[1:2,2],
        [
          MarkerElement(color = Makie.wong_colors()[i], marker=:circle, markersize = 15)
          for i in 1:length(rvs.instrument_names)
        ],
        rvs.instrument_names,
        "instrument",
        valign=:top,
        halign=:left,
        width=Relative(1),
        # height=Relative(1),
    )
    Legend(
        gs[3,2],
        [
            LineElement(color = :blue,linewidth=4,),
            MarkerElement(color = :red, strokecolor=:black, strokewidth=2, marker=:circle, markersize = 15),
        ],
        [
            Makie.rich("maximum\n", Makie.rich("a posteriori", font=:italic), "\nmodel"),
            "binned",
        ],
        valign=:top,
        halign=:left,
        # width=Relative(1),
        # height=Relative(1),
    )
    Makie.rowsize!(gs, 1, Auto(2))
    Makie.rowsize!(gs, 2, Auto(1))
    Makie.rowsize!(gs, 3, Auto(2))

end


end