


##################################################
# RV vs time plot
# Note that the "full" plot requires to use
# rvpostplot(). It shows the jitter, residuals,
# and phase folded curves properly for a single sample.
# This just gives a global view.


function rvtimeplot(
    model,
    results,
    fname="$(model.system.name)-rvplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700,600),
        figure...
    )
    rvtimeplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function rvtimeplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    # If less than 500 samples, just show all of them
    N  = min(size(results, 1)*size(results, 3), 500),
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii = (
        N == size(results, 1)*size(results, 3) ? 
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3),N)
    ),
    axis=(;),
    colormap=:plasma,
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    kwargs...
)
    gs = gridspec_or_fig


    # Find the minimum time step size
    periods = Float64[]
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, :)
        append!(periods, period.(orbs))
    end
    min_period = minimum(periods)
    med_period = median(periods)

    # Always try and show at least one full period
    # TODO: decide t_start and t_end based on epochs of data if astrometry 
    # available, otherwise other data, otherwise arbitrary
    t_start = years2mjd(1990)
    t_stop = years2mjd(2020)
    t_cur = t_stop-t_start
    if t_cur < med_period
        t_extend = med_period - t_cur
        t_start -= t_extend/2
        t_stop += t_extend/2
    end

    ts = range(t_start, t_stop, step=min_period / 150)
    if length(ts) > 1500
        # Things can get very slow if we have some very short period 
        # solutions
        ts = range(t_start, t_stop, length=1500)
    end
    # Make sure we include the approx. data epochs
    ts = sort([ts; years2mjd(1990.25); years2mjd((1990.25 + 2016) / 2); years2mjd(2016)])
    
    date_pos, date_strs = _date_ticks(ts)
    ax= Axis(
        gs[1, 1];
        ylabel="rv [m/s]",
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=top_time_axis,
        xticklabelsvisible=top_time_axis,
        axis...
    )
    ax_secondary = Axis(
        gs[1, 1];
        xlabel="MJD",
        xgridvisible=false,
        ygridvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        ylabelvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    linkxaxes!(ax, ax_secondary)
    xlims!(ax, extrema(ts))
    xlims!(ax_secondary, extrema(ts))



    rv_model_t = zeros(length(ii), length(ts))
    color_model_t = zeros(length(ii), length(ts))
    for (i_planet, planet_key) in enumerate(keys(model.system.planets))
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        # Draws from the posterior
        mass = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol

        # Now time-series
        sols = orbitsolve.(orbs, ts')
  
        rv_model_t .+= radvel.(sols, mass) 

              
        # We subtract off the "average" instrument RV offset here.
        # This varies by instrument, and isn't centred around a net hierarchical value.
        # For now just take the average of the different zero points per row.
        if i_planet == 1 && haskey(results, :rv0_1)
            rv0s = []
            for i in 1:10
                k = Symbol("rv0_$i")
                if haskey(results, k)
                    push!(rv0s, results[k][ii])
                end
            end
            ave_rv0 = mean(stack(rv0s),dims=2)
            rv_model_t .-= ave_rv0
        end
            

        color_model_t .+= rem2pi.(
            meananom.(sols), RoundDown) .+ 0 .* ii
    end

    lines!(ax,
        concat_with_nan(ts' .+ 0 .* rv_model_t),
        concat_with_nan(rv_model_t);
        alpha=min.(1, 100 / length(ii)),
        color=concat_with_nan(color_model_t),
        colorrange=(0,2pi),
        colormap
    )



    # # Now overplot the data points, if any.
    # for like_obj in model.system.observations
    #     if nameof(typeof(like_obj)) != :StarAbsoluteRVLikelihood
    #         continue
    #     end
    #     epoch = vec(like_obj.table.epoch)
    #     rv = vec(like_obj.table.rv)
    #     σ_rv = vec(like_obj.table.σ_rv)
    #     Makie.errorbars!(
    #         ax, epoch, rv, σ_rv;
    #         color=:black,
    #         linewidth=3,
    #     )
    #     Makie.scatter!(
    #         ax, epoch, rv;
    #         color=:white,
    #         strokewidth=0.5,
    #         strokecolor=:black,
    #         markersize=2,
    #     )
    # end

    if colorbar
        Colorbar(
            gs[1,2];
            colormap,
            label="orbit fraction past periastron",
            colorrange=(0,1)
        )
    end
end
