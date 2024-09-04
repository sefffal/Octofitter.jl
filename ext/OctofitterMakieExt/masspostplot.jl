


##################################################
# Mass vs. semi-major axis plot
function Octofitter.masspostplot(
    model,
    results,
    fname="$(model.system.name)-masspostplot!.png",
    args...;
    figure=(;),
    # massloglog=nothing,
    kwargs...
)

    # # auto-determine if we should use a log-log plot
    # if isnothing(massloglog)

    # end

    fig = Figure(;
        size=(500,300),
        figure...,
    )
    Octofitter.masspostplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end

const mass_mjup_label = rich("mass [M", subscript("jup"), "]")
function Octofitter.masspostplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    axis=(;),
    kwargs...
)
    gs = gridspec_or_fig

    ax_hist = Axis(gs[1:2,1];
        # ylabel=rich("mass [M", subscript("jup"), "]"),
        xlabel=mass_mjup_label,
        xgridvisible=false,
        ygridvisible=false,
    )
    ylims!(ax_hist, low=0)
    ax_scat_sma = Axis(gs[1,2];
        ylabel=mass_mjup_label,
        xlabel="sma [AU]",
        # xscale=log10,
        # xticks=2 .^ (0:12),
        # yticks=2 .^ (0:12),
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )
    ax_scat_ecc = Axis(gs[2,2];
        ylabel=mass_mjup_label,
        xlabel="eccentricity",
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )
    xlims!(ax_scat_ecc, 0, 1)
    cred_intervals = []
    for planet_key in keys(model.system.planets)
        mk = Symbol("$(planet_key)_mass")
        if !haskey(results, mk)
            continue
        end
        sma = vec(results["$(planet_key)_a"])
        mass = vec(results[mk])
        ecc = vec(results["$(planet_key)_e"])
        stephist!(
            ax_hist,
            mass,
            linewidth=3
        )
        scatter!(ax_scat_sma, sma, mass;
            markersize=2,
            rasterize=4,
        )
        scatter!(ax_scat_ecc, ecc, mass;
            markersize=2,
            rasterize=4,
        )
        low,mid,high = quantile(mass, (0.16, 0.5, 0.84))
        label = margin_confidence_default_formatter(mid-low,mid,high-mid)
        push!(
            cred_intervals,
            Makie.latexstring("$(planet_key)_{mass} = "*label)
        )
    end
    if !isempty(cred_intervals)
        Label(
            gs[0,1],
            reduce(*, cred_intervals),
            tellwidth=false
        )
    end


end
