
##################################################
# HGCA Plot
const pmra_label = rich("μ", subscript("α*"), " [mas/yr]")
const pmdec_label = rich("μ", subscript("δ"), " [mas/yr]")
function Octofitter.hgcaplot(
    model,
    results,
    fname="$(model.system.name)-hgcaplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700, 600),
        figure...
    )
    Octofitter.hgcaplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function Octofitter.hgcaplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    # If less than 500 samples, just show all of them
    # N=min(size(results, 1) * size(results, 3), 500),
    ts,
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii=(
        N == size(results, 1) * size(results, 3) ?
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3), N)
    ),
    axis=(;),
    colormap=:plasma,
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    alpha=min.(1, 100 / length(ii)),
    kwargs...
)
    gs = gridspec_or_fig

    date_pos, date_strs, xminorticks = _date_ticks(ts)
    ax_velra = Axis(
        gs[1, 1:3];
        ylabel=pmra_label,
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=top_time_axis,
        xticklabelsvisible=top_time_axis,
        xminorticks,
        xminorticksvisible=top_time_axis,
        axis...
    )
    ax_veldec = Axis(
        gs[2, 1:3];
        xlabel="MJD",
        ylabel=pmdec_label,
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    linkxaxes!(ax_velra, ax_veldec)

    ax_dat1 = Axis(
        gs[3, 1],
        xlabel=pmra_label,
        ylabel=pmdec_label,
        autolimitaspect=1.0,
        title="H",
        # titlecolor=Makie.wong_colors()[1],
        xticklabelrotation=pi / 4,
        xgridvisible=false,
        ygridvisible=false,
    )
    ax_dat2 = Axis(
        gs[3, 2],
        xlabel=pmra_label,
        ylabel=pmdec_label,
        autolimitaspect=1.0,
        title="G-H",
        # titlecolor=Makie.wong_colors()[2],
        xticklabelrotation=pi / 4,
        xgridvisible=false,
        ygridvisible=false,
        ylabelvisible=false,
    )
    ax_dat3 = Axis(
        gs[3, 3],
        xlabel=pmra_label,
        ylabel=pmdec_label,
        autolimitaspect=1.0,
        title="G",
        # titlecolor=Makie.wong_colors()[3],
        xticklabelrotation=pi / 4,
        xgridvisible=false,
        ygridvisible=false,
        ylabelvisible=false,
    )
    # linkxaxes!(ax_dat1,ax_dat2,ax_dat3)
    # linkyaxes!(ax_dat1,ax_dat2,ax_dat3)


    xlims!(ax_velra, extrema(ts))
    xlims!(ax_veldec, extrema(ts))

    pmra_model_t = zeros(length(ii), length(ts))
    pmdec_model_t = zeros(length(ii), length(ts))
    color_model_t = zeros(length(ii), length(ts))
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        # Draws from the posterior
        mass = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol

        # Now time-series
        sols = orbitsolve.(orbs, ts')
        # Can we use the existing simulator for this please?
        pmra_model_t .+= pmra.(sols, mass)
        pmdec_model_t .+= pmdec.(sols, mass)
        color_model_t .= rem2pi.(
            meananom.(sols), RoundDown) .+ 0 .* ii
    end
    if haskey(results, :pmra)
        pmra_model_t .+= results[:pmra][ii]
    end
    if haskey(results, :pmdec)
        pmdec_model_t .+= results[:pmdec][ii]
    end
    if colorbar
        Colorbar(
            gs[1:2, 4];
            colormap,
            label="mean anomaly",
            colorrange=(0,2pi),
            ticks=(
                [0,pi/2,pi,3pi/2,2pi],
                ["0", "π/2", "π", "3π/2", "2π"]
            )
        )
    end
    lines!(ax_velra,
        concat_with_nan(ts' .+ 0 .* pmra_model_t),
        concat_with_nan(pmra_model_t);
        alpha,
        color=concat_with_nan(color_model_t),
        colorrange=(0, 2pi),
        colormap
    )
    lines!(ax_veldec,
        concat_with_nan(ts' .+ 0 .* pmdec_model_t),
        concat_with_nan(pmdec_model_t);
        alpha,
        color=concat_with_nan(color_model_t),
        colorrange=(0, 2pi),
        colormap
    )

    # Now over plot any astrometry
    like_objs = filter(model.system.observations) do like_obj
        nameof(typeof(like_obj)) == :HGCALikelihood
    end
    if isempty(like_objs)
        return
    end
    hgca_like = only(like_objs)

    tx = [
        hgca_like.table.epoch_ra_hip_mjd
        (hgca_like.table.epoch_ra_hip_mjd + hgca_like.table.epoch_ra_gaia_mjd) / 2
        hgca_like.table.epoch_ra_gaia_mjd
    ]
    ty = [
        hgca_like.table.epoch_dec_hip_mjd
        (hgca_like.table.epoch_dec_hip_mjd + hgca_like.table.epoch_dec_gaia_mjd) / 2
        hgca_like.table.epoch_dec_gaia_mjd
    ]
    x = [
        hgca_like.table.pmra_hip
        hgca_like.table.pmra_hg
        hgca_like.table.pmra_gaia
    ]
    y = [
        hgca_like.table.pmdec_hip
        hgca_like.table.pmdec_hg
        hgca_like.table.pmdec_gaia
    ]

    cor = [
        hgca_like.table.pmra_pmdec_hip
        hgca_like.table.pmra_pmdec_hg
        hgca_like.table.pmra_pmdec_gaia
    ]

    σ₁ = [
        hgca_like.table.pmra_hip_error
        hgca_like.table.pmra_hg_error
        hgca_like.table.pmra_gaia_error
    ]
    σ₂ = [
        hgca_like.table.pmdec_hip_error
        hgca_like.table.pmdec_hg_error
        hgca_like.table.pmdec_gaia_error
    ]

    # # One more model-plot: add scatter points to existing lines at the 
    # # data epochs, colored correctly
    # scatx = vec(stack(map(tx) do tx_i
    #     pmra_model_t[:, argmin(abs.(ts .- tx_i))]
    # end))
    # scaty = vec(stack(map(ty) do ty_i
    #     pmdec_model_t[:, argmin(abs.(ts .- ty_i))]
    # end))
    # scatter!(
    #     ax_vel2d,
    #     scatx,
    #     scaty,
    #     color=repeat(Makie.wong_colors()[1:length(x)], outer=length(ii)),
    #     markersize=5,
    # )

    error_ellipses = broadcast(x, y, σ₁, cor, σ₂) do x, y, σ₁, cor, σ₂
        Σ = [
            σ₁^2 cor*σ₁*σ₂
            cor*σ₁*σ₂ σ₂^2
        ]
        vals, vecs = eigen(Σ) # should be real and sorted by real eigenvalue
        length_major = sqrt(vals[2])
        length_minor = sqrt(vals[1])
        λ = vecs[:, 2]
        α = atan(λ[2], λ[1])

        xvals = [
            # Major axis
            x - length_major * cos(α),
            x + length_major * cos(α),
            NaN,
            # Minor axis
            x - length_minor * cos(α + π / 2),
            x + length_minor * cos(α + π / 2),
            NaN,
        ]
        yvals = [
            # Major axis
            y - length_major * sin(α),
            y + length_major * sin(α),
            NaN,
            # Minor axis
            y - length_minor * sin(α + π / 2),
            y + length_minor * sin(α + π / 2),
            NaN,
        ]
        xvals, yvals
    end

    

    colsize!(gs, 1, Auto(1 // 3))
    colsize!(gs, 2, Auto(1 // 3))
    colsize!(gs, 3, Auto(1 // 3))
    rowsize!(gs, 3, Aspect(3, 1.0))

    ## Model
    ## Compute these for all results, not just `ii`
    θ_systems_from_chain = Octofitter.mcmcchain2result(model, results)
    # Display all points, unless there are more than 10k 
    jj = 1:size(results,1)*size(results,3)
    if size(results,1)*size(results,3) > 5_000
        jj = ii
    end
    for (θ_system, i) in zip(θ_systems_from_chain, jj)
        orbits = map(keys(model.system.planets)) do planet_key
            Octofitter.construct_elements(results, planet_key, i)
        end
        hgca_like_sim = Octofitter.generate_from_params(hgca_like, θ_system, orbits)
        # HIP Epoch
        Makie.scatter!(
            ax_dat1,
            [only(hgca_like_sim.table.pmra_hip)],
            [only(hgca_like_sim.table.pmdec_hip)],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [(only(hgca_like_sim.table.epoch_ra_hip_mjd))],
            [only(hgca_like_sim.table.pmra_hip)],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [(only(hgca_like_sim.table.epoch_dec_hip_mjd))],
            [only(hgca_like_sim.table.pmdec_hip)],
            color=:black,
            markersize=3,
        )
        # HG Epoch
        Makie.scatter!(
            ax_dat2,
            [only(hgca_like_sim.table.pmra_hg)],
            [only(hgca_like_sim.table.pmdec_hg)],
            color=:black,
            markersize=2,
        )
        hg_ra_epoch = (only(hgca_like_sim.table.epoch_ra_hip_mjd) +
                       only(hgca_like_sim.table.epoch_ra_gaia_mjd)) / 2
        hg_dec_epoch = (only(hgca_like_sim.table.epoch_dec_hip_mjd) +
                        only(hgca_like_sim.table.epoch_dec_gaia_mjd)) / 2
        Makie.scatter!(
            ax_velra,
            [(hg_ra_epoch)],
            [only(hgca_like_sim.table.pmra_hg)],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [(hg_dec_epoch)],
            [only(hgca_like_sim.table.pmdec_hg)],
            color=:black,
            markersize=3,
        )
        # GAIA Epoch
        Makie.scatter!(
            ax_dat3,
            [only(hgca_like_sim.table.pmra_gaia)],
            [only(hgca_like_sim.table.pmdec_gaia)],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [(only(hgca_like_sim.table.epoch_ra_gaia_mjd))],
            [only(hgca_like_sim.table.pmra_gaia)],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [(only(hgca_like_sim.table.epoch_dec_gaia_mjd))],
            [only(hgca_like_sim.table.pmdec_gaia)],
            color=:black,
            markersize=3,
        )
    end

    ## Data 
    Makie.lines!(
        ax_dat1,
        error_ellipses[1][1],
        error_ellipses[1][2],
        color=Makie.wong_colors()[1],
        linewidth=3.5,
    )
    Makie.scatter!(
        ax_dat1, x[1], y[1],
        color=Makie.wong_colors()[1],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
    Makie.lines!(
        ax_dat2,
        error_ellipses[2][1],
        error_ellipses[2][2],
        color=Makie.wong_colors()[2],
        linewidth=3.5,
    )
    Makie.scatter!(
        ax_dat2, x[2], y[2],
        color=Makie.wong_colors()[2],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
    Makie.lines!(
        ax_dat3,
        error_ellipses[3][1],
        error_ellipses[3][2],
        color=Makie.wong_colors()[3],
        linewidth=3.5,
    )
    Makie.scatter!(
        ax_dat3, x[3], y[3],
        color=Makie.wong_colors()[3],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )


    # 1D plots: stroke twice for contrast
    Makie.errorbars!(
        ax_velra,
        tx,
        x,
        σ₁,
        color=Makie.wong_colors()[[1,2,3]],
    )
    Makie.errorbars!(
        ax_veldec,
        ty,
        y,
        σ₂,
        color=Makie.wong_colors()[[1,2,3]],
    )
    Makie.scatter!(
        ax_velra,
        tx,
        x,
        σ₁,
        color=Makie.wong_colors()[1:length(x)],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
    Makie.scatter!(
        ax_veldec,
        ty,
        y,
        σ₂,
        color=Makie.wong_colors()[1:length(x)],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )

    return [ax_velra, ax_veldec]
end