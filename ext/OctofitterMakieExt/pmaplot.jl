
##################################################
# HGCA Plot
function pmaplot(
    model,
    results,
    fname="$(model.system.name)-pmaplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700, 600),
        figure...
    )
    pmaplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function pmaplot!(
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
    gs_row = 0
    ax_velra = Axis(
        gs[gs_row += 1, 1:3];
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
        gs[gs_row += 1, 1:3];
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

    # linkxaxes!(ax_dat1,ax_dat2,ax_dat3)
    # linkyaxes!(ax_dat1,ax_dat2,ax_dat3)


    xlims!(ax_velra, extrema(ts))
    xlims!(ax_veldec, extrema(ts))

    pmra_model_t = zeros(length(ii), length(ts))
    pmdec_model_t = zeros(length(ii), length(ts))
    color_model_t = zeros(length(ii), length(ts))
    absolute_orbits = false
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)
        absolute_orbits |= first(orbs) isa AbsoluteVisual
        # Draws from the posterior
        mass = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol

        # Now time-series
        sols = orbitsolve.(orbs, ts')
        try
            pmra(first(sols), first(mass))
        catch
            continue
        end
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
        colormap,
        rasterize=4,
    )
    lines!(ax_veldec,
        concat_with_nan(ts' .+ 0 .* pmdec_model_t),
        concat_with_nan(pmdec_model_t);
        alpha,
        color=concat_with_nan(color_model_t),
        colorrange=(0, 2pi),
        colormap,
        rasterize=4,
    )




    cor_axes = Axis[]



    #########################################
    # HGCA

    # Now over plot any astrometry
    like_objs = filter(model.system.observations) do like_obj
        nameof(typeof(like_obj)) == :HGCAObs ||
        nameof(typeof(like_obj)) == :HGCAInstantaneousObs ||
    end
    if !isempty(like_objs)
        hgca_like = only(like_objs)



        ax_dat1 = Axis(
            gs[gs_row+=1, 1],
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
            gs[gs_row, 2],
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
            gs[gs_row, 3],
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

        # The HGCA catalog values have an non-linearity correction added.
        # If we are doing our own rigorous propagation we don't need this
        # correction. We could subtract it from the measurements, but 
        # here we just add it to our model so that they match
        if absolute_orbits
            hg_nonlinear_dpmra = hgca_like.hgca.nonlinear_dpmra[1]
            hg_nonlinear_dpmdec = hgca_like.hgca.nonlinear_dpmdec[1]
            hip_nonlinear_dpmra = 2hgca_like.hgca.nonlinear_dpmra[1]
            hip_nonlinear_dpmdec = 2hgca_like.hgca.nonlinear_dpmdec[1]
        else
            hg_nonlinear_dpmra = 
            hg_nonlinear_dpmdec = 
            hip_nonlinear_dpmra = 
            hip_nonlinear_dpmdec = zero(hgca_like.hgca.nonlinear_dpmra[1])
        end


        tx = [
            hgca_like.hgca.epoch_ra_hip_mjd
            (hgca_like.hgca.epoch_ra_hip_mjd + hgca_like.hgca.epoch_ra_gaia_mjd) / 2
            hgca_like.hgca.epoch_ra_gaia_mjd
        ]
        ty = [
            hgca_like.hgca.epoch_dec_hip_mjd
            (hgca_like.hgca.epoch_dec_hip_mjd + hgca_like.hgca.epoch_dec_gaia_mjd) / 2
            hgca_like.hgca.epoch_dec_gaia_mjd
        ]
        x = [
            hgca_like.hgca.pmra_hip -  hip_nonlinear_dpmra
            hgca_like.hgca.pmra_hg - hg_nonlinear_dpmra
            hgca_like.hgca.pmra_gaia
        ]
        y = [
            hgca_like.hgca.pmdec_hip -  hip_nonlinear_dpmdec
            hgca_like.hgca.pmdec_hg -  hg_nonlinear_dpmdec
            hgca_like.hgca.pmdec_gaia
        ]

        cor = [
            hgca_like.hgca.pmra_pmdec_hip
            hgca_like.hgca.pmra_pmdec_hg
            hgca_like.hgca.pmra_pmdec_gaia
        ]

        σ₁ = [
            hgca_like.hgca.pmra_hip_error
            hgca_like.hgca.pmra_hg_error
            hgca_like.hgca.pmra_gaia_error
        ]
        σ₂ = [
            hgca_like.hgca.pmdec_hip_error
            hgca_like.hgca.pmdec_hg_error
            hgca_like.hgca.pmdec_gaia_error
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



        ## Model
        ## Compute these for all results, not just `ii`
        θ_systems_from_chain = Octofitter.mcmcchain2result(model, results)
        # Display all points, unless there are more than 10k 
        jj = 1:size(results,1)*size(results,3)
        if size(results,1)*size(results,3) > 5_000
            jj = ii
        end
        sims = []
        for (θ_system, i) in zip(θ_systems_from_chain, jj)
            orbits = map(keys(model.system.planets)) do planet_key
                Octofitter.construct_elements(model, results, planet_key, i)
            end
            θ_obs = (;)
            name = Octofitter.normalizename(likelihoodname(hgca_like))
            if hasproperty(θ_system, :observations) && hasproperty(θ_system.observations, name)
                θ_obs = θ_system.observations[name]
            end
            if hasproperty(hgca_like, :table)
                solutions = map(orbits) do orbit
                    return orbitsolve.(orbit, hgca_like.table.epoch)
                end
                sim = Octofitter.simulate(hgca_like, θ_system, θ_obs, orbits, solutions, 0)
            else
                solutions = [() for _ in length(model.system.planets)]
                sim = Octofitter.simulate(hgca_like, θ_system, θ_obs, orbits, solutions, -1)
            end
            push!(sims, sim)
        end
        # HIP Epoch
        Makie.scatter!(
            ax_dat1,
            [only(sim.pmra_hip_model) for sim in sims ],
            [only(sim.pmdec_hip_model) for sim in sims],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [only(hgca_like.hgca.epoch_ra_hip_mjd) for sim in sims],
            [only(sim.pmra_hip_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [only(hgca_like.hgca.epoch_dec_hip_mjd) for sim in sims],
            [only(sim.pmdec_hip_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        # HG Epoch
        Makie.scatter!(
            ax_dat2,
            [only(sim.pmra_hg_model) for sim in sims],
            [only(sim.pmdec_hg_model) for sim in sims],
            color=:black,
            markersize=2,
        )
        hg_ra_epoch = [(only(hgca_like.hgca.epoch_ra_hip_mjd) +
                        only(hgca_like.hgca.epoch_ra_gaia_mjd)) / 2 for sim in sims]
        hg_dec_epoch = [(only(hgca_like.hgca.epoch_dec_hip_mjd) +
                        only(hgca_like.hgca.epoch_dec_gaia_mjd)) / 2 for sim in sims]
        Makie.scatter!(
            ax_velra,
            hg_ra_epoch,
            [only(sim.pmra_hg_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            hg_dec_epoch,
            [only(sim.pmdec_hg_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        # GAIA Epoch
        Makie.scatter!(
            ax_dat3,
            [only(sim.pmra_gaia_model) for sim in sims],
            [only(sim.pmdec_gaia_model) for sim in sims],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [only(hgca_like.hgca.epoch_ra_gaia_mjd) for sim in sims],
            [only(sim.pmra_gaia_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [only(hgca_like.hgca.epoch_dec_gaia_mjd) for sim in sims],
            [only(sim.pmdec_gaia_model) for sim in sims],
            color=:black,
            markersize=3,
        )

        ## Data 
        Makie.lines!(
            ax_dat1,
            error_ellipses[1][1],
            error_ellipses[1][2],
            color=Makie.wong_colors()[1],
            linewidth=3.5
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
            strokecolor=Makie.wong_colors()[1:length(x)],
            color=:transparent,
            markersize=10,
            strokewidth=1.5,
        )
        Makie.scatter!(
            ax_veldec,
            ty,
            y,
            σ₂,
            strokecolor=Makie.wong_colors()[1:length(x)],
            color=:transparent,
            markersize=10,
            strokewidth=1.5,
        )


        append!(cor_axes, [ax_dat1, ax_dat2, ax_dat3])
    end





    #######################
    # G DR3-2 difference

    # Now over plot any astrometry
    like_objs = filter(model.system.observations) do like_obj
        nameof(typeof(like_obj)) == :GaiaDifferenceLike
    end
    if !isempty(like_objs)


        # Remove G if redundent 
        if !isempty(cor_axes)
            delete!(cor_axes[end])
        end
            

        ax_dat1 = Axis(
            gs[gs_row+=1, 1],
            xlabel=pmra_label,
            ylabel=pmdec_label,
            autolimitaspect=1.0,
            title="DR2",
            # titlecolor=Makie.wong_colors()[4],
            xticklabelrotation=pi / 4,
            xgridvisible=false,
            ygridvisible=false,
        )
        ax_dat2 = Axis(
            gs[gs_row, 2],
            xlabel=pmra_label,
            ylabel=pmdec_label,
            autolimitaspect=1.0,
            title="DR 2-3",
            # titlecolor=Makie.wong_colors()[5],
            xticklabelrotation=pi / 4,
            xgridvisible=false,
            ygridvisible=false,
            ylabelvisible=false,
        )
        ax_dat3 = Axis(
            gs[gs_row, 3],
            xlabel=pmra_label,
            ylabel=pmdec_label,
            autolimitaspect=1.0,
            title="DR3",
            # titlecolor=Makie.wong_colors()[6],
            xticklabelrotation=pi / 4,
            xgridvisible=false,
            ygridvisible=false,
            ylabelvisible=false,
        )

        gaialike = only(like_objs)

        ## Model
        ## Compute these for all results, not just `ii`
        θ_systems_from_chain = Octofitter.mcmcchain2result(model, results)
        # Display all points, unless there are more than 10k 
        jj = 1:size(results,1)*size(results,3)
        if size(results,1)*size(results,3) > 5_000
            jj = ii
        end
        sims = []
        name = Octofitter.normalizename(likelihoodname(gaialike))
        for (θ_system, i) in zip(θ_systems_from_chain, jj)
            orbits = map(keys(model.system.planets)) do planet_key
                Octofitter.construct_elements(model, results, planet_key, i)
            end
            solutions = map(orbits) do orbit
                return orbitsolve.(orbit, gaialike.table.epoch)
            end

            θ_obs = (;)
            if hasproperty(θ_system, :observations) && hasproperty(θ_system.observations, name)
                θ_obs = θ_system.observations[name]
            end

            sim = Octofitter.simulate(gaialike, θ_system,θ_obs, orbits, solutions, 0)
            push!(sims, sim[2])
        end


        # Add these to the catalog DR2
        if hasproperty(θ_systems_from_chain[1], :dr2_systematic_Δra) &&
            hasproperty(θ_systems_from_chain[1], :dr2_systematic_Δdec)
            dr2_systematic_Δra = median([θnt.dr2_systematic_Δra for θnt in θ_systems_from_chain])
            dr2_systematic_Δdec = median([θnt.dr2_systematic_Δdec for θnt in θ_systems_from_chain])
        else
            dr2_systematic_Δra = 0.0
            dr2_systematic_Δdec = 0.0
        end
        dr2_systematic_Δμ_ra = median([θnt.dr2_systematic_Δμ_ra for θnt in θ_systems_from_chain])
        dr2_systematic_Δμ_dec = median([θnt.dr2_systematic_Δμ_dec for θnt in θ_systems_from_chain])
    
        σ_dr2_systematic_Δμ_ra = std([θnt.dr2_systematic_Δμ_ra for θnt in θ_systems_from_chain])
        σ_dr2_systematic_Δμ_dec = std([θnt.dr2_systematic_Δμ_dec for θnt in θ_systems_from_chain])
        if !isfinite(σ_dr2_systematic_Δμ_ra)
            σ_dr2_systematic_Δμ_ra = 0
        end
        if !isfinite(σ_dr2_systematic_Δμ_dec)
            σ_dr2_systematic_Δμ_dec = 0
        end

        tx = [
            Octofitter.meta_gaia_DR2.ref_epoch_mjd
            (Octofitter.meta_gaia_DR2.ref_epoch_mjd+Octofitter.meta_gaia_DR3.ref_epoch_mjd)/2
            Octofitter.meta_gaia_DR3.ref_epoch_mjd
        ]
        ty = [
            Octofitter.meta_gaia_DR2.ref_epoch_mjd
            (Octofitter.meta_gaia_DR2.ref_epoch_mjd+Octofitter.meta_gaia_DR3.ref_epoch_mjd)/2
            Octofitter.meta_gaia_DR3.ref_epoch_mjd
        ]
        # TODO: have to add median calibration offset
        Δt = (Octofitter.meta_gaia_DR3.ref_epoch_mjd - Octofitter.meta_gaia_DR2.ref_epoch_mjd)/Octofitter.julian_year
        x = [
            gaialike.dr2.pmra + dr2_systematic_Δμ_ra
            ((gaialike.dr3.ra-gaialike.dr2.ra)*60*60*1000*cosd((gaialike.dr3.dec+gaialike.dr2.dec)/2) - dr2_systematic_Δra)/Δt
            gaialike.dr3.pmra
        ]
        y = [
            gaialike.dr2.pmdec + dr2_systematic_Δμ_dec
            ((gaialike.dr3.dec-gaialike.dr2.dec)*60*60*1000 - dr2_systematic_Δdec)/Δt
            gaialike.dr3.pmdec
        ]

        # Calculate uncertainties in the differences
        σ_Δra = sqrt(gaialike.dr3.ra_error^2 + gaialike.dr2.ra_error^2)
        σ_Δdec = sqrt(gaialike.dr3.dec_error^2 + gaialike.dr2.dec_error^2)
        
        # Calculate covariances for each DR
        cov_dr3 = gaialike.dr3.ra_dec_corr * gaialike.dr3.ra_error * gaialike.dr3.dec_error
        cov_dr2 = gaialike.dr2.ra_dec_corr * gaialike.dr2.ra_error * gaialike.dr2.dec_error
        
        # Calculate total covariance for the differences
        cov_diff = cov_dr3 + cov_dr2
        
        # Calculate correlation coefficient for the differences
        ρ_diff = cov_diff / (σ_Δra * σ_Δdec)

        cor = [
            gaialike.dr2.pmra_pmdec_corr
            ρ_diff
            gaialike.dr3.pmra_pmdec_corr
        ]
        σ₁ = [
            hypot(gaialike.dr2.pmra_error, σ_dr2_systematic_Δμ_ra)
            σ_Δra
            gaialike.dr3.pmra_error
        ]
        σ₂ = [
            hypot(gaialike.dr2.pmdec_error, σ_dr2_systematic_Δμ_dec)
            σ_Δdec
            gaialike.dr3.pmdec_error
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
        # colsize!(gs, 2, Auto(1 // 3))
        colsize!(gs, 3, Auto(1 // 3))
        rowsize!(gs, gs_row, Aspect(3, 1.0))
    
    
        # HIP Epoch
        Makie.scatter!(
            ax_dat1,
            [only(sim.pmra_dr2_model) for sim in sims],
            [only(sim.pmdec_dr2_model) for sim in sims],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [tx[1] for _ in sims],
            [only(sim.pmra_dr2_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [tx[1] for _ in sims],
            [only(sim.pmdec_dr2_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        # # HG Epoch
        Makie.scatter!(
            ax_dat2,
            [only(sim.pmra_dr32_model) for sim in sims],
            [only(sim.pmdec_dr32_model) for sim in sims],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [tx[2] for _ in sims],
            [only(sim.pmra_dr32_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [tx[2] for _ in sims],
            [only(sim.pmdec_dr32_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        # GAIA Epoch
        Makie.scatter!(
            ax_dat3,
            [only(sim.pmra_dr3_model) for sim in sims],
            [only(sim.pmdec_dr3_model) for sim in sims],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [tx[3] for _ in sims],
            [only(sim.pmra_dr3_model) for sim in sims],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [tx[3] for _ in sims],
            [only(sim.pmdec_dr3_model) for sim in sims],
            color=:black,
            markersize=3,
        )

        ## Data 
        Makie.lines!(
            ax_dat1,
            error_ellipses[1][1],
            error_ellipses[1][2],
            color=Makie.wong_colors()[4],
            linewidth=3.5
        )
        Makie.scatter!(
            ax_dat1, x[1], y[1],
            color=Makie.wong_colors()[4],
            markersize=10,
            strokewidth=1.5,
            strokecolor=:black
        )
        Makie.lines!(
            ax_dat2,
            error_ellipses[2][1],
            error_ellipses[2][2],
            color=Makie.wong_colors()[5],
            linewidth=3.5,
        )
        Makie.scatter!(
            ax_dat2, x[2], y[2],
            color=Makie.wong_colors()[5],
            markersize=10,
            strokewidth=1.5,
            strokecolor=:black
        )
        Makie.lines!(
            ax_dat3,
            error_ellipses[3][1],
            error_ellipses[3][2],
            color=Makie.wong_colors()[6],
            linewidth=3.5,
        )
        Makie.scatter!(
            ax_dat3, x[3], y[3],
            color=Makie.wong_colors()[6],
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
            color=Makie.wong_colors()[[4,5,6]],
        )
        Makie.errorbars!(
            ax_veldec,
            ty,
            y,
            σ₂,
            color=Makie.wong_colors()[[4,5,6]],
        )
        Makie.scatter!(
            ax_velra,
            tx,
            x,
            σ₁,
            strokecolor=Makie.wong_colors()[[4,5,6]],
            color=:transparent,
            markersize=10,
            strokewidth=1.5,
        )
        Makie.scatter!(
            ax_veldec,
            ty,
            y,
            σ₂,
            strokecolor=Makie.wong_colors()[[4,5,6]],
            color=:transparent,
            markersize=10,
            strokewidth=1.5,
        )


        for ax in cor_axes 
            ax.xlabelvisible = false
        end
        append!(cor_axes, [ax_dat1, ax_dat3])
    end


    if !isempty(cor_axes)
        xspace = maximum(Makie.tight_xticklabel_spacing!, cor_axes)
        for ax in cor_axes
            ax.xticklabelspace = xspace + 20
        end


        if gs.size[1] == 3
            colsize!(gs, 1, Auto(1 // 3))
            colsize!(gs, 2, Auto(1 // 3))
            colsize!(gs, 3, Auto(1 // 3))
        else
            colsize!(gs, 1, Auto(1 // 2))
            colsize!(gs, 2, Auto(1 // 2))
        end
        rowsize!(gs, gs_row, Aspect(3, 1.0))
    end

    return [ax_velra, ax_veldec]
end