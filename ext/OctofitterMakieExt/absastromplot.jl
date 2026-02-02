
##################################################
# HGCA Plot
function absastromplot(
    model,
    results,
    fname="$(model.system.name)-absastromplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700, 600),
        figure...
    )
    absastromplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function absastromplot!(
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

    # Add in a non-linear effect due to light travel time changes
    els = Octofitter.construct_elements(model, results, first(keys(model.system.planets)), ii)
    for (j, orb) in enumerate(els)
        for (i,t1) in enumerate(ts)
            sol = orbitsolve(orb, t1)
            Δt = 100
            t2 = t1 + Δt
            sol′ = orbitsolve(orb,t2)
            diff_lt_app_pmra = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmra2
            diff_lt_app_pmdec = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmdec2
 
            pmra_model_t[j,i] += results[:pmra][ii[j]] +diff_lt_app_pmra
            pmdec_model_t[j,i] += results[:pmdec][ii[j]] +diff_lt_app_pmdec
        end
    end

    
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)
    
        # Draws from the posterior
        mass = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol

        # Now time-series
        sols = orbitsolve.(orbs, ts')
        try
            pmra(first(sols), first(mass))
        catch
            continue
        end
        # TODO: Can we use the existing simulator for this please?
        pmra_model_t .+= pmra.(sols, mass)
        pmdec_model_t .+= pmdec.(sols, mass)
        
        color_model_t .= rem2pi.(
            meananom.(sols), RoundDown) .+ 0 .* ii
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
        nameof(typeof(like_obj)) == :G23HObs
    end
    if !isempty(like_objs)
        absastrom = only(like_objs)

        # The HGCA catalog values have an non-linearity correction added.
        # If we are doing our own rigorous propagation we don't need this
        # correction. We could subtract it from the measurements, but 
        # here we just add it to our model so that they match
        # if absolute_orbits
            hg_nonlinear_dpmra = absastrom.catalog.nonlinear_dpmra[1]
            hg_nonlinear_dpmdec = absastrom.catalog.nonlinear_dpmdec[1]
            hip_nonlinear_dpmra = 2absastrom.catalog.nonlinear_dpmra[1]
            hip_nonlinear_dpmdec = 2absastrom.catalog.nonlinear_dpmdec[1]
        # else
        #     hg_nonlinear_dpmra = 
        #     hg_nonlinear_dpmdec = 
        #     hip_nonlinear_dpmra = 
        #     hip_nonlinear_dpmdec = zero(absastrom.catalog.nonlinear_dpmra[1])
        # end


        θ_systems_from_chain = Octofitter.mcmcchain2result(model, results)
        # Display all points, unless there are more than 5k 
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
            name = Octofitter.normalizename(likelihoodname(absastrom))
            θ_obs = θ_system.observations[name]
            if hasproperty(absastrom, :table)
                solutions = map(orbits) do orbit
                    return orbitsolve.(orbit, absastrom.table.epoch)
                end
                sim = Octofitter.simulate(absastrom, θ_system, θ_obs, orbits, solutions, 0)
            else
                solutions = [() for _ in length(model.system.planets)]
                sim = Octofitter.simulate(absastrom, θ_system, θ_obs, orbits, solutions, -1)
            end
            push!(sims, sim)
        end


        component_flags =  [
            :ra_hip, :dec_hip,
            :ra_hg, :dec_hg,
            :ra_dr2, :dec_dr2,
            :ra_dr32, :dec_dr32,
            :ra_dr3, :dec_dr3,
            :ueva_dr3,
        ]
        sim_mask = [flag ∈ absastrom.table.kind for flag in component_flags]

        mask = findall([contains(string(k),"ra") for k in absastrom.table.kind[:]])
        sim_mask_ra = findall([contains(string(k),"ra") && k ∈ absastrom.table.kind for k in component_flags])
        x = vec(absastrom.table.epoch)
        y = vec(absastrom.table.pm)
        σ = vec(absastrom.table.σ_pm)
        if !isempty(mask)
            # Apply nonlinear correction to get true proper motion for plotting
            y_corrected = copy(y)
            # if absolute_orbits
                # Remove nonlinear corrections from Hipparcos and HGCA measurements 
                # to show true proper motion consistently with model lines
                hip_ra_mask = [k == :ra_hip for k in absastrom.table.kind]
                hg_ra_mask = [k == :ra_hg for k in absastrom.table.kind]
                y_corrected[hip_ra_mask] .-= hip_nonlinear_dpmra
                y_corrected[hg_ra_mask] .-= hg_nonlinear_dpmra
            # end
            
            scatter!(
                ax_velra,
                concat_with_nan(stack(x[mask] for _ in 1:length(sims))),
                concat_with_nan(stack(map(sim->sim.μ[sim_mask_ra], sims))),
                color=:grey,
                markersize=2,
            )
            errorbars!(
                ax_velra,
                x[mask], y_corrected[mask],
                x[mask] .- (absastrom.table.start_epoch[mask]),
                (absastrom.table.stop_epoch[mask]) .- x[mask],
                color=:black,
                direction=:x, linewidth=0.7
            )
            errorbars!(ax_velra, x[mask], y_corrected[mask], σ[mask], color=:black, linewidth=0.7)
            scatter!(ax_velra, x[mask], y_corrected[mask], color=:black)
            ymin = minimum(y_corrected[mask] .-  10σ[mask])
            ymax = maximum(y_corrected[mask] .+  10σ[mask])
            ylims!(ax_velra, ymin, ymax)
        end

        mask = findall([contains(string(k),"dec") for k in absastrom.table.kind[:]])
        sim_mask_dec = findall([contains(string(k),"dec") && k ∈ absastrom.table.kind for k in component_flags])
        x = vec(absastrom.table.epoch)
        y = vec(absastrom.table.pm)
        σ = vec(absastrom.table.σ_pm) 
        if !isempty(mask)
            # Apply nonlinear correction to get true proper motion for plotting
            y_corrected = copy(y)
            # if absolute_orbits
                # Remove nonlinear corrections from Hipparcos and HGCA measurements 
                # to show true proper motion consistently with model lines
                hip_dec_mask = [k == :dec_hip for k in absastrom.table.kind]
                hg_dec_mask = [k == :dec_hg for k in absastrom.table.kind]
                y_corrected[hip_dec_mask] .-= hip_nonlinear_dpmdec
                y_corrected[hg_dec_mask] .-= hg_nonlinear_dpmdec
            # end
            
            scatter!(
                ax_veldec,
                concat_with_nan(stack(x[mask] for _ in 1:length(sims))),
                concat_with_nan(stack(map(sim->sim.μ[sim_mask_dec], sims))),
                color=:grey,
                markersize=2,
            )
            errorbars!(
                ax_veldec,
                x[mask], y_corrected[mask],
                x[mask] .- (absastrom.table.start_epoch[mask]),
                (absastrom.table.stop_epoch[mask]) .- x[mask],
                color=:black,
                direction=:x, linewidth=0.7
            )
            errorbars!(ax_veldec, x[mask], y_corrected[mask], σ[mask], color=:black, linewidth=0.7)
            scatter!(ax_veldec, x[mask], y_corrected[mask], color=:black)
            ymin = minimum(y_corrected[mask] .-  10σ[mask])
            ymax = maximum(y_corrected[mask] .+  10σ[mask])
            ylims!(ax_veldec, ymin, ymax)
        end


        jj = findall(sim_mask)
        mask = [flag ∈ component_flags for flag in vec(absastrom.table.kind)]
        ax_Z = Axis(
            gs[gs_row += 1, 1:3];
            ylabel="Z-scores",
            xgridvisible=false,
            ygridvisible=false,
            xticksvisible=bottom_time_axis,
            xticklabelsvisible=bottom_time_axis,
            xlabelvisible=bottom_time_axis,
            xticks = (
                1:length(jj),
                replace.(string.(vec(absastrom.table.kind[mask])),"_" => " ", "ra"=>"μα*", "dec"=>"μδ")
            ),
            xticklabelrotation=pi/2,
            axis...
        )

        # plot all of ra,dec, ueva on axis in order of epochs
        # plot iad in another panel

        μ_h_cat, Σ_h = isnothing(absastrom.catalog.dist_hip) ? ([0.,0.], zeros(2,2)) : params(absastrom.catalog.dist_hip) 
        μ_hg_cat, Σ_hg = isnothing(absastrom.catalog.dist_hg) ? ([0.,0.], zeros(2,2)) : params(absastrom.catalog.dist_hg) 
        μ_dr2_cat, Σ_dr2 = params(absastrom.catalog.dist_dr2)
        μ_dr32_cat, Σ_dr32 = params(absastrom.catalog.dist_dr32)
        μ_dr3_cat, Σ_dr3 = params(absastrom.catalog.dist_dr3)
        Σ_h = Matrix(Σ_h)
        Σ_hg = Matrix(Σ_hg)
        Σ_dr2 = Matrix(Σ_dr2)
        Σ_dr32 = Matrix(Σ_dr32)
        Σ_dr3 = Matrix(Σ_dr3)

        z_scores = map(sims) do sim
            y = collect(vec(absastrom.table.pm))
            
            # Apply the same nonlinear corrections as in the plots above for consistency
            # if absolute_orbits
                # Remove nonlinear corrections from Hipparcos and HGCA measurements 
                # to show true proper motion consistently
            hip_ra_idx = findfirst(==(:ra_hip), absastrom.table.kind)
            hip_dec_idx = findfirst(==(:dec_hip), absastrom.table.kind) 
            hg_ra_idx = findfirst(==(:ra_hg), absastrom.table.kind)
            hg_dec_idx = findfirst(==(:dec_hg), absastrom.table.kind)
            
            if !isnothing(hip_ra_idx)
                y[hip_ra_idx] -= hip_nonlinear_dpmra
            end
            if !isnothing(hip_dec_idx)
                y[hip_dec_idx] -= hip_nonlinear_dpmdec
            end
            if !isnothing(hg_ra_idx)
                y[hg_ra_idx] -= hg_nonlinear_dpmra
            end
            if !isnothing(hg_dec_idx)
                y[hg_dec_idx] -= hg_nonlinear_dpmdec
            end
            # end
            
            idx = findfirst(==(:ueva_dr3), absastrom.table.kind)
            if !isnothing(idx)
                y[idx] = sim.μ_1_3
            end
           
            resids = sim.μ[jj] .- y[mask]

            sigmas = [
                sqrt.(diag(Σ_h))
                sqrt.(diag(Σ_hg))
                sqrt.(diag(Σ_dr2))
                sqrt.(diag(Σ_dr32))
                sqrt.(diag(Σ_dr3))
                sim.UEVA_unc
            ][jj]
            return resids ./ sigmas
        end

        hlines!(ax_Z, [0], linewidth=3, color=:black)
        hlines!(ax_Z, [-1,1], linewidth=0.5, color=:black)
        boxplot!(
            ax_Z,
            vec(stack(1:length(jj) for _ in 1:length(sims))),
            vec(stack(z_scores)),
            markersize=2,
            color=:grey
        )
        # scatter!(
        #     ax_Z,
        #     concat_with_nan(stack(1:length(jj) for _ in 1:length(sims))),
        #     concat_with_nan(stack(z_scores)),
        # )


    end




    return [ax_velra, ax_veldec,]
end