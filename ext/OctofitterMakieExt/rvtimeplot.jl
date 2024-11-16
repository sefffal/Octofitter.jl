


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
    ts,
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
    alpha=min.(1, 100 / length(ii)),
    kwargs...
)
    gs = gridspec_or_fig
    nt_format = vec(Octofitter.mcmcchain2result(model, results))

    # Detect if should use km/s
    use_kms = false
    kms_mult = 1.0
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        sols = orbitsolve.(orbs, ts')
        if !hasproperty(first(nt_format[ii]).planets[planet_key], :mass)
            continue
        end
        planet_mass = map(nt->nt.planets[planet_key].mass, nt_format[ii])
        rv_model_t = radvel.(sols, planet_mass.*Octofitter.mjup2msol)
        maximum(abs, rv_model_t)
        if maximum(abs, rv_model_t) > 1000
            use_kms = true
            kms_mult = 1e-3
        end
    end

    date_pos, date_strs, xminorticks = _date_ticks(ts)

    # The only safe way to display this data without misinterpretation is
    # to plot it separately by instrument. 
    # For an all up plot, users can see a single draw from `rvpostplot`.
    all_axes = Axis[]
    like_objs = filter(model.system.observations) do like_obj
        nameof(typeof(like_obj)) ∈ (
            :MarginalizedStarAbsoluteRVLikelihood,
            :StarAbsoluteRVLikelihood
        )
    end


    # CASE: no data, just plot model.
    if isempty(like_objs)
        @info "is empyty"
        ax = Axis(
            gs[1, 1];
            ylabel= use_kms ? "rv [km/s]" : "rv [m/s]",
            xaxisposition=:top,
            xticks=(date_pos, date_strs),
            xgridvisible=false,
            ygridvisible=false,
            xminorticks,
            xminorticksvisible=top_time_axis,
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
        push!(all_axes, ax)
        push!(all_axes, ax_secondary)
        
        rv_star_model_t = zeros(length(ii), length(ts))
        color_model_t = zeros(length(ii), length(ts))
        for (i_planet, planet_key) in enumerate(keys(model.system.planets))
        
            orbs = Octofitter.construct_elements(results, planet_key, ii)
            sols = orbitsolve.(orbs, ts')

            # Draws from the posterior
            key = Symbol("$(planet_key)_mass")
            if haskey(results, key)
                M_planet = results[key][ii] .* Octofitter.mjup2msol
            else
                M_planet = zeros(length(ii))
            end
            M_tot = results["M"][ii]

            ### Star RV    
            rv_star_model_t .+= radvel.(sols, M_planet)

            color_model_t .= rem2pi.(
                meananom.(sols), RoundDown) .+ 0 .* ii
        end    

        lines!(ax,
            concat_with_nan(ts' .+ 0 .* rv_star_model_t),
            concat_with_nan(rv_star_model_t) .* kms_mult;
            alpha,
            color = concat_with_nan(color_model_t),
            colorrange=(0,2pi),
            colormap,
            label="A",
            rasterize=4,
        )
    end

    # CASE: one or more RV likelihood objects.
    # The only safe way to display these is with a separate panel per instrument
    gs_row = 0
    for (i_like, like_obj) in enumerate(like_objs)
        gs_row += 1
        ax = Axis(
            gs[gs_row, 1];
            ylabel= use_kms ? "rv [km/s]" : "rv [m/s]",
            xaxisposition=:top,
            xticks=(date_pos, date_strs),
            xgridvisible=false,
            ygridvisible=false,
            xminorticks,
            xminorticksvisible=top_time_axis && i_like==1,
            xticksvisible=top_time_axis && i_like==1,
            xticklabelsvisible=top_time_axis && i_like==1,
            axis...
        )
        ax_secondary = Axis(
            gs[gs_row, 1];
            xlabel="MJD",
            xgridvisible=false,
            ygridvisible=false,
            yticksvisible=false,
            yticklabelsvisible=false,
            ylabelvisible=false,
            xticksvisible=bottom_time_axis && i_like == length(like_objs),
            xticklabelsvisible=bottom_time_axis && i_like == length(like_objs),
            xlabelvisible=bottom_time_axis && i_like == length(like_objs),
            axis...
        )
        linkxaxes!(ax, ax_secondary)
        xlims!(ax, extrema(ts))
        xlims!(ax_secondary, extrema(ts))
        push!(all_axes, ax)
        push!(all_axes, ax_secondary)
        

        # For each instrument's data, subtract the mean rv0 across all draws
        # Then for each orbit draw, add the difference between the mean rv0 across all instruments, and the above

        orbits = map(keys(model.system.planets)) do planet_key
            Octofitter.construct_elements(results, planet_key, ii)
        end
        epoch = vec(like_obj.table.epoch)
        rv = collect(vec(like_obj.table.rv))
        σ_rv = vec(like_obj.table.σ_rv)

        # instead of each 
        rv0 = map(enumerate(ii)) do (i_sol, i)
            θ_system = nt_format[i]

            # StarAbsoluteRVLikelihood have an RV0 parameter in the model
            if hasproperty(like_obj,:offset_symbol)
                return θ_system[like_obj.offset_symbol]
            end

            # MarginalizedStarAbsoluteRVLikelihood do not, we can compute it in
            # closed form (just like in the likelihood evaluation)
            planet_orbits_this_sample = getindex.(orbits, i_sol)
            return _find_rv_zero_point_maxlike(like_obj, θ_system, planet_orbits_this_sample)
        end'

        rv_star_model_t = zeros(length(ii), length(ts))
        color_model_t = zeros(length(ii), length(ts))
        for (i_planet, planet_key) in enumerate(keys(model.system.planets))
        
            orbs = Octofitter.construct_elements(results, planet_key, ii)
            sols = orbitsolve.(orbs, ts')

            # Draws from the posterior
            key = Symbol("$(planet_key)_mass")
            if haskey(results, key)
                M_planet = results[key][ii] .* Octofitter.mjup2msol
            else
                M_planet = zeros(length(ii))
            end
            M_tot = results["M"][ii]
            M_star = M_tot - M_planet

            # Star RV  influence from this planet   
            rv_star_model_t .+= radvel.(sols, M_planet)

            color_model_t .= rem2pi.(
                meananom.(sols), RoundDown) .+ 0 .* ii

        end    

        # Add RV0 offset to make model match data, but subtract
        # mean offset to keep it near zero
        rv_star_model_t .+= rv0'
        rv_star_model_t .-= mean(rv0)

        lines!(ax,
            concat_with_nan(ts' .+ 0 .* rv_star_model_t),
            concat_with_nan(rv_star_model_t) .* kms_mult;
            alpha,
            color = concat_with_nan(color_model_t),
            colorrange=(0,2pi),
            colormap,
            label="A",
            rasterize=4,
        )


        epoch = vec(like_obj.table.epoch)
        rv = collect(vec(like_obj.table.rv))
        σ_rv = vec(like_obj.table.σ_rv)

        jitter_symbol = like_obj.jitter_symbol
        jitter = map(sample->sample[jitter_symbol], nt_format)
        σ_tot = sqrt.(σ_rv .^2 .+ mean(jitter) .^2)

        # Subtract off mean of all draws RV0 to keep near zero
        rv .-= mean(rv0)
    
        Makie.errorbars!(
            ax, epoch, rv .* kms_mult, σ_tot .* kms_mult;
            color = :grey,
            linewidth=1,
        )
        Makie.errorbars!(
            ax, epoch, rv .* kms_mult, σ_rv .* kms_mult;
            color = :black,
            linewidth=2,
        )
        Makie.scatter!(
            ax, epoch, rv .* kms_mult;
            color = :white,
            strokewidth=2,
            strokecolor=:black,
            markersize=8, #1.5,
        )

        # Label inside top right corner
        Makie.text!(ax,
            1.0, 1.0,
            text=string(like_obj.instrument_name),
            align=(:right,:top),
            space=:relative,
            offset = (-4, -4),
        );
        
    end

    if colorbar
        Colorbar(
            gs[1,2];
            colormap,
            label="mean anomaly",
            colorrange=(0,2pi),
            ticks=(
                [0,pi/2,pi,3pi/2,2pi],
                ["0", "π/2", "π", "3π/2", "2π"]
            )
        )
    end

    return all_axes
end




function rvtimeplot_relative(
    model,
    results,
    fname="$(model.system.name)-rvplot-relative.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700,600),
        figure...
    )
    rvtimeplot_relative!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function rvtimeplot_relative!(
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
    ts,
    axis=(;),
    colormap=:plasma,
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    alpha=min.(1, 100 / length(ii)),
    kwargs...
)
    gs = gridspec_or_fig
    nt_format = vec(Octofitter.mcmcchain2result(model, results))

    date_pos, date_strs, xminorticks = _date_ticks(ts)

    # Detect if should use km/s
    use_kms = false
    kms_mult = 1.0
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        sols = orbitsolve.(orbs, ts')
        rv_model_t = radvel.(sols)
        if median(abs.(vec(rv_model_t))) > 1000
            use_kms = true
            kms_mult = 1e-3
        end
    end

    ax= Axis(
        gs[1, 1];
        ylabel= use_kms ? "relative rv [km/s]" : "relative rv [m/s]",
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
        xminorticks,
        xminorticksvisible=top_time_axis,
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

    for (i_planet, planet_key) in enumerate(keys(model.system.planets))
        # Plot each relative RV separately
        orbs = Octofitter.construct_elements(results, planet_key, ii)

        sols = orbitsolve.(orbs, ts')
        
        rv_model_t = radvel.(sols)
        color_model_t = rem2pi.(meananom.(sols), RoundDown) .+ 0 .* ii

        lines!(
            ax,
            concat_with_nan(ts' .+ 0 .* ii),
            concat_with_nan(rv_model_t) .* kms_mult;
            color=concat_with_nan(color_model_t),
            colormap,
            alpha,
            rasterize=4,
        )
    end

    # Now overplot the data points, if any.
    # We can do this for relative RV since there is no zero point offset
    for planet_key in keys(model.system.planets)
        planet = getproperty(model.system.planets, planet_key)
        for like_obj in planet.observations
            if nameof(typeof(like_obj)) != :PlanetRelativeRVLikelihood
                continue
            end
            epoch = vec(like_obj.table.epoch)
            rv = vec(like_obj.table.rv)
            σ_rv = vec(like_obj.table.σ_rv)
            jitter = map(nt_format) do θ_system
                θ_system.planets[planet_key][like_obj.jitter_symbol]
            end
            σ_tot = median(sqrt.(σ_rv .^2 .+ jitter' .^2),dims=2)[:]
            Makie.errorbars!(
                ax, epoch, rv .* kms_mult, σ_tot .* kms_mult;
                color =  :grey,
                linewidth=1,
            )
            Makie.errorbars!(
                ax, epoch, rv.*kms_mult, σ_rv .* kms_mult;
                color =  :black,
                linewidth=2,
            )
            Makie.scatter!(
                ax, epoch, rv.*kms_mult;
                color=:white,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )
        end
    end

    if colorbar
        Colorbar(
            gs[1,2];
            colormap,
            label="mean anomaly",
            colorrange=(0,2pi),
            ticks=(
                [0,pi/2,pi,3pi/2,2pi],
                ["0", "π/2", "π", "3π/2", "2π"]
            )
        )
    end
    return [ax]
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
        if hasproperty(θ_system.planets[planet_i], :mass)
            planet_mass = θ_system.planets[planet_i].mass
        else
            planet_mass = 0.0
        end

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
