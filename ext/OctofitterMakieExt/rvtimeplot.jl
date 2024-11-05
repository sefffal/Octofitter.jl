


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
    planet_rv=nothing,
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
        planet_mass = map(nt->nt.planets[planet_key].mass, nt_format[ii])
        rv_model_t = radvel.(sols, planet_mass.*Octofitter.mjup2msol)
        maximum(abs, rv_model_t)
        if maximum(abs, rv_model_t) > 1000
            use_kms = true
            kms_mult = 1e-3
        end
    end

    date_pos, date_strs, xminorticks = _date_ticks(ts)
    ax= Axis(
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

    rv_star_model_t = zeros(length(ii), length(ts))
    color_model_t = zeros(length(ii), length(ts))
    for (i_planet, planet_key) in enumerate(keys(model.system.planets))
      
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        sols = orbitsolve.(orbs, ts')

        # Draws from the posterior
        M_planet = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol
        # M_planet = (results["$(planet_key)_mass"][ii].*0 .+ 10) .* Octofitter.mjup2msol
        M_tot = results["M"][ii]
        M_star = M_tot - M_planet

        ### Star RV    
        rv_star_model_t .+= radvel.(sols, M_planet)

        color_model_t .= rem2pi.(
            meananom.(sols), RoundDown) .+ 0 .* ii

        ### Planet RV
        if planet_rv

            # rv_star = -radvel(planet_orbit, epochs[epoch_i])*planet_mass*Octofitter.mjup2msol/M_tot
            # rv_planet_vs_star = radvel(planet_orbit, epochs[epoch_i])
            # rv_planet = rv_star + rv_planet_vs_star
            # rv_planet_buf[epoch_i] -= rv_planet
            # rv_planet_model_t = radvel.(sols) .* M_star ./ M_tot
            
            rv_planet_model_t = rv_star_model_t  .+ radvel.(sols)

            if haskey(results, Symbol("$(planet_key)_rv0_1"))
                rv0s = []
                for i in 1:10
                    k = Symbol("$(planet_key)_rv0_$i")
                    if haskey(results, k)
                        push!(rv0s, results[k][ii])
                    end
                end
                ave_rv0 = mean(stack(rv0s),dims=2)
                rv_planet_model_t .+= ave_rv0
            end
            lines!(ax,
                concat_with_nan(ts' .+ 0 .* rv_planet_model_t),
                concat_with_nan(rv_planet_model_t) .* kms_mult;
                alpha,
                # color=concat_with_nan(color_planet_model_t),
                # colorrange=(0,2pi),
                # colormap
                color=Makie.wong_colors()[1+i_planet],
                label=string(planet_key),
                rasterize=4,
            )

        end


    end

    lines!(ax,
        concat_with_nan(ts' .+ 0 .* rv_star_model_t),
        concat_with_nan(rv_star_model_t) .* kms_mult;
        alpha,
        color = planet_rv ? Makie.wong_colors()[1] : concat_with_nan(color_model_t),
        colorrange=(0,2pi),
        colormap,
        label="A",
        rasterize=4,
    )

    if planet_rv
        Makie.axislegend(ax)
    end


    
    
    i_like = 0
    orbits = map(keys(model.system.planets)) do planet_key
        Octofitter.construct_elements(results, planet_key, ii)
    end
    for like_obj in model.system.observations
        if nameof(typeof(like_obj)) != :MarginalizedStarAbsoluteRVLikelihood &&
            nameof(typeof(like_obj)) != :StarAbsoluteRVLikelihood
            continue
        end
        i_like += 1
        epoch = vec(like_obj.table.epoch)
        rv = collect(vec(like_obj.table.rv))
        σ_rv = vec(like_obj.table.σ_rv)

        jitter_symbol = like_obj.jitter_symbol

        jitter = map(sample->sample[jitter_symbol], nt_format)
        # rv0 = map(sample->sample.rv0[inst_idx], nt_format)'
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

        σ_tot = sqrt.(σ_rv .^2 .+ mean(jitter) .^2 .+ var(rv0))
        rv .-= mean(rv0,dims=2)
        Makie.errorbars!(
            ax, epoch, rv .* kms_mult, σ_tot .* kms_mult;
            color = planet_rv ? Makie.wong_colors()[1] : :grey,
            linewidth=1,
        )
        Makie.errorbars!(
            ax, epoch, rv .* kms_mult, σ_rv .* kms_mult;
            color = Makie.wong_colors()[i_like],
            linewidth=2,
        )
        Makie.scatter!(
            ax, epoch, rv .* kms_mult;
            color = Makie.wong_colors()[i_like],
            strokewidth=2,
            strokecolor=:black,
            markersize=8, #1.5,
        )
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
