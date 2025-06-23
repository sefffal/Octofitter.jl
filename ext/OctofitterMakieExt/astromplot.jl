


#################################################
# Astrometry plot

function Octofitter.astromplot(
    model,
    results,
    fname="$(model.system.name)-astromplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(; size=(700,600), figure...)
    update(fig) do fig
        Octofitter.astromplot!(fig.layout, model, results, args...; kwargs...)
    end
    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function Octofitter.astromplot!(
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
    colormap=Makie.cgrad([Makie.wong_colors()[1], "#DDDDDD"]),
    colormap_instruments=Makie.cgrad(:Egypt,categorical=true),
    colormap_epochs=Makie.cgrad(:Lakota,categorical=true),
    colorbar=true,
    mark_epochs_mjd=Float64[],
    alpha=min.(1, 100 / length(ii)),
    show_post_pred_legend=true,
    show_instrument_names=true,
    use_arcsec=nothing,
    ts,
    kwargs...
)
    gs = gridspec_or_fig

    if length(model.system.planets) == 1
        colormaps = Dict(
            first(keys(model.system.planets)) => colormap
        )
    else
        colormaps = Dict(
            begin
                c = Makie.wong_colors()[i]
                planet_key => Makie.cgrad([c, "#DDDDDD"])
            end
            for (i,planet_key) in enumerate(keys(model.system.planets))
        )
    end

    EAs = range(0, 2pi, length=150)


    # Detect if should use arcseconds instead of mas for plotting
    if isnothing(use_arcsec)
        use_arcsec = false
        for planet_key in keys(model.system.planets)
            orbs = Octofitter.construct_elements(model, results, planet_key, ii)

            # Draws from the posterior
            sols = orbitsolve_eccanom.(orbs, EAs')
            # sols_0 = orbitsolve.(orbs, epoch_0)
            try
                raoff(first(sols))
            catch
                continue
            end
            s = maximum(projectedseparation, sols)
            if s > 1500
                use_arcsec = true
            end
        end
    end
    axis_mult = use_arcsec ? 1e-3 : 1.0 

    ax = Axis(
        gs[1, 1];
        autolimitaspect=1,
        xreversed=true,
        xlabel= use_arcsec ? "Δα* [as]" : "Δα* [mas]",
        ylabel= use_arcsec ? "Δδ [as]" : "Δδ [mas]",
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )

    # Start by plotting the orbits

    # # Find earliest epoch
    # epoch_0 = mjd("2020")
    # for planet in model.system.planets
    #     for like_obj in planet.observations
    #         if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood
    #             epoch_0 = min(epoch_0, minimum(like_obj.table.epoch))
    #         end
    #     end
    # end


    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)

        # Draws from the posterior
        # Normally we want to sample in equal steps of eccentric anomaly to get
        # nice smooth curves with few points
        # if first(orbs) isa AbsoluteVisual
            # ...but this orbit type have certain non-linear effects which make the orbits
            # not repeat each orbit period, and we need to step in time.
            # Solve from an initial epoch, then in equal steps of mean anomaly (best we can do, ideally we'd do steps of eccentic anomaly)
            ts_prime = first(ts) .+ range(0, stop=1, length=150)' .* period.(orbs)

            sols = orbitsolve.(orbs, ts_prime)
        # else
        #     # Orbit is perfectly periodic, so take equal steps in 
        #     sols = orbitsolve_eccanom.(orbs, EAs')
        # end
        

        # Skip orbits that don't have enough information to plot in sky coordinates
        try
            raoff(first(sols))
        catch
            continue
        end


        ra_host_perturbation = zeros(size(ts_prime))
        dec_host_perturbation = zeros(size(ts_prime))
        for planet_key′ in keys(model.system.planets)
            if !haskey(results, Symbol("$(planet_key′)_mass"))
                continue
            end

            other_planet_mass = results["$(planet_key′)_mass"][ii]
            orbit_other = Octofitter.construct_elements(model, results, planet_key′, ii)
            # Only account for interior planets
            mask = semimajoraxis.(orbit_other) .< semimajoraxis.(orbs)

            sols′ = orbitsolve.(orbit_other, ts_prime)

            ra_host_perturbation .+= mask .* raoff.(sols′, other_planet_mass.*Octofitter.mjup2msol)
            dec_host_perturbation .+= mask .* decoff.(sols′, other_planet_mass.*Octofitter.mjup2msol)
        end

        ra_model = (raoff.(sols) .- ra_host_perturbation)
        dec_model = (decoff.(sols) .- dec_host_perturbation)
        
        lines!(ax,
            concat_with_nan(ra_model).*axis_mult,
            concat_with_nan(dec_model).*axis_mult;
            color=concat_with_nan(
                # rem2pi.(EAs' .- eccanom.(sols_0), RoundDown) .+ 0 .* ii
                rem2pi.(meananom.(sols), RoundDown) .+ 0 .* ii
            ),
            alpha=alpha,
            # transparency=true,
            colormap=colormaps[planet_key],
            rasterize=4,
        )
    end

    # Now over plot any astrometry
    like_objs = []
    for planet in model.system.planets
        append!(like_objs, planet.observations)
    end
    append!(like_objs, model.system.observations)

    # Colour data based on the instrument name
    rel_astrom_likes = filter(like_objs) do like_obj
        nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood ||
        (nameof(typeof(like_obj)) == :ObsPriorAstromONeil2019  && nameof(typeof(like_obj.wrapped_like)) == :PlanetRelAstromLikelihood)
    end
    rel_astrom_names = sort(unique(likelihoodname.(rel_astrom_likes)))
    n_rel_astrom = length(rel_astrom_names)

    i_like_obj = 0
    for like_obj in like_objs
        if  nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood ||
            (nameof(typeof(like_obj)) == :ObsPriorAstromONeil2019  && nameof(typeof(like_obj.wrapped_like)) == :PlanetRelAstromLikelihood)

            i_like_obj = findfirst(==(likelihoodname(like_obj)), rel_astrom_names)
            x = Float64[]
            y = Float64[]
            xs = Float64[]
            ys = Float64[]

            if hasproperty(like_obj.table, :sep)
                # Plot astrometry point
                x = @. like_obj.table.sep * sin(like_obj.table.pa)
                y = @. like_obj.table.sep * cos(like_obj.table.pa)
                # And overplot with uncertainty ellipse
                if hasproperty(like_obj.table, :cor)
                    cor = like_obj.table.cor
                else
                    cor = 0
                end
                x = like_obj.table.sep
                y = pi / 2 .- like_obj.table.pa
                σ₁ = like_obj.table.σ_sep
                σ₂ = like_obj.table.σ_pa
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
                        x .+ range(
                            start=- length_major,
                            stop =+ length_major,
                            length=9
                        ).*cos(α)
                        NaN
                        # Minor axis
                        x .+ range(
                            start=- length_minor,
                            stop =+ length_minor,
                            length=9
                        ).* cos(α + π / 2)
                        NaN
                    ]
                    yvals = [
                        y .+ range(
                            start= - length_major,
                            stop = + length_major,
                            length=9
                        ).*sin(α)
                        NaN
                        # Minor axis
                        y .+ range(
                            start= - length_minor,
                            stop = + length_minor,
                            length=9
                        ).* sin(α + π / 2)
                        NaN
                    ]
                    xvals, yvals
                end
                rs = Base.vcat(getindex.(error_ellipses, 1)...)
                thetas = Base.vcat(getindex.(error_ellipses, 2)...)
                xs = rs .* cos.(thetas)
                ys = rs .* sin.(thetas)
                x, y = x .* cos.(y), x .* sin.(y)
            elseif hasproperty(like_obj.table, :ra)
                x = like_obj.table.ra
                y = like_obj.table.dec

                if hasproperty(like_obj.table, :cor)
                    cor = like_obj.table.cor
                else
                    cor = 0
                end

                σ₁ = like_obj.table.σ_ra
                σ₂ = like_obj.table.σ_dec

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
                xs = Base.vcat(getindex.(error_ellipses, 1)...)
                ys = Base.vcat(getindex.(error_ellipses, 2)...)
            end
            Makie.lines!(ax, xs.*axis_mult, ys.*axis_mult, color=:white, linewidth=3,rasterize=4,)
            Makie.lines!(ax, xs.*axis_mult, ys.*axis_mult, color=:black, linewidth=2,rasterize=4,)
            if n_rel_astrom == 1
                color = :white
            else
                color = colormap_instruments[mod1(i_like_obj,end)]
            end
            Makie.scatter!(
                ax,
                vec(x).*axis_mult,
                vec(y).*axis_mult;
                color,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )

        # If the model, instead/in addition to astrometry, includes one of the following 
        # "position-like" observations, add scatter points at the posterior projected locatinos.
        elseif nameof(typeof(like_obj)) in (:ImageLikelihood, :LogLikelihoodMap, :InterferometryLikelihood, :GRAVITYWideCPLikelihood)
            
            for planet_key in keys(model.system.planets)
                orbs = Octofitter.construct_elements(model, results, planet_key, ii)
                # Now time-series
                sols = orbitsolve.(orbs, like_obj.table.epoch')
                Makie.scatter!(
                    ax,
                    vec(raoff.(sols)).*axis_mult,
                    vec(decoff.(sols)).*axis_mult,
                    color=:black,
                    markersize=4,
                )
            end
        end
    end

    row_i = 1
    if show_instrument_names && n_rel_astrom > 1
        row_i += 1
        elements = [
            MarkerElement(marker=:circle ,color=colormap_instruments[mod1(i,end)],strokewidth=1,strokecolor=:black)
            for i in 1:n_rel_astrom
        ]
        Legend(
            gs[row_i,1:2],
            elements,
            rel_astrom_names,
            "Instrument",
            position=:rb,
            # backgroundcolor=(:white, 0.65),
            tellwidth=false,
            tellheight=true,
            width=Relative(1.0)
        )
    end

    # The user can ask us to plot the position at a particular date
    if !isempty(mark_epochs_mjd) 
        i = 0
        labels = map(mark_epochs_mjd) do epoch_mjd
            replace(string(mjd2date(epoch_mjd)), "T"=>" ") # TODO: better to use a format string
        end
        # Find N last characers in common for all labels and chop them off.
        local label_last_shared_index
        for j in reverse(1:minimum(length.(labels)))
            label_last_shared_index = j
            if !all(s-> s∈(':', '0', ' '), getindex.(labels, j))
                break
            end
        end
        label_last_shared_index = max(10, label_last_shared_index)
        labels = map(labels) do label
            label[1:label_last_shared_index]
        end
        for i in eachindex(mark_epochs_mjd)
            epoch_mjd = mark_epochs_mjd[i]
            for (i_planet, planet_key) in enumerate(keys(model.system.planets))
                orbs = Octofitter.construct_elements(model, results, planet_key, ii)
                color = colormap_epochs[mod1(i,end)]
                sols = orbitsolve.(orbs, epoch_mjd)
                kwargs = (;)
                if i_planet == 1
                    kwargs = (;kwargs..., label = labels[i])
                end
                Makie.scatter!(
                    ax,
                    vec(raoff.(sols)).*axis_mult,
                    vec(decoff.(sols)).*axis_mult;
                    color,
                    alpha,
                    marker=:rect,
                    markersize=6,
                    strokewidth=1,
                    strokecolor=(:black,alpha),
                    kwargs...
                )

                println("planet $planet_key astometry\t $epoch_mjd [MJD]")
                @printf("\traoff:\t%.4g ± %.4g [mas]\n", mean(vec(raoff.(sols))), std(vec(raoff.(sols))))
                @printf("\tdecoff:\t%.4g ± %.4g [mas]\n", mean(vec(decoff.(sols))), std(vec(decoff.(sols))))
                @printf("\tsep:\t%.4g ± %.4g [mas]\n", mean(vec(projectedseparation.(sols))), std(vec(projectedseparation.(sols))))
                @printf("\tpa:\t%.4g ± %.4g [deg]\n", mean(vec(rad2deg.(posangle.(sols)))), std(vec(rad2deg.(posangle.(sols)))))
                println()
            end
        end
       
        if show_post_pred_legend
            row_i += 1
            Legend(
                gs[row_i,1:2],
                ax,
                "Posterior Predictions",
                position=:rb,
                # backgroundcolor=(:white, 0.65),
                tellwidth=false,
                tellheight=true,
                width=Relative(1.0)
            )
        end
    end

    if colorbar
        if length(colormaps) == 1
            Colorbar(
                gs[1,2];
                colormap,
                label="mean anomaly →",
                colorrange=(0,2pi),
                ticks=(
                    [0,pi/2,pi,3pi/2,2pi],
                    ["0", "π/2", "π", "3π/2", "2π"]
                )
            )
        else
            col = 1
            for planet_key in keys(colormaps)
                col += 1
                Colorbar(
                    gs[1,col];
                    colormap=colormaps[planet_key],
                    label="mean anomaly →",
                    colorrange=(0,2pi),
                    ticks=(
                        [0,pi/2,pi,3pi/2,2pi],
                        ["0", "π/2", "π", "3π/2", "2π"]
                    ),
                    ticksvisible=col - 1 == length(colormaps),
                    ticklabelsvisible=col - 1 == length(colormaps),
                    labelvisible=col - 1 == length(colormaps)
                )
                Label(gs[1,col,Top()],string(planet_key))
                Makie.colgap!(gs, col-1, 4)
            end
        end
    end

    scatter!(ax, [0],[0],marker='★', markersize=20, color=:white, strokecolor=:black, strokewidth=1.5)
end




function physorbplot!(
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
    colormap=Makie.cgrad([Makie.wong_colors()[1], "#DDDDDD"]),
    colormap_epochs=Makie.cgrad(:Lakota,categorical=true),
    colorbar=true,
    mark_epochs_mjd=Float64[],
    alpha=min.(1, 100 / length(ii)),
    ts,
    kwargs...
)
    gs = gridspec_or_fig

    if length(model.system.planets) == 1
        colormaps = Dict(
            first(keys(model.system.planets)) => colormap
        )
    else
        colormaps = Dict(
            begin
                c = Makie.wong_colors()[i]
                planet_key => Makie.cgrad([c, "#FAFAFA"])
            end
            for (i,planet_key) in enumerate(keys(model.system.planets))
        )
    end

    EAs = range(0, 2pi, length=150)


    ax = Axis(
        gs[1, 1];
        autolimitaspect=1,
        xreversed=true,
        xlabel= "Δx [au]",
        ylabel= "Δy [au]",
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )

    # Start by plotting the orbits

    # # Find earliest epoch
    # epoch_0 = mjd("2020")
    # for planet in model.system.planets
    #     for like_obj in planet.observations
    #         if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood
    #             epoch_0 = min(epoch_0, minimum(like_obj.table.epoch))
    #         end
    #     end
    # end


    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)

        # Draws from the posterior
        if first(orbs) isa AbsoluteVisual
            # Solve from an initial epoch, then in equal steps of mean anomaly (best we can do, ideally we'd do steps of eccentic anomaly)
            MAs = range(0, 2pi, length=150)
            ts_prime = first(ts) .+ range(0, 1 + 1/150, length=150)' .* period.(orbs)
            sols = orbitsolve.(orbs, ts_prime)
        else
            # Orbit is perfectly periodic, so take equal steps in 
            sols = orbitsolve_eccanom.(orbs, EAs')
        end

        try
            posx(first(sols))
        catch
            continue
        end

        lines!(ax,
            concat_with_nan(posx.(sols)),
            concat_with_nan(posy.(sols)),
            color=concat_with_nan(
                # rem2pi.(EAs' .- eccanom.(sols_0), RoundDown) .+ 0 .* ii
                rem2pi.(meananom.(sols), RoundDown) .+ 0 .* ii
            ),
            alpha=alpha,
            transparency=true,
            colormap=colormaps[planet_key],
            rasterize=4,
        )
    end


    # The user can ask us to plot the position at a particular date
    if !isempty(mark_epochs_mjd)
        i = 0
        labels = map(mark_epochs_mjd) do epoch_mjd
            replace(string(mjd2date(epoch_mjd)), "T"=>" ") # TODO: better to use a format string
        end
        # Find N last characers in common for all labels and chop them off.
        local label_last_shared_index
        for j in reverse(1:minimum(length.(labels)))
            label_last_shared_index = j
            if !all(s-> s∈(':', '0', ' '), getindex.(labels, j))
                break
            end
        end
        label_last_shared_index = max(10, label_last_shared_index)
        labels = map(labels) do label
            label[1:label_last_shared_index]
        end
        for i in eachindex(mark_epochs_mjd)
            epoch_mjd = mark_epochs_mjd[i]
            for planet_key in keys(model.system.planets)
                orbs = Octofitter.construct_elements(model, results, planet_key, ii)                
                color = colormap_epochs[mod1(i,end)]
                sols = orbitsolve.(orbs, epoch_mjd)
                Makie.scatter!(
                    ax,
                    vec(posx.(sols)),
                    vec(posy.(sols));
                    color,
                    markersize=6,
                    strokewidth=1,
                    strokecolor=:black,
                    label = labels[i],
                )
            end
        end
        Legend(
            gs[2,1:2],
            ax,
            "Posterior Predictions",
            position=:rb,
            # backgroundcolor=(:white, 0.65),
            tellwidth=false,
            tellheight=true,
            width=Relative(1.0)
        )
    end

    if colorbar
        if length(colormaps) == 1
            Colorbar(
                gs[1,2];
                colormap,
                label="mean anomaly →",
                colorrange=(0,2pi),
                ticks=(
                    [0,pi/2,pi,3pi/2,2pi],
                    ["0", "π/2", "π", "3π/2", "2π"]
                )
            )
        else
            col = 1
            for planet_key in keys(colormaps)
                col += 1
                Colorbar(
                    gs[1,col];
                    colormap=colormaps[planet_key],
                    label="mean anomaly ($planet_key) →",
                    colorrange=(0,2pi),
                    ticks=(
                        [0,pi/2,pi,3pi/2,2pi],
                        ["0", "π/2", "π", "3π/2", "2π"]
                    ),
                    ticksvisible=col - 1 == length(colormaps),
                    ticklabelsvisible=col - 1 == length(colormaps),
                )
            end
        end
    end

    scatter!(ax, [0],[0],marker='★', markersize=20, color=:white, strokecolor=:black, strokewidth=1.5)
end
