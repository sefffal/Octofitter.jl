


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
    Octofitter.astromplot!(fig.layout, model, results, args...; kwargs...)

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
    colormap=:plasma,
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


    # Detect if should use arcseconds instead of mas for plotting
    use_arcsec = false
    axis_mult = 1.0

    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)

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
            axis_mult = 1e-3
            use_arcsec = true
        end
    end

    ax = Axis(
        gs[1, 1];
        autolimitaspect=1,
        xreversed=true,
        xlabel= use_arcsec ? "Δra [as]" : "Δra [mas]",
        ylabel= use_arcsec ? "Δdec [as]" : "Δdec [mas]",
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
        orbs = Octofitter.construct_elements(results, planet_key, ii)

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
            raoff(first(sols))
        catch
            continue
        end

        lines!(ax,
            concat_with_nan(raoff.(sols)).*axis_mult,
            concat_with_nan(decoff.(sols)).*axis_mult;
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

    # Now over plot any astrometry
    like_objs = []
    for planet in model.system.planets
        append!(like_objs, planet.observations)
    end
    append!(like_objs, model.system.observations)

    for like_obj in like_objs
        if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood

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
                        x - length_major * cos(α),
                        x,
                        x + length_major * cos(α),
                        NaN,
                        # Minor axis
                        x - length_minor * cos(α + π / 2),
                        x,
                        x + length_minor * cos(α + π / 2),
                        NaN,
                    ]
                    yvals = [
                        # Major axis
                        y - length_major * sin(α),
                        y,
                        y + length_major * sin(α),
                        NaN,
                        # Minor axis
                        y - length_minor * sin(α + π / 2),
                        y,
                        y + length_minor * sin(α + π / 2),
                        NaN,
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
            Makie.scatter!(
                ax,
                vec(x).*axis_mult,
                vec(y).*axis_mult,
                color=:white,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )
        # If the model, instead/in addition to astrometry, includes one of the following 
        # "position-like" observations, add scatter points at the posterior projected locatinos.
        elseif nameof(typeof(like_obj)) in (:ImageLikelihood, :LogLikelihoodMap, :InterferometryLikelihood, :GRAVITYWideCPLikelihood)
            
            for planet_key in keys(model.system.planets)
                orbs = Octofitter.construct_elements(results, planet_key, ii)
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
                orbs = Octofitter.construct_elements(results, planet_key, ii)
                color = Makie.wong_colors()[mod1(i,end)]
                sols = orbitsolve.(orbs, epoch_mjd)
                Makie.scatter!(
                    ax,
                    vec(raoff.(sols)).*axis_mult,
                    vec(decoff.(sols)).*axis_mult;
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
    colormap=:plasma,
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
        orbs = Octofitter.construct_elements(results, planet_key, ii)

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
                orbs = Octofitter.construct_elements(results, planet_key, ii)
                color = Makie.wong_colors()[mod1(i,end)]
                sols = orbitsolve.(orbs, epoch_mjd)
                Makie.scatter!(
                    ax,
                    vec(posx.(sols)),
                    vec(posy.(sols)),
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
