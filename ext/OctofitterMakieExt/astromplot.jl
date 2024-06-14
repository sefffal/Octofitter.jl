


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
    kwargs...
)
    gs = gridspec_or_fig

    ax = Axis(
        gs[1, 1];
        autolimitaspect=1,
        xreversed=true,
        xlabel="Δra [mas]",
        ylabel="Δdec [mas]",
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )

    # Start by plotting the orbits
    EAs = range(0, 2pi, length=150)

    # Find earliest epoch
    epoch_0 = mjd("2020")
    for planet in model.system.planets
        for like_obj in planet.observations
            x = Float64[]
            y = Float64[]
            xs = Float64[]
            ys = Float64[]
            if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood
                epoch_0 = min(epoch_0, minimum(like_obj.table.epoch))
            end
        end
    end

    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)

        # Draws from the posterior
        sols = orbitsolve_eccanom.(orbs, EAs')
        # sols_0 = orbitsolve.(orbs, epoch_0)

        lines!(ax,
            concat_with_nan(raoff.(sols)),
            concat_with_nan(decoff.(sols));
            color=concat_with_nan(
                # rem2pi.(EAs' .- eccanom.(sols_0), RoundDown) .+ 0 .* ii
                rem2pi.(meananom.(sols), RoundDown) .+ 0 .* ii
            ),
            alpha=min.(1, 100 / length(ii)),
            transparency=true,
            colormap
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
            Makie.lines!(ax, xs, ys, color=:white, linewidth=3)
            Makie.lines!(ax, xs, ys, color=:black, linewidth=2)
            Makie.scatter!(
                ax,
                vec(x),
                vec(y),
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
                    vec(raoff.(sols)),
                    vec(decoff.(sols)),
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
                    vec(raoff.(sols)),
                    vec(decoff.(sols));
                    color,
                    markersize=6,
                    strokewidth=1,
                    strokecolor=:black,
                    label = labels[i]
                )
            end
        end
        axislegend(ax, "Posterior Predictions", position=:lt)
    end

    if colorbar 
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
    end

    scatter!(ax, [0],[0],marker='⭐', markersize=30, color=:black)
end
