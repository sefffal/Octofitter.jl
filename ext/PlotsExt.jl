module PlotsExt
using Octofitter
using Plots
using RecipesBase

"""
    plotchains(
        chain, planet_key;
        N=1500,
        ii = rand(1:size(chain,1)*size(chain,3), N),
        color=length(model.system.planets) == 0 || !haskey(chain, string(planet_key)*"_a") ? nothing : string(planet_key)*"_a",
        colorbartitle=color,
        clims=nothing,
        cmap=:plasma,
        alpha=30/length(ii),
        attime=nothing,
        kwargs...,
    )

Draw samples from a posterior chain for a given planet given by name `planet_key` and visualize them in some way.
Use `kind` to control what plot is made. A few options: :astrometry, :radvel, :trueanom, :meananom, :eccanom, :x, :y, :z, (:x, :y), :raoff, :decoff, :pmra, :pmdec, :accra, :accdec, :radvel, :posangle, :projectedseparation.
See PlanetOrbits documentation for more details.

Inputs:
* chain                   The chain to draw from
* planet_key              Planet name in the model (symbol)
* N=1500                  Number of samples to draw for the plot 
* kind=nothing            Specify what kind of plot to make. 
* ii=...                  Specific row numbers to use, if you want to e.g. plot the same 100 samples in a few different plots
* color="\$planet_key_a"   Column name to to map colors to. Semi-major axis by default but can be any column or an arbitrary array.
* colorbartitle=color     Name for colourbar
* clims=nothing           Tuple of colour limits (min and max)
* cmap=:plasma            Colormap
* alpha=...               Transparency of the lines
"""
function Octofitter.plotchains(args...;kwargs...)
    p = Plots.plot()
    plotchains!(p, args...;kwargs...)
end

"""
    plotchains!(
        plot, ...;
        kwargs...,
    )

See plotchains
"""
function Octofitter.plotchains!(
    p, chain, planet_key;
    N=1500,
    ii = rand(1:size(chain,1)*size(chain,3), N),
    color=!haskey(chain, string(planet_key)*"_a") ? nothing : string(planet_key)*"_a",
    colorbartitle=color,
    clims=nothing,
    cmap=:plasma,
    alpha=30/length(ii),
    attime=nothing,
    mass=nothing,
    kwargs...,
)
    # Construct orbits from chain
    orbits = Octofitter.construct_elements(chain, planet_key, ii)
    if !isnothing(attime)
        orbits = orbitsolve.(orbits, attime)
    end

    # Handle the colours of each element
    if isnothing(color)
        colours = nothing
    elseif typeof(color) <: Union{Symbol,AbstractString}
        if haskey(chain, Symbol(color))
            colours = chain[color][ii]
        else
            colours = nothing
            kwargs = (;kwargs..., color)
        end
    elseif typeof(color) <: AbstractArray # They passed in the values directly
        colours = color[ii]
    else
        error("unkown color option")
    end

    if !isnothing(colours)
        if isnothing(clims)
            clims=quantile(Iterators.flatten(
                extrema(colours)
            ), (0.01,0.99))
        end
        minc = first(clims)
        dc = last(clims) - minc
        clims = (minc, dc+minc)
    end


    cmap= Plots.cgrad(cmap)
    for i in eachindex(orbits)
        if !isnothing(colours)
            c = colours[i]
            kwargs_loc = (;color=cmap[(c-minc)/dc], kwargs...)
        else
            kwargs_loc = kwargs
        end
        m = nothing
        if !isnothing(mass)
            if length(mass) > 1
                m = mass[ii][i]*mjup2msol
            else
                m = mass*mjup2msol
            end
        end
        Plots.plot!(p, orbits[i]; label="", alpha, mass=m, colorbar=true, kwargs_loc...)
    end
    # Colorbar
    # Plots.scatter!(p, [NaN],[NaN]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)
    Plots.scatter!(p, [0],[0]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)
    Plots.scatter!(p, [0],[0], marker=(:star, :black, 5), markerstrokewidth=1, markerstrokecolor=:white,  label="")
    return p
    # # Star at centre

end



"""
Plot a parameter against time.

Example:
```
plot(
    timeplot(chains_pma, :b, :mass, :ra),
    timeplot(chains_pma, :b, :mass, :dec),
    timeplot(chains_pma, :b, :mass, :rv),
    timeplot(chains_pma, :b, :mass, :pmra),
    timeplot(chains_pma, :b, :mass, :pmdec),
    layout = @layout([
        A
        B
        C
        D E
    ]),
    framestyle=:box,
    size=(500,700)
)
```
"""
function timeplot(args...; kwargs...)
    Plots.plot()
    timeplot!(args...; kwargs...)
end
function timeplot!(
    model,
    chain,
    planet_key,
    color,
    prop,
    ;
    N = 1500,
    ii = rand(1:size(chain,1)*size(chain,3), N),
    alpha = 0.05,
    cmap=:plasma,
    clims = quantile(vec(chain[color]),(0.01, 0.99)),
    kwargs...
)
    # Handle the colours of each element
    if typeof(color) <: Union{Symbol,AbstractString}
        color = chain[color][ii]
    else # They passed in the values directly
        color = color[ii]
    end
    planet = model.system.planets[planet_key]
    if prop == :ra
        ylabel = "RA (mas)"
    elseif prop == :dec
        ylabel = "DEC (mas)"
    elseif prop == :sep
        ylabel = "SEP (mas)"
    elseif prop == :pa
        ylabel = "PA (°)"
    elseif prop == :rv
        ylabel = "RV m/s"
    elseif prop == :pmra
        ylabel = "∂RA mas/yr"
    elseif prop == :pmdec
        ylabel = "∂DEC mas/yr"
    else
        error("Unsupported property. Choose :ra, :dec, :sep, :pa, :rv, :pmra, or :pmdec")
    end
    p1 = Plots.plot!(;
        ylabel,
        legend=:none
    )
    xerr = nothing
    if planet_key isa String || planet_key isa Symbol
        planet_keys = tuple(planet_key)
    else
        planet_keys = planet_key
    end
    all_epochs_planet = mapreduce(vcat, keys(model.system.planets)) do planet_key
        planet = model.system.planets[planet_key]
        reduce(vcat, Any[hasproperty(t.table, :epoch) ? t.table.epoch : t.table.ra_epoch for t in planet.observations if hasproperty(t.table, :epoch) || hasproperty(t.table, :ra_epoch)], init=[])
    end
    all_epochs_star = mapreduce(vcat, model.system.observations, init=Float64[]) do obs
        if hasproperty(obs.table, :epoch)
            return obs.table.epoch
        else
            return Float64[]
        end
    end
    all_epochs = [all_epochs_planet; all_epochs_star]
    if isempty(all_epochs) 
        all_epochs = mjd() .+ [-365*2, +365*2]
    end
    t = range((extrema(all_epochs) .+ [-365, 365])..., length=100)
    # t = range((extrema(all_epochs) .+ [-365, 365])..., length=500)
    y = nothing
    ribbon = nothing
    zcolor = nothing
    if prop == :ra
        if !isnothing(astrometry(planet))
            y = astrometry(planet).table.ra
            yerr = astrometry(planet).table.σ_ra
            x = astrometry(planet).table.epoch
        end
        elements = Octofitter.construct_elements(chain, planet_key, ii)
        fit = raoff.(elements, t')
    elseif prop == :dec
        if !isnothing(astrometry(planet))
            y = astrometry(planet).table.dec
            yerr = astrometry(planet).table.σ_dec
            x = astrometry(planet).table.epoch
        end
        elements = Octofitter.construct_elements(chain, planet_key, ii)
        fit = decoff.(elements, t')
    elseif prop == :sep
        if !isnothing(astrometry(planet))
            if hasproperty(astrometry(planet).table, :sep)
                y = astrometry(planet).table.sep
                yerr = astrometry(planet).table.σ_sep
            else
                xx = astrometry(planet).table.ra
                yy = astrometry(planet).table.dec
                xxerr = astrometry(planet).table.σ_ra
                yyerr = astrometry(planet).table.σ_dec
                y_unc = sqrt.((xx .± xxerr).^2 .+ (yy .± yyerr).^2)
                y = Measurements.value.(y_unc)
                yerr = Measurements.uncertainty.(y_unc)
            end
            x = astrometry(planet).table.epoch
        end
        elements = Octofitter.construct_elements(chain, planet_key, ii)
        fit = projectedseparation.(elements, t')
        if Symbol("$(planet_key)_mass") in keys(chain)
            fit  = fit# .- projectedseparation.(elements, t', chain["$(planet_key)_mass"][ii].*Octofitter.mjup2msol)
        end
    elseif prop == :pa
        if !isnothing(astrometry(planet))
            if hasproperty(astrometry(planet).table, :pa)
                y = rad2deg.(rem2pi.(astrometry(planet).table.pa, RoundNearest))
                yerr = rad2deg.(astrometry(planet).table.σ_pa)
            else
                xx = astrometry(planet).table.ra
                yy = astrometry(planet).table.dec
                xxerr = astrometry(planet).table.σ_ra
                yyerr = astrometry(planet).table.σ_dec
                y_unc = atand.((xx .± xxerr), (yy .± yyerr))
                y = Measurements.value.(y_unc)
                yerr = Measurements.uncertainty.(y_unc)
            end
            x = astrometry(planet).table.epoch
        end
        elements = Octofitter.construct_elements(chain, planet_key, ii)
        fit = rad2deg.(posangle.(elements, t'))
        # Avoid ugly wrapping
        diffs = diff(fit, dims=1)
        fit[findall(abs.(diffs) .> 100)] .= NaN
    elseif prop == :pmra
        if !isnothing(propermotionanom(model.system))
            hgca = propermotionanom(model.system).table[1]
            t = range(years2mjd(hgca.epoch_ra_hip)-365*5, years2mjd(hgca.epoch_ra_gaia)+365*5, length=100)
            y = [hgca.pmra_hip, hgca.pmra_hg, hgca.pmra_gaia]
            yerr = [hgca.pmra_hip_error, hgca.pmra_hg_error, hgca.pmra_gaia_error]
            x = years2mjd.([hgca.epoch_ra_hip, (hgca.epoch_ra_hip + hgca.epoch_ra_gaia)/2, hgca.epoch_ra_gaia])
            xerr = [4*365/2, 25*365/2, 3*365/2]
        end
        fit = 0
        for planet_key in planet_keys
            elements = Octofitter.construct_elements(chain, planet_key, ii)
            fit = fit .+ pmra.(elements, t', collect(chain["$(planet_key)_mass"][ii]).*Octofitter.mjup2msol)
        end
        if :pmra in keys(chain)
            fit = fit .+ chain["pmra"][ii]
        end
    elseif prop == :pmdec
        if !isnothing(propermotionanom(model.system))
            hgca = propermotionanom(model.system).table[1]
            t = range(years2mjd(hgca.epoch_dec_hip)-365*5, years2mjd(hgca.epoch_dec_gaia)+365*5, length=100)
            y = [hgca.pmdec_hip, hgca.pmdec_hg, hgca.pmdec_gaia]
            yerr = [hgca.pmdec_hip_error, hgca.pmdec_hg_error, hgca.pmdec_gaia_error]
            x = years2mjd.([hgca.epoch_dec_hip, (hgca.epoch_dec_hip + hgca.epoch_dec_gaia)/2, hgca.epoch_dec_gaia])
            xerr = [4*365/2, 25*365/2, 3*365/2]
        end
        fit = 0
        for planet_key in planet_keys
            elements = Octofitter.construct_elements(chain, planet_key, ii)
            fit = fit .+ pmdec.(elements, t', collect(chain["$(planet_key)_mass"][ii]).*Octofitter.mjup2msol)
        end
        if :pmdec in keys(chain)
            fit = fit .+ chain["pmdec"][ii]
        end
    elseif prop == :rv
        # ribbon = chain[jitter_inst][ii]
        fit = 0
        for planet_key in planet_keys
            elements = Octofitter.construct_elements(chain, planet_key, ii)
            fit = fit .+ radvel.(elements, t', collect(chain["$(planet_key)_mass"][ii]).*mjup2msol)
        end
        for obs in model.system.observations
            # TODO: make this pluggable instead of this hacky workaround
            if startswith(string(typeof(obs)), "RadialVelocity")
                if haskey(chain,:rv0)
                    barycentric_rv_inst_1 = median(vec(chain["rv0"]))
                    jitter = barycentric_rv_inst_1 = median(vec(chain["jitter"]))
                end
                if haskey(chain,:rv0_1)
                    barycentric_rv_inst_1 = median(vec(chain["rv0_1"]))
                    jitter = barycentric_rv_inst_1 = median(vec(chain["jitter_1"]))
                end
                if haskey(chain,:rv0_2)
                    barycentric_rv_inst_2 = median(vec(chain["rv0_2"]))
                    jitter = barycentric_rv_inst_2 = median(vec(chain["jitter_2"]))
                end
                if haskey(chain,:rv0_3)
                    barycentric_rv_inst_3 = median(vec(chain["rv0_3"]))
                    jitter = barycentric_rv_inst_3 = median(vec(chain["jitter_3"]))
                end
                if haskey(chain,:rv0_4)
                    barycentric_rv_inst_4 = median(vec(chain["rv0_4"]))
                    jitter = barycentric_rv_inst_4 = median(vec(chain["jitter_4"]))
                end

                if hasproperty(obs.table, :inst_idx)
                    idxes = length(unique(obs.table.inst_idx))
                else
                    idxes = 1
                end
                y = [[] for _ in 1:idxes]
                x = [[] for _ in 1:idxes]
                yerr = [[] for _ in 1:idxes]
                for row in obs.table
                    if !hasproperty(row,:inst_idx) || row.inst_idx == 1
                        inst_idx = 1
                        barycentric_rv_inst = barycentric_rv_inst_1
                    elseif row.inst_idx == 2
                        inst_idx = 2
                        barycentric_rv_inst = barycentric_rv_inst_2
                    elseif row.inst_idx == 3
                        inst_idx = 3
                        barycentric_rv_inst = barycentric_rv_inst_3
                    elseif row.inst_idx == 4
                        inst_idx = 4
                        barycentric_rv_inst = barycentric_rv_inst_4
                    end
                    push!(x[inst_idx], row.epoch)
                    push!(y[inst_idx], row.rv - barycentric_rv_inst)
                    push!(yerr[inst_idx], row.σ_rv + jitter)
                end
            end
        end
    end
    for i in axes(fit, 1)
        c = (color[i] - first(clims))/(last(clims)-first(clims))
        Plots.plot!(
            mjd2date.(t), view(fit, i, :);
            # line_z=chain[color][ii],
            color=Plots.cgrad(cmap)[c],
            alpha,
            fillalpha=alpha/2,
            ribbon=isnothing(ribbon) ? nothing : ribbon[i],
            kwargs...,
        )
    end
    if !isnothing(y)
        if eltype(y) <: AbstractArray
            for j in eachindex(y)
                Plots.scatter!(
                    p1,
                    mjd2date.(x[j]),
                    y[j];
                    yerr=!isnothing(yerr) ? yerr[j] : nothing,
                    xerr=!isnothing(xerr) ? xerr[j] : nothing,
                    markerstrokewidth=1.5,
                    markerstrokecolor=:auto,
                    markersize=2.5,
                    color=j
                )
            end
        else
            Plots.scatter!(
                p1,
                mjd2date.(x),
                y; yerr, xerr,
                color=1,
                markerstrokewidth=1.5,
                markerstrokecolor=:auto,
                markersize=2.5,
            )
        end
    end
    p1
end

function Octofitter.timeplotgrid(
    model,
    chains;
    color = "$(first(keys(model.system.planets)))_e",
    clims = quantile(vec(chains[color]),(0.01, 0.99)),
    cmap = :plasma,
    alpha=0.1,
    N=1500,
    ii = rand(1:size(chains,1)*size(chains,3), N),
    kwargs...
)
    kwargs = (; 
        ii,
        clims,
        cmap,
        alpha,
        kwargs...
    )
    planet_keys = keys(model.system.planets)

    ppost = Plots.plot()
    for (i,planet_key) in enumerate(planet_keys)
        # plotchains!(ppost, chains, planet_key; color, body=:primary, mass=chains["$(planet_key)_mass"], rev=false, colorbar=nothing, kwargs...)
        plotchains!(ppost, chains, planet_key; color, rev=false, colorbar=nothing, kwargs...)
        astrom = astrometry(model.system.planets[planet_key])
        if !isnothing(astrom)
            Plots.plot!(ppost, astrom, linecolor=1)
        end
        Plots.scatter!([0],[0],marker=(:star, :white, :black, 5),label="")
    end
    psep = Plots.plot()
    for planet_key in planet_keys
        timeplot!(model, chains, planet_key, color, :sep; clims, kwargs...)
    end
    ppa = Plots.plot()
    for planet_key in planet_keys
        timeplot!(model, chains, planet_key, color, :pa; clims, kwargs...)
    end
    
    ppmra = Plots.plot()
    timeplot!(model, chains, planet_keys, color, :pmra;clims,  kwargs...)
    
    ppmdec = Plots.plot()
    timeplot!(model, chains, planet_keys, color, :pmdec;clims,  kwargs...)
    
    prv = Plots.plot()
    timeplot!(model, chains, planet_keys, color, :rv; clims, kwargs...)
    
    pcolormap = Plots.scatter([NaN,NaN],[NaN,NaN],framestyle=:none,foreground=:transparent, markersize=0, marker_z=collect(clims); colorbar=true,color=cmap, colorbartitle=string(color), clims)
    layout = Plots.@layout [
        [
            A B C 
            D E F
        ] G{0.025w}
    ]

    # Horizontal
    Plots.plot(
        ppost,
        psep,
        ppa,
        prv,
        ppmra,
        ppmdec,
        pcolormap,
        layout = layout,
        framestyle=:box,
        grid=false,
        size=(1400,900),
        margin=4Plots.mm
    )
end

# function pmaplot(model, chains; kwargs...)
#     pmaplot(chains, propermotionanom(model.system); kwargs...)
# end
# function pmaplot(model, chains, pma::ProperMotionAnomHGCA; color, N=1500, ii =rand(1:size(chains,1)*size(chains,3), N), kwargs...)
    
#     hgca = pma.table

#     # Roughly over what time period were the observations made?
#     dt_gaia = 1038 # EDR3: days between  Date("2017-05-28") - Date("2014-07-25")
#     dt_hip = 4*365
#     # How many points over Δt should we average the proper motion and stellar position
#     # at each epoch? This is because the PM is not an instantaneous measurement.
#     N_ave = 25

#     ra_gaia_model = zeros(size(ii))
#     dec_gaia_model = zeros(size(ii))
#     pmra_gaia_model = zeros(size(ii))
#     pmdec_gaia_model = zeros(size(ii))
#     ra_hip_model = zeros(size(ii))
#     dec_hip_model = zeros(size(ii))
#     pmra_hip_model = zeros(size(ii))
#     pmdec_hip_model = zeros(size(ii))
#     pmra_hg_model = zeros(size(ii))
#     pmdec_hg_model = zeros(size(ii))

#     planets = keys(chains.info.model.system.planets)

#     for (j,jj) in enumerate(ii)
#         elements = [
#             Octofitter.construct_elements(chains, p, jj)
#             for p in planets
#         ]

#         # First epoch: Hipparcos
#         for i in eachindex(elements)
#             orbit = elements[i]
#             p = planets[i]
#             # Average multiple observations over a timescale +- dt/2
#             # to approximate what HIPPARCOS would have measured.
#             for δt = range(-dt_hip/2, dt_hip/2, N_ave)
#                 # RA and dec epochs are usually slightly different
#                 # Note the unit conversion here from jupiter masses to solar masses to 
#                 # make it the same unit as the stellar mass (element.mu)
#                 # TODO: we can't yet use the orbitsolve interface here for the pmra calls,
#                 # meaning we calculate the orbit 2x as much as we need.
#                 o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_hip[1])+δt)
#                 o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_hip[1])+δt)
#                 ra_hip_model[j] += -raoff(o_ra) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
#                 dec_hip_model[j] += -decoff(o_dec) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
#                 pmra_hip_model[j] += pmra(o_ra, chains["$(p)_mass"][jj]*mjup2msol)
#                 pmdec_hip_model[j] += pmdec(o_dec, chains["$(p)_mass"][jj]*mjup2msol)
#             end
#         end
#         ra_hip_model[j]/=N_ave
#         dec_hip_model[j]/=N_ave
#         pmra_hip_model[j]/=N_ave
#         pmdec_hip_model[j]/=N_ave
    
#         # First epoch: Hipparcos
#         for i in eachindex(elements)
#             orbit = elements[i]
#             p = planets[i]
#             # Average multiple observations over a timescale +- dt/2
#             # to approximate what HIPPARCOS would have measured.
#             for δt = range(-dt_gaia/2, dt_gaia/2, N_ave)
#                 # RA and dec epochs are usually slightly different
#                 # Note the unit conversion here from jupiter masses to solar masses to 
#                 # make it the same unit as the stellar mass (element.mu)
#                 # TODO: we can't yet use the orbitsolve interface here for the pmra calls,
#                 # meaning we calculate the orbit 2x as much as we need.
#                 o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_gaia[1])+δt)
#                 o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_gaia[1])+δt)
#                 ra_gaia_model[j] += -raoff(o_ra) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
#                 dec_gaia_model[j] += -decoff(o_dec) * chains["$(p)_mass"][jj]*mjup2msol/orbit.M
#                 pmra_gaia_model[j] += pmra(o_ra, chains["$(p)_mass"][jj]*mjup2msol)
#                 pmdec_gaia_model[j] += pmdec(o_dec, chains["$(p)_mass"][jj]*mjup2msol)
#             end
#         end
#         ra_gaia_model[j]/=N_ave
#         dec_gaia_model[j]/=N_ave
#         pmra_gaia_model[j]/=N_ave
#         pmdec_gaia_model[j]/=N_ave
    
    
#         # Model the GAIA-Hipparcos delta-position velocity
#         pmra_hg_model[j] = (ra_gaia_model[j] - ra_hip_model[j])/(years2mjd(hgca.epoch_ra_gaia[1]) - years2mjd(hgca.epoch_ra_hip[1]))
#         pmdec_hg_model[j] = (dec_gaia_model[j] - dec_hip_model[j])/(years2mjd(hgca.epoch_dec_gaia[1]) - years2mjd(hgca.epoch_dec_hip[1]))
#     end
#     # # Compute the likelihood at all three epochs (Hipparcos, GAIA-Hip, GAIA)
#     # pmra_model = (pmra_hip_model, pmra_hg_model, pmra_gaia_model)
#     # pmdec_model = (pmdec_hip_model, pmdec_hg_model, pmdec_gaia_model)
#     # pmra_meas = (hgca.pmra_hip[1], hgca.pmra_hg[1], hgca.pmra_gaia[1])
#     # pmdec_meas = (hgca.pmdec_hip[1], hgca.pmdec_hg[1], hgca.pmdec_gaia[1])
#     # σ_pmra = (hgca.pmra_hip_error[1], hgca.pmra_hg_error[1], hgca.pmra_gaia_error[1])
#     # σ_pmdec = (hgca.pmdec_hip_error[1], hgca.pmdec_hg_error[1], hgca.pmdec_gaia_error[1])
    
#     # return Plots.plot(
#     #     Plots.histogram2d(pmra_hip_model, pmdec_hip_model),
#     #     Plots.histogram2d(pmra_hg_model, pmdec_hg_model),
#     #     Plots.histogram2d(pmra_gaia_model, pmdec_gaia_model),
#     # )

#     pma_extrema_x = extrema(vcat(pmra_hip_model, pmra_hg_model, pmra_gaia_model, ))
#     pma_extrema_y = extrema(vcat(pmdec_hip_model, pmdec_hg_model, pmdec_gaia_model, ))

#     if !isnothing(color)
#         subplots = [
#             Plots.scatter(
#                 pmra_hip_model .+ chains["pmra"][ii],
#                 pmdec_hip_model .+ chains["pmdec"][ii],
#                 marker_z=chains[color][ii],
#                 title="Hipparcos",
#             ),
#             Plots.scatter(
#                 pmra_hg_model .+ chains["pmra"][ii],
#                 pmdec_hg_model .+ chains["pmdec"][ii],
#                 marker_z=chains[color][ii],
#                 title="Gaia-Hipparcos",
#             ),
#             Plots.scatter(
#                 pmra_gaia_model .+ chains["pmra"][ii],
#                 pmdec_gaia_model .+ chains["pmdec"][ii],
#                 marker_z=chains[color][ii],
#                 title="Gaia",
#             ),
#         ]
#     else
#         error("Not implemented")
#     end

#     for p in subplots
#         # Plots.scatter!(p, [0], [0], marker=(5, :circle, :red),label="")
#         # Plots.hline!(p, [0], color=:black, label="", linewidth=1)
#         # Plots.vline!(p, [0], color=:black, label="", linewidth=1)
#         Plots.xlabel!(p, "Δμ ra - mas/yr")
#     end
#     Plots.ylabel!(subplots[1], "Δμ dec - mas/yr")


#     pmra_meas = (hgca.pmra_hip[1], hgca.pmra_hg[1], hgca.pmra_gaia[1])
#     pmdec_meas = (hgca.pmdec_hip[1], hgca.pmdec_hg[1], hgca.pmdec_gaia[1])
#     σ_pmra = (hgca.pmra_hip_error[1], hgca.pmra_hg_error[1], hgca.pmra_gaia_error[1])
#     σ_pmdec = (hgca.pmdec_hip_error[1], hgca.pmdec_hg_error[1], hgca.pmdec_gaia_error[1])
#     for i in 1:3
#         Plots.scatter!(
#             subplots[i],
#             [pmra_meas[i]], [pmdec_meas[i]],
#             xerr=[σ_pmra[i]], yerr=[σ_pmdec[i]],
#             markersize=6,
#             markerstrokewidth=3,
#             color=:white
#         )
#         Plots.scatter!(
#             subplots[i],
#             [pmra_meas[i]], [pmdec_meas[i]],
#             xerr=[σ_pmra[i]], yerr=[σ_pmdec[i]],
#             markersize=6,
#             markerstrokewidth=1,
#             color=:black
#         )
#     end

#     pma_plot =  Plots.plot(
#         subplots...,
#         alpha=1,
#         legend=false,
#         colorbar=false,
#         label="",
#         markerstrokewidth=0,
#         ms=1,
#         aspectratio=1,
#         layout=(1,3),
#         framestyle=:box,
#         minorticks=true,
#         grid=false,
#         # xlims=pma_extrema_x,
#         # ylims=pma_extrema_y,
#         # clims=clims_pma
#     )

#     return pma_plot


# end
end
