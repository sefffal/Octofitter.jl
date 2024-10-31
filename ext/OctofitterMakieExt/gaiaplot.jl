#=

Want to display, for a given sample, the posterior
from GAIA and the fitted posterior.
=#
using DataFrames
samples_dr3 = DataFrame(
    rand(gaialike.dist_dr3, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr3.ra .-= gaialike.dr3.ra
samples_dr3.dec .-= gaialike.dr3.dec

samples_dr2 = DataFrame(
    rand(gaialike.dist_dr2, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr2.ra .-= gaialike.dr2.ra
samples_dr2.dec .-= gaialike.dr2.dec


# θ_systems_from_chain = Octofitter.mcmcchain2result(model, chain)
# i = argmax(chain[:logpost][:])
# i = rand(1:size(chain,1))
# orbits = map(keys(model.system.planets)) do planet_key
#     Octofitter.construct_elements(chain, planet_key, i)
# end
# P = Octofitter.simulate(gaialike, θ_systems_from_chain[i], orbits)

post_dr3 = DataFrame(
    rand(out.P_dr3, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr3.ra .-= gaialike.dr3.ra
post_dr3.dec .-= gaialike.dr3.dec

post_dr2 = DataFrame(
    rand(out.P_dr2, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr2.ra .-= gaialike.dr2.ra
post_dr2.dec .-= gaialike.dr2.dec


fig = Figure(size=(700,700))

pairplot(
    fig[1:5,1:5],
    PairPlots.Series(samples_dr3,label="DR3",color=:blue)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(post_dr3,label="DR3 Model",color=:blue,linestyle=:dash)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(samples_dr2,label="DR2",color=:red)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5,),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(post_dr2,label="DR2 Model",color=:red,linestyle=:dash)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    )
)
fig
##
ax = Axis(
    fig[1:2,4:5],
    xlabel="[mas]",
    ylabel="[mas]",
    backgroundcolor="#EEE",
    autolimitaspect=1
)

# scatterlines!(ax,gaialike.table.epoch, resids)
N = length(gaialike.table.epoch)
resids = out.sol.u[1:N]
uncs = fill(out.sol.u[N+1],N)
for i in eachindex(resids)
    s,c = sincos(gaialike.table.var" scanAngle[rad]"[i])
    lines!(ax,
        [resids[i] .* c - uncs[i].*c, resids[i] .* c + uncs[i].*c],
        [resids[i] .* s - uncs[i].*s, resids[i] .* s + uncs[i].*s],
        color=:black
        # 
    )
end
p = scatter!(ax,
    vec(resids .* cos.(gaialike.table.var" scanAngle[rad]")),
    vec(resids .* sin.(gaialike.table.var" scanAngle[rad]")),
    color=vec(gaialike.table.epoch),
    colormap=:turbo
)
Colorbar(fig[1,3], p, label="MJD")
fig

##
fig = Figure()
axl = Axis(
    fig[1,1],
    xlabel="MJD",ylabel="along scan residual [mas]",
    # yticklabelcolor=Makie.wong_colors()[1],
    # ylabelcolor=Makie.wong_colors()[1],
    # leftspinecolor=Makie.wong_colors()[1],
)
axr = Axis(
    fig[2,1],
    xlabel="MJD",
    ylabel="scan angle [deg]",
    # yticklabelcolor=Makie.wong_colors()[2],
    # ylabelcolor=Makie.wong_colors()[2],
    # rightspinecolor=Makie.wong_colors()[2],
)
N = length(gaialike.table.epoch)
resids = out.sol.u[1:N]
# uncs = out.sol.u[N+1:2N]
# astrometric_excess_noise_dr3 = out.sol.u[2N+1]
# astrometric_excess_noise_dr2 = out.sol.u[2N+2]
# mask = abs.(uncs) .> 10
# uncs[mask] .= NaN
# resids[mask] .= NaN
uncs = fill(out.sol.u[N+1],N)
astrometric_excess_noise_dr3 = gaialike.dr3.astrometric_excess_noise
astrometric_excess_noise_dr2 = gaialike.dr2.astrometric_excess_noise
scatterlines!(axl,gaialike.table.epoch, resids)
# ylims!(axl,-100,100)
iend = findlast(<=(Octofitter.meta_gaia_DR2.stop_mjd), gaialike.table.epoch)
errorbars!(axl, gaialike.table.epoch[begin:iend], resids[begin:iend], sqrt.(uncs.^2 .+ astrometric_excess_noise_dr2.^2)[begin:iend], color=:grey, linewidth=4)
errorbars!(axl, gaialike.table.epoch, resids, sqrt.(uncs.^2 .+ astrometric_excess_noise_dr3.^2), color=:black, linewidth=1)
# scatterlines!(axr,gaialike.table.epoch, rad2deg.(gaialike.table.var" scanAngle[rad]"),color=Makie.wong_colors()[2])
scatterlines!(axr,gaialike.table.epoch, rad2deg.(rem2pi.(gaialike.table.var" scanAngle[rad]",RoundDown)),color=Makie.wong_colors()[2])
linkxaxes!(axl,axr)
Makie.rowsize!(fig.layout, 2, Auto(0.5))
fig
##

fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="[mas]",
    ylabel="[mas]",
)

# scatterlines!(ax,gaialike.table.epoch, resids)
for i in eachindex(resids)
    s,c = sincos(gaialike.table.var" scanAngle[rad]"[i])
    lines!(ax,
        [resids[i] .* c - uncs[i].*c, resids[i] .* c + uncs[i].*c],
        [resids[i] .* s - uncs[i].*s, resids[i] .* s + uncs[i].*s],
        color=:black
        # 
    )
end
p = scatter!(ax,
    resids .* cos.(gaialike.table.var" scanAngle[rad]"),
    resids .* sin.(gaialike.table.var" scanAngle[rad]"),
    color=gaialike.table.epoch,
    colormap=:turbo
)
Colorbar(fig[1,2], p, label="MJD")
fig

## Al
fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="MJD",
    ylabel="along scan residual [mas]"
)
s = sin.(gaialike.table.var" scanAngle[rad]")
c = cos.(gaialike.table.var" scanAngle[rad]")
scatterlines!(ax, gaialike.table.epoch, resids .* c)
errorbars!(ax, gaialike.table.epoch, resids .* c, uncs .* c)
scatterlines!(ax, gaialike.table.epoch, resids .* s)
errorbars!(ax, gaialike.table.epoch, resids .* s, uncs .* s)

fig
# Next up need to add in astrometric excess noise separately per mission,
# and something to regularize it

##
xs,ys,_= Octofitter._simulate_skypath_observations(gaialike, orbits[1], chain[i,:"b_mass",1])
ax = Axis(
    fig[1,3],
    xreversed=true,
    autolimitaspect=1
)
scatterlines!(ax,xs,ys)
fig
##
octoplot(model,chain[i,:,:])