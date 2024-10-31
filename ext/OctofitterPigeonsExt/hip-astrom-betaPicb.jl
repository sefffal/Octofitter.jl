using Octofitter
using Distributions
using Pigeons
using CairoMakie
##
hip_like = Octofitter.HipparcosIADLikelihood(;hip_id=27321)
##
astrom_like1,astrom_like2 = Octofitter.Whereistheplanet_astrom("betpic",object=1)

# Visualize the astrometry data
# fig,ax,pl=scatter(astrom_like2.table.ra,astrom_like2.table.dec)
# x = astrom_like1.table.sep .* sin.(astrom_like1.table.pa)
# y = astrom_like1.table.sep .* cos.(astrom_like1.table.pa)
# scatter!(ax, x,y)
# fig
##

@planet B AbsoluteVisual{KepOrbit} begin
    a ~ LogUniform(5, 20)
    e ~ Uniform(0,0.9)
    ω ~ Uniform(0,2pi)
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)# = 0 #~ UniformCircular()
    mass = system.M_sec
    θ ~ Uniform(0,2pi)#UniformCircular()
    tp = θ_at_epoch_to_tperi(system,B,57423.0) # epoch of GAIA measurement
end astrom_like1 astrom_like2
const bp_radeg = hip_like.hip_sol.radeg
const bp_dedeg = hip_like.hip_sol.dedeg
@system BetaPic begin
    M_pri ~ truncated(Normal(1.75,0.05), lower=0.03) # Msol
    M_sec ~ Uniform(0, 100) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol

    # Position and projected velocities of barycentre 
    # at a given reference epoch 
    # at a given reference epoch 
    plx ~ truncated(Normal(50.9307, 5 ), lower=10)
    ref_epoch = years2mjd(1991.25)
    # ref_epoch = years2mjd(2016.25)
    rv = 20.0e3
    ra ~ Normal(bp_radeg, 0.1)
    dec ~ Normal(bp_dedeg, 0.1)
    # ref_epoch = years2mjd(2016.15)
            
    # Priors on the center of mass proper motion
    pmra ~ Normal(4.65, 20.11)
    pmdec ~ Normal(83.1,  20.15)#~ Normal(83.1,  0.15)
end hip_like B
model = Octofitter.LogDensityModel(BetaPic)

##
chain,pt = octofit_pigeons(model,n_rounds=7, n_chains=30,n_chains_variational=20)
##
octoplot(model,chain)
## HMC actually performs a lot better on this example
chainhmc = octofit(model)
##
octoplot(model,chainhmc)
##

fig = Figure()

ax = Axis(fig[1,1:2],
    xlabel="α* [mas]",
    ylabel="δ [mas]",
)


scatterlines!(hip_like.table.Δα✱, hip_like.table.Δδ,label="Hipparcos model", markersize=2)

ax2 = Axis(
    fig[2,1:2],
    xlabel="epoch [MJD]",
    ylabel="along-scan residual [mas]"
)
errorbars!(ax2, hip_like.table.epoch,  zeros(size(hip_like.table.epoch)), hip_like.table.sres_renorm, color=:black)

nts = Octofitter.mcmcchain2result(model, chain,)

# ii = rand(axes(chain),1)
_, ii = findmax(chain[:logpost][:])
for i in ii
    orbC = Octofitter.construct_elements(chain, :B, i)
    sim = Octofitter.simulate(hip_like, nts[i,], [orbC])
    # Model
    scatterlines!(ax,
        sim.α✱_model_with_perturbation[:,1],
        sim.δ_model_with_perturbation[:,1],
        # label="Our Model",
        color=Makie.wong_colors()[2],
        alpha=1
    )

    # Data
    for i in axes(hip_like.table.α✱ₘ,1)
        # @show point1 point2
        lines!(ax, hip_like.table.α✱ₘ[i][1:2], hip_like.table.δₘ[i][1:2],color=:black)
    end

    resid = map(eachindex(hip_like.table.epoch)) do i
        # point = α✱ₘ
        point =  [
            sim.α✱_model_with_perturbation[i],
            sim.δ_model_with_perturbation[i]
        ]
        line_point_1 =  [hip_like.table.α✱ₘ[i][1], hip_like.table.δₘ[i][1]]
        line_point_2 =  [hip_like.table.α✱ₘ[i][2], hip_like.table.δₘ[i][2]]
        Octofitter.distance_point_to_line(point, line_point_1, line_point_2)
    end

  
    scatter!(ax2, hip_like.table.epoch, resid, markersize=3, color=Makie.wong_colors()[2], alpha=1)
end
Legend(fig[1,3],ax)
ylims!(ax2, low=0)
fig
##
using PairPlots
pairplot((;mass=chain[:B_mass][:]))