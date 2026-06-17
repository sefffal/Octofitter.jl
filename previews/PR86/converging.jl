
using Octofitter
using OctofitterRadialVelocity
using CairoMakie
using PairPlots
using Distributions
using PlanetOrbits

planet_sim_mass = 0.001 # solar masses here


orb_template = orbit(
    a = 1.0,
    e = 0.4,
    # i= pi/4, # You can remove I think
    # Ω = 0.1, # You can remove I think
    ω = 1π/4, # radians
    M = 1.0, # Total mass, not stellar mass FYI
    plx=100.0,
    tp =58829 # Epoch of periastron passage. 
)
# Makie.lines(orb_template)

posts = []
models = []
##
for n_epochs = 2:20
    epochs = 58849 .+ range(0,length=n_epochs, step=50)
    rvlike = MarginalizedStarAbsoluteRVLikelihood(
        Table(
            epoch=epochs,
            rv=radvel.(orb_template, epochs, planet_sim_mass),
            σ_rv=fill(5.0, size(epochs)),
        ),
        jitter=:j1
    )
    # scatter(rvlike.table.epoch, rvlike.table.rv)|>display

    first_epoch_for_tp_ref = first(epochs)
    @planet b RadialVelocityOrbit begin
        e ~ Uniform(0,0.9)
        a ~ Uniform(0, 5)
        mass ~ Uniform(0, 10)

        ω ~ Uniform(0,2pi)

        # τ ~ UniformCircular(1.0)
        τ ~ Uniform(0.0, 1.0)
        tp =  b.τ*√(b.a^3/system.M)*365.25 + $first_epoch_for_tp_ref # reference epoch for τ. Choose to be near data
    end 

    @system SimualtedSystem begin
        M ~ truncated(Normal(1, 0.04),lower=0.3, upper=3) # (Baines & Armstrong 2011).
        plx = 100.0
        j1 ~ LogUniform(0.1, 100)
    end rvlike b

    model = Octofitter.LogDensityModel(SimualtedSystem)

    # chain,pt = octofit_pigeons(model, explorer=SliceSampler(), n_rounds=9, variational=nothing, n_chains_variational=0, n_chains=8)
    chain,pt = octofit_pigeons(model, explorer=SliceSampler(), n_rounds=9,  n_chains_variational=8, n_chains=8)
    push!(posts, chain)
    push!(models, model)
    # display(chain)
end
##
N_epochs = 2:20#(50-length(posts)+1):50
meds = map(posts) do chn
    quantile(vec(chn[:b_mass]), 0.5)
end
low = map(posts) do chn
    quantile(vec(chn[:b_mass]), 0.14)
end
high = map(posts) do chn
    quantile(vec(chn[:b_mass]), 0.86)
end
fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="N epochs",
    ylabel="Mass"
)

band!(ax, N_epochs, low, high, color=(:gray, 0.5))
scatterlines!(ax, N_epochs, meds, color=:black)
# hlines!(ax, eccentricity(orb_template))
hlines!(ax, planet_sim_mass/Octofitter.mjup2msol)

fig