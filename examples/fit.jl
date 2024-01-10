using Octofitter
using CairoMakie
using PairPlots
using Distributions
using Plots:Plots
using PlanetOrbits
##
orb_template = orbit(
    a = 1.0,
    e = 0.1,
    i= pi/4,
    Ω = 0.1,
    ω = 1π/4,
    M = 1.0,
    plx=100.0,
    m=0,
    tp =58849
)
Plots.plot(orb_template)
##

astrom = PlanetRelAstromLikelihood(
    [(; epoch=epoch+58849,
    ra=raoff(orb_template, epoch) ,
    dec=decoff(orb_template, epoch) ,
    σ_ra = 0.2,
    σ_dec = 0.2,
    )
    for epoch in [20:8:45;]]...
)
Plots.plot(orb_template, aspectratio=1, lw=0, label="")
Plots.plot!(astrom, aspectratio=1, framestyle=:box)
##
@planet b Visual{KepOrbit} begin
    e ~ Uniform(0,0.9)
    a ~ LogUniform(0.01, 50)
    mass ~ truncated(Normal(1, 1), lower=0)
    i ~ Sine()
    Ω ~ UniformCircular()
    ω ~ UniformCircular()

    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P +  58849
end astrom


# @planet b Visual{CartesianOrbit} begin
#     x ~ Normal(0, 1)
#     y ~ Normal(0, 1)
#     z ~ Normal(0, 1)
#     vx ~ Normal(0, 10)
#     vy ~ Normal(0, 10)
#     vz ~ Normal(0, 10)
# end astrom


@system test begin
    M ~ truncated(Normal(1, 0.04),lower=0) # (Baines & Armstrong 2011).
    plx = 100.0
end b

model = Octofitter.LogDensityModel(test; autodiff=:ForwardDiff,verbosity=4)
# model = Octofitter.LogDensityModel(test; autodiff=:FiniteDiff,verbosity=4)

##
using Random
Random.seed!(2)

results = octofit(model)#, adaptation=1500, iterations=50000)
plotchains(results, :b, color=:b_e); Plots.plot!(astrom)
##

octocorner(model, results,)
octoplot(model, results)


##



# Fit 
@time results = octofit(model)
display(results)

# Plot
p = Octofitter.plotchains(results, :b, color=:b_e)
Plots.plot!(p, astrom)
display(p)



# Plots.plot(results[:b_a])



tbl = (;
    M=totalmass.(els),
    b_i=inclination.(els),
    b_e=eccentricity.(els),
    b_ω=results["b_ω"][:],
    b_Ω=results["b_Ω"][:],
    b_τ=results["b_τ"][:]
)
corner_theme = 
    PairPlots.Hist(bins=15),
    PairPlots.Scatter(markersize=3,filtersigma=3,),
    PairPlots.MarginStepHist(bins=15),
    PairPlots.Contour(linewidth=2, sigmas=1:3, color=:black),
    PairPlots.MarginConfidenceLimits()

pairplot(
    tbl =>corner_theme,
    PairPlots.Truth((;
        M=totalmass(orb_template),
        b_i=inclination(orb_template),
        b_e=eccentricity(orb_template),
        b_ω=orb_template.parent.ω,
        b_Ω=orb_template.parent.Ω,    
    ))
)