# Required packages:
# Octofitter, Distributions, Plots, CairoMakie, PairPlots, Random
# See Project and Manifest files.

using Octofitter, Distributions, PlanetOrbits
using OctofitterWhereistheplanet
cd(@__DIR__)

## Define model

# Note: for more information, see docs https://sefffal.github.io/Octofitter.jl/dev/


# Reference epoch for tau (orbital position) in convention used by Orbitize.
const tref = 58849.0
const theta_epoch = mean(OctofitterWhereistheplanet.astrom("hr8799";object=4)[2].table.epoch)


using Octofitter.StaticArrays
function τ_from_θ(system,planet)
    (;plx,M) = system
    (;B,G,A,F,e,θ) = planet
    # Thiele-Innes transformation matrix
    T = @SArray [
        A F
        B G
    ]
    x_over_r, y_over_r = T \ @SArray[
        cos(θ)
        sin(θ)
    ]
    # Note: this is missing the radius factor but we don't need it for true anomaly.
    # r = sqrt((sol.x*sol.elem.B+sol.y*sol.elem.G)^2 + (sol.x*sol.elem.A + sol.y*sol.elem.F)^2 )

    # True anomaly is now straightforward to calculate in the deprojected coordinate space
    ν = atan(y_over_r,x_over_r)

    # Mean anomaly (see Wikipedia page)
    MA = atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π-e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν))

    # Math in order to get semi-major axis -> mean motion and period (already in TI constructor)
    u = (A^2 + B^2 + F^2 + G^2)/2
    v = A*G - B * F
    α = sqrt(u + sqrt((u+v)*(u-v)))
    a = α/plx
    periodyrs = √(a^3/M)
    perioddays = periodyrs * PlanetOrbits.year2day # period [days]
    n = 2π/√(a^3/M) 
    # Get epoch of periastron passage
    tₚ = theta_epoch - MA/n*PlanetOrbits.year2day - tref
    # Tau: periastron passage epoch / orbital period
    τ = rem(tₚ/perioddays,1,RoundDown)
    return τ
end

struct SmaLogUni <: Octofitter.AbstractLikelihood end
function Octofitter.ln_like(::SmaLogUni, θ_planet, orbit) 
    # Convert back to Campbell elements to maintain a prior on semi-major axis.
    campbell_orbit = Visual{KepOrbit}(orbit)
    return logpdf(LogUniform(0.001, 1000), campbell_orbit.a)
end

@planet b ThieleInnesOrbit begin
    e ~ Uniform(0.0, 1.0)

    # The adaptation phase will evidently perform better 
    # if working on parameter values inside [-2,2]
    A_scaled ~ Normal(0, 1)
    B_scaled ~ Normal(0, 1)
    F_scaled ~ Normal(0, 1)
    G_scaled ~ Normal(0, 1)
    A = 3000 * b.A_scaled
    B = 3000 * b.B_scaled
    F = 3000 * b.F_scaled
    G = 3000 * b.G_scaled
    tref = tref

    θ ~ UniformCircular()
    τ = τ_from_θ(system, b)

    mass = 5

end OctofitterWhereistheplanet.astrom("hr8799";object=1)[2]


@planet c ThieleInnesOrbit begin
    e ~ Uniform(0.0, 1.0)

    # The adaptation phase will evidently perform better 
    # if working on parameter values inside [-2,2]
    A_scaled ~ Normal(0, 1)
    B_scaled ~ Normal(0, 1)
    F_scaled ~ Normal(0, 1)
    G_scaled ~ Normal(0, 1)
    A = 3000 * c.A_scaled
    B = 3000 * c.B_scaled
    F = 3000 * c.F_scaled
    G = 3000 * c.G_scaled
    tref = tref

    θ ~ UniformCircular()
    τ = τ_from_θ(system, c)

    mass = 7

end OctofitterWhereistheplanet.astrom("hr8799";object=2)[2]


@planet d ThieleInnesOrbit begin
    e ~ Uniform(0.0, 1.0)

    # The adaptation phase will evidently perform better 
    # if working on parameter values inside [-2,2]
    A_scaled ~ Normal(0, 1)
    B_scaled ~ Normal(0, 1)
    F_scaled ~ Normal(0, 1)
    G_scaled ~ Normal(0, 1)
    A = 3000 * d.A_scaled
    B = 3000 * d.B_scaled
    F = 3000 * d.F_scaled
    G = 3000 * d.G_scaled
    tref = tref

    θ ~ UniformCircular()
    τ = τ_from_θ(system, d)
    mass = 7


end OctofitterWhereistheplanet.astrom("hr8799";object=3)[2]


@planet e ThieleInnesOrbit begin
    e ~ Uniform(0.0, 1.0)

    # The adaptation phase will evidently perform better 
    # if working on parameter values inside [-2,2]
    # A_scaled ~ Normal(0, 1)
    # B_scaled ~ Normal(0, 1)
    # F_scaled ~ Normal(0, 1)
    # G_scaled ~ Normal(0, 1)

    # A = 3000 * e.A_scaled
    # B = 3000 * e.B_scaled
    # F = 3000 * e.F_scaled
    # G = 3000 * e.G_scaled

    AFr ~ LogUniform(10, 10_000)
    AFt ~ UniformCircular()
    BGr ~ LogUniform(10, 10_000)
    BGt ~ UniformCircular()
    A = e.AFr * cos(e.AFt)
    F = e.AFr * sin(e.AFt)
    B = e.BGr * cos(e.BGt)
    G = e.BGr * sin(e.BGt)
    
    tref = tref

    θ ~ UniformCircular()
    τ = τ_from_θ(system, e)

    mass = 7

end OctofitterWhereistheplanet.astrom("hr8799";object=4)[2]



# Simple model for the star (with one orbitting planet, b).
# Note this is not a nested model, just a nice way of structuring multi-planet 
# models in the code.
@system HR8799 begin
    # Parallax distance to the star (milliarcseconds of apparent motion per Earth year)
    # Sets overall scale of the model.
    plx ~ gaia_plx(;gaia_id=2832463659640297472)
    # Mass of the star in solar masses: determines orbital speed.
    M   ~ truncated(Normal(1.5, 0.2),lower=0)

    # True proper motion of the system barycentre.
    # The variance on this will be quite low, but let the mean be determined
    # entirely from the data.
    pmra  = 0 #~ Normal(0, 500)
    pmdec = 0 #~ Normal(0, 500)
# end b c d e HGCALikelihood(;gaia_id=2832463659640297472 )
end e


## Generate custom LogDensityModel
model = Octofitter.LogDensityModel(HR8799; autodiff=:ForwardDiff, verbosity=4)


## Sample using AdvancedHMC
using Random
rng = Random.Xoshiro(1)


# Code is in sampling.jl beginning at line 629.
chain = Octofitter.advancedhmc(
    rng,
    model, 
    0.85; # Target acceptance. Cranking this up to 0.95+ will reduce the number of divergences to a narrow region.
    num_chains = 1,
    # Number of adaptation steps, and iterations after adaption.
    # Would recommend 20,000 but 2,000 are enough for quick tests.
    adaptation = 1_000, #, 20_000, 
    iterations = 1_000, #, 20_000, 
    verbosity = 4,
    tree_depth = 14,

    # Can pass initial parameters for testing but the default initialization 
    # proceedure is quite reliable.
    # initial_parameters
)


##
using Plots: Plots
# Have to change parameter names to the name of the planet(s) used
p = Plots.plot()
N = 500
# Octofitter.plotchains!(p, chain, :b; kind=:astrometry, color=:b_e, N, framestyle=:box)
# Octofitter.plotchains!(p, chain, :c; kind=:astrometry, color=:c_e, N, framestyle=:box)
# Octofitter.plotchains!(p, chain, :d; kind=:astrometry, color=:d_e, N, framestyle=:box)
Octofitter.plotchains!(p, chain, :e; kind=:astrometry, color=:e_e, N, framestyle=:box)
Plots.plot!(model.system.planets.e.observations[4],xlims=:symmetric,ylims=:symmetric)
##
octoplot(model, chain)
## Corner plot
using CairoMakie
using PairPlots
fig = pairplot(chain => (PairPlots.Scatter(),))
save(joinpath(@__DIR__,"octofitter-pairplot.png"), fig)


## Diagnostics
num_errs = map(samp->samp.stat.numerical_error, chain.info.states[1])

# Corner plot with divergences highlighted in red
using CairoMakie
using PairPlots

# marker size: tweak smaller for runs with many iterations
ms = 1.5

c = Makie.wong_colors()
colnames = [
    :plx, :M, :b_e, :b_A, :b_B, :b_F, :b_G, :b_θ, :b_θx, :b_θy
]
selectcols(chain) = (;
    iter=collect(chain.value.axes[1]),
    M=vec(chain[:M]),
    b_A=vec(chain[:e_A]),
    b_B=vec(chain[:e_B]),
    b_F=vec(chain[:e_F]),
    b_G=vec(chain[:e_G]),
    b_θ=vec(chain[:e_θ]),
)
fig = pairplot(
    PairPlots.Series(selectcols(chain[.! num_errs,:,:]),label="sample") => (PairPlots.Scatter(color=c[1],markersize=ms),PairPlots.MarginDensity(color=(c[1],0.4),linewidth=1)),
    PairPlots.Series(selectcols(chain[num_errs,:,:]),label="divergence") => (PairPlots.Scatter(color=:red,markersize=ms),PairPlots.MarginDensity(color=(:red,0.4),linewidth=1)),
)
save(joinpath(@__DIR__, "pairplot-divergences-2.png"), fig, px_per_unit=2)



## Diagnostics
num_errs = map(samp->samp.stat.numerical_error, chain.info.states[1])

# Corner plot with divergences highlighted in red
using CairoMakie
using PairPlots

# marker size: tweak smaller for runs with many iterations
ms = 1.5

c = Makie.wong_colors()
colnames = [
    :plx, :M, :b_e, :b_A, :b_B, :b_F, :b_G, :b_θ, :b_θx, :b_θy
]
selectcols(chain) = (;
    iter=collect(chain.value.axes[1]),
    M=vec(chain[:M]),
    b_e=vec(chain[:e_e]),
    b_A=vec(chain[:e_A]),
    b_B=vec(chain[:e_B]),
    b_F=vec(chain[:e_F]),
    b_G=vec(chain[:e_G]),
    b_θ=vec(chain[:e_θ]),
)
fig = pairplot(
    PairPlots.Series(selectcols(chain[.! num_errs,:,:]),label="sample") => (PairPlots.Scatter(color=c[1],markersize=ms),PairPlots.MarginDensity(color=(c[1],0.4),linewidth=1)),
    PairPlots.Series(selectcols(chain[num_errs,:,:]),label="divergence") => (PairPlots.Scatter(color=:red,markersize=ms),PairPlots.MarginDensity(color=(:red,0.4),linewidth=1)),
)
save(joinpath(@__DIR__, "pairplot-divergences-2.png"), fig, px_per_unit=2)

