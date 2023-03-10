# Required packages:
# Octofitter, Distributions, Plots, CairoMakie, PairPlots, Random
# See Project and Manifest files.

using Octofitter, Distributions, PlanetOrbits
using OctofitterWhereistheplanet
cd(@__DIR__)

## Define model

# Note: for more information, see docs https://sefffal.github.io/Octofitter.jl/dev/


# Reference epoch for orbital position.
const tref2 = 55500.0


using Octofitter.StaticArrays
function τ_from_θ(system,planet)
    (;plx,M) = system
    (;B_scaled,G_scaled,A_scaled,F_scaled,e,θx,θy) = planet
    B = 2000*B_scaled
    G = 2000*G_scaled
    A = 2000*A_scaled
    F = 2000*F_scaled
    θ = atan(θy,θx)
    t=tref=tref2
    # Thiele-Innes transformation matrix
    T = @SArray [
        A F
        B G
    ]
    xx, yy = T \ @SArray[
        cos(θ)
        sin(θ)
    ] # Note: this is missing the radius factor but we don't need it for true anomaly.
    # r = sqrt((sol.x*sol.elem.B+sol.y*sol.elem.G)^2 + (sol.x*sol.elem.A + sol.y*sol.elem.F)^2 )

    # True anomaly is now straightforward to calculate in the deprojected coordinate space
    ν = atan(yy,xx)

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
    tₚ = t - MA/n*PlanetOrbits.year2day - tref
    # Tau: periastron passage epoch / orbital period
    τ = rem(tₚ/perioddays,1,RoundDown)
    return τ
end

struct CustomLikelihood1 <: Octofitter.AbstractLikelihood end
function Octofitter.ln_like(::CustomLikelihood1, θ_planet, orbit) 
    # Convert back to Campbell elements to maintain a prior on semi-major axis.
    campbell_orbit = VisualOrbit(orbit)
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
    tref = tref1

    θ ~ UniformCircular()
    τ = τ_from_θ(system, b)

end CustomLikelihood1() OctofitterWhereistheplanet.astrom("hr8799";object=1)[2]


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
    tref = tref1

    θ ~ UniformCircular()
    τ = τ_from_θ(system, c)

end CustomLikelihood1() OctofitterWhereistheplanet.astrom("hr8799";object=2)[2]


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
    tref = tref1

    θ ~ UniformCircular()
    τ = τ_from_θ(system, d)

end CustomLikelihood1() OctofitterWhereistheplanet.astrom("hr8799";object=3)[2]


@planet e ThieleInnesOrbit begin
    e ~ Uniform(0.0, 1.0)

    # The adaptation phase will evidently perform better 
    # if working on parameter values inside [-2,2]
    A_scaled ~ Normal(0, 1)
    B_scaled ~ Normal(0, 1)
    F_scaled ~ Normal(0, 1)
    G_scaled ~ Normal(0, 1)
    A = 3000 * e.A_scaled
    B = 3000 * e.B_scaled
    F = 3000 * e.F_scaled
    G = 3000 * e.G_scaled
    tref = tref1

    θ ~ UniformCircular()
    τ = τ_from_θ(system, e)

end CustomLikelihood1() OctofitterWhereistheplanet.astrom("hr8799";object=4)[2]



# Simple model for the star (with one orbitting planet, b).
# Note this is not a nested model, just a nice way of structuring multi-planet 
# models in the code.
@system HR8799 begin
    # Parallax distance to the star (milliarcseconds of apparent motion per Earth year)
    # Sets overall scale of the model.
    plx_scale ~ Normal(0, 1.0)
    plx = 24.461952 + 0.045453273 * system.plx_scale
    # Mass of the star in solar masses: determines orbital speed.
    M_scale   ~ Normal(0, 1)
    M   = 1.5 + 0.2system.M_scale
end b


## Generate custom LogDensityModel
model = Octofitter.LogDensityModel(HR8799; autodiff=:ForwardDiff, verbosity=4)


## Sample using AdvancedHMC
using Random
rng = Random.Xoshiro(1)


# Code is in sampling.jl beginning at line 629.
chain = Octofitter.advancedhmc(
    rng,
    model, 
    0.9; # Target acceptance. Cranking this up to 0.95+ will reduce the number of divergences to a narrow region.
    num_chains = 1,
    # Number of adaptation steps, and iterations after adaption.
    # Would recommend 20,000 but 2,000 are enough for quick tests.
    adaptation = 2_000, #, 20_000, 
    iterations = 2_000, #, 20_000, 
    verbosity = 4,
    tree_depth = 14,

    # Can pass initial parameters for testing but the default initialization 
    # proceedure is quite reliable.
    # initial_parameters
)



# Optional: save and restore chains
# Octofitter.savechain(joinpath(@__DIR__, "chains.fits"), chain);
# Octofitter.loadchain(joinpath(@__DIR__, "chains.fits"), chain);

##  Standard model output: plot orbits against data
using Plots: Plots
octoplot(model, chain)


## Corner plot
using CairoMakie
using PairPlots
fig = pairplot(chain)
save(joinpath(@__DIR__,"octofitter-pairplot.png"), fig)


## Diagnostics
tree_depths = map(samp->samp.stat.tree_depth, chain.info.states[1])

## Corner plot with divergences highlighted in red
using CairoMakie
using PairPlots

# marker size: tweak smaller for runs with many iterations
ms = 1.5

c = Makie.wong_colors()
colnames = [
    :plx, :M, :b_e, :b_A, :b_B, :b_F, :b_G, :b_τ, :b_τx, :b_τy
]
fig = pairplot(
    PairPlots.Series(chain[.! num_err,colnames,:],label="sample") => (PairPlots.Scatter(color=c[1],markersize=ms),PairPlots.MarginDensity(color=(c[1],0.4),linewidth=1)),
    PairPlots.Series(chain[num_err,colnames,:],label="divergence") => (PairPlots.Scatter(color=:red,markersize=ms),PairPlots.MarginDensity(color=(:red,0.4),linewidth=1)),
)
save(joinpath(@__DIR__, "pairplot-divergences.png"), fig, px_per_unit=2)

