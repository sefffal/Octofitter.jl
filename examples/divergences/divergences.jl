# Required packages:
# Octofitter, Distributions, Plots, CairoMakie, PairPlots, Random
# See Project and Manifest files.

using Octofitter, Distributions
cd(@__DIR__)

## Define model

# Note: for more information, see docs https://sefffal.github.io/Octofitter.jl/dev/

# Data:
# Position of the planet relative to the star measured at 27 epochs.
# This is a simulated but representative dataset.

# Octofitter also supports radial velocity of the star, stellar motion, etc.
# but relative astrometry is the simplest case.
astrom = Astrometry(
    (epoch=58849.0,  pa=-2.33643,  sep=615.194, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=58879.0,  pa=-2.30229,  sep=606.432, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=58909.0,  pa=-2.26984,  sep=663.252, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=58939.0,  pa=-2.23892,  sep=635.889, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=58969.0,  pa=-2.20938,  sep=610.268, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=58999.0,  pa=-2.1811 ,  sep=669.674, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59029.0,  pa=-2.15397,  sep=666.777, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59059.0,  pa=-2.12789,  sep=654.323, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59089.0,  pa=-2.10278,  sep=713.583, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59215.0,  pa=-2.00625,  sep=747.252, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59245.0,  pa=-1.98508,  sep=736.368, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59275.0,  pa=-1.96451,  sep=709.685, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59305.0,  pa=-1.9445 ,  sep=791.394, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59335.0,  pa=-1.92502,  sep=777.824, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59365.0,  pa=-1.90602,  sep=773.808, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59395.0,  pa=-1.88748,  sep=866.536, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59425.0,  pa=-1.86938,  sep=789.465, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59455.0,  pa=-1.85167,  sep=839.135, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59945.0,  pa=-1.60397,  sep=986.338, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=59975.0,  pa=-1.59064,  sep=941.645, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60005.0,  pa=-1.57747,  sep=959.192, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60035.0,  pa=-1.56444,  sep=928.867, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60065.0,  pa=-1.55155,  sep=952.933, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60095.0,  pa=-1.53879,  sep=977.518, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60125.0,  pa=-1.52616,  sep=950.104, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60155.0,  pa=-1.51365,  sep=953.884, σ_sep=30.0, σ_pa=0.0139626, cor=0),
    (epoch=60185.0,  pa=-1.50125,  sep=985.232, σ_sep=30.0, σ_pa=0.0139626, cor=0),
)

# Reference epoch for orbital position.
const tref1 = mean(astrom.table.epoch)


# Define the planet `b` model
@planet b ThieleInnesOrbit astrom begin
    e ~ Uniform(0.0, 1.0)

    # A, B, F, and G are Thiele-Innes constants that represent
    # a Keplerian orbit. They are in units of milliarcseconds, or
    # effectively the apparent position.
    A ~ Normal(0, 2000)
    B ~ Normal(0, 2000)
    F ~ Normal(0, 2000)
    G ~ Normal(0, 2000)

    # This is a shortcut for creating a variable defined on a circular domain
    # It creates:
    # τx ~ Normal(0, 1)
    # τy ~ Normal(0, 1)
    #  τ = atan(b.τy, b.τx)/2pi and a LogNormal(log(1.0), 0.02) prior on sqrt(b.τx^2 + b.τy^2)
    τ ~ UniformCircular(1.0)

    # τ represents the planet's position along its orbit in units of fractional orbital period
    # at the epoch `tref`.
    # I believe it is the most problematic parameter.

    # Fixed reference epoch.
    tref = tref1
end

# Alternative, traditional parameterization using Campbell angles (semi-major axis, inclination, etc.)
# These lead to even more divergences.
# @planet b VisualOrbit astrom begin
#     a ~ LogUniform(0.001, 10000.0)
#     e ~ Uniform(0.0, 0.99)
#     i ~ Sine()
#     ω ~ UniformCircular()
#     Ω ~ UniformCircular()
#     τ ~ UniformCircular(1.0)
# end


# Simple model for the star (with one orbitting planet, b).
# Note this is not a nested model, just a nice way of structuring multi-planet 
# models in the code.
@system HD1234 b begin
    # Parallax distance to the star (milliarcseconds of apparent motion per Earth year)
    # Sets overall scale of the model.
    plx ~ truncated(Normal(100, 5.0), lower=0)
    # Mass of the star in solar masses: determines orbital speed.
    M   ~ truncated(Normal(1.5, 0.2), lower=0)
end


## Generate custom LogDensityModel
# Bijectors.jl are used to remap variables from their restricted ranges to
# the full real line and corrects the prior density.
model = Octofitter.LogDensityModel(HD1234; autodiff=:ForwardDiff, verbosity=4)
# :FiniteDiff also supported if loaded.

## Test model with various adverserial inputs
Octofitter.checkmodel(model)

## Sample using AdvancedHMC
using Random
rng = Random.Xoshiro(1)


# Code is in sampling.jl beginning at line 629.
chain = Octofitter.advancedhmc(
    rng,
    model, 
    0.65; # Target acceptance. Cranking this up to 0.95+ will reduce the number of divergences to a narrow region.
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
num_err = map(samp->samp.stat.numerical_error, chain.info.states[1])

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
save(joinpath(@__DIR__, "pairplot-divergences.png"), fig, px_per_unit=4)


## Plot tree depth vs all parameters
using CairoMakie
using PairPlots
c = Makie.cgrad(:turbo, rev=true)
fig = pairplot(
    PairPlots.Series(chain[tree_depths .== 1,colnames,:], label="TD=1", color=c[1/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 2,colnames,:], label="TD=2", color=c[2/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 3,colnames,:], label="TD=3", color=c[3/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 4,colnames,:], label="TD=4", color=c[4/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 5,colnames,:], label="TD=5", color=c[5/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 6,colnames,:], label="TD=6", color=c[6/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 7,colnames,:], label="TD=7", color=c[7/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 8,colnames,:], label="TD=8", color=c[8/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 9,colnames,:], label="TD=9", color=c[9/13])    => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 10,colnames,:], label="TD=10", color=c[10/13]) => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 11,colnames,:], label="TD=11", color=c[11/13]) => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 12,colnames,:], label="TD=12", color=c[12/13]) => (PairPlots.Scatter(markersize=ms),),
    PairPlots.Series(chain[tree_depths .== 13,colnames,:], label="TD=13", color=c[13/13]) => (PairPlots.Scatter(markersize=ms),),
)

save(joinpath(@__DIR__, "pairplot-tree-depth.png"), fig, px_per_unit=4)
