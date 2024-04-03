# Required packages:
# Octofitter, Distributions, Plots, CairoMakie, PairPlots, Random
# See Project and Manifest files.

using Octofitter, Distributions, PlanetOrbits
using OctofitterWhereistheplanet
cd(@__DIR__)

sbc_index = parse(Int, ARGS[1])
@info "Running SBC trial $sbc_index"
outname = @sprintf("sbctest-%04d", sbc_index)

# Reference epoch for tau (orbital position) in convention used by Orbitize.
const tref = 58849.0
const theta_epoch = mean(Octofitter.Whereistheplanet_astrom("hr8799";object=4)[2].table.epoch)


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

    AFr ~ LogUniform(10, 10_000)
    AFt ~ UniformCircular()
    BGr ~ LogUniform(10, 10_000)
    BGt ~ UniformCircular()
    A = b.AFr * cos(b.AFt)
    F = b.AFr * sin(b.AFt)
    B = b.BGr * cos(b.BGt)
    G = b.BGr * sin(b.BGt)
    
    tref = tref

    θ ~ UniformCircular()
    τ = τ_from_θ(system, b)

    mass = 7

end Octofitter.Whereistheplanet_astrom("hr8799";object=4)[2]



# Simple model for the star (with one orbitting planet, b).
# Note this is not a nested model, just a nice way of structuring multi-planet 
# models in the code.
@system system begin
    # Parallax distance to the star (milliarcseconds of apparent motion per Earth year)
    # Sets overall scale of the model.
    plx ~ gaia_plx(;gaia_id=2832463659640297472)
    # Mass of the star in solar masses: determines orbital speed.
    M   ~ truncated(Normal(1.5, 0.2),lower=0)

    # True proper motion of the system barycentre.
    # The variance on this will be quite low, but let the mean be determined
    # entirely from the data.
    # pmra  = 0 #~ Normal(0, 500)
    # pmdec = 0 #~ Normal(0, 500)
end b


## Generate custom LogDensityModel
model = Octofitter.LogDensityModel(system; autodiff=:ForwardDiff, verbosity=4)


## Sample using AdvancedHMC
using Random
rng = Random.Xoshiro(1)

##
settings = Dict(
    :target_accept=>0.85,
    :num_chains=>1,
    :adaptation=>1000,
    :iterations=>1000,
    :thinning=>50,
    :tree_depth=>14,
    :verbosity=>4
)
Octofitter.calibrate(HR8799, settings, outname)
