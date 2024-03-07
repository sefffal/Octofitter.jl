

using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using CairoMakie
using PairPlots
using CSV
using DataFrames
using Distributions

rv_file = download("https://raw.githubusercontent.com/California-Planet-Search/radvel/master/example_data/k2-131.txt")
rv_dat_raw = CSV.read(rv_file, DataFrame, delim=' ')
rv_dat = DataFrame();
rv_dat.epoch = jd2mjd.(rv_dat_raw.time)
rv_dat.rv = rv_dat_raw.mnvel
rv_dat.σ_rv = rv_dat_raw.errvel
tels = sort(unique(rv_dat_raw.tel))
rv_dat.inst_idx = map(rv_dat_raw.tel) do tel
    findfirst(==(tel), tels)
end
rvlike = StarAbsoluteRVLikelihood(
    rv_dat,
    instrument_names=["harps-n", "psf"],
)



@planet b RadialVelocityOrbit begin
    e = 0
    ω = 0.0

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ truncated(Normal(0.3693038/365.256360417, 0.0000091/365.256360417),lower=0.00001)
    a = cbrt(system.M * b.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp = b.τ*b.P*365.256360417 + 57782 # reference epoch for τ. Choose an MJD date near your data.
    
    # minimum planet mass [jupiter masses]. really m*sin(i)
    mass ~ LogUniform(0.001, 10)
end


@system k2_132 begin
    # total mass [solar masses]
    M ~ truncated(Normal(0.82, 0.02),lower=0) # (Baines & Armstrong 2011).

    # HARPS-N
    rv0_1 ~ Normal(0,10000) # m/s
    jitter_1 ~ LogUniform(0.01,100) # m/s

    # FPS
    rv0_2 ~ Normal(0,10000) # m/s
    jitter_2 ~ LogUniform(0.01,100) # m/s
end rvlike b


model = Octofitter.LogDensityModel(k2_132)


# in a script.jl:
struct MyTargetFlag end 
using Pigeons
Pigeons.instantiate_target(flag::MyTargetFlag) = model
