
###
using Octofitter, Distributions, PlanetOrbits, Plots

## Generate synthetic transit event
u = [0.1, 0.4] # quad limb dark
ld = QuadLimbDark(u)

o = orbit(M=1,a=0.1,i=π/2,Ω=π/2,e=0,ω=0,τ=0)
ts =  0:0.01:60
tds = DirectDetections.transit_depth.(o, ts, 0.1, 696340e3, ld)
i = argmin(tds)

synth_obs = Table([
    # (epoch=t, phot= 4 < t < 6 ? DirectDetections.transit_depth(o, t, 0.1, 696340e3, ld)+0.0001randn() : 1 + 0.0001randn(), σ_phot=0.0001)
    (epoch=t, phot=DirectDetections.transit_depth(o, t, 0.1, 696340e3, ld)+0.0001randn(), σ_phot=0.0001)
    # for t in ts[i-40:i+40]
    # for t in ts[i-240:i+100]
    for t in ts[i-140:i+100]
])

p=plot()
scatter!(p, synth_obs.epoch, synth_obs.phot, yerror=synth_obs.σ_phot, xlims=(extrema(synth_obs.epoch) .+ (0,0)), color=:black, label="simulated data", xlabel="days", ylabel="transit depth")
p


## Construct system model

# We can add one or more named planets parameterized by a given orbit type
@named b = Planet{KepOrbit}(
    # Our variables block. Can be fixed, or distributions as priors.
    Variables(
        e = 0,           #
        τ = UniformCircular(1.0),   # 
        a = TruncatedNormal(0.1, 0.05, 0, Inf),
        i = π/2,
        ω = 0.0,
        Ω = π/2,
        r = Uniform(0.00, 0.3)
    ),
)

lightcurvedata = LightCurve4(
    # Limb darkening prescription from Transits.jl
    Transits.QuadLimbDark,
    # Observation table
    synth_obs
)

@named system = System(
    Variables(
        M = 0.78,
        R = 696340e3,

        # Limb darkening parameterization (u1 through u4)
        u1 = 0.1,#TruncatedNormal(0.1, 0.01, 0, 1),
        u2 = 0.4,#TruncatedNormal(0.4, 0.07, 0, 1)

    ),
    # Add observations (comment out to sample from priors only)
    lightcurvedata,
    # And planets
    b
)

## Sample from chains

results = DirectDetections.hmc(
    system, 0.65;


    # It appears we must take more initial gueses to find a good starting point
    # with transit data. 
    # If our initial conditions put a completely flat part of the transit model
    # overtop the observations, warm up will fail.
    initial_samples=250_000,

    adaptation =  2000,
    iterations =  2000,
    tree_depth = 9,

    verbosity = 4,
)

##
function transitplot(chain, planetkey)
    ii = rand(1:size(chain,1)*size(chain,3), 150)
    els = DirectDetections.construct_elements(chain, planetkey, ii)

    lightcurvedata = first(filter(obs->isa(obs, LightCurve4), chain.info.model.observations))

    tspan = range(minimum(lightcurvedata.table.epoch)-5, maximum(lightcurvedata.table.epoch)+5, step=0.02)
    # tspan = 0:0.005:60
    return tspan,  map(els, ii) do orb, i
        params = tuple()
        if haskey(chain, :u1)
            params = tuple(chain["u1"][i])
            if haskey(chain, :u2)
                params = tuple(chain["u1"][i], chain["u2"][i])
                if haskey(chain, :u3)
                    params = tuple(chain["u1"][i], chain["u2"][i], chain["u3"][i])
                    if haskey(chain, :u4)
                        params = tuple(chain["u1"][i], chain["u2"][i], chain["u3"][i], chain["u4"][i])
                    end
                end
            end
        end
        ld = DirectDetections.limbdarkfunc(lightcurvedata)(SVector(params))
        Rₛₜₐᵣ = chain["R"][i]
        r = chain["$planetkey.r"][i]

        phot_star =  DirectDetections.transit_depth.(orb, tspan, r, Rₛₜₐᵣ, ld)

    end


end
tss, tps = transitplot(results, :b);
p=plot()
plot!.(Ref(tss), tps,label="",color=:grey)


scatter!(p, lightcurvedata.table.epoch, lightcurvedata.table.phot, yerror=lightcurvedata.table.σ_phot, color=:black, label="simulated data", xlabel="days", ylabel="transit depth")
println("Done.")
p
