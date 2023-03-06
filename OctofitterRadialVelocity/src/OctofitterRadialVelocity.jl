module OctofitterRadialVelocity

using Octofitter
using PlanetOrbits
using Tables, TypedTables
using Distributions
using DataDeps
using LoopVectorization
using StrideArrays

# Radial Velocity data type
const rv_cols = (:epoch, :rv, :σ_rv)
struct RadialVelocityLikelihood{TTable<:Table} <: Octofitter.AbstractLikelihood
    table::TTable
    function RadialVelocityLikelihood(observations...)
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Ecpected columns $rv_cols")
        end
        return new{typeof(table)}(table)
    end
end
RadialVelocityLikelihood(observations::NamedTuple...) = RadialVelocityLikelihood(observations)
export RadialVelocityLikelihood

"""
Radial velocity likelihood.
"""
function Octofitter.ln_like(rv::RadialVelocityLikelihood, θ_system, elements, num_epochs::Val{L}=Val(length(rv.table))) where L
    T = typeof(θ_system.M)
    ll = zero(T)

    epochs = rv.table.epoch
    σ_rvs = rv.table.σ_rv
    inst_idxs = rv.table.inst_idx
    rvs = rv.table.rv

    # # # TODO: This is a debug override. This is forcing omega to be shifted by pi
    # # # to compare with orvara.
    # nt = θ_system.planets.B
    # nt = merge(nt,(;ω=nt.ω+π))
    # elements = (Octofitter.construct_elements(VisualOrbit, θ_system, nt),)

    # single_instrument_mode = !hasproperty(rv.table, :inst_idx)
    barycentric_rv_inst = (
        hasproperty(rv.table, :rv0_1) ? rv.table.rv0_1 : zero(T),
        hasproperty(rv.table, :rv0_2) ? rv.table.rv0_2 : zero(T),
        hasproperty(rv.table, :rv0_3) ? rv.table.rv0_3 : zero(T),
        hasproperty(rv.table, :rv0_4) ? rv.table.rv0_4 : zero(T),
    )
    jitter_inst = (
        hasproperty(rv.table, :jitter_1) ? rv.table.jitter_1 : zero(T),
        hasproperty(rv.table, :jitter_2) ? rv.table.jitter_2 : zero(T),
        hasproperty(rv.table, :jitter_3) ? rv.table.jitter_3 : zero(T),
        hasproperty(rv.table, :jitter_4) ? rv.table.jitter_4 : zero(T),
    )
    # Vector of radial velocity of the star at each epoch. Go through and sum up the influence of
    # each planet and put it into here. 
    # Then loop through and get likelihood.
    # Hopefully this is more efficient than looping over each planet at each epoch and adding up the likelihood.
    rv_star = StrideArray{T}(undef, (StaticInt(num_epochs),))
    fill!(rv_star, 0)
    for planet_i in eachindex(elements)
        orbit = elements[planet_i]
        # Need to structarrays orbit???
        planet_mass = θ_system.planets[planet_i].mass
        @turbo for epoch_i in eachindex(epochs)
            rv_star[epoch_i] += radvel(orbit, epochs[epoch_i], planet_mass*Octofitter.mjup2msol)
        end
    end
    @turbo for i in eachindex(epochs)
        # Each measurement is tagged with a jitter and rv zero point variable.
        # We then query the system variables for them.
        # A previous implementation used symbols instead of indices but it was too slow.
        inst_idx = inst_idxs[i]
        resid = rv_star[i] + barycentric_rv_inst[inst_idx] - rvs[i]
        σ² = σ_rvs[i]^2 + jitter_inst[inst_idx]^2
        χ² = -0.5resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²

        # Leveraging Distributions.jl to make this clearer:
        # ll += logpdf(Normal(0, sqrt(σ²)), resid)
    end

    return ll
end




# Generate new radial velocity observations for a planet
function Octofitter.genobs(like::RadialVelocityLikelihood, elem::VisualOrbit, θ_planet)

    # Get epochs and uncertainties from observations
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocityLikelihood(radvel_table)
end




# Generate new radial velocity observations for a star
function Octofitter.genobs(like::RadialVelocityLikelihood, elems::Vector{<:VisualOrbit}, θ_system)

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = radvel.(reshape(elems, :, 1), epochs, transpose(planet_masses))
    noise = randn(length(epochs)) .* θ_system.jitter
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocityLikelihood(radvel_table)
end

mjd2jd(mjd) = mjd - 2400000.5
jd2mjd(jd) = jd + 2400000.5



include("harps.jl")
include("hires.jl")


# Plot recipe for astrometry data
# @recipe function f(rv::RadialVelocityLikelihood)

# end




function __init__()

    register(DataDep("HARPS_RVBank",
        """
        Dataset:     A public HARPS radial velocity database corrected for systematic errors
        Author:      Trifonov et al.
        License:     CC0-1.0
        Publication: https://www.aanda.org/articles/aa/full_html/2020/04/aa36686-19/aa36686-19.html
        Website:     https://www2.mpia-hd.mpg.de/homes/trifonov/HARPS_RVBank.html

        A public HARPS radial velocity database corrected for systematic errors. (2020)
        
        File size: 132MiB
        """,
        "17b2a7f47569de11ff1747a96997203431c81586ffcf08212ddaa250bb879a40",
        "https://www2.mpia-hd.mpg.de/homes/trifonov/HARPS_RVBank_v1.csv",
    ))

    register(DataDep("HIRES_rvs",
        """
        Dataset:     A public HIRES radial velocity database corrected for systematic errors
        Author:      Butler et al.
        License:     
        Publication: https://ui.adsabs.harvard.edu/abs/2017yCat..51530208B/abstract
        Website:     https://ebps.carnegiescience.edu/data/hireskeck-data

        File size: 3.7MiB
        """,
        "https://drive.google.com/uc?id=10xCy8UIH8wUAnNJ8zCN8kfFzazWnw-f_&export=download",
        "ad68c2edb69150318e8d47e34189fe104f2a5194a4fcd363c78c741755893251",
        post_fetch_method=unpack
    ))




    return
end

end
