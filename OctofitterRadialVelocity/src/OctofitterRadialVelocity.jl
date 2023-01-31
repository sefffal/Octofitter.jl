module OctofitterRadialVelocity

using Octofitter
using PlanetOrbits
using Tables, TypedTables
using Distributions
using DataDeps

# Radial Velocity data type
const rv_cols = (:epoch, :rv, :σ_rv)
struct RadialVelocity{TTable<:Table} <: Octofitter.AbstractObs
    table::TTable
    function RadialVelocity(observations...)
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Ecpected columns $rv_cols")
        end
        return new{typeof(table)}(table)
    end
end
RadialVelocity(observations::NamedTuple...) = RadialVelocity(observations)
export RadialVelocity

"""
Radial velocity likelihood.
"""
function Octofitter.ln_like(rv::RadialVelocity, θ_system, elements)
    T = Float64
    ll = zero(T)

    for i in eachindex(rv.table.epoch)
        rv_star = zero(T)
        for j in eachindex(elements)
            θ_planet = θ_system.planets[j]
            orbit = elements[j]

            if θ_planet.mass < 0
                return -Inf 
            end

            rv_star += radvel(orbit, rv.table.epoch[i], θ_planet.mass*Octofitter.mjup2msol)
        end
        # Each measurement is tagged with a jitter and rv zero point variable.
        # We then query the system variables for them.
        # A previous implementation used symbols instead of indices but it was too slow.
        # barycentric_rv_inst = getproperty(θ_system, rv.table.rv0[i])::Float64
        # jitter_inst = getproperty(θ_system, rv.table.jitter[i])::Float64
        if !hasproperty(rv.table, :inst_idx)
            barycentric_rv_inst = θ_system.rv0
            jitter_inst = θ_system.jitter
        else
            inst_idx = rv.table.inst_idx[i]
            if inst_idx == 1
                barycentric_rv_inst = θ_system.rv0_1
                jitter_inst = θ_system.jitter_1
            elseif inst_idx == 2
                barycentric_rv_inst = θ_system.rv0_2
                jitter_inst = θ_system.jitter_2
            elseif inst_idx == 3
                barycentric_rv_inst = θ_system.rv0_3
                jitter_inst = θ_system.jitter_3
            elseif inst_idx == 4
                barycentric_rv_inst = θ_system.rv0_4
                jitter_inst = θ_system.jitter_4
            else
                error("More than four radial velocity instruments are not yet supported for performance reasons")
            end
        end
        resid = rv_star + barycentric_rv_inst - rv.table.rv[i]
        σ² = rv.table.σ_rv[i]^2 + jitter_inst^2
        χ² = -0.5resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end

    return ll
end




# Generate new radial velocity observations for a planet
function Octofitter.genobs(obs::RadialVelocity, elem::VisualOrbit, θ_planet)

    # Get epochs and uncertainties from observations
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocity(radvel_table)
end




# Generate new radial velocity observations for a star
function Octofitter.genobs(obs::RadialVelocity, elems::Vector{<:VisualOrbit}, θ_system)

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = radvel.(reshape(elems, :, 1), epochs, transpose(planet_masses))
    noise = randn(length(epochs)) .* θ_system.jitter
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return RadialVelocity(radvel_table)
end





function HARPS_search(target, catalog=datadep"HARPS_RVBank")

    
    rvbank = CSV.read(catalog, Table)

    # GJ436

    # open(joinpath(catalog, "list.dat"), read=true) do io
    #     map(readlines) 
    return filter(row->row.target==target, eachrow(rvbank))
end



# Plot recipe for astrometry data
# @recipe function f(rv::RadialVelocity)

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
        # "42916f0008b8dd9ead18f8e793c009800b2dcf02b523d219463e1ecf021af6d1",
        "https://www2.mpia-hd.mpg.de/homes/trifonov/HARPS_RVBank_v1.csv",
    ))
    return
end

end