module OctofitterVisibilities

using Octofitter
using PlanetOrbits
using Tables, TypedTables

const vis_cols = (:epoch, :col1, :col2)
struct VisibiltiesLikelihood{TTable<:Table} <: Octofitter.AbstractLikelihood
    table::TTable
    function VisibiltiesLikelihood(observations...)
        table = Table(observations...)
        if !issubset(vis_cols, Tables.columnnames(table))
            error("Ecpected columns $vis_cols")
        end
        return new{typeof(table)}(table)
    end
end
VisibiltiesLikelihood(observations::NamedTuple...) = VisibiltiesLikelihood(observations)
export VisibiltiesLikelihood

"""
Visibliitiy modelling likelihood for point sources.
"""
function Octofitter.ln_like(vis::VisibiltiesLikelihood, θ_system, orbits, num_epochs::Val{L}=Val(length(rv.table))) where L
    T = typeof(θ_system.M)
    ll = zero(T)

    # Access the data here: 
    epochs = vis.table.epoch
    bands = vis.table.band
    # and observations too

    # `θ_system` contains all parameter values, e.g. photometry of each planet
    
    # `orbits` creates pre-constructed orbit objects (one per planet) though these could be created from `θ_system`
    # We can use `orbits` to calculate the position of a planet on a given epoch.

    # Loop through planets
    for i in eachindex(orbits)
        # All parameters relevant to this planet
        θ_planet = θ_system.planets[i]
        # orbit object pre-created from above parameters (shared between all likelihood functions)
        orbit = orbits[i]

        # Get model flux parameter in this band (band provided as a symbol, e.g. :L along with data in table row.)
        f_band = getproperty(θ_planet, band)

        
        for j in eachindex(epochs)
            epoch = epochs[j]
            band = bands[j]
            
            sol = orbitsolve(orbit, epoch)
            Δra  = raoff(sol)  # in mas
            Δdec = decoff(sol) # in mas

            # Accumulate into likelihood
            # ll += ...
        end
    end


    return ll
end




# # Generate new observations for a planet (I don't think this is relevant)
# function Octofitter.generate_from_params(like::VisibiltiesLikelihood, orbit::PlanetOrbits.AbstractOrbit, θ_planet)

# end




# Generate new observations for a system of possibly multiple planets
function Octofitter.generate_from_params(like::VisibiltiesLikelihood, orbits::Vector{<:VisualOrbit}, θ_system)

    # # Get epochs, uncertainties, and planet masses from observations and parameters
    # epochs = like.table.epoch 
    # σ_rvs = like.table.σ_rv 
    # planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # # Generate new star radial velocity data
    # rvs = radvel.(reshape(orbits, :, 1), epochs, transpose(planet_masses))
    # noise = randn(length(epochs)) .* θ_system.jitter
    # rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    # radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    # return with same number of rows: band, epoch
    # position(s) of point sources according to orbits, θ_system

    return VisibiltiesLikelihood(new_vis_like_table)
end




end
