"""
    PlanetRelativeRVLikelihood(
        (;inst_idx=1, epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;inst_idx=1, epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;inst_idx=1, epoch=5100.2,  rv=7.90,  σ_rv=.11),
    )

Represents a likelihood function of relative astometry between a host star and a secondary body.
`:epoch` (mjd), `:rv` (m/s), and `:σ_rv` (m/s), and `:inst_idx` are all required.

`:inst_idx` is used to distinguish RV time series between insturments so that they may optionally
be fit with different zero points and jitters.
In addition to the example above, any Tables.jl compatible source can be provided.
"""
struct PlanetRelativeRVLikelihood{TTable<:Table,GP} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_names::Vector{String}
    gaussian_process::GP
    function PlanetRelativeRVLikelihood(observations...; instrument_names=nothing, gaussian_process=nothing)
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
       
        # We sort the data first by instrument index then by epoch to make some later code faster
        m = maximum(table.epoch)
        ii = sortperm(table.epoch)
        table = table[ii,:]
        if isnothing(instrument_names)
            instrument_names = ["1"]
        end

        rows = map(eachrow(table)) do row′
            row = (;row′[1]..., inst_idx=1)
            return row
        end
        table = Table(rows)

        return new{typeof(table),typeof(gaussian_process)}(table, instrument_names, gaussian_process)
    end
end
PlanetRelativeRVLikelihood(observations::NamedTuple...;kwargs...) = PlanetRelativeRVLikelihood(observations; kwargs...)
export PlanetRelativeRVLikelihood

"""
Radial velocity likelihood.
"""
function Octofitter.ln_like(
    rvlike::PlanetRelativeRVLikelihood,
    θ_planet,
    planet_orbit::AbstractOrbit,
    num_epochs::Val{L}=Val(length(rvlike.table))
) where L
    T = typeof(first(θ_planet))
    ll = zero(T)

    # We don't support multiple instruments for relative RVs.
    # This isn't essential because there should be no zero point offsets that
    # differ by instrument.
    # The only thing that could differ would be the jitter

    # Data for this instrument:
    epochs = vec(rvlike.table.epoch)
    σ_rvs = vec(rvlike.table.σ_rv)
    rvs = vec(rvlike.table.rv)
    jitter = θ_planet.jitter
    noise_var = zeros(T, length(rvs))

    # RV "data" calculation: measured RV + our barycentric rv calculation
    rv_star_buf = zeros(T, length(rvs))
    rv_star_buf .+= rvs
    
    ## Zygote version:
    # rv_star_buf = @view(rv_star_buf_all_inst[istart:iend])

    # Go through all planets and subtract their modelled influence on the RV signal:
    # You could consider `rv_star` as the residuals after subtracting these.
    orbit = planet_orbit
    # Threads.@threads 
    for epoch_i in eachindex(epochs)
        rv_star_buf[epoch_i] -= radvel(orbit, epochs[epoch_i])
    end

    # The noise variance per observation is the measurement noise and the jitter added
    # in quadrature
    noise_var .= σ_rvs.^2 .+ jitter^2

    # Two code paths, depending on if we are modelling the residuals by 
    # a Gaussian process or not.
    if isnothing(rvlike.gaussian_process)
        # Don't fit a GP
        fx = MvNormal(Diagonal((noise_var)))
    else
        # Fit a GP
        local gp
        try
            gp = @inline rvlike.gaussian_process(θ_system)
        catch err
            if err isa PosDefException
                # @warn "err" exception=(err, catch_backtrace()) maxlog=1
                return -Inf
            elseif err isa ArgumentError
                # @warn "err" exception=(err, catch_backtrace()) maxlog=1
                return -Inf
            else
                rethrow(err)
            end
        end                

        fx = gp(epochs, noise_var)
        try
            ll += logpdf(fx, rv_star_buf)
        catch err
            if err isa PosDefException || err isa DomainError
                return -Inf
            else
                rethrow(err)
            end
        end
    end

    ll += logpdf(fx, rv_star_buf)

    return ll
end


# Generate new radial velocity observations for a planet
function Octofitter.generate_from_params(like::PlanetRelativeRVLikelihood, θ_planet, elem::PlanetOrbits.AbstractOrbit)

    # Get epochs and uncertainties from observations
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return PlanetRelativeRVLikelihood(radvel_table)
end


# Generate new radial velocity observations for a star
function Octofitter.generate_from_params(like::PlanetRelativeRVLikelihood, θ_system, orbits::Vector{<:Visual{KepOrbit}})

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 

    # Generate new star radial velocity data
    rvs = radvel.(reshape(orbits, :, 1), epochs)
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return PlanetRelativeRVLikelihood(radvel_table)
end



# # Plot recipe for astrometry data
# @recipe function f(rv::PlanetRelativeRVLikelihood)
   
#     xguide --> "time (mjd)"
#     yguide --> "radvel (m/s)"

#     multiple_instruments = hasproperty(rv.table,:inst_idx) && 
#                            length(unique(rv.table.inst_idx)) > 1
#     if !multiple_instruments
#         @series begin
#             color --> :black
#             label := nothing
#             seriestype := :scatter
#             markersize--> 0
#             yerr := rv.table.σ_rv
#             rv.table.epoch, rv.table.rv
#         end
#     else
#         for inst_idx in sort(unique(rv.table.inst_idx))
#             @series begin
#                 label := nothing
#                 seriestype := :scatter
#                 markersize--> 0
#                 color-->inst_idx
#                 markerstrokecolor-->inst_idx
#                 yerr := rv.table.σ_rv[rv.table.inst_idx.==inst_idx]
#                 rv.table.epoch[rv.table.inst_idx.==inst_idx], rv.table.rv[rv.table.inst_idx.==inst_idx]
#             end
#         end
#     end


# end