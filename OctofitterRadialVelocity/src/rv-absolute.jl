
# Radial Velocity data type
const rv_cols = (:epoch, :rv, :σ_rv)

"""
    StarAbsoluteRVLikelihood(
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
struct StarAbsoluteRVLikelihood{TTable<:Table,GP,TF,offset_symbol,jitter_symbol} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_name::String
    gaussian_process::GP
    trend_function::TF
    offset_symbol::Symbol
    jitter_symbol::Symbol
    function StarAbsoluteRVLikelihood(
        observations...;
        offset,
        jitter,
        trend_function=(θ_system, epoch)->zero(Octofitter._system_number_type(θ_system)),
        instrument_name="",
        gaussian_process=nothing
    )
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
        if hasproperty(table, :inst_idx) && length(unique(table.inst_idx)) > 1
            error("Deprecated: data from separate RV instruments should now be placed into different StarAbsoluteRVLikelihood likelihood objects, rather than specified by an inst_idx parameter.")
        end
        rows = map(eachrow(table)) do row′
            row = (;row′[1]..., rv=float(row′[1].rv[1]))
            return row
        end
        # We sort the data first by instrument index then by epoch to make some later code faster
        m = maximum(r->r.epoch, rows)
        ii = sortperm([r.epoch for r in rows])
        table = Table(rows[ii])
        if isnothing(instrument_name)
            instrument_name = string.(1:maximum(table.inst_idx))
        end
        # Pre-allocate symbols for offset and jitter term access
        offset_symbols = Tuple(Symbol("rv0_", i) for i in eachindex(instrument_name))
        jitter_symbols = Tuple(Symbol("jitter_", i) for i in eachindex(instrument_name))

        return new{typeof(table),typeof(gaussian_process),typeof(trend_function),offset, jitter, }(table, instrument_name, gaussian_process, trend_function, offset, jitter, )
    end
end
StarAbsoluteRVLikelihood(observations::NamedTuple...;kwargs...) = StarAbsoluteRVLikelihood(observations; kwargs...)
function Octofitter.likeobj_from_epoch_subset(obs::StarAbsoluteRVLikelihood, obs_inds)
    return StarAbsoluteRVLikelihood(
        obs.table[obs_inds,:,1]...;
        offset_symbol=obs.offset_symbol,
        jitter_symbol=obs.jitter_symbol,
        trend_function=obs.trend_function,
        instrument_name=obs.instrument_name,
        gaussian_process=obs.gaussian_process,
    )
end
export StarAbsoluteRVLikelihood


function _getparams(::StarAbsoluteRVLikelihood{TTable,GP,TF,offset_symbol,jitter_symbol}, θ_system) where {TTable,GP,TF,offset_symbol,jitter_symbol}
    offset = getproperty(θ_system, offset_symbol)
    jitter = getproperty(θ_system, jitter_symbol)
    return (;offset, jitter)
end

"""
Absolute radial velocity likelihood (for a star).
"""
function Octofitter.ln_like(
    rvlike::StarAbsoluteRVLikelihood,
    θ_system,
    planet_orbits::Tuple,
    orbit_solutions,
    orbit_solutions_i_epoch_start
)
    L = length(rvlike.table.epoch)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    (;offset, jitter) = _getparams(rvlike, θ_system)
    

    @no_escape begin

        # Each RV instrument index can have it's own barycentric RV offset and jitter.
        # Grab the offsets/jitters and store in a tuple. 
        # Then we can index by inst_idxs to look up the right offset and jitter.

        # Vector of radial velocity of the star at each epoch. Go through and sum up the influence of
        # each planet and put it into here. 
        rv_buf =  @alloc(T, L) # @SArray
        rv_var_buf =   @alloc(T, L) # @SArray

        # Data for this instrument:
        # epochs = rvlike.table.epoch
        # σ_rvs = rvlike.table.σ_rv
        # rvs = rvlike.table.rv

        # RV "data" calculation: measured RV + our barycentric rv calculation
        rv_buf .= rvlike.table.rv .- offset .- rvlike.trend_function(θ_system, rvlike.table.epoch)

        # Go through all planets and subtract their modelled influence on the RV signal:
        # You could consider `rv_star` as the residuals after subtracting these.
        
        for planet_i in eachindex(planet_orbits)
            orbit = planet_orbits[planet_i]
            planet_mass = θ_system.planets[planet_i].mass
            for epoch_i in eachindex(rvlike.table.epoch)
                rv_buf[epoch_i] -= radvel(
                    orbit_solutions[planet_i][epoch_i+orbit_solutions_i_epoch_start],
                    planet_mass*Octofitter.mjup2msol
                )
            end
            # Zygote version:
            # rv_buf -= radvel.(orbit, epochs, planet_mass*Octofitter.mjup2msol)
        end        

        # The noise variance per observation is the measurement noise and the jitter added
        # in quadrature
        rv_var_buf .= rvlike.table.σ_rv.^2 .+ jitter^2

        # Two code paths, depending on if we are modelling the residuals by 
        # a Gaussian process or not.
        if isnothing(rvlike.gaussian_process)
            # Don't fit a GP
            fx = MvNormal(Diagonal((rv_var_buf)))
            ll += logpdf(fx, rv_buf)
        else
            # Fit a GP
            local gp
            try
                gp = @inline rvlike.gaussian_process(θ_system)
            catch err
                if err isa PosDefException
                    @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    ll = convert(T,-Inf)
                elseif err isa ArgumentError
                    @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    ll = convert(T,-Inf)
                else
                    rethrow(err)
                end
            end             
            
            if isfinite(ll)
                fx = gp(rvlike.table.epoch, rv_var_buf)
                try
                    ll += logpdf(fx, rv_buf)
                catch err
                    if err isa PosDefException || err isa DomainError
                        @warn "err" exception=(err, catch_backtrace()) θ_system
                        ll = convert(T,-Inf)
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end
    return ll
end




# Generate new radial velocity observations for a star
function Octofitter.generate_from_params(like::StarAbsoluteRVLikelihood, θ_system, orbits::Vector{<:Visual{KepOrbit}})

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = radvel.(reshape(orbits, :, 1), epochs, transpose(planet_masses))
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return StarAbsoluteRVLikelihood(radvel_table)
end



# # Plot recipe for astrometry data
# @recipe function f(rv::StarAbsoluteRVLikelihood)
   
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