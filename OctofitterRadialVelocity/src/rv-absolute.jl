
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
struct StarAbsoluteRVLikelihood{TTable<:Table,GP,TOffsets,TJitters} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_names::Vector{String}
    gaussian_process::GP
    offset_symbols::TOffsets
    jitter_symbols::TJitters
    function StarAbsoluteRVLikelihood(observations...; instrument_names=nothing, gaussian_process=nothing)
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
        rows = map(eachrow(table)) do row′
            row = (;inst_idx=1, row′[1]..., rv=float(row′[1].rv[1]))
            return row
        end
        # We sort the data first by instrument index then by epoch to make some later code faster
        m = maximum(r->r.epoch, rows)
        ii = sortperm([
            r.inst_idx .* (10m) .+ r.epoch
            for r in rows
        ])
        table = Table(rows[ii])
        # Check instrument indexes are contiguous
        if length(unique(table.inst_idx)) != maximum(table.inst_idx)
            error("instrument indexes must run from 1:N without gaps")
        end
        if isnothing(instrument_names)
            instrument_names = string.(1:maximum(table.inst_idx))
        end
        # Pre-allocate symbols for offset and jitter term access
        offset_symbols = Tuple(Symbol("rv0_", i) for i in eachindex(instrument_names))
        jitter_symbols = Tuple(Symbol("jitter_", i) for i in eachindex(instrument_names))

        return new{typeof(table),typeof(gaussian_process), typeof(offset_symbols), typeof(jitter_symbols)}(table, instrument_names, gaussian_process, offset_symbols, jitter_symbols)
    end
end
StarAbsoluteRVLikelihood(observations::NamedTuple...;kwargs...) = StarAbsoluteRVLikelihood(observations; kwargs...)
export StarAbsoluteRVLikelihood

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

    # Each RV instrument index can have it's own barycentric RV offset and jitter.
    # Grab the offsets/jitters and store in a tuple. 
    # Then we can index by inst_idxs to look up the right offset and jitter.
    barycentric_rv_inst = θ_system.rv0
    # barycentric_rv_inst = SVector{length(rvlike.offset_symbols), T}(getproperty(θ_system, symbol) for symbol in rvlike.offset_symbols)
    jitter_inst =  θ_system.jitter #SVector{length(rvlike.jitter_symbols), T}(getproperty(θ_system, symbol) for symbol in rvlike.jitter_symbols)
    max_inst_idx = maximum(rvlike.table.inst_idx)

    # Vector of radial velocity of the star at each epoch. Go through and sum up the influence of
    # each planet and put it into here. 
    rv_star_buf_all_inst =  zeros(T, L) # @SArray
    noise_var_all_inst =   zeros(T, L) # @SArray

    # Work through RV data and model one instrument at a time
    for inst_idx in 1:max_inst_idx
        if isnan(jitter_inst[inst_idx]) || isnan(barycentric_rv_inst[inst_idx])
            # @warn("`jitter_$inst_idx` and `rv0_$inst_idx` must be provided (were NaN)", maxlog=5)
        end

        # Find the range in the table occupied by data from this instrument

        istart = only(findfirst(==(inst_idx), vec(rvlike.table.inst_idx)))
        iend = only(findlast(==(inst_idx), vec(rvlike.table.inst_idx)))
        # istart = searchsortedfirst(vec(rvlike.table.inst_idx), inst_idx)
        # iend = searchsortedlast(vec(rvlike.table.inst_idx), inst_idx)
        # # Make Enzyme happy:
        # istart = 1
        # for i in eachindex(rvlike.table.inst_idx)
        #     if rvlike.table.inst_idx[i] == inst_idx
        #         istart = i
        #         break
        #     end
        # end
        # iend = L
        # for i in reverse(eachindex(rvlike.table.inst_idx))
        #     if rvlike.table.inst_idx[i] == inst_idx
        #         iend = i
        #         break
        #     end
        # end

        # Data for this instrument:
        epochs = @view rvlike.table.epoch[istart:iend]
        σ_rvs = @view rvlike.table.σ_rv[istart:iend]
        rvs = @view rvlike.table.rv[istart:iend]
        noise_var = @view noise_var_all_inst[istart:iend]

        # RV "data" calculation: measured RV + our barycentric rv calculation
        rv_star_buf = @view rv_star_buf_all_inst[istart:iend] 
        rv_star_buf .+= rvs
        rv_star_buf .-= barycentric_rv_inst[inst_idx] 

        # # Zygote version:
        # rv_star_buf = @view(rv_star_buf_all_inst[istart:iend]) .+ barycentric_rv_inst[inst_idx]

        # rv_star_buf .+= barycentric_rv_inst[inst_idx] 
        # rv_star_buf += barycentric_rv_inst[inst_idx] 

        # Go through all planets and subtract their modelled influence on the RV signal:
        # You could consider `rv_star` as the residuals after subtracting these.
        
        for planet_i in eachindex(planet_orbits)
            orbit = planet_orbits[planet_i]
            planet_mass = θ_system.planets[planet_i].mass
            for epoch_i_inst in eachindex(epochs)
                epoch_i = epoch_i_inst + istart
                rv_star_buf[epoch_i_inst] -= radvel(
                    orbit_solutions[planet_i][epoch_i+orbit_solutions_i_epoch_start],
                    planet_mass*Octofitter.mjup2msol
                )
            end
            # Zygote version:
            # rv_star_buf -= radvel.(orbit, epochs, planet_mass*Octofitter.mjup2msol)
        end        

        # The noise variance per observation is the measurement noise and the jitter added
        # in quadrature
        noise_var .= σ_rvs.^2 .+ jitter_inst[inst_idx]^2

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
                    @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    return convert(T,-Inf)
                elseif err isa ArgumentError
                    @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    return convert(T,-Inf)
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

    end
    return ll
end



# Generate new radial velocity observations for a planet
function Octofitter.generate_from_params(like::StarAbsoluteRVLikelihood, θ_planet, elem::PlanetOrbits.AbstractOrbit)

    # Get epochs and uncertainties from observations
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 

    # Generate new planet radial velocity data
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return StarAbsoluteRVLikelihood(radvel_table)
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