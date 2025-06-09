"""
    PlanetRelativeRVLikelihood(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);
        
        instrument_name="inst name",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100.0)  # RV jitter (m/s)
        end
    )

    # Example with Gaussian Process:
    PlanetRelativeRVLikelihood(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);
        
        instrument_name="inst name",
        gaussian_process = θ_obs -> GP(θ_obs.gp_η₁^2 * SqExponentialKernel() ∘ ScaleTransform(1/θ_obs.gp_η₂)),
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100.0)     # RV jitter (m/s)
            gp_η₁ ~ LogUniform(1.0, 100.0)      # GP amplitude
            gp_η₂ ~ LogUniform(1.0, 100.0)      # GP length scale
        end
    )

Represents a likelihood function of relative radial velocity between a host star and a secondary body.
`:epoch` (mjd), `:rv` (m/s), and `:σ_rv` (m/s) are all required.

In addition to the example above, any Tables.jl compatible source can be provided.

The `jitter` variable should be defined in the variables block and represents additional 
uncertainty to be added in quadrature to the formal measurement errors.

When using a Gaussian process, the `gaussian_process` parameter should be a function that takes
`θ_obs` (observation parameters) and returns a GP kernel. GP hyperparameters should be defined
in the variables block and accessed via `θ_obs.parameter_name`.
"""
struct PlanetRelativeRVLikelihood{TTable<:Table,GP} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_name::String
    gaussian_process::GP
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    function PlanetRelativeRVLikelihood(observations...; instrument_name="1", gaussian_process=nothing, variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin;end))
        (priors,derived)=variables
        table = Table(observations...)
        if !Octofitter.equal_length_cols(table)
            error("The columns in the input data do not all have the same length")
        end
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
       
        # We sort the data first by instrument index then by epoch to make some later code faster
        ii = sortperm(table.epoch)
        table = table[ii,:]

        rows = map(eachrow(table)) do row′
            row = (;row′[1]..., rv=float(row′[1].rv[1]))
            return row
        end
        table = Table(rows)
        ii = sortperm(table.epoch)
        table = table[ii]

        return new{typeof(table),typeof(gaussian_process)}(table, instrument_name, gaussian_process, priors, derived)
    end
end
PlanetRelativeRVLikelihood(observations::NamedTuple...;kwargs...) = PlanetRelativeRVLikelihood(observations; kwargs...)
function Octofitter.likeobj_from_epoch_subset(obs::PlanetRelativeRVLikelihood, obs_inds)
    return PlanetRelativeRVLikelihood(
        obs.table[obs_inds,:,1]...;
        instrument_name=instrument_name(obs),
        gaussian_process=obs.gaussian_process,
        variables=(obs.priors, obs.derived)
    )
end
export PlanetRelativeRVLikelihood



"""
Radial velocity likelihood.
"""
function Octofitter.ln_like(rvlike::PlanetRelativeRVLikelihood, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)


    T = Octofitter._system_number_type(θ_planet)
    ll = zero(T)

    this_orbit = orbits[i_planet]

    # We don't support multiple instruments for relative RVs.
    # This isn't essential because there should be no zero point offsets that
    # differ by instrument.
    # The only thing that could differ would be the jitter

    # Data for this instrument:
    epochs = vec(rvlike.table.epoch)
    σ_rvs = vec(rvlike.table.σ_rv)
    rvs = vec(rvlike.table.rv)
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)

    @no_escape begin
        noise_var = @alloc(T, length(rvs))
        fill!(noise_var, 0)

        # RV "data" calculation: measured RV + our barycentric rv calculation
        rv_star_buf = @alloc(T, length(rvs))
        fill!(rv_star_buf, 0)
        rv_star_buf .+= rvs
        
        ## Zygote version:
        # rv_star_buf = @view(rv_star_buf_all_inst[istart:iend])

        # Go through all planets and subtract their modelled influence on the RV signal:
        # You could consider `rv_star` as the residuals after subtracting these.
        # Threads.@threads 
        for i_epoch in eachindex(epochs)
            sol = orbit_solutions[i_planet][i_epoch+orbit_solutions_i_epoch_start]
            @assert isapprox(rvlike.table.epoch[i_epoch], PlanetOrbits.soltime(sol), rtol=1e-2)
            # Relative RV due to planet
            rv_star_buf[i_epoch] -= radvel(sol)

            for (i_other_planet, key) in enumerate(keys(θ_system.planets))
                orbit_other = orbits[i_other_planet]
                # Account for perturbation due to any inner planets with non-zero mass
                if semimajoraxis(orbit_other) < semimajoraxis(this_orbit)
                    θ_planet′ = θ_system.planets[key]
                    if !hasproperty(θ_planet′, :mass)
                        continue
                    end
                    mass_other = θ_planet′.mass*Octofitter.mjup2msol
                    sol′ = orbit_solutions[i_other_planet][i_epoch + orbit_solutions_i_epoch_start]
                    
                    rv_star_buf[i_epoch] -= radvel(sol′, mass_other)
                    
                    @assert isapprox(rvlike.table.epoch[i_epoch], PlanetOrbits.soltime(sol′), rtol=1e-2)
                end
            end
        end

        # The noise variance per observation is the measurement noise and the jitter added
        # in quadrature
        noise_var .= σ_rvs.^2 .+ jitter^2

        # Two code paths, depending on if we are modelling the residuals by 
        # a Gaussian process or not.
        if isnothing(rvlike.gaussian_process)
            # Don't fit a GP
            fx = MvNormal(Diagonal((noise_var)))
            ll += logpdf(fx, rv_star_buf)
        else
            # Fit a GP
            local gp
            try
                gp = @inline rvlike.gaussian_process(θ_obs)
                fx = gp(epochs, noise_var)
                ll += logpdf(fx, rv_star_buf)
            catch err
                if err isa PosDefException
                    # @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    ll = convert(T, -Inf)
                elseif err isa ArgumentError
                    # @warn "err" exception=(err, catch_backtrace()) maxlog=1
                    ll = convert(T, -Inf)
                else
                    rethrow(err)
                end
            end                
        end
    end

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

    return PlanetRelativeRVLikelihood(
        radvel_table;
        instrument_name=instrument_name(like),
        gaussian_process=like.gaussian_process,
        variables=(like.priors, like.derived)
    )
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

    return PlanetRelativeRVLikelihood(
        radvel_table;
        instrument_name=instrument_name(like),
        gaussian_process=like.gaussian_process,
        variables=(like.priors, like.derived)
    )
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