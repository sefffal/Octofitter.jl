"""
    PlanetRelativeRVObs(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);

        name="inst name",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100.0)  # RV jitter (m/s)
        end
    )

    # Example with Gaussian Process:
    PlanetRelativeRVObs(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);

        name="inst name",
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
struct PlanetRelativeRVObs{TTable<:Table,GP} <: Octofitter.AbstractObs
    table::TTable
    name::String
    gaussian_process::GP
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    function PlanetRelativeRVObs(
        observations;
        name::String,
        gaussian_process=nothing,
        variables::Union{Nothing,Tuple{Octofitter.Priors,Octofitter.Derived}}=nothing,
    )    
        if isnothing(variables)
            variables = @variables begin
                jitter ~ LogUniform(0.001, 100)
            end
            @info "Added the following observation variables:"
            display(variables[1])
            display(variables[2])
        end
        (priors,derived)=variables
        table = Table(observations)
        if !Octofitter.equal_length_cols(table)
            error("The columns in the input data do not all have the same length")
        end
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
       
        # We sort the data first by instrument index then by epoch to make some later code faster
        ii = sortperm(table.epoch)
        table = table[ii,:]

        if any(>=(mjd("2050")),  table.epoch) || any(<=(mjd("1950")),  table.epoch)
            @warn "The data you entered fell outside the range year 1950 to year 2050. The expected input format is MJD (modified julian date). We suggest you double check your input data!"
        end

        rows = map(eachrow(table)) do row′
            row = (;row′[1]..., rv=float(row′[1].rv[1]))
            return row
        end
        table = Table(rows)
        ii = sortperm(table.epoch)
        table = table[ii]

        return new{typeof(table),typeof(gaussian_process)}(table, name, gaussian_process, priors, derived)
    end
end
PlanetRelativeRVObs(observations::NamedTuple...;kwargs...) = PlanetRelativeRVObs(observations; kwargs...)
function Octofitter.likeobj_from_epoch_subset(obs::PlanetRelativeRVObs, obs_inds)
    return PlanetRelativeRVObs(
        obs.table[obs_inds,:,1]...;
        name=likelihoodname(obs),
        gaussian_process=obs.gaussian_process,
        variables=(obs.priors, obs.derived)
    )
end

# Backwards compatibility
const PlanetRelativeRVLikelihood = PlanetRelativeRVObs

export PlanetRelativeRVObs, PlanetRelativeRVLikelihood


# In-place simulation logic for PlanetRelativeRVObs (performance-critical)
function Octofitter.simulate!(rv_model_buf, rvlike::PlanetRelativeRVObs, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    this_orbit = orbits[i_planet]
    
    # Data for this instrument:
    epochs = vec(rvlike.table.epoch)
    
    # Compute the model RV values (what we expect to observe)
    fill!(rv_model_buf, 0)
    
    # Add RV contribution from this planet and any inner planets:
    for i_epoch in eachindex(epochs)
        sol = orbit_solutions[i_planet][i_epoch+orbit_solutions_i_epoch_start]
        @assert isapprox(rvlike.table.epoch[i_epoch], PlanetOrbits.soltime(sol), rtol=1e-2)
        # Relative RV due to planet
        rv_model_buf[i_epoch] += radvel(sol)

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
                
                rv_model_buf[i_epoch] += radvel(sol′, mass_other)
                
                @assert isapprox(rvlike.table.epoch[i_epoch], PlanetOrbits.soltime(sol′), rtol=1e-2)
            end
        end
    end
    
    return (rv_model = rv_model_buf, epochs = epochs)
end

# Allocating simulation logic for PlanetRelativeRVObs (convenience method)
function Octofitter.simulate(rvlike::PlanetRelativeRVObs, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_planet)
    rv_model_buf = Vector{T}(undef, length(rvlike.table.epoch))
    return Octofitter.simulate!(rv_model_buf, rvlike, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
end


"""
Radial velocity likelihood.
"""
function Octofitter.ln_like(rvlike::PlanetRelativeRVObs, ctx::Octofitter.PlanetObservationContext)
    (; θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start) = ctx
    T = Octofitter._system_number_type(θ_planet)
    ll = zero(T)
    
    # Data for this instrument:
    epochs = vec(rvlike.table.epoch)
    σ_rvs = vec(rvlike.table.σ_rv)
    rvs = vec(rvlike.table.rv)
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)

    @no_escape begin
        # Allocate buffers using bump allocator
        rv_model_buf = @alloc(T, length(epochs))
        rv_residuals = @alloc(T, length(epochs))
        noise_var = @alloc(T, length(epochs))
        
        # Use in-place simulation method to get model values
        sim = Octofitter.simulate!(rv_model_buf, rvlike, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)

        # Compute residuals: observed - model
        rv_residuals .= rvs .- rv_model_buf

        # The noise variance per observation is the measurement noise and the jitter added
        # in quadrature
        noise_var .= σ_rvs.^2 .+ jitter^2

        # Two code paths, depending on if we are modelling the residuals by 
        # a Gaussian process or not.
        if isnothing(rvlike.gaussian_process)
            # Don't fit a GP
            fx = MvNormal(Diagonal((noise_var)))
            ll += logpdf(fx, rv_residuals)
        else
            # Fit a GP
            local gp
            try
                gp = @inline rvlike.gaussian_process(θ_obs)
                fx = gp(epochs, noise_var)
                ll += logpdf(fx, rv_residuals)
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
function Octofitter.generate_from_params(like::PlanetRelativeRVObs, θ_planet, elem::PlanetOrbits.AbstractOrbit; add_noise)

    # Get epochs and uncertainties from observations
    epochs = like.table.epoch
    σ_rvs = like.table.σ_rv

    # Generate new planet radial velocity data
    rvs = radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    if add_noise
        radvel_table.rv .+= randn.() .* σ_rvs
    end

    return PlanetRelativeRVObs(
        radvel_table;
        name=likelihoodname(like),
        gaussian_process=like.gaussian_process,
        variables=(like.priors, like.derived)
    )
end


# Generate new radial velocity observations for a star
function generate_from_params(like::PlanetRelativeRVObs, θ_system,  θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start; add_noise)

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch
    σ_rvs = like.table.σ_rv

    # Generate new star radial velocity data
    rvs = radvel.(reshape(orbits, :, 1), epochs)
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    if add_noise
        jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
        radvel_table.rv .+= randn.() .* hypot.(σ_rvs, jitter)
    end

    return PlanetRelativeRVObs(
        radvel_table;
        name=likelihoodname(like),
        gaussian_process=like.gaussian_process,
        variables=(like.priors, like.derived)
    )
end



# # Plot recipe for astrometry data
# @recipe function f(rv::PlanetRelativeRVObs)
   
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