using LinearAlgebra
using Interpolations
#=
This file implements a variant of InterferometryLikelihood
that accounts for Fiber positioning throughput loss.
=#
const required_cols_fiber = (required_cols..., )

struct FiberInterferometryLikelihood{TTable<:Table,TInterp} <: AbstractInterferometryLikelihood
    table::TTable
    fiber_coupling_interpolator::TInterp
end
function FiberInterferometryLikelihood(observations...)
    input_table = Table(observations...)
    if :filename ∈ Tables.columnnames(input_table)
        rows = map(_prepare_input_row, eachrow(input_table))
        table = Table(rows)
    else
        table = input_table
    end

    if !issubset(required_cols_fiber, Tables.columnnames(table))
        error("Expected columns $vis_cols")
    end


    # Create an interpolator object that maps separation to fiber coupling efficiency
    @info "Pre-calculating fiber coupling efficiency over grid"
    sep_mas = 0:2:100
    λs = range(extrema(vec(table.eff_wave[1]))..., length=15)
    fiber_coupling = stack([
        fiber_coupling_fraction(sep_mas, λ)
        for λ in λs
    ])
    coupling_interp = LinearInterpolation((sep_mas,λs), fiber_coupling, extrapolation_bc=0.0)

    return FiberInterferometryLikelihood(table, coupling_interp)
end
FiberInterferometryLikelihood(observations::NamedTuple...) = FiberInterferometryLikelihood(observations)
export FiberInterferometryLikelihood

# theta = range(0, 100, length=250)
# Credit: W. Balmer, D. Bakely, and others.
function fiber_coupling_fraction(theta, lambda_w=2.2e-6)
    D = 8
    x = range(-D * 2, D * 2, length=500)#0)
    y = range(-D * 2, D * 2, length=500)#0)
    r = LinearAlgebra.norm.(x, y')
    m = r .< D / 2
    # arcseconds
    phase = reshape(x, :, 1, 1) ./ lambda_w .* reshape(theta, 1, 1, :) * 1e-3 / (180 / pi * 3600) * 2pi
    w_0 = 0.32D
    field_pup = @. m * exp(1im * phase)
    field_fiber = @. exp(-1 * r^2 / (2 * w_0^2))
    Inj = abs2.(sum(field_pup .* field_fiber, dims=(1, 2))) / abs(sum(m .* field_fiber))^2
    return Inj[1,1,:]
end



"""
Visibliitiy modelling likelihood for point sources.
"""
function Octofitter.ln_like(vis::FiberInterferometryLikelihood, θ_system, orbits, num_epochs::Val{L}=Val(length(vis.table))) where {L}

    T = typeof(θ_system.M)
    ll = zero(T)

    # Access the data here: 
    epochs = vis.table.epoch
    band = vis.table.band

    # Add an extra optional uncertainty
    # in quadrature
    σ_cp_jitter = hasproperty(θ_system, :σ_cp_jitter) ? θ_system.σ_cp_jitter : zero(T)

    # Loop through epochs
    for i_epoch in eachindex(epochs)

        epoch = epochs[i_epoch]
        this_band = band[i_epoch]

        index_cps1 = vis.table.index_cps1[i_epoch]
        index_cps2 = vis.table.index_cps2[i_epoch]
        index_cps3 = vis.table.index_cps3[i_epoch]
        use_vis2 = vis.table.use_vis2[i_epoch]

        cps_model = zeros(T, size(vis.table.cps_data[i_epoch][:, 1]))
        cvis_model = zeros(complex(T), size(vis.table.u[i_epoch][:, 1]))

        contrasts = T[getproperty(θ_planet, this_band) for θ_planet in θ_system.planets]
        sols = [orbitsolve(orbits[i_planet], epoch) for i_planet in 1:length(θ_system.planets)]

        throughputs = broadcast(sols, contrasts,  vis.table.eff_wave[i_epoch]') do sol, flux_ratio, wavelength_m
            # Model the fiber as placed at the photocentre of the two bodies
            secondary_offset_mas = projectedseparation(sol)
            # Now calculate throughput loss on the secondary due to it being offset wrt. the 
            # fiber (assumed to be at photocentre)
            fiber_offset_mas  = (flux_ratio*secondary_offset_mas)/(1.0 + flux_ratio)
            coupling =  vis.fiber_coupling_interpolator(fiber_offset_mas, wavelength_m)
            return coupling
        end

        # Loop through wavelengths
        for i_wave in axes(vis.table.u[i_epoch], 2)
            u = @views vis.table.u[i_epoch][:, i_wave]
            v = @views vis.table.v[i_epoch][:, i_wave]
            cps_data = @views vis.table.cps_data[i_epoch][:, i_wave]
            σ_cp = @views vis.table.dcps[i_epoch][:, i_wave]
            vis2_data = @views vis.table.vis2_data[i_epoch][:, i_wave]
            dvis2 = @views vis.table.dvis2[i_epoch][:, i_wave]

            # to normalize complex visibilities 
            cvis_model .= 0
            cps_model .= 0
            norm_factor_model = zero(T)

            # Consider all planets
            for i_planet in eachindex(orbits)
                # All parameters relevant to this planet
                # Get model contrast parameter in this band (band provided as a symbol, e.g. :L along with data in table row.)
                contrast = contrasts[i_planet]
                throughput = throughputs[i_planet] 
                Δra = raoff(sols[i_planet])  # in mas
                Δdec = decoff(sols[i_planet]) # in mas

                # add complex visibilities from all planets at a single epoch, for this wavelength
                cvis_bin!(cvis_model; Δdec, Δra, contrast=contrast*throughput, u, v)
                norm_factor_model += contrast
            end
            cvis_model .+= 1.0 #add contribution from the primary primary
            cvis_model .*= 1.0 / (1.0 + norm_factor_model)
            # Compute closure phases
            closurephase!(cps_model; vis=cvis_model, index_cps1, index_cps2, index_cps3)
            if use_vis2
                #compute squared visibilities
                # TODO: optimize
                vis2 = abs.(cvis_model) .^ 2
                const_v2 = -sum(log.(2π * (dvis2 .^ 2))) / 2
                #calculate vis2 ln likelihood
                lnlike_v2 = lnlike_v2 .+ -0.5 * sum((vis2_data .- vis2) .^ 2 ./ dvis2 .^ 2) .+ const_v2
                ll += lnlike_v2
            end

            # Calculate cp ln likelihood
            const_cp = zero(eltype(σ_cp))
            for I in eachindex(σ_cp)
                const_cp -= log(2π * (σ_cp[I] ^ 2 + σ_cp_jitter^2))
            end
            const_cp /= 2
            lnlike_cp = zero(T)
            for I in eachindex(σ_cp, cps_data, cps_model)
                σ² = σ_cp_jitter^2 + σ_cp[I]^2
                lnlike_cp -= 0.5 * (cps_data[I] - cps_model[I])^2 / σ²
            end
            lnlike_cp += const_cp


            # Accumulate into likelihood
            ll += lnlike_cp
        end
    end


    return ll
end




# Generate new observations for a system of possibly multiple planets
function Octofitter.generate_from_params(vis_input::FiberInterferometryLikelihood, θ_system, orbits::Vector{<:AbstractOrbit})

    vis = deepcopy(vis_input)

    T = typeof(θ_system.M)
    ll = zero(T)

    # Access the data here: 
    epochs = vis.table.epoch
    band = vis.table.band

    # Loop through epochs
    for i_epoch in eachindex(epochs)

        epoch = epochs[i_epoch]
        this_band = band[i_epoch]

        index_cps1 = vis.table.index_cps1[i_epoch]
        index_cps2 = vis.table.index_cps2[i_epoch]
        index_cps3 = vis.table.index_cps3[i_epoch]
        use_vis2 = vis.table.use_vis2[i_epoch]

        cps_model = zeros(T, size(vis.table.cps_data[i_epoch][:, 1]))
        cvis_model = zeros(complex(T), size(vis.table.u[i_epoch][:, 1]))

        contrasts = T[getproperty(θ_planet, this_band) for θ_planet in θ_system.planets]
        sols = [orbitsolve(orbits[i_planet], epoch) for i_planet in 1:length(θ_system.planets)]

        throughputs = broadcast(sols, contrasts,  vis.table.eff_wave[i_epoch]') do sol, flux_ratio, wavelength_m
            # Model the fiber as placed at the photocentre of the two bodies
            secondary_offset_mas = projectedseparation(sol)
            # Now calculate throughput loss on the secondary due to it being offset wrt. the 
            # fiber (assumed to be at photocentre)
            fiber_offset_mas  = (flux_ratio*secondary_offset_mas)/(1.0 + flux_ratio)
            coupling =  vis.fiber_coupling_interpolator(fiber_offset_mas, wavelength_m)
            return coupling
        end

        # Loop through wavelengths
        for i_wave in axes(vis.table.u[i_epoch], 2)
            u = @views vis.table.u[i_epoch][:, i_wave]
            v = @views vis.table.v[i_epoch][:, i_wave]
            cps_data = @views vis.table.cps_data[i_epoch][:, i_wave]
            σ_cp = @views vis.table.dcps[i_epoch][:, i_wave]
            vis2_data = @views vis.table.vis2_data[i_epoch][:, i_wave]
            dvis2 = @views vis.table.dvis2[i_epoch][:, i_wave]

            # to normalize complex visibilities 
            cvis_model .= 0
            cps_model .= 0
            norm_factor_model = zero(T)

            # Consider all planets
            for i_planet in eachindex(orbits)
                # All parameters relevant to this planet
                # Get model contrast parameter in this band (band provided as a symbol, e.g. :L along with data in table row.)
                contrast = contrasts[i_planet]
                throughput = throughputs[i_planet] 
                Δra = raoff(sols[i_planet])  # in mas
                Δdec = decoff(sols[i_planet]) # in mas

                # add complex visibilities from all planets at a single epoch, for this wavelength
                cvis_bin!(cvis_model; Δdec, Δra, contrast=contrast*throughput, u, v)
                norm_factor_model += contrast
            end
            cvis_model .+= 1.0 #add contribution from the primary primary
            cvis_model .*= 1.0 / (1.0 + norm_factor_model)

            # Compute closure phases
            closurephase!(cps_model; vis=cvis_model, index_cps1, index_cps2, index_cps3)

            cps_data .= cps_model
            vis2_data .= abs.(cvis_model) .^ 2

        end
    end
    return vis
end