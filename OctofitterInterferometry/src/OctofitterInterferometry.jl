module OctofitterInterferometry

using Octofitter
using PlanetOrbits
using Tables, TypedTables

using FITSIO, OIFITS

abstract type AbstractInterferometryLikelihood <: Octofitter.AbstractLikelihood end

const required_cols = (:epoch, :u, :v, :cps_data, :dcps, :vis2_data, :dvis2, :index_cps1, :index_cps2, :index_cps3, :spectrum_var, :use_vis2)
struct InterferometryLikelihood{TTable<:Table} <: AbstractInterferometryLikelihood
    table::TTable
end
function InterferometryLikelihood(
    observations...;
)
    input_table = Table(observations...)
    if :filename ∈ Tables.columnnames(input_table)
        rows = map(eachrow(input_table)) do row
            _prepare_input_row(row)
        end
        table = Table(rows)
    else
        table = input_table
    end

    if !issubset(required_cols, Tables.columnnames(table))
        error("Expected columns $vis_cols")
    end
    return InterferometryLikelihood{typeof(table)}(table)
end
InterferometryLikelihood(observations::NamedTuple...) = InterferometryLikelihood(observations)
export InterferometryLikelihood

# Prepare closure phases etc
function _prepare_input_row(row)
    row = only(row)
    (;wavelength_min_meters, wavelength_max_meters) = (;wavelength_min_meters=-Inf, wavelength_max_meters=Inf, row...)
    FITS(row.filename, "r") do f
        local wavs, vis2s, cps
        try
            wavs = read(OIDataBlock, f["OI_WAVELENGTH", 10])
            vis2s = read(OIDataBlock, f["OI_VIS2", 10])
            cps = read(OIDataBlock, f["OI_T3", 10])
        catch
            try
                wavs = read(OIDataBlock, f["OI_WAVELENGTH"])
                vis2s = read(OIDataBlock, f["OI_VIS2"])
                cps = read(OIDataBlock, f["OI_T3"])
            catch
                throw(KeyError("Could not find keys OI_WAVELENGTH, OI_VIS2, and OI_T3 in $(row.filename)"))
            end
        end
        # read data
        eff_wave = wavs.eff_wave
        vis2 = vis2s.vis2data
        vis2_err = vis2s.vis2err
        ut = vis2s.ucoord
        vt = vis2s.vcoord
        vis2_index = vis2s.sta_index
        cp = cps.t3phi
        cp_err = cps.t3phierr
        cp_index = cps.sta_index
        #convert u,v to units of wavelength
        u = ut ./ eff_wave' # Units of inverse wavelength
        v = vt ./ eff_wave' # Units of inverse wavelength

        # Clamp CP err to a minimum of 2 degrees
        if any(==(0), cp_err)
            @warn "Some closure phase errors are exactly 0. This will lead to numerical issues. Either verify the data, or provide a non-zero `σ_cp_jitter` variable when sampling."
            @warn "clamping uncertainties to at least 2 degrees"
            cp_err .= max.(2, cp_err)
        end

        mask = trues(length(eff_wave))
        mask .= wavelength_min_meters .< eff_wave .< wavelength_max_meters

        # These say what baseline (cp1) should be added to (cp2) and then subtract (cp3)
        # to get a closure phase in our modelling.
        cp_inds1, cp_inds2, cp_inds3 = cp_indices(; vis2_index, cp_index)
        # @warn "WARN!! Inflating errors for debug 25deg"
        return (;
            row...,
            row.epoch,
            row.spectrum_var,
            use_vis2 = hasproperty(row, :use_vis2) ? row.use_vis2 : nothing,
            u=u[:, mask],
            v=v[:, mask],
            eff_wave=eff_wave[mask],
            cps_data=transpose(cp)[:, mask],
            dcps=transpose(cp_err)[:, mask],
            vis2_data=transpose(vis2)[:, mask],
            dvis2=transpose(vis2_err)[:, mask],
            index_cps1=cp_inds1,
            index_cps2=cp_inds2,
            index_cps3=cp_inds3,
        )
    end
end

"""
Visibliitiy modelling likelihood for point sources.
"""
function Octofitter.ln_like(vis::InterferometryLikelihood, θ_system, orbits,     orbit_solutions,
    orbit_solutions_i_epoch_start)

    T = typeof(θ_system.M)
    ll = zero(T)

    # Access the data here: 
    epochs = vis.table.epoch
    spectrum_var = vis.table.spectrum_var

    # Add an extra optional uncertainty
    # in quadrature
    σ_cp_jitter = hasproperty(θ_system, :σ_cp_jitter) ? θ_system.σ_cp_jitter : zero(T)

    # Loop through epochs
    for i_epoch in eachindex(epochs)

        epoch = epochs[i_epoch]
        this_spectrum_var = spectrum_var[i_epoch]

        index_cps1 = vis.table.index_cps1[i_epoch]
        index_cps2 = vis.table.index_cps2[i_epoch]
        index_cps3 = vis.table.index_cps3[i_epoch]
        use_vis2 = vis.table.use_vis2[i_epoch]

        cps_model = zeros(T, size(vis.table.cps_data[i_epoch][:, 1]))
        cvis_model = zeros(complex(T), size(vis.table.u[i_epoch][:, 1]))

        contrasts = T[getproperty(θ_planet, this_spectrum_var) for θ_planet in θ_system.planets]

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
                Δra = raoff(orbit_solutions[i_planet][i_epoch+orbit_solutions_i_epoch_start])  # in mas
                Δdec = decoff(orbit_solutions[i_planet][i_epoch+orbit_solutions_i_epoch_start]) # in mas

                # add complex visibilities from all planets at a single epoch, for this wavelength
                cvis_bin!(cvis_model; Δdec, Δra, contrast, u, v)
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



function cvis_bin!(cvis; Δdec, Δra, contrast, u, v)
    #u,v: baselines [wavelengths]
    #Δdec: dec offset [mas]
    #Δra: ra offset [mas]
    #contrast: secondary/primary  contrast 
    ##################################
    #returns complex visibilities of a point source at position Δdec,Δra

    l2 = contrast

    # Threads.@threads 
    for I in eachindex(cvis, u, v)
        arg = -2π * (u[I] * Δra + v[I] * Δdec) * π / (180 * 3600 * 1000)
        if !isfinite(arg)
            # @warn "non finite complex vis (maxlog=5)" l2  Δra Δdec maxlog=5
            cvis[I] = NaN
            continue
        end
        s, c = sincos(arg)
        cvis[I] += l2 * (c + s * im)
    end

    return cvis
end


function closurephase!(cp_out; vis::AbstractArray, index_cps1::AbstractArray, index_cps2::AbstractArray, index_cps3::AbstractArray)
    #vis: complex visibilities
    #i_cps1,i_cps2,i_cps3: closure triangle indices 
    ##################################
    #returns closure phases [degrees]

    # visphi = rad2deg.(atan.(imag.(vis), real.(vis))) #convert to degrees
    # visphi = mod.(visphi .+ 10980.0, 360.0) .- 180.0 #phase in range (-180,180]

    # visphi = rad2deg.(rem2pi.(atan.(imag.(vis), real.(vis)), RoundNearest)) #convert to degrees
    # cp_out .= @views visphi[index_cps1, :] .+ visphi[index_cps2, :] .- visphi[index_cps3, :]

    visphi(vis) = rad2deg.(rem2pi.(atan.(imag.(vis), real.(vis)), RoundNearest))

    for i_cp in eachindex(index_cps1, index_cps2, index_cps3)
        for i_obs in axes(cp_out, 2)
            visphi1 = visphi(vis[index_cps1[i_cp], i_obs])
            visphi2 = visphi(vis[index_cps2[i_cp], i_obs])
            visphi3 = visphi(vis[index_cps3[i_cp], i_obs])
            cp_out[i_cp, i_obs] = visphi1 + visphi2 - visphi3
        end
    end

    # cp_out .= @views visphi[index_cps1, :] .+ visphi[index_cps2, :] .- visphi[index_cps3, :]
    return cp_out
end

"""
Extracts indices for calculating closure phases from visibility and closure phase station indices
"""
function cp_indices(; vis2_index::Matrix{<:Int64}, cp_index::Matrix{<:Int64})
    i_cps1 = zeros(Int64, size(cp_index)[2])
    i_cps2 = zeros(Int64, size(cp_index)[2])
    i_cps3 = zeros(Int64, size(cp_index)[2])

    nh = maximum(vis2_index) #number of stations making up your interferometer
    nb = Int64(nh * (nh - 1) / 2) #number of baselines
    ncp = Int64(nh * (nh - 1) * (nh - 2) / 6) #total number of closure phases

    for i in range(1, size(cp_index)[2]), j in range(1, size(vis2_index)[2])
        if cp_index[1, i] == vis2_index[1, j] && cp_index[2, i] == vis2_index[2, j]
            if floor((j - 1) / nb) == floor((i - 1) / ncp)
                i_cps1[i] = j
            end
        end
        if cp_index[2, i] == vis2_index[1, j] && cp_index[3, i] == vis2_index[2, j]
            if floor((j - 1) / nb) == floor((i - 1) / ncp)
                i_cps2[i] = j
            end
        end
        if cp_index[1, i] == vis2_index[1, j] && cp_index[3, i] == vis2_index[2, j]
            if floor((j - 1) / nb) == floor((i - 1) / ncp)
                i_cps3[i] = j
            end
        end
    end
    return i_cps1, i_cps2, i_cps3
end

# Generate new observations for a system of possibly multiple planets
function Octofitter.generate_from_params(like::InterferometryLikelihood, θ_system, orbits::Vector{<:AbstractOrbit})

    # # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch
    bands = like.table.band
    u = like.table.u
    v = like.table.v
    cps_data = like.table.cps_data
    dcps = like.table.dcps
    vis2_data = like.table.vis2_data
    dvis2 = like.table.dvis2
    i_cps1 = like.table.index_cps1
    i_cps2 = like.table.index_cps2
    i_cps3 = like.table.index_cps3
    use_vis2 = like.table.use_vis2
    cp_all = []
    vis2_all = []
    for j in eachindex(epochs)
        band = bands[j]
        epoch = epochs[j]
        complexvis_model = zeros(Complex{Float64}, length(u[j]))
        norm_factor = 0.0 # to normalize complex visibilities 
        for i in eachindex(orbits)
            # All parameters relevant to this planet
            θ_planet = θ_system.planets[i]

            # orbit object pre-created from above parameters (shared between all likelihood functions)
            orbit = orbits[i]
            contrast = getproperty(θ_planet, band) #secondary/primary
            sol = orbitsolve(orbit, epoch)
            Δra = raoff(sol)  # in mas
            Δdec = decoff(sol) # in mas
            # Accumulate into cvis
            cvis_bin!(complexvis_model; Δdec, Δra, contrast, u=u[j], v=v[j])
            norm_factor += contrast
        end
        complexvis_model .+= 1.0 #add contribution from the primary 
        complexvis_model ./= 1.0 + norm_factor
        # compute closure phase
        cp_model = closurephase(vis=complexvis_model, index_cps1=i_cps1[j], index_cps2=i_cps2[j], index_cps3=i_cps3[j])
        #compute squared visibilities
        vis2_model = abs.(complexvis_model) .^ 2

        append!(cp_all, [cp_model])
        append!(vis2_all, [vis2_model])
    end
    cp_all = identity.(cp_all)
    vis2_all = identity.(vis2_all)
    new_vis_like_table = Table(epoch=epochs, u=u, v=v, cps_data=cp_all, dcps=dcps,
        vis2_data=vis2_all, dvis2=dvis2, index_cps1=i_cps1,
        index_cps2=i_cps2, index_cps3=i_cps3, band=bands, use_vis2=use_vis2)


    # return with same number of rows: band, epoch
    # position(s) of point sources according to orbits, θ_system
    return InterferometryLikelihood(new_vis_like_table)
end

include("GRAVITY.jl")
end