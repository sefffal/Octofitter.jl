# ---------------------------------------------------
# OI-FITS input
#
# This reads the instrument's own data products and knows nothing about
# orbits, so the de-primarying work does not touch it. Both the plain
# closure-phase and the kernel-phase models consume the same rows, which is one
# of the two reasons the two v1 types merged (the other being the position/flux
# front-end).
#
# What did change: v1 read the binary tables through OIFITS.jl's
# `read(OIDataBlock, ::FITSIO.HDU)`. OIFITS 2.0 moved to AstroFITS and dropped
# that method, so with the OIFITS version this package resolves to today, *every*
# file failed to load. The columns are read straight from FITSIO here instead,
# which also preserves v1's data-block *selection* — prefer `EXTVER = 10`, the
# GRAVITY science-combiner convention, and fall back to the first block of each
# kind — something the replacement `OIDataSet` API cannot express.
# ---------------------------------------------------

"""
    _prepare_input_row(row)

Read one OI-FITS file into the row layout the likelihood expects: `u`/`v` in
inverse wavelengths, closure phases and squared visibilities laid out as
(observable × wavelength), and the three index vectors that say which
baselines form each closure triangle.

Optional fields on the input row:

  - `wavelength_min_meters`, `wavelength_max_meters` — restrict the spectral
    channels kept.
  - `oi_extver` — which OI-FITS data-block version to prefer (default 10,
    GRAVITY's science-combiner convention). Falls back to the first block of
    each kind.
"""
function _prepare_input_row(row)
    row = only(row)
    (; wavelength_min_meters, wavelength_max_meters, oi_extver) =
        (; wavelength_min_meters=-Inf, wavelength_max_meters=Inf, oi_extver=10, row...)
    FITS(row.filename, "r") do f
        wavs = _oi_hdu(f, "OI_WAVELENGTH", oi_extver, row.filename)
        vis2s = _oi_hdu(f, "OI_VIS2", oi_extver, row.filename)
        cps = _oi_hdu(f, "OI_T3", oi_extver, row.filename)

        # read data
        eff_wave = collect(Float64, read(wavs, "EFF_WAVE"))
        Λ = length(eff_wave)
        vis2 = _channelwise(read(vis2s, "VIS2DATA"), Λ, "VIS2DATA")
        vis2_err = _channelwise(read(vis2s, "VIS2ERR"), Λ, "VIS2ERR")
        ut = read(vis2s, "UCOORD")
        vt = read(vis2s, "VCOORD")
        vis2_index = Int64.(read(vis2s, "STA_INDEX"))
        cp = _channelwise(read(cps, "T3PHI"), Λ, "T3PHI")
        cp_err = _channelwise(read(cps, "T3PHIERR"), Λ, "T3PHIERR")
        cp_index = Int64.(read(cps, "STA_INDEX"))
        # convert u,v to units of wavelength
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
        return (;
            row...,
            row.epoch,
            use_vis2=hasproperty(row, :use_vis2) ? row.use_vis2 : false,
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

# v1's data-block preference, kept: GRAVITY writes its science-combiner blocks
# with EXTVER = 10, and picking the first HDU instead would silently select a
# polarization channel.
function _oi_hdu(f, name, extver, filename)
    try
        return f[name, extver]
    catch err
        err isa InterruptException && rethrow()
    end
    try
        return f[name]
    catch err
        err isa InterruptException && rethrow()
        throw(KeyError("Could not find an $name extension in $filename"))
    end
end

# FITSIO returns a column with one element per row as a `Vector` and a column
# with a per-row array as a `Matrix`, so a single-channel file and a
# multi-channel one come back with different ranks. Normalize to
# (channel × row), which is what the rest of this file assumes.
_channelwise(x::AbstractMatrix, Λ, name) = size(x, 1) == Λ ? x : _err_shape(x, Λ, name)
_channelwise(x::AbstractVector, Λ, name) = Λ == 1 ? reshape(x, 1, :) : _err_shape(x, Λ, name)
@noinline _err_shape(x, Λ, name) = error(
    "OI-FITS column $name has shape $(size(x)), which is not $Λ spectral channels " *
    "by however many rows. The OI_WAVELENGTH block and the data block probably " *
    "belong to different instruments; select the right one with `oi_extver`.")

"""
    cp_indices(; vis2_index, cp_index)

Extract, for each closure triangle, the three baseline indices whose phases
combine as `+1 +2 -3` to give that closure phase.
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
