# This code provides helpers for loading data from the HIPPARCOS GAIA Catalog of Accelerations.
# 

using FITSIO
using Tables
using Measurements

"""
    ProperMotionAnomHGCA(;gaia_id=1234)

Load proper motion anomaly data from the HIPPARCOS-GAIA Catalog of Accelerations (Brandt et al)
for a star with catalog id `gaia_id`.
The resulting velocities are in mas/yr and have the long term trend between HIPPARCOS and GAIA
already subtracted out. e.g. we would expect 0 pma if there is no companion.
"""
function ProperMotionAnomHGCA(;gaia_id,catalog=(datadep"HGCA_eDR3")*"/HGCA_vEDR3.fits")

    ## Load the Hipparcos-GAIA catalog of accelerations
    hgca = FITS(catalog) do hdu
        Tables.columntable(hdu[2])
    end

    idx = findfirst(==(gaia_id), hgca.gaia_source_id)

    # Proper motion anomaly
    # The difference between the ~instant proper motion measured by GAIA compared to the 
    # long term trend between Hipparcos and GAIA
    Δμ_gaia_ra = (hgca.pmra_gaia[idx] ± hgca.pmra_gaia_error[idx]) - (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx])
    Δμ_gaia_dec = (hgca.pmdec_gaia[idx] ± hgca.pmdec_gaia_error[idx]) - (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx])

    Δμ_hip_ra = (hgca.pmra_hip[idx] ± hgca.pmra_hip_error[idx]) - (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx])
    Δμ_hip_dec = (hgca.pmdec_hip[idx] ± hgca.pmdec_hip_error[idx]) - (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx])


    return ProperMotionAnom(
        # Hipparcos epoch
        (;
            ra_epoch=years2mjd(hgca.epoch_ra_hip[idx]),
            dec_epoch=years2mjd(hgca.epoch_dec_hip[idx]),

            pm_ra=Measurements.value(Δμ_hip_ra),
            σ_pm_ra=Measurements.uncertainty(Δμ_hip_ra),

            pm_dec=Measurements.value(Δμ_hip_dec),
            σ_pm_dec=Measurements.uncertainty(Δμ_hip_dec),  
        ),
        # GAIA epoch
        (;
            ra_epoch=years2mjd(hgca.epoch_ra_gaia[idx]),
            dec_epoch=years2mjd(hgca.epoch_dec_gaia[idx]),

            pm_ra=Measurements.value(Δμ_gaia_ra),
            σ_pm_ra=Measurements.uncertainty(Δμ_gaia_ra),

            pm_dec=Measurements.value(Δμ_gaia_dec),
            σ_pm_dec=Measurements.uncertainty(Δμ_gaia_dec),
            
        ),
    )
end
export ProperMotionAnomHGCA


"""
    gaia_plx(gaia_id=12123)

Get a distribution (TruncatedNormal) of parallax distance in mas of a source with 
GAIA catalog id `gaia_id`.
"""
function gaia_plx(;gaia_id,catalog=(datadep"HGCA_eDR3")*"/HGCA_vEDR3.fits") 
    
    # Load the Hipparcos-GAIA catalog of accelerations as a basic column table
    hgca = FITS(catalog) do hdu
        Tables.columntable(hdu[2])
    end

    idx = findfirst(==(gaia_id), hgca.gaia_source_id)
    return TruncatedNormal(hgca.parallax_gaia[idx,], hgca.parallax_gaia_error[idx,], 0, Inf)
end
export gaia_plx
