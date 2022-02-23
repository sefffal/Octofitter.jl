# This code provides helpers for loading data from the HIPPARCOS GAIA Catalog of Accelerations.
# 

using FITSIO
using Measurements

"""
    ProperMotionAnomHGCA(;gaia_id=1234)

Load proper motion anomaly data from the HIPPARCOS-GAIA Catalog of Accelerations (Brandt et al)
for a star with catalog id `gaia_id`.
The resulting velocities are in mas/yr and have the long term trend between HIPPARCOS and GAIA
already subtracted out. e.g. we would expect 0 pma if there is no companion.
"""
function ProperMotionAnomHGCA2(;gaia_id,catalog=(datadep"HGCA_eDR3")*"/HGCA_vEDR3.fits")

    ## Load the Hipparcos-GAIA catalog of accelerations
    hgca = FITS(catalog) do hdu
        Table(hdu[2])
    end

    # Find the row with a 
    idx = findfirst(==(gaia_id), hgca.gaia_source_id)

    # # Proper motion anomaly
    # # The difference between the ~instant proper motion measured by GAIA compared to the 
    # # long term trend between Hipparcos and GAIA
    # Δμ_gaia_ra = (hgca.pmra_gaia[idx] ± hgca.pmra_gaia_error[idx]) - (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx])
    # Δμ_gaia_dec = (hgca.pmdec_gaia[idx] ± hgca.pmdec_gaia_error[idx]) - (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx])

    # Δμ_hip_ra = (hgca.pmra_hip[idx] ± hgca.pmra_hip_error[idx]) - (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx])
    # Δμ_hip_dec = (hgca.pmdec_hip[idx] ± hgca.pmdec_hip_error[idx]) - (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx])

    return ProperMotionAnomHGCA2(hgca[idx])


    return ProperMotionAnomHGCA2(
        # Hipparcos epoch
        (;
            ra_epoch=years2mjd(hgca.epoch_ra_hip[idx]),
            dec_epoch=years2mjd(hgca.epoch_dec_hip[idx]),
            dt=4*365,
            
            pm_ra=hgca.pmra_hip[idx],
            σ_pm_ra=hgca.pmra_hip_error[idx],

            pm_dec=hgca.pmdec_hip[idx],
            σ_pm_dec=hgca.pmdec_hip_error[idx],
        ),
        # Hipparcos - GAIA
        (;
            ra_epoch=(years2mjd(hgca.epoch_ra_hip[idx]) + years2mjd(hgca.epoch_ra_gaia[idx]))/2,
            dec_epoch=(years2mjd(hgca.epoch_dec_hip[idx]) + years2mjd(hgca.epoch_dec_gaia[idx]))/2,
            dt=years2mjd(hgca.epoch_ra_gaia[idx]) - years2mjd(hgca.epoch_ra_hip[idx]),

            pm_ra=hgca.pmra_hg[idx],
            σ_pm_ra=hgca.pmra_hg_error[idx],

            pm_dec=hgca.pmdec_hg[idx],
            σ_pm_dec=hgca.pmdec_hg_error[idx],
        ),
        # GAIA epoch
        (;
            ra_epoch=years2mjd(hgca.epoch_ra_gaia[idx]),
            dec_epoch=years2mjd(hgca.epoch_dec_gaia[idx]),
            dt=3*365,
            pm_ra=hgca.pmra_gaia[idx],
            σ_pm_ra=hgca.pmra_gaia_error[idx],
            pm_dec=hgca.pmdec_gaia[idx],
            σ_pm_dec=hgca.pmdec_gaia_error[idx],
        ),
    )
    # ProperMotionAnomHG(
    #     # Hipparcos - GAIA
    #     (;
    #         ra_epoch_1=years2mjd(hgca.epoch_ra_hip[idx]),
    #         ra_epoch_2=years2mjd(hgca.epoch_ra_gaia[idx]),
    #         dec_epoch_1=years2mjd(hgca.epoch_dec_hip[idx]),
    #         dec_epoch_2=years2mjd(hgca.epoch_dec_gaia[idx]),
    #         dt_1=4*365,
    #         dt_2=3*365,

    #         pm_ra=hgca.pmra_hg[idx],
    #         σ_pm_ra=hgca.pmra_hg_error[idx],

    #         pm_dec=hgca.pmdec_hg[idx],
    #         σ_pm_dec=hgca.pmdec_hg_error[idx],
    #     ),
    # )

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


# function ghca_pmra(;gaia_id)