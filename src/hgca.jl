

using FITSIO
using Tables
using Measurements




function ProperMotionAnomHGCA(;gaia_id,catalog=(datadep"HGCA_eDR3")*"/HGCA_vEDR3.fits")

    ## Load the Hipparcos-GAIA catalog of accelerations
    hgca = FITS(catalog) do hdu
        Tables.columntable(hdu[2])
    end

    idx = findfirst(==(gaia_id), hgca.gaia_source_id)

    # Proper motion anomaly
    # The difference between the ~instant proper motion measured by GAIA compared to the 
    # long term trend between Hipparcos and GAIA
    # Δμ_gaia_ra = (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx]) - (hgca.pmra_gaia[idx] ± hgca.pmra_gaia_error[idx])
    # Δμ_gaia_dec = (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx]) - (hgca.pmdec_gaia[idx] ± hgca.pmdec_gaia_error[idx])

    # Δμ_hip_ra = (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx]) - (hgca.pmra_hip[idx] ± hgca.pmra_hip_error[idx])
    # Δμ_hip_dec = (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx]) - (hgca.pmdec_hip[idx] ± hgca.pmdec_hip_error[idx])


    # The above was backwards. We want to subtract the long term trend from each of the "instant" measurements
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


function gaia_plx(;gaia_id,catalog=(datadep"HGCA_eDR3")*"/HGCA_vEDR3.fits") 
    
    ## Load the Hipparcos-GAIA catalog of accelerations
    hgca = FITS(catalog) do hdu
        Tables.columntable(hdu[2])
    end

    idx = findfirst(==(gaia_id), hgca.gaia_source_id)
    return Normal(hgca.parallax_gaia[idx,], hgca.parallax_gaia_error[idx,])
end
export gaia_plx
