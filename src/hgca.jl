
using DataDeps
using FITSIO
using Tables
using Measurements

function init_datadeps()
    register(DataDep("HGCA_eDR3",
        """
        Dataset: Hipparcos Gaia Catalog of Accelerations
        Author: Brandt et al
        License: arXiv.org Non-exclusive license to distribute
        Website: https://arxiv.org/abs/2105.11662

        A catalog by Brandt et al containing cross calibrated proper motions
        between Hipparcos and GAIA eDR3.

        File size: 19MiB
        """,
        "http://physics.ucsb.edu/~tbrandt/HGCA_vEDR3.fits",
        "23684d583baaa236775108b360c650e79770a695e16914b1201f290c1826065c"
    ))
end

# ## Load the Hipparcos-GAIA catalog of accelerations
# hgca = FITS(raw"C:\Users\William\Downloads\HGCA_vEDR3.fits") do hdu
#     DataFrame(hdu[2])
# end

# # sort!(hgca, :chisq)

# ## Include all stars in the catalog closer than this paralax cut:
# # # Projection
# # xy = aitoff.(hgca.gaia_ra, hgca.gaia_dec)
# # x,y = eachrow(reinterpret(reshape, Float64, xy))

# # a_star = @. hgca.gaia_ra * cosd(hgca.gaia_dec)

# ## Proper motion anomaly
# # The difference between the ~instant proper motion measured by GAIA compared to the 
# # long term trend between Hipparcos and GAIA
# Δμ_gaia_ra = @. (hgca.pmra_hg ± hgca.pmra_hg_error) - (hgca.pmra_gaia ± hgca.pmra_gaia_error)
# Δμ_gaia_dec = @. (hgca.pmdec_hg ± hgca.pmdec_hg_error) - (hgca.pmdec_gaia ± hgca.pmdec_gaia_error)
# Δμ_gaia = @. √(Δμ_gaia_ra^2 + Δμ_gaia_dec^2)

# Δμ_hip_ra = @. (hgca.pmra_hg ± hgca.pmra_hg_error) - (hgca.pmra_hip ± hgca.pmra_hip_error)
# Δμ_hip_dec = @. (hgca.pmdec_hg ± hgca.pmdec_hg_error) - (hgca.pmdec_hip ± hgca.pmdec_hip_error)
# Δμ_hip = @. √(Δμ_hip_ra^2 + Δμ_hip_dec^2)


function ProperMotionAnomHGCA(;gaia_id,catalog=(datadep"HGCA_eDR3")*"/HGCA_vEDR3.fits")

    ## Load the Hipparcos-GAIA catalog of accelerations
    hgca = FITS(catalog) do hdu
        Tables.columntable(hdu[2])
    end

    idx = findfirst(==(gaia_id), hgca.gaia_source_id)

    # Proper motion anomaly
    # The difference between the ~instant proper motion measured by GAIA compared to the 
    # long term trend between Hipparcos and GAIA
    Δμ_gaia_ra = (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx]) - (hgca.pmra_gaia[idx] ± hgca.pmra_gaia_error[idx])
    Δμ_gaia_dec = (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx]) - (hgca.pmdec_gaia[idx] ± hgca.pmdec_gaia_error[idx])

    Δμ_hip_ra = (hgca.pmra_hg[idx] ± hgca.pmra_hg_error[idx]) - (hgca.pmra_hip[idx] ± hgca.pmra_hip_error[idx])
    Δμ_hip_dec = (hgca.pmdec_hg[idx] ± hgca.pmdec_hg_error[idx]) - (hgca.pmdec_hip[idx] ± hgca.pmdec_hip_error[idx])

    return ProperMotionAnom(
        # GAIA epoch
        (;
            ra_epoch=years2mjd(hgca.epoch_ra_gaia[idx]),
            dec_epoch=years2mjd(hgca.epoch_dec_gaia[idx]),

            pm_ra=Measurements.value(Δμ_gaia_ra),
            σ_pm_ra=Measurements.uncertainty(Δμ_gaia_ra),

            pm_dec=Measurements.value(Δμ_gaia_dec),
            σ_pm_dec=Measurements.uncertainty(Δμ_gaia_dec),
            
        ),
        # Hipparcos epoch
        (;
            ra_epoch=years2mjd(hgca.epoch_ra_hip[idx]),
            dec_epoch=years2mjd(hgca.epoch_dec_hip[idx]),

            pm_ra=Measurements.value(Δμ_hip_ra),
            σ_pm_ra=Measurements.uncertainty(Δμ_hip_ra),

            pm_dec=Measurements.value(Δμ_hip_dec),
            σ_pm_dec=Measurements.uncertainty(Δμ_hip_dec),  
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
