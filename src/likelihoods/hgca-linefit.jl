using HTTP: HTTP

struct HGCALinearFitLikelihood2{TTable1<:Table,TTable2<:Table,TTable3<:Table} <: AbstractLikelihood
    # Table of our three HGCA observations
    table::TTable1
    # Table of Hipparcos scans retrieved from IAD
    hip_table::TTable2
    # Table of GAIA scans retrieved from GHOST
    gaia_table::TTable3
end
export HGCALinearFitLikelihood2


"""
    HGCALinearFitLikelihood2(;gaia_id=1234)

Load proper motion anomaly data from the HIPPARCOS-GAIA Catalog of Accelerations (Brandt et al)
for a star with catalog id `gaia_id`.
The resulting velocities are in mas/yr and have the long term trend between HIPPARCOS and GAIA
already subtracted out. e.g. we would expect 0 pma if there is no companion.
"""
function HGCALinearFitLikelihood2(;
    gaia_id,
    hgca_catalog=(datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits",
    hip_catalog=(datadep"Hipparcos_IAD"),
    expand_uncertainties_factor=1)

    ###############################################
    # Load HGCA for this target

    # Load the Hipparcos-GAIA catalog of accelerations (downloaded automatically with datadeps)
    hgca_all = FITS(hgca_catalog, "r") do fits
        Table(fits[2])
    end

    # Available columns (for reference)
    # chisq            crosscal_pmdec_hg  crosscal_pmdec_hip   crosscal_pmra_hg   crosscal_pmra_hip  epoch_dec_gaia          epoch_dec_hip
    # epoch_ra_gaia    epoch_ra_hip       gaia_dec             gaia_ra            gaia_source_id     hip_id                  nonlinear_dpmdec
    # nonlinear_dpmra  parallax_gaia      parallax_gaia_error  pmdec_gaia         pmdec_gaia_error   pmdec_hg                pmdec_hg_error
    # pmdec_hip        pmdec_hip_error    pmra_gaia            pmra_gaia_error    pmra_hg            pmra_hg_error           pmra_hip
    # pmra_hip_error   pmra_pmdec_gaia    pmra_pmdec_hg        pmra_pmdec_hip     radial_velocity    radial_velocity_error   radial_velocity_source

    # Find the row with a GAIA source id match
    idx = findfirst(==(gaia_id), hgca_all.gaia_source_id)
    hgca_all.pmra_hip_error[idx] *= expand_uncertainties_factor
    hgca_all.pmdec_hip_error[idx] *= expand_uncertainties_factor
    hgca_all.pmra_hg_error[idx] *= expand_uncertainties_factor
    hgca_all.pmdec_hg_error[idx] *= expand_uncertainties_factor
    hgca_all.pmra_gaia_error[idx] *= expand_uncertainties_factor
    hgca_all.pmdec_gaia_error[idx] *= expand_uncertainties_factor

    # Convert measurement epochs to MJD.
    # The HGCA doesn't say, but we assume these are actually Julian years and not decimal years.
    J2000_mjd = 51544.5 # year J2000 in MJD
    epoch_ra_hip_mjd = (hgca_all.epoch_ra_hip[idx] - 2000)*julian_year + J2000_mjd
    epoch_dec_hip_mjd = (hgca_all.epoch_dec_hip[idx] - 2000)*julian_year + J2000_mjd
    epoch_ra_gaia_mjd = (hgca_all.epoch_ra_gaia[idx] - 2000)*julian_year + J2000_mjd
    epoch_dec_gaia_mjd = (hgca_all.epoch_dec_gaia[idx] - 2000)*julian_year + J2000_mjd

    hgca = (;
        NamedTuple(hgca_all[idx])...,
        epoch_ra_hip_mjd,
        epoch_dec_hip_mjd,
        epoch_ra_gaia_mjd,
        epoch_dec_gaia_mjd,
    )

    # Hipparcos epoch
    c = hgca.pmra_pmdec_hip[1] * hgca.pmra_hip_error[1] * hgca.pmdec_hip_error[1]
    dist_hip = MvNormal(@SArray[
        hgca.pmra_hip_error[1]^2 c
        c hgca.pmdec_hip_error[1]^2
    ])
    # Hipparcos - GAIA epoch
    c = hgca.pmra_pmdec_hg[1] * hgca.pmra_hg_error[1] * hgca.pmdec_hg_error[1]
    dist_hg = MvNormal(@SArray[
        hgca.pmra_hg_error[1]^2 c
        c hgca.pmdec_hg_error[1]^2
    ])
    # GAIA epoch
    c = hgca.pmra_pmdec_gaia[1] * hgca.pmra_gaia_error[1] * hgca.pmdec_gaia_error[1]
    dist_gaia = MvNormal(@SArray[
        hgca.pmra_gaia_error[1]^2 c
        c hgca.pmdec_gaia_error[1]^2
    ])
    hgca_table = Table([(;hgca...,dist_hip,dist_hg,dist_gaia)])

    ###############################################
    # Load the GHOST data for this target
    gaia_table = FlexTable( GHOST_forecast(;gaia_id, catalog=hgca_catalog) )
    gaia_table.epoch = jd2mjd.(
        gaia_table.var" ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]"
    )
    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    gaia_table.cpsi = cos.(π/2 .+ gaia_table.var" scanAngle[rad]")
    gaia_table.spsi = sin.(π/2 .+ gaia_table.var" scanAngle[rad]")

    ###############################################
    hip_id = hgca.hip_id[]
    # Load the Hipparcos data for this target
    # TODO: factor out shared functionality from hipparcos.jl
    file = @sprintf("H%06d.d", hip_id)
    fname = joinpath(hip_catalog, "ResRec_JavaTool_2014", file[1:4], file)

    lines = readlines(fname)

    # # See table 1 of https://arxiv.org/pdf/1108.4971 for units

    # HIP    MCE    NRES NC isol_n SCE  F2     F1
    # 27321  27251  111  1  5      0    -1.63  0 
    hip, mce, nres, nc, isol_n, sce, f2, f1 = parse.(Float64, split(lines[7])[2:end])

    # # Hp      B-V    VarAnn NOB NR
    # # 3.9077  0.171  0      111 0  
    # hp, b_m_v, varann, nob, nr = parse.(Float64, split(lines[9])[2:end])

    # (
    #     radeg,
    #     dedeg,
    #     plx,
    #     pm_ra,
    #     pm_de,
    #     e_ra,
    #     e_de,
    #     e_plx,
    #     e_pmra,
    #     e_pmde,
    #     dpmra,
    #     dpmde,
    #     e_dpmra,
    #     e_dpmde,
    #     ddpmra,
    #     ddpmde,
    #     e_ddpmra,
    #     e_ddpmde,
    #     upsra,
    #     upsde,
    #     e_upsra,
    #     e_upsde,
    #     var
    # ) = tryparse.(Float64, split(lines[11])[2:end])
    # hip_sol = (;
    #     hip, mce, nres, nc, isol_n, sce, f2, f1,
    #     hp, b_m_v, varann, nob, nr,
    #     radeg, dedeg, plx, pm_ra, pm_de, e_ra, e_de, e_plx,
    #     e_pmra, e_pmde, dpmra, dpmde, e_dpmra, e_dpmde, ddpmra, ddpmde,
    #     e_ddpmra, e_ddpmde, upsra, upsde, e_upsra, e_upsde, var
    # )

    if isol_n != 5
        error("Only stars with solution types 5 are currently supported. If you need solution type 1, please open an issue on GitHub and we would be happy to add it.")
    end

    iad_table_rows = NamedTuple[]
    for line in lines[13:end]
        if startswith(line, '#')
            continue
        end
        iorb, epoch, parf, cpsi, spsi, res, sres = split(line)
        push!(iad_table_rows, (;
            iorb=parse(Int, iorb),
            epoch_yrs=parse(Float64, epoch),
            parf=parse(Float64, parf),
            cpsi=parse(Float64, cpsi),
            spsi=parse(Float64, spsi),
            res=parse(Float64, res),
            sres=parse(Float64, sres),
        ))
    end

    iad_table = FlexTable(iad_table_rows)
    # transform epoch to MJD
    # iad_table.epoch = years2mjd.(hipparcos_catalog_epoch_decimal_year .+ iad_table.epoch_yrs)
    iad_table.epoch = hipparcos_catalog_epoch_mjd .+ iad_table.epoch_yrs*julian_year

    # Remove rejected scans, if any
    iad_table.reject = iad_table.sres .<= 0
    if any(iad_table.reject)
        @warn "rejected scans present" count(iad_table.reject)
    end



    return HGCALinearFitLikelihood2(hgca_table, Table(iad_table), Table(gaia_table),)

end
export HGCALinearFitLikelihood2


"""
Specific HGCA proper motion modelling. Model the GAIA-Hipparcos/Δt proper motion
using 5 position measurements averaged at each of their epochs.
"""
function ln_like(hgca_like::HGCALinearFitLikelihood2, θ_system, elements, _L=1) #=length of observations: we know this is one=#
    ll = 0.0

    (;
        pmra_hip_model,
        pmdec_hip_model,
        pmra_gaia_model,
        pmdec_gaia_model,
        pmra_hg_model,
        pmdec_hg_model,
    ) = _simulate_hgca(hgca_like, θ_system, elements)

    # Hipparcos epoch
    resids_hip = @SArray[
        pmra_hip_model - hgca_like.table.pmra_hip[1],
        pmdec_hip_model - hgca_like.table.pmdec_hip[1]
    ]
    ll += logpdf(hgca_like.table.dist_hip[1], resids_hip)

    # Hipparcos - GAIA epoch

    # TODO: We have to undo the spherical coordinate correction that was done in the HGCA catalog
    # since our calculations use real, not tangent plane, coordinates
    resids_hg = @SArray[
        pmra_hg_model - hgca_like.table.pmra_hg[1]
        pmdec_hg_model - hgca_like.table.pmdec_hg[1]
    ]
    ll += logpdf(hgca_like.table.dist_hg[1], resids_hg)

    # GAIA epoch
    resids_gaia = @SArray[
        pmra_gaia_model - hgca_like.table.pmra_gaia[1],
        pmdec_gaia_model - hgca_like.table.pmdec_gaia[1]
    ]
    ll += logpdf(hgca_like.table.dist_gaia[1], resids_gaia)

    return ll
end


function _simulate_hgca_linefit(pma, θ_system, orbits)

    #=
    Think through how I really want to do this and don't just
    follow the crowd.

    Variables and derived quantities:
    * corrected RA & DEC at any epoch

    Data:
    * visit times and scan angles from Hipparcos and from GAIA
    * the Hipparcos and GAIA PMs, calibrated by the HGCA, corrected for non-linearity and RV
    * the GAIA-Hipparcos total PM, calibrated by the HGCA, corrected for non-linearity and RV

    The best I can do:
    * Model what Hipparcos would have seen, fit to IAD using along-scan residuals for likelihood.
    * Model what GAIA will have seen, compute the best-fitting sky-path, and compare to catalog values
        * We in theory have uncertainties on our fit, and uncertainties on the catalog values. Hmm.
          Maybe it's okay to use compare our best fit to the catalog values using catalog uncertainties.
        * This will in general require optimization as an inner function. That's probably okay but will
          take some care
        * Proceedure:
            1. Take true orbit solution with impacts from planets etc, and calculate the perturbed sky path.
            2. Save positions along the sky path
            3. Against these positions, optimize a model of a zero planet AbsoluteVisual orbit (AbsoluteVisual{FixedPos}?)
            4. compare these optimized values against the catalog values and uncertainties.

        * To do the model I can do nested minimization. Fit a 0 planet AbsoluteVisual model get a sky path
        * What's more, I can do this separately for GAIA DR2 and GAIA DR3! That gets an extra point.
    * Long term trend gets handled automatically, right? Since the H and the G epochs will be pinned to 
      physical positions by our model and by the data.
    * What do we get from the HGCA? 
      * HGCA gives a reference frame correction factor. We could use this to bump the hipparcos data.
      * HGCA gives us an error inflation factor. I should read more about this. 
        But I could handle this also by adding my own jitter terms to the fit.




    How can I structure this?
    * I already have Hipparcos done.
    * I can implement a GAIA likelihood following the Hipparcos approach, but instead of comparing to residuals we compute and compare catalog values
        * This will be easier to test independently
    
    * What do we
        * Implement an HGCA_Detailed?_Likelihod object that can hold 
    =#

    return (;
        pmra_hip_model,
        pmdec_hip_model,
        pmra_gaia_model,
        pmdec_gaia_model,
        pmra_hg_model,
        pmdec_hg_model,
    )
end

