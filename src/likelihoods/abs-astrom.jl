

struct GaiaHipparcosUEVAJointLikelihood_v1{TTableH,TTableG,TCat,fluxratio_var} <: AbstractLikelihood
    hip_table::TTableH
    gaia_table::TTableG
    catalog::TCat
    A_prepared_5_hip::Matrix{Float64}
    A_prepared_5_dr2::Matrix{Float64}
    A_prepared_5_dr3::Matrix{Float64}
    fluxratio_var::Symbol
    include_iad::Bool
    ueva_mode::Symbol
end
GaiaHipparcosUEVAJointLikelihood = GaiaHipparcosUEVAJointLikelihood_v1


function _getparams(::GaiaHipparcosUEVAJointLikelihood_v1{THGCA,THip,TGaia,fluxratio_var}, θ_planet) where {THGCA,THip,TGaia,fluxratio_var}
    if fluxratio_var == :__NA
        return (;fluxratio=zero(Octofitter._system_number_type(θ_planet)))
    end
    fluxratio = getproperty(θ_planet, fluxratio_var)
    return (;fluxratio)
end

function GaiaHipparcosUEVAJointLikelihood_v1(;
        gaia_id,
        fluxratio_var=nothing,
        scanlaw_table=nothing,
        catalog,
        include_iad=false,
        ueva_mode::Symbol=:RUWE,
    )

    # Load the HCGA
    catalog = FITS(catalog, "r") do fits
        t = Table(fits[2])
        idx = findfirst(==(gaia_id), t.gaia_source_id)
        return NamedTuple(t[idx])
    end

    # Convert measurement epochs to MJD.
    # Careful: these are Julian years, not decimal years (T. Brant., private communications)
    J2000_mjd = 51544.5 # year J2000 in MJD
    catalog = (;
        catalog...,
        epoch_ra_hip_mjd=(catalog.epoch_ra_hip - 2000) * julian_year + J2000_mjd,
        epoch_dec_hip_mjd=(catalog.epoch_dec_hip - 2000) * julian_year + J2000_mjd,
        epoch_ra_dr2_mjd=(catalog.epoch_ra_dr2 - 2000) * julian_year + J2000_mjd,
        epoch_dec_dr2_mjd=(catalog.epoch_dec_dr2 - 2000) * julian_year + J2000_mjd,
        epoch_ra_dr3_mjd=(catalog.epoch_ra_dr3 - 2000) * julian_year + J2000_mjd,
        epoch_dec_dr3_mjd=(catalog.epoch_dec_dr3 - 2000) * julian_year + J2000_mjd,
    )


    @warn "TODO: make sure column makes it into final catalog, loading from Gaia for now"
    dr3 = Octofitter._query_gaia_dr3(;gaia_id)
    catalog = (;catalog..., astrometric_chi2_al_dr3=dr3.astrometric_chi2_al)

    # Load the Hipparcos IAD data for epochs and scan angles
    hip_like = HipparcosIADLikelihood(;
        catalog.hip_id,
        ref_epoch_ra=catalog.epoch_ra_hip_mjd,
        ref_epoch_dec=catalog.epoch_dec_hip_mjd
    )
    A_prepared_5_hip = hip_like.A_prepared_5
    hip_table = hip_like.table

    # Load the Gaia scanlaw etc
    # gaia_like = GaiaCatalogFitLikelihood(; gaia_id_dr3=gaia_id)

    # Besides epoch and catalog, I'm not sure we will really use this data table
    # except maybe for plotting
    # table = Table(;
    #     epoch=[hipparcos_catalog_epoch_mjd, meta_gaia_DR2.ref_epoch_mjd],
    #     catalog=[:hipparcos, :gaia],
    #     ra=[hip_like.hip_sol.radeg, gaia_like.gaia_sol.ra],
    #     dec=[hip_like.hip_sol.dedeg, gaia_like.gaia_sol.dec],
    #     plx=[hip_like.hip_sol.plx, gaia_like.gaia_sol.parallax],
    #     pmra=[hip_like.hip_sol.pm_ra, gaia_like.gaia_sol.pmra],
    #     pmdec=[hip_like.hip_sol.pm_de, gaia_like.gaia_sol.pmdec],
    # )

    # Precompute MvNormal distributions for correlation between ra and dec
    # Hipparcos epoch
    c = catalog.pmra_pmdec_hip[1] * catalog.pmra_hip_error[1] * catalog.pmdec_hip_error[1]
    dist_hip = MvNormal(
        @SVector([catalog.pmra_hip, catalog.pmdec_hip]),
        @SArray[
            catalog.pmra_hip_error[1]^2 c
            c catalog.pmdec_hip_error[1]^2
        ]
    )
    # Hipparcos - GAIA epoch
    c = catalog.pmra_pmdec_hg[1] * catalog.pmra_hg_error[1] * catalog.pmdec_hg_error[1]
    dist_hg = MvNormal(
        @SVector([catalog.pmra_hg, catalog.pmdec_hg]),
        @SArray [
            catalog.pmra_hg_error[1]^2 c
            c catalog.pmdec_hg_error[1]^2
        ]
    )

    # GAIA DR2 epoch
    c = catalog.pmra_pmdec_dr2[1] * catalog.pmra_dr2_error[1] * catalog.pmdec_dr2_error[1]
    dist_dr2 = MvNormal(
        @SVector([catalog.pmra_dr2, catalog.pmdec_dr2]),
        @SArray [
            catalog.pmra_dr2_error[1]^2 c
            c catalog.pmdec_dr2_error[1]^2
        ]
    )


    # GAIA DR3-DR2 epoch
    c = catalog.pmra_pmdec_dr32[1] * catalog.pmra_dr32_error[1] * catalog.pmdec_dr32_error[1]
    dist_dr32 = MvNormal(
        @SVector([catalog.pmra_dr32, catalog.pmdec_dr32]),
        @SArray [
            catalog.pmra_dr3_error[1]^2 c
            c catalog.pmdec_dr3_error[1]^2
        ]
    )

    # GAIA DR3 epoch
    c = catalog.pmra_pmdec_dr3[1] * catalog.pmra_dr3_error[1] * catalog.pmdec_dr3_error[1]
    dist_dr3 = MvNormal(
        @SVector([catalog.pmra_dr3, catalog.pmdec_dr3]),
        @SArray [
            catalog.pmra_dr3_error[1]^2 c
            c catalog.pmdec_dr3_error[1]^2
        ]
    )

    catalog = (; catalog..., dist_hip, dist_hg, dist_dr2, dist_dr32, dist_dr3)

    if isnothing(fluxratio_var)
        fluxratio_var = :__NA
    end




    if isnothing(scanlaw_table)
        # @warn "No scan law table provided. We will fetch an approximate solution from the GHOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GHOST_forecast(catalog.ra_dr3,catalog.dec_dr3))
        forecast_table.epoch = jd2mjd.(forecast_table.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        forecast_table.scanAngle_rad = forecast_table.scanAngle_rad_
    else
        @info "Scanlaw table from the `scanninglaw` python package was provided, will not use GHOST."
        forecast_table = FlexTable(scanlaw_table)
        forecast_table.epoch = tcb_at_gaia_2mjd.(forecast_table.times)
        forecast_table.scanAngle_rad = deg2rad.(forecast_table.angles)
    end
    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.scanAngle_rad)

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    gaia_table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)



    # Now remove any known gaps -- data sourced from HTOF.py; authors G.M. Brandt et al
    gaps_dr2 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiadr2_08252020.csv"), FlexTable)
    gaps_edr23 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiaedr3_12232020.csv"), FlexTable)
    gaps = Table(
        start_mjd=obmt2mjd.(vcat(gaps_dr2.start,gaps_edr23.start)),
        stop_mjd=obmt2mjd.(vcat(gaps_dr2.end,gaps_edr23.end)),
        note=[gaps_dr2.comment; gaps_edr23.description]
    )
    gaia_table = filter(eachrow(gaia_table)) do row
        row = row[]
        for gap in eachrow(gaps)
            gap = gap[]
            if gap.start_mjd <= row.epoch <= gap.stop_mjd
                @info "Detected known gap in Gaia scans; skipping." window=row.epoch note=gap.note
                return false
            end
        end
        return true
    end
    gaia_table = Table(map(dat->dat[], gaia_table))


    # Determine fraction of epochs in DR2 that overlap with DR3

    istart_dr2 = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaia_table.epoch))
    iend_dr2 = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr2)
        istart_dr2 = 1
    end
    if isnothing(iend_dr2)
        iend_dr2 = length(gaia_table.epoch)
    end

    istart_dr3 = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end

    min_epoch = +Inf
    max_epoch = -Inf
    min_epoch = min(min_epoch,meta_gaia_DR2.start_mjd)
    max_epoch = max(max_epoch,meta_gaia_DR2.stop_mjd)
    min_epoch = min(min_epoch,meta_gaia_DR3.start_mjd)
    max_epoch = max(max_epoch,meta_gaia_DR3.stop_mjd)
    gaia_table = Table(gaia_table[min_epoch .<= gaia_table.epoch .<= max_epoch,:])

    # DR2
    A_prepared_5_dr2 = prepare_A_5param(gaia_table, catalog.epoch_ra_dr2_mjd,  catalog.epoch_dec_dr2_mjd)
    
    # DR3
    A_prepared_5_dr3 = prepare_A_5param(gaia_table, catalog.epoch_ra_dr3_mjd,  catalog.epoch_dec_dr3_mjd)
    
    return GaiaHipparcosUEVAJointLikelihood_v1{
        typeof(hip_table),
        typeof(gaia_table),
        typeof(catalog),
        fluxratio_var,
    }(
        hip_table,
        gaia_table,
        catalog,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        fluxratio_var,
        include_iad,
        ueva_mode,
    )

end





function ln_like(like::GaiaHipparcosUEVAJointLikelihood_v1, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    sim = simulate(like, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

    if isnothing(sim)
        return convert(T,-Inf)
    end

    (; μ_h, μ_hg, μ_dr2, μ_dr32, μ_dr3, UEVA_model, UEVA_unc, μ_1_3) = sim       

    absolute_orbits = false
    for orbit in orbits
        absolute_orbits |= orbit isa AbsoluteVisual
        # TODO: could check in a more user-friendly way
        # that we don't have a mismatch of different orbit types
        # for different planets?
    end

    # If we have propagated the barycentric motion ourselves, we want to remove the
    # nonlinear correction already applied to the HGCA by Tim Brandt (private communications)/

    # Rather than subtract it from the HGCA observed values (which are here, already
    # baked into the pre-computed MvNormal distributions), just add them to the model
    # values
    μ_hg += @SVector [
        like.catalog.nonlinear_dpmra,
        like.catalog.nonlinear_dpmdec,
    ]

    # also have to remove the HGCA's nonlinear_dpmra/dec from the hipparcos epoch
    # Note: factor of two needed since dpmra is defined to the HG epoch, so H epoch
    # is twice as much. (T. Brandt, private communications).
    μ_h += @SVector [
        2like.catalog.nonlinear_dpmra,
        2like.catalog.nonlinear_dpmdec,
    ]

    ll += logpdf(like.catalog.dist_hip, μ_h)
    ll += logpdf(like.catalog.dist_hg, μ_hg)
    ll += logpdf(like.catalog.dist_dr3, μ_dr2)
    ll += logpdf(like.catalog.dist_dr3, μ_dr32)
    ll += logpdf(like.catalog.dist_dr3, μ_dr3)

    # UEVA: EAN/RUWE
    if !isfinite(UEVA_unc) || UEVA_unc <= eps()
        ll = -Inf
    else
        try
            ll += logpdf(Normal(μ_1_3, UEVA_unc), UEVA_model)
        catch
            @error "Invalid" μ_1_3 UEVA_unc UEVA_model
            error("invalid values encountered")
        end
    end


    return ll
end


function simulate(like::GaiaHipparcosUEVAJointLikelihood_v1, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    T = _system_number_type(θ_system)


    # Generate simulated observations from this sample draw
    # (;missed_transits) = θ_system
    (;σ_att, σ_AL, σ_calib, gaia_n_dof, missed_transits) = θ_system
    σ_formal = sqrt(σ_att^2 + σ_AL^2)

    if eltype(missed_transits) <: AbstractFloat
        missed_transits = Int.(missed_transits)
    end
    if length(unique(missed_transits)) < length(missed_transits)
        return nothing
    end
    ii = sort(setdiff(1:length(like.gaia_table.epoch), missed_transits))
    gaia_table = like.gaia_table[ii,:]
    A_prepared_5_dr3 = like.A_prepared_5_dr3[ii,:]
    A_prepared_5_dr2 = like.A_prepared_5_dr2[ii,:]

    # Now we fit a no-planet (zero mass planet) sky path model to this data.
    # These should be fit using the appropriate catalog reference epoch so 
    # that they can be compared correctly.

    absolute_orbits = false
    for orbit in orbits
        absolute_orbits |= orbit isa AbsoluteVisual
        # TODO: could check in a more user-friendly way
        # that we don't have a mismatch of different orbit types
        # for different planets?
    end



    # I guess we add that delta PM to our propagated PM, and compare vs the catalog.

    # Helper functions to either get the static pmra from the orbital elements,
    # or, if using an AbsoluteVisualOrbit, get the propagated pmra at the
    # current epoch accounting for barycentric motion.
    function propagate_astrom(orbit::PlanetOrbits.AbsoluteVisualOrbit, epoch_ra, epoch_dec)
        sol_ra = orbitsolve(orbit, epoch_ra)
        cmp_ra = sol_ra.compensated
        sol_dec = orbitsolve(orbit, epoch_dec)
        cmp_dec = sol_dec.compensated
        # Account for the instantaneous differential light travel time apparent acceleration.
        # Treat as linear for the duration of Gaia or Hipparcos
        t1 = max(epoch_ra, epoch_dec)
        Δt = 100
        t2 = t1 + Δt
        sol = epoch_ra >= epoch_dec ? sol_ra : sol_dec
        sol′ = orbitsolve(orbit,t2)
        # This isn't right! This is double counting the proper motion which already goes into ra/dec
        # Take change in delta_time and multiply it by pmra/pmdec
        diff_lt_app_pmra = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmra2
        diff_lt_app_pmdec = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmdec2
        return cmp_ra.ra2, cmp_dec.dec2, cmp_ra.pmra2+diff_lt_app_pmra, cmp_dec.pmdec2+diff_lt_app_pmdec
    end
    function propagate_astrom(orbit::Any, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end

    @no_escape begin


    ################################
    # Hipparcos
    Δα_mas = @alloc(T, size(like.hip_table,1))
    fill!(Δα_mas, 0)
    Δδ_mas = @alloc(T, size(like.hip_table,1))
    fill!(Δδ_mas, 0)

    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        (;fluxratio) = _getparams(like, θ_planet)
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            like.hip_table, orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions[i_planet],
            orbit_solutions_i_epoch_start[i_planet], T
        )
    end
    out = fit_5param_prepared(like.A_prepared_5_hip, like.hip_table, Δα_mas, Δδ_mas, like.hip_table.res, like.hip_table.sres)
    out = fit_5param_prepared(like.A_prepared_5_hip, like.hip_table, Δα_mas, Δδ_mas)
    if like.include_iad
        out = fit_5param_prepared(like.A_prepared_5_hip, like.hip_table, Δα_mas, Δδ_mas, like.hip_table.res, like.hip_table.sres)
    else
        out = fit_5param_prepared(like.A_prepared_5_hip, like.hip_table, Δα_mas, Δδ_mas)
    end
    Δα_h, Δδ_h, Δpmra_h, Δpmdec_h = out.parameters
    α_h₀, δ_h₀, pmra_h₀, pmdec_h₀ = propagate_astrom(first(orbits), like.catalog.epoch_ra_hip_mjd, like.catalog.epoch_dec_hip_mjd)
    μ_h = @SVector [pmra_h₀ + Δpmra_h, pmdec_h₀ + Δpmdec_h]


    ################################
    # DR2
    istart = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaia_table.epoch))
    iend = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart)
        istart = 1
    end
    if isnothing(iend)
        iend = length(gaia_table.epoch)
    end
    Δα_mas = @alloc(T, iend-istart+1); fill!(Δα_mas, 0)
    Δδ_mas = @alloc(T, iend-istart+1); fill!(Δδ_mas, 0)
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        (;fluxratio) = _getparams(like, θ_planet)
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            gaia_table[istart:iend], orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions[i_planet],
            orbit_solutions_i_epoch_start[i_planet], T
        )
    end
    out = fit_5param_prepared(A_prepared_5_dr2[istart:iend,:], gaia_table[istart:iend], Δα_mas, Δδ_mas)
    # out = fit_4param_prepared(hgca_like.gaialike.A_prepared_4, gaia_table, Δα_mas, Δδ_mas)
    Δα_dr2, Δδ_dr2, Δpmra_dr2, Δpmdec_dr2 = out.parameters
    # Rigorously propagate the linear proper motion component in spherical coordinates
    # Account for within-gaia differential light travel time 
    α_dr2₀, δ_dr2₀, pmra_dr2₀, pmdec_dr2₀ = propagate_astrom(first(orbits), like.catalog.epoch_ra_dr2_mjd, like.catalog.epoch_dec_dr2_mjd)
    μ_dr2 = @SVector [pmra_dr2₀ + Δpmra_dr2, pmdec_dr2₀ + Δpmdec_dr2]

    ################################
    # DR3
    istart = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaia_table.epoch))
    iend = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart)
        istart = 1
    end
    if isnothing(iend)
        iend = length(gaia_table.epoch)
    end
    Δα_mas = @alloc(T, iend-istart+1); fill!(Δα_mas, 0)
    Δδ_mas = @alloc(T, iend-istart+1); fill!(Δδ_mas, 0)
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        (;fluxratio) = _getparams(like, θ_planet)
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            gaia_table[istart:iend], orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions[i_planet],
            orbit_solutions_i_epoch_start[i_planet], T,
        )
    end

    out_dr3 = fit_5param_prepared(A_prepared_5_dr3[istart:iend,:], gaia_table[istart:iend], Δα_mas, Δδ_mas, 0.0, σ_formal; include_chi2=Val(true))
    Δα_dr3, Δδ_dr3, Δpmra_dr3, Δpmdec_dr3 = out_dr3.parameters
    # Rigorously propagate the linear proper motion component in spherical coordinates
    # Account for within-gaia differential light travel time 
    α_dr3₀, δ_dr3₀, pmra_dr3₀, pmdec_dr3₀ = propagate_astrom(first(orbits), like.catalog.epoch_ra_dr3_mjd, like.catalog.epoch_dec_dr3_mjd)
    μ_dr3 = @SVector [pmra_dr3₀ + Δpmra_dr3, pmdec_dr3₀ + Δpmdec_dr3]


    ########################


    # Simple linear approximation: don't deal with curvature & secular acceleration directly
    if absolute_orbits

        # HG
        Δα_hg_prop = (α_dr3₀ - α_h₀)*60*60*1000*cosd((δ_dr3₀ + δ_h₀)/2)
        Δδ_hg_prop = (δ_dr3₀ - δ_h₀)*60*60*1000
        pmra_hg_model = (Δα_dr3 - Δα_h + Δα_hg_prop) / (
            like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
        )*julian_year
        pmdec_hg_model = (Δδ_dr3 - Δδ_h + Δδ_hg_prop) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
        )*julian_year


        # DR3-DR2
        Δα_dr32_prop = (α_dr3₀ - α_dr2₀)*60*60*1000*cosd((δ_dr3₀ + δ_dr2₀)/2)
        Δδ_dr32_prop = (δ_dr3₀ - δ_dr2₀)*60*60*1000
        pmra_dr32_model = (Δα_dr3 - Δα_dr2 + Δα_dr32_prop) / (
            like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd
        )*julian_year
        pmdec_dr32_model = (Δδ_dr3 - Δδ_dr2 + Δδ_dr32_prop) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd
        )*julian_year

    else
        pmra_hg_model = (Δα_dr3 - Δα_h) / (
                like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
        )*julian_year + θ_system.pmra
        pmdec_hg_model = (Δδ_dr3 - Δδ_h) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
        )*julian_year + θ_system.pmdec

        pmra_dr32_model = (Δα_dr3 - Δα_dr2) / (
                like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd
        )*julian_year + θ_system.pmra
        pmdec_dr32_model = (Δδ_dr3 - Δδ_dr2) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd
        )*julian_year + θ_system.pmdec


    end


    μ_hg = @SVector [pmra_hg_model, pmdec_hg_model]
    μ_dr32 = @SVector [pmra_dr32_model, pmdec_dr32_model]



    ##############################
    # DR3 UEVA calculation
    # From Gaia catalog:
    (;
        astrometric_chi2_al_dr3,         # Chi squared of along-scan measurements
        astrometric_n_good_obs_al_dr3,   # Number of good AL observations (N)  
        astrometric_matched_transits_dr3,# Number of field of view transits (N_FoV)
        # phot_g_mean_mag,             # G magnitude
        # bp_rp,                       # BP-RP color
        # ra, dec,                     # Position
        astrometric_excess_noise_dr3,
        ruwe_dr3,
    ) = like.catalog

    # Observed UEVA
    if like.ueva_mode == :EAN
        UEVA = astrometric_excess_noise_dr3^2 + σ_att^2 + σ_AL^2
    elseif like.ueva_mode == :RUWE
        # normalization factor for that G mag & BP-RP
        # eqn. (3) Kiefer et al 2024
        # TODO: should this 5 change to be 6 if more dof?
        if gaia_n_dof == 6
            @warn "RUWE calculation: TODO: should the 5 become a 6 when doing 6-dof fits?" maxlog=10
        end
        u0 = 1/ruwe_dr3*sqrt(astrometric_chi2_al_dr3/(astrometric_n_good_obs_al_dr3-5))
        UEVA = (ruwe_dr3 * u0)^2 * σ_formal^2
    else
        error("Unsupported mode (should be :EAN or :RUWE, was $(like.ueva_mode)")
    end

    N = astrometric_n_good_obs_al_dr3
    N_FoV = astrometric_matched_transits_dr3
    N_AL = N/N_FoV # TODO: why isn't this an integer always?
    
    chi2_astro_scaled = out_dr3.chi_squared_astro * N_AL

    # Expected UEVA for a single star
    μ_UEVA_single = (N_AL/(N_AL*N_FoV - gaia_n_dof)) * 
    ((N_FoV - gaia_n_dof)*σ_calib^2 + N_FoV*σ_AL^2)

    # And its variance
    σ_UEVA_single = sqrt(2*N_AL/(N_AL*N_FoV - gaia_n_dof)^2 * 
    (
        N_AL*(N_FoV - gaia_n_dof)*σ_calib^4 + 
        N_FoV*σ_AL^4 + 2*N_FoV*σ_AL^2*σ_calib^2
    ))

    UEVA_model_1 = (chi2_astro_scaled * σ_formal^2) / (N_AL * (N_FoV - gaia_n_dof))
    
    # Compare to expected single-star distribution
    μ_1_3 = UEVA^(1/3) 
    UEVA_unc = σ_UEVA_single * UEVA^(-2/3) /3  # divide by 3 due to cube root transformation
    UEVA_model = cbrt(UEVA_model_1 + μ_UEVA_single)

    end

    
    return (;

        # UEVA: EAN/RUWE
        UEVA_model,
        UEVA_unc,
        μ_1_3,

        # Packaged up nicely
        μ_h,
        μ_hg,
        μ_dr2,
        μ_dr32,
        μ_dr3,

        # Individual
        pmra_hip_model=μ_h[1],
        pmdec_hip_model=μ_h[2],
        pmra_hg_model=μ_hg[1],
        pmdec_hg_model=μ_hg[2],
        pmra_dr2_model=μ_dr2[1],
        pmdec_dr2_model=μ_dr2[2],
        pmra_dr32_model=μ_dr32[1],
        pmdec_dr32_model=μ_dr32[2],
        pmra_dr3_model=μ_dr3[1],
        pmdec_dr3_model=μ_dr3[2],

        # modelled_gaia_parameters_dr3=modelled_gaia_parameters_dr3,
        # modelled_gaia_parameters_dr2=modelled_gaia_parameters_dr2,#.-correction,
        
        # pmra_dr3_model = modelled_gaia_parameters_dr3[3],
        # pmdec_dr3_model = modelled_gaia_parameters_dr3[4],

        # # pmra_dr2_model = modelled_gaia_parameters_dr2[3] - θ_system.dr2_systematic_Δμ_ra,
        # # pmdec_dr2_model = modelled_gaia_parameters_dr2[4] - θ_system.dr2_systematic_Δμ_dec,
        # # pmra_dr32_model=((
        # #     modelled_gaia_parameters_dr3[1]-modelled_gaia_parameters_dr2[1]
        # # )*60*60*1000*cosd((like.dr3.ra-like.dr2.ra)/2) + θ_system.dr2_systematic_Δra)/Δt,
        # # pmdec_dr32_model=((
        # #     modelled_gaia_parameters_dr3[2]-modelled_gaia_parameters_dr2[2]
        # # )*60*60*1000 + θ_system.dr2_systematic_Δdec)/Δt,

        # pmra_dr2_model = modelled_gaia_parameters_dr2[3],
        # pmdec_dr2_model = modelled_gaia_parameters_dr2[4],
        # pmra_dr32_model=((
        #     modelled_gaia_parameters_dr3[1]-modelled_gaia_parameters_dr2[1]
        # )*60*60*1000*cosd((like.dr3.ra-like.dr2.ra)/2))/Δt,
        # pmdec_dr32_model=((
        #     modelled_gaia_parameters_dr3[2]-modelled_gaia_parameters_dr2[2]
        # )*60*60*1000)/Δt,
        # correction,

        # residual_pmra_dr2_model = μ_dr2_dr3[3] - μ_dr2_dr3_modelled[3],
        # residual_pmdec_dr2_model =  μ_dr2_dr3[4] - μ_dr2_dr3_modelled[4],
        # residual_pmra_dr3_model = μ_dr2_dr3[7] - μ_dr2_dr3_modelled[7],
        # residual_pmdec_dr3_model =  μ_dr2_dr3[8] - μ_dr2_dr3_modelled[8],
        # residual_pmra_dr32_model= (μ_dr2_dr3[1] - μ_dr2_dr3_modelled[1])*60*60*1000*cosd((like.dr3.ra-like.dr2.ra)/2)/Δt,
        # residual_pmdec_dr32_model= (μ_dr2_dr3[2] - μ_dr2_dr3_modelled[2])*60*60*1000/Δt,

    )
end