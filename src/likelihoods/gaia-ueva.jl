
"""

To provide the scanlaw, run in Python:
```
import astropy.coordinates, scanninglaw, pandas
o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg')
t = scanninglaw.times.Times(version='dr3_nominal')
t.query(o,return_angles=True)
```
"""
struct GaiaUEVALikelihood{TCat,TTable,TDist} <: AbstractLikelihood
    # Source ID from each given catalog, if available
    gaia_id::Int
    mode::Symbol
    # Catalog values of key parameters, if available
    dr3::TCat
    # Scan law table
    table::TTable
    # Precomputed MVNormal distribution from Gaia catalog
    dist::TDist
end
function GaiaUEVALikelihood(;
    gaia_id,
    scanlaw_table,
    mode::Symbol
)
    
    # DR3
    # Query Gaia archive for DR3 solution
    dr3 = Octofitter._query_gaia_dr3(;gaia_id=gaia_id)
    μ_dr3 = [
        dr3.parallax,
        dr3.ra, # dec
        dr3.dec, # dec
        dr3.pmra, 
        dr3.pmdec,
    ]
    # TODO: /cosd() or *cosd().-- I think it's / cosd()
    σ_dr3 = [
        dr3.parallax_error,
        dr3.ra_error, # mas
        dr3.dec_error/cosd(dr3.dec), # mas
        dr3.pmra_error ,
        dr3.pmdec_error,
    ]
    C_dr3 = [
        # plx                   ra                      dec                     pmra                    pmdec
        1                       dr3.ra_parallax_corr    dr3.dec_parallax_corr   dr3.parallax_pmra_corr  dr3.parallax_pmdec_corr
        dr3.ra_parallax_corr    1                       dr3.ra_dec_corr         dr3.ra_pmra_corr        dr3.ra_pmdec_corr
        dr3.dec_parallax_corr   dr3.ra_dec_corr         1                       dr3.dec_pmra_corr       dr3.dec_pmdec_corr
        dr3.parallax_pmra_corr  dr3.ra_pmra_corr        dr3.dec_pmra_corr       1                       dr3.pmra_pmdec_corr
        dr3.parallax_pmdec_corr dr3.ra_pmdec_corr       dr3.dec_pmdec_corr      dr3.pmra_pmdec_corr     1
    ]
    Σ_dr3 = Diagonal(σ_dr3) * C_dr3 * Diagonal(σ_dr3)
    dist_dr3 = MvNormal(μ_dr3, Hermitian(Σ_dr3))

    # Compute parallax factors
    forecast_table = FlexTable(scanlaw_table)
    forecast_table.epoch = tcb_at_gaia_2mjd.(forecast_table.times)
    forecast_table.scanAngle_rad = deg2rad.(forecast_table.angles)

    # Get the Earth's position at those epochs
    earth_pos_vel = FlexTable(geocentre_position_query.(forecast_table.epoch))

    f = @. earth_pos_vel.x * sind(dr3.ra)-earth_pos_vel.y*cosd(dr3.ra)
    g = @. earth_pos_vel.x * cosd(dr3.ra) * sind(dr3.dec) + 
           earth_pos_vel.y * sind(dr3.ra) * sind(dr3.dec) -
           earth_pos_vel.z * cosd(dr3.dec)
    forecast_table.parallaxFactorAlongScan = @. f*sin(forecast_table.scanAngle_rad) + g*cos(forecast_table.scanAngle_rad)

    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.scanAngle_rad)

    # merge the Gaia scan prediction and geocentre position results into one table
    table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)
    table = Table(table)

    return GaiaUEVALikelihood{typeof(dr3),typeof(table),typeof(dist_dr3)}(gaia_id, mode, dr3, table, dist_dr3)
end

function ln_like(gaialike::GaiaUEVALikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 
    ll, _ = simulate(gaialike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    return ll
end

function simulate(gaialike::GaiaUEVALikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 



    # Modelled UEVA
    T = _system_number_type(θ_system)
    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol


    # (;σ_att, σ_AL, σ_calib, gaia_n_dof) = θ_system
    (;σ_att, σ_AL, σ_calib, gaia_n_dof) = θ_system

    σ_formal = sqrt(σ_att^2 + σ_AL^2)


    # Generate simulated observations from this sample draw
    # These are generated using whatever reference epoch the user has specified that corresponds
    # TODO: in this specific case I should make sure these are generated e.g. using the same RV that gaia used,
    # not whatever the user specified in the model
    sample_orbit = only(orbits)
    # (;α_model,δ_model,αₘ,δₘ) = _simulate_skypath_observations(gaialike, sample_orbit, planet_mass_msol, orbit_solutions, orbit_solutions_i_epoch_start, T)

    # @show α_model δ_model

    #TODO: next problem. Planet mass doesn't seem to be impacting model

    ll = zero(T)

    # TODO: need to deterministically remove the right epochs,
    # or need to marginalize over these
    # ii_skip = rand(1:length(gaialike.table.epoch), 3)
    # ii_skip = 1:8:40
    # ii_skip = 1:18:50
    ii_skip = 1:-1
    ii = sort(setdiff(1:length(gaialike.table.epoch), ii_skip))
    # ii = 1:length(gaialike.table.epoch)

    # out = fit_5param(
    #     # epochs and scan angle
    #     gaialike.table[ii,:],
    #     # simulated observations in ra and dec
    #     α_model[ii,:] .- gaialike.dr3.ra,
    #     δ_model[ii,:] .- gaialike.dr3.dec,
    #     # catalog ra and dec -- we fit in this tangent plane. TODO: might need to iterate for high PM stars.
    #     meta_gaia_DR3.ref_epoch_mjd,
    #     meta_gaia_DR3.ref_epoch_mjd;
    #     # uncertainty to use for chi^2 calculation
    #     σ_formal,
    #     include_chi2=true,
    #     # n_noise_realizations=1000
    # )
    # α₀, δ₀, ϖ, μα, μδ = out.parameters
    # # @show α₀, δ₀, ϖ, μα, μδ
    # modelled_gaia_parameters = [
    #     ϖ,
    #     gaialike.dr3.ra + α₀/60/60/1000/cosd(gaialike.dr3.dec),
    #     gaialike.dr3.dec + δ₀/60/60/1000, # cos delta?
    #     μα,
    #     μδ
    # ]


    T = _system_number_type(θ_system)

    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass * Octofitter.mjup2msol

    sample_orbit = only(orbits)

    ll = zero(T)

    Δα_mas = zeros(T, size(gaialike.table,1))
    Δδ_mas = zeros(T, size(gaialike.table,1))

    # Δα_mas = (α_model .- gaialike.dr3.ra).*60*60*1000*cosd(gaialike.dr3.dec)
    # Δδ_mas = (δ_model .- gaialike.dr3.dec).*60*60*1000
    # α_modelδ_model
    
    # δ_model[ii,:] .- gaialike.dr3.dec

    # @show α_model δ_model
    for (orbit, θ_planet) in zip(orbits, θ_system.planets)
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        # (;fluxratio) = _getparams(gaialike, θ_planet)
        fluxratio =0.0# TODO
        # planet_mass_msol=0.0
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            gaialike.table, orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions,
            orbit_solutions_i_epoch_start, T
        )
    end

    # function propagate_astrom(orbit::PlanetOrbits.AbsoluteVisualOrbit, ref_epoch)
    #     sol = orbitsolve(orbit, ref_epoch)
    #     cmp = sol.compensated
    #     return cmp.ra2, cmp.dec2, cmp.pmra2, cmp.pmdec2
    # end
    # function propagate_astrom(orbit::Any, _, _)
    #     return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    # end


    # ref_epoch = years2mjd(gaialike.gaia_sol.ref_epoch)
    out = fit_5param(
        gaialike.table[ii],
        Δα_mas[ii],
        Δδ_mas[ii],
        meta_gaia_DR3.ref_epoch_mjd,
        meta_gaia_DR3.ref_epoch_mjd;
        # uncertainty to use for chi^2 calculation
        σ_formal,
        include_chi2=true,
    )
    Δα_g, Δδ_g, Δplx, Δpmra_g, Δpmdec_g = out.parameters
    # @show  gaialike.table[ii].epoch Δα_mas[ii] Δδ_mas[ii] Δα_g Δδ_g Δplx Δpmra_g Δpmdec_g

    # α_g₀, δ_g₀, pmra_g₀, pmdec_g₀ = propagate_astrom(first(orbits), ref_epoch)
    # modelled_gaia_parameters_dr3 = @SVector [
    #     gaialike.gaia_sol.parallax + Δplx,
    #     α_g₀ + Δα_g/60/60/1000/cosd(δ_g₀),
    #     δ_g₀ + Δδ_g/60/60/1000,
    #     pmra_g₀ + Δpmra_g,
    #     pmdec_g₀ + Δpmdec_g
    # ]



    # First part of the likelihood -- do we match the catalog plx, ra, dec, pmra, pmdec parameters?
    # ll += logpdf(gaialike.dist, modelled_gaia_parameters)


    # Next part of the likelihood -- do the distribances introduced by the planet lead to fitting noise
    # that matches what Gaia reported, through the UEVA model.

    # From Gaia catalog:
    (;
        astrometric_chi2_al,         # Chi squared of along-scan measurements
        astrometric_n_good_obs_al,   # Number of good AL observations (N)  
        astrometric_matched_transits,# Number of field of view transits (N_FoV)
        phot_g_mean_mag,             # G magnitude
        bp_rp,                       # BP-RP color
        ra, dec,                     # Position
        astrometric_excess_noise,
        ruwe,
    ) = gaialike.dr3



    # σ_formal = hypot(sigma_att, sigma_al)

    # Observed UEVA
    if gaialike.mode == :EAN
        UEVA = astrometric_excess_noise^2 + σ_att^2 + σ_AL^2
    elseif gaialike.mode == :RUWE
        # normalization factor for that G mag & BP-RP
        # eqn. (3) Kiefer et al 2024
        u0 = 1/gaialike.dr3.ruwe*sqrt(gaialike.dr3.astrometric_chi2_al/(gaialike.dr3.astrometric_n_good_obs_al-5))
        UEVA = (ruwe * u0)^2 * σ_formal^2
    else
        error("Unsupported mode (should be :EAN or :RUWE)")
    end


    N = astrometric_n_good_obs_al
    N_FoV = astrometric_matched_transits
    N_AL = N/N_FoV # 
    # N_AL = floor(Int, N/N_FoV)
    # # TODO: problem of dropped FoV transits -- need to marginalize over it.

    chi2_astro_scaled = out.chi_squared_astro * N_AL


    # Expected UEVA for a single star
    μ_UEVA_single = (N_AL/(N_AL*N_FoV - gaia_n_dof)) * 
    ((N_FoV - gaia_n_dof)*σ_calib^2 + N_FoV*σ_AL^2)

    # And its variance
    σ_UEVA_single = sqrt(2*N_AL/(N_AL*N_FoV - gaia_n_dof)^2 * 
    (
        N_AL*(N_FoV - gaia_n_dof)*σ_calib^4 + 
        N_FoV*σ_AL^4 + 2*N_FoV*σ_AL^2*σ_calib^2
    ))

    UEVA_model = (chi2_astro_scaled * σ_formal^2) / (N_AL * (N_FoV - gaia_n_dof))
    
    # Compare to expected single-star distribution
    resid = UEVA^(1/3) - (UEVA_model + μ_UEVA_single)^(1/3)
    UEVA_unc = σ_UEVA_single/3  # divide by 3 due to cube root transformation
    ll += logpdf(Normal(0, UEVA_unc), resid)

    # @show θ_system.planets.b.mass
    # @show UEVA
    # @show UEVA_model

    # @show gaia_n_dof
    # @show σ_formal
    # @show σ_att
    # @show σ_AL
    # @show astrometric_n_good_obs_al
    # @show astrometric_matched_transits
    # @show N_AL
    # @show chi2_astro_scaled
    # @show μ_UEVA_single
    # @show σ_UEVA_single

    # # Photocenter scaling
    # @show l = 0
    # @show q = planet_mass_msol / (θ_system.M - planet_mass_msol)
    # @show photocenter_scalar = abs(q - l)/((1 + l) * (1 + q))

    # # For a dark companion (l=0):
    # @show photocenter_scalar = q/(1 + q)  # where q = m2/m1

    # # Semi-major axis in mas
    # @show  a = 1000 * ((θ_system.M)*θ_system.planets.b.P^2)^(1/3) / (1000/θ_system.plx)  # dist in pc

    # # Photocenter semi-major axis
    # @show  a_phot = a * photocenter_scalar

    return ll, (;
        # modelled_gaia_parameters,
        UEVA,
        UEVA_model
    )
end