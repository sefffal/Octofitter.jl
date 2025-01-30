
# https://gea.esac.esa.int/archive/documentation/GDR3/Introduction/chap_cu0int/cu0int_sec_release_framework/cu0int_ssec_time_coverage.html
# The approximate relation – with an accuracy of order 1 second – between OBMT (in revolutions) and TCB (in Julian years) at the position of Gaia is:
# # TCB≃J2015.0+(OBMT−1717.6256rev)/(1461rev yr−1).
# # This relation is valid for OBMT>500 rev (3 March 2014).
# # TCB≃J2015.0+(OBMT−1717.6256rev)/(1461rev yr−1).
# function gaia_obmt_rev_2mjd(obmt::Number)
#     if obmt <= 500
#         @warn "Gaia OBMT to TCB relation not valid for OBMT<500 rev"
#     end
#     return 57023.25 + (obmt - 1717.6265)/1461*julian_year
# end
# function gaia_obmt_days_2mjd(obmt::Number)
#     if obmt <= 500
#         @warn "Gaia OBMT to TCB relation not valid for OBMT<500 rev"
#     end
#     return jd2mjd(2457023.75 + obmt - 1717.6256/4)
# end

struct GaiaDivergenceLikelihood_3{TCat3,TCat2,TCat1,TTable,TDist3,TDist2,TDist1} <: AbstractLikelihood
    # Source ID from each given catalog, if available
    source_id_dr3::Union{Nothing,Int}
    source_id_dr2::Union{Nothing,Int}
    source_id_dr1::Union{Nothing,Int}
    # Catalog values of key parameters, if available
    dr3::TCat3
    dr2::TCat2
    dr1::TCat1
    dist_dr3::TDist3
    dist_dr2::TDist2
    dist_dr1::TDist1
    # GHOST predicted observations
    table::TTable
end
function GaiaDivergenceLikelihood_3(;
    source_id_dr1=nothing,
    source_id_dr2=nothing,
    source_id_dr3=nothing,
    scanlaw_table=nothing,

    model_parallax=true
)
    if all(isnothing, (source_id_dr1,source_id_dr2,source_id_dr3))
        throw(ArgumentError("Please provide at least one of `source_id_dr1`, `source_id_dr2`, or `source_id_dr3`"))
    end

    min_epoch = +Inf
    max_epoch = -Inf
    if isnothing(source_id_dr1)
        dr1 = nothing
        dist_dr1 = nothing
    else
        min_epoch = min(min_epoch,meta_gaia_DR1.start_mjd)
        max_epoch = max(max_epoch,meta_gaia_DR1.stop_mjd)
        # ...
        dr1 = Octofitter._query_gaia_dr1(;gaia_id=source_id_dr1)
        ra_deg = dr1.ra
        dec_deg = dr1.dec
        μ = [
            dr1.parallax,
            dr1.ra,# deg
            dr1.dec,# deg
            dr1.pmra, 
            dr1.pmdec,
        ]
        σ = [
            dr1.parallax_error,
            dr1.ra_error/ 60/60/1000,
            dr1.dec_error / 60/60/1000,
            dr1.pmra_error ,
            dr1.pmdec_error,
        ]
        C = [
            # plx                   ra                      dec                     pmra                    pmdec
            1                       dr1.ra_parallax_corr    dr1.dec_parallax_corr   dr1.parallax_pmra_corr  dr1.parallax_pmdec_corr
            dr1.ra_parallax_corr    1                       dr1.ra_dec_corr         dr1.ra_pmra_corr        dr1.ra_pmdec_corr
            dr1.dec_parallax_corr   dr1.ra_dec_corr         1                       dr1.dec_pmra_corr       dr1.dec_pmdec_corr
            dr1.parallax_pmra_corr  dr1.ra_pmra_corr        dr1.dec_pmra_corr       1                       dr1.pmra_pmdec_corr
            dr1.parallax_pmdec_corr dr1.ra_pmdec_corr       dr1.dec_pmdec_corr      dr1.pmra_pmdec_corr     1
        ]
        Σ = Diagonal(σ) * C * Diagonal(σ)
        if !model_parallax
            μ = μ[2:end]
            Σ = Σ[2:end,2:end]
        end
        dist_dr1 = MvNormal(μ,Hermitian(Σ))
    end

    if isnothing(source_id_dr2)
        dr2 = nothing
        dist_dr2 = nothing
    else
        min_epoch = min(min_epoch,meta_gaia_DR2.start_mjd)
        max_epoch = max(max_epoch,meta_gaia_DR2.stop_mjd)
        # Query Gaia archive for DR3 solution
        dr2 = Octofitter._query_gaia_dr2(;gaia_id=source_id_dr2)
        ra_deg = dr2.ra
        dec_deg = dr2.dec
        μ = [
            dr2.parallax,
            dr2.ra,# deg
            dr2.dec,# deg
            dr2.pmra, 
            dr2.pmdec,
        ]
        σ = [
            dr2.parallax_error,
            dr2.ra_error/ 60/60/1000,
            dr2.dec_error / 60/60/1000,
            dr2.pmra_error ,
            dr2.pmdec_error,
        ]
        C = [
            # plx                   ra                      dec                     pmra                    pmdec
            1                       dr2.ra_parallax_corr    dr2.dec_parallax_corr   dr2.parallax_pmra_corr  dr2.parallax_pmdec_corr
            dr2.ra_parallax_corr    1                       dr2.ra_dec_corr         dr2.ra_pmra_corr        dr2.ra_pmdec_corr
            dr2.dec_parallax_corr   dr2.ra_dec_corr         1                       dr2.dec_pmra_corr       dr2.dec_pmdec_corr
            dr2.parallax_pmra_corr  dr2.ra_pmra_corr        dr2.dec_pmra_corr       1                       dr2.pmra_pmdec_corr
            dr2.parallax_pmdec_corr dr2.ra_pmdec_corr       dr2.dec_pmdec_corr      dr2.pmra_pmdec_corr     1
        ]
        Σ = Diagonal(σ) * C * Diagonal(σ)
        if !model_parallax
            μ = μ[2:end]
            Σ = Σ[2:end,2:end]
        end
        dist_dr2 = MvNormal(μ,Hermitian(Σ))
    end

    if isnothing(source_id_dr3)
        dr3 = nothing
        dist_dr3 = nothing
    else
        min_epoch = min(min_epoch,meta_gaia_DR3.start_mjd)
        max_epoch = max(max_epoch,meta_gaia_DR3.stop_mjd)
        # Query Gaia archive for DR3 solution
        dr3 = Octofitter._query_gaia_dr3(;gaia_id=source_id_dr3)
        ra_deg = dr3.ra
        dec_deg = dr3.dec
        μ = [
            dr3.parallax,
            dr3.ra,# deg
            dr3.dec,# deg
            dr3.pmra, 
            dr3.pmdec,
        ]
        σ = [
            dr3.parallax_error,
            dr3.ra_error/ 60/60/1000,
            dr3.dec_error / 60/60/1000,
            dr3.pmra_error ,
            dr3.pmdec_error,
        ]
        C = [
            # plx                   ra                      dec                     pmra                    pmdec
            1                       dr3.ra_parallax_corr    dr3.dec_parallax_corr   dr3.parallax_pmra_corr  dr3.parallax_pmdec_corr
            dr3.ra_parallax_corr    1                       dr3.ra_dec_corr         dr3.ra_pmra_corr        dr3.ra_pmdec_corr
            dr3.dec_parallax_corr   dr3.ra_dec_corr         1                       dr3.dec_pmra_corr       dr3.dec_pmdec_corr
            dr3.parallax_pmra_corr  dr3.ra_pmra_corr        dr3.dec_pmra_corr       1                       dr3.pmra_pmdec_corr
            dr3.parallax_pmdec_corr dr3.ra_pmdec_corr       dr3.dec_pmdec_corr      dr3.pmra_pmdec_corr     1
        ]
        Σ = Diagonal(σ) * C * Diagonal(σ)
        if !model_parallax
            μ = μ[2:end]
            Σ = Σ[2:end,2:end]
        end
        dist_dr3 = MvNormal(μ,Hermitian(Σ))
    end

    if isnothing(scanlaw_table)
        @warn "No scan law table provided. We will fetch an approximate solution from the GHOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GHOST_forecast(ra_deg,dec_deg))
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
    table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)
    # table.derived_along_scan_resid = zeros(length(table.epoch))
    # table.derived_initial_α = zeros(length(table.epoch))
    # table.derived_initial_δ = zeros(length(table.epoch))
    # table.derived_αₘ = fill(SVector{2,Float64}(0,0), length(table.epoch))
    # table.derived_δₘ = fill(SVector{2,Float64}(0,0), length(table.epoch))
    table = Table(table[min_epoch .<= table.epoch .<= max_epoch,:])

    return GaiaDivergenceLikelihood_3(
        source_id_dr3,
        source_id_dr2,
        source_id_dr1,
        dr3,
        dr2,
        dr1,
        dist_dr3,
        dist_dr2,
        dist_dr1,
        table,
    )
end

function getrv(dr)
    if !hasproperty(dr, :radial_velocity)
        0.
    elseif isnothing(dr.radial_velocity)
        0.
    else
        dr.radial_velocity
    end
end







function ln_like(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 
    ll, _ = simulate(gaialike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    return ll
end


function simulate(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    T = _system_number_type(θ_system)

    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol

    sample_orbit = only(orbits)

    # Re-create orbit object with offsets and frame rotation applied
    # TODO: hacky
    # combined = (;
    #     θ_system...,
    #     θ_system.planets.b...,
    #     ra = θ_system.ra+θ_system.Δra/60/60/1000,
    #     dec = θ_system.dec+θ_system.Δdec/60/60/1000,
    #     pmra = θ_system.pmra+θ_system.Δμ_ra,
    #     pmdec = θ_system.pmra+θ_system.Δμ_dec,
    # )
    # offset_orbit = orbit(;combined...)
    offset_orbit = sample_orbit

    # Generate simulated observations from this sample draw
    # These are generated using whatever reference epoch the user has specified that corresponds
    # to their ra, dec, plx, etc variables
    (;α_model,δ_model,αₘ,δₘ) = _simulate_skypath_observations(gaialike, offset_orbit::AbsoluteVisual, planet_mass_msol, orbit_solutions, orbit_solutions_i_epoch_start, T)


    ll = zero(T)

    # Now we fit a no-planet (zero mass planet) sky path model to this data.
    # These should be fit using the appropriate catalog reference epoch so 
    # that they can be compared correctly.
    modelled_gaia_parameters_dr3 = modelled_gaia_parameters_dr2 = nothing
    ## Model 1: DR3
    if !isnothing(gaialike.dr3)
        istart = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaialike.table.epoch))
        iend = findlast(<=(meta_gaia_DR3.start_mjd), vec(gaialike.table.epoch))
        if isnothing(istart)
            istart = 1
        end
        if isnothing(iend)
            iend = length(gaialike.table.epoch)
        end
        data = (;
            αₘ=αₘ[istart:iend],
            δₘ=δₘ[istart:iend],
            table=gaialike.table[istart:iend],
            ref_epoch=meta_gaia_DR3.ref_epoch_mjd,
            along_scan_uncertainty_mas=θ_system.excess_noise_dr3,
            catalog_rv_m_s=getrv(gaialike.dr3)*1e3,
            full3D=false,
        )
        modelled_gaia_parameters_dr3, loglikedr3 = _find_maxlike_gaia_model_soln(data)

        # # OPTION 3
        # # Basic constants
        # σ_att = 0.076  # mas, default attitude noise
        # u₀ = compute_u0(phot_g_mean_mag, bp_rp)  # need to implement this lookup
        # σ_AL = compute_sigma_AL(phot_g_mean_mag, bp_rp)  # need to implement this lookup
        # RUWE = (1/u₀) × √(χ²_astro/(N-5))

        
        # OPTION2
        if size(gaialike.dist_dr3) == (5,)
            ll += logpdf(gaialike.dist_dr3, modelled_gaia_parameters_dr3)
        else
            # Asked to ignore parallax
            ll += logpdf(gaialike.dist_dr3, view(modelled_gaia_parameters_dr3, 2:5))
        end

        # OPTION 1
        # Compute the posterior distribution Gaia would have reported from this data
    #     P, sol_u = _posterior_dist_from_optimized_gaia_model(guess, data)
    #     if isnothing(P)
    #         # Optimization can in some cases fail or fail to report covariances
    #         # This should be improved to happen less often.
    #         return -Inf
    #     end
    #     Q = gaialike.dist_dr3
    #     kl = Distributions.kldivergence(P,Q)
    #     ll += -kl




    end


    ## Model 2: DR2
    if !isnothing(gaialike.dr2)
        istart = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaialike.table.epoch))
        iend = findlast(<=(meta_gaia_DR2.start_mjd), vec(gaialike.table.epoch))
        if isnothing(istart)
            istart = 1
        end
        if isnothing(iend)
            iend = length(gaialike.table.epoch)
        end
        data = (;
            αₘ=αₘ[istart:iend],
            δₘ=δₘ[istart:iend],
            table=gaialike.table[istart:iend],
            ref_epoch=meta_gaia_DR2.ref_epoch_mjd,
            along_scan_uncertainty_mas=θ_system.excess_noise_dr2,
            catalog_rv_m_s=getrv(gaialike.dr2)*1e3,
            full3D=false,
        )
        modelled_gaia_parameters_dr2, loglikedr2 = _find_maxlike_gaia_model_soln(data)
        if size(gaialike.dist_dr2) == 5
            ll += logpdf(gaialike.dist_dr2, modelled_gaia_parameters_dr2)
        else
            # Ignore parallax
            ll += logpdf(gaialike.dist_dr2, view(modelled_gaia_parameters_dr2, 2:5))
        end

        # Compute the posterior distribution Gaia would have reported from this data
        # P, sol_u = _posterior_dist_from_optimized_gaia_model(guess, data)
        # if isnothing(P)
        #     # Optimization can in some cases fail or fail to report covariances
        #     # This should be improved to happen less often.
        #     return -Inf
        # end
        # Q = gaialike.dist_dr2
        # kl = Distributions.kldivergence(P,Q)
        # ll += -kl

    end

    return ll, (modelled_gaia_parameters_dr3, modelled_gaia_parameters_dr2)
end



"""
Internal, nested likelihood model. 
This is our best reconstruction of the likelihood function the Gaia team would have tried
to fit (for the 5-parameter solution). Notably it includes no effects of secondary perturbers.
It does include the effect of perspective acceleration due to radial velocity, which should 
be set to the catalog values (not eg the user's model inputs).

We use this likelihood to reproduce their results.
"""
function _gaia_skypath_loglikelihood(args,(;
    αₘ,δₘ,
    table,
    ref_epoch,
    along_scan_uncertainty_mas,
    full3D,
    catalog_rv_m_s,
))
    (
        plx,
        ra,
        dec,
        pmra,
        pmdec,
    ) = args

    ll = zero(eltype(args))


    plx0 = plx
    dist0 = 1000/plx0
    δdist_pc_δt_sec = catalog_rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
    for i in eachindex(table.epoch)
        Δdist_pc = δdist_pc_δt_sec * (table.epoch[i] - ref_epoch)
        dist1 = dist0 + Δdist_pc
        plx_at_time = 1000 / dist1
        # Calculate sky coordinate Delta at each scan epoch from the catalog position
        # using the eath's motion ("ephemeris"), x,y,z in AU.
        # TODO: should I rather be using the catalog Ra and Dec in here? Unclear
        ra_apparent_Δmas = plx_at_time * (
            table.x[i] * sind(ra) -
            table.y[i] * cosd(ra)
        ) + (table.epoch[i] - ref_epoch)/julian_year * pmra
        dec_apparent_Δmas = plx_at_time * (
            table.x[i] * cosd(ra) * sind(dec) +
            table.y[i] * sind(ra) * sind(dec) -
            table.z[i] * cosd(dec)
        ) + (table.epoch[i] - ref_epoch)/julian_year * pmdec
        ra_apparent_deg = ra + ra_apparent_Δmas/(60*60*1000)/cosd(dec)
        dec_apparent_deg = dec + dec_apparent_Δmas/(60*60*1000)
        point = @SVector [ra_apparent_deg, dec_apparent_deg] 
        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
        line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
        # TODO: distance point to line should account for spherical coordinates!
        resid = distance_point_to_line(point, line_point_1, line_point_2) # mas
        # uncertainty_mas = sqrt(along_scan_uncertainty_mas[i]^2 + astrometric_excess_noise^2)
        # ll += logpdf(Normal(0,uncertainty_mas/60/60/1000), resid)
        ll += logpdf(Normal(0,along_scan_uncertainty_mas/60/60/1000), resid)
    end
    return ll
end

# """
# Fit the Gaia 5-parameter model to some data, and report the best fit means and covariances
# as a multivariate normal distribution.
# Uncertainties are calculated using the inverse hessian of the likelihood at the best-fit location.
# """
# function _find_maxlike_gaia_model_soln(guess, data)
#     _cost_func_1 = OptimizationFunction(
#         (args...)->-_gaia_skypath_loglikelihood(args...),
#         # Optimization.OptimizationBase.DifferentiationInterface.SecondOrder(AutoEnzyme(), AutoForwardDiff())
#         Optimization.OptimizationBase.DifferentiationInterface.SecondOrder(AutoForwardDiff(), AutoForwardDiff())
#     )
#     # _optim_alg_1 = LBFGS(
#     #     m=4,
#     #     # linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
#     #     # alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
#     # )
#     # # _optim_alg_1 = NelderMead()
#     _optim_alg_1 = NewtonTrustRegion()
#     # OptimizationProblem(_cost_func_1, guess, data)
#     prob = OptimizationProblem{true}(_cost_func_1, guess, data)
#     sol = solve(prob, _optim_alg_1,
#         # NewtonTrustRegion(), # Is a good choice too, but about 2x as slow.
#         # x_abstol=NaN,
#         # x_reltol=NaN,
#         # f_abstol=NaN,
#         # f_reltol=NaN,
#         g_abstol= 1e-2,
#         g_reltol= 1e-2,
#         # g_tol=1e-15, f_tol=NaN,  x_tol=NaN, x_reltol=NaN, 
#         allow_f_increases=true,
#         iterations=10000,
#     )
#     # @show sol.u
#     out = sol.u::Vector{eltype(guess)}
#     return out
# end
# function _find_maxlike_gaia_model_soln(data)
#     # Unpack data
#     (; αₘ, δₘ, table, ref_epoch, along_scan_uncertainty_mas, catalog_rv_m_s, plx₀) = data
#     n_obs = length(table.epoch)
    
    
#     # Estimate initial position from mean of scan endpoints (in degrees)
#     α₀ = mean(mean([p[1] for p in αₘ]))
#     δ₀ = mean(mean([p[2] for p in δₘ]))
    
#     # Initialize design matrix A and observation vector b
#     A = zeros(n_obs, 5)
#     b = zeros(n_obs)
#     W = Diagonal(fill(1/along_scan_uncertainty_mas^2, n_obs))
    
#     # Constants for RV correction
#     δdist_pc_δt_sec = catalog_rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
    
#     for i in 1:n_obs
#         # Time difference in years
#         Δt = (table.epoch[i] - ref_epoch)/365.25
        
#         # RV correction to parallax
#         Δdist_pc = δdist_pc_δt_sec * (table.epoch[i] - ref_epoch)
#         plx_factor = 1.0 / (1.0 + Δdist_pc * plx₀/1000.0)
        
#         # Design matrix components
#         # Column 1: Parallax effect
#         # The parallax factor needs to be computed based on the solar system barycenter's position
#         # relative to the observation direction at each epoch
#         π_factor = plx_factor * (
#             table.cosϕ[i] * table.parallax_factors_ra[i] + 
#             table.sinϕ[i] * table.parallax_factors_dec[i]
#         )
#         A[i, 1] = π_factor
        
#         # Columns 2-3: Position (RA, Dec)
#         # Project the position offset onto the scanning direction
#         A[i, 2] = table.cosϕ[i]
#         A[i, 3] = table.sinϕ[i]
        
#         # Columns 4-5: Proper motions (pmRA, pmDec)
#         # Time evolution of position projected onto scanning direction
#         A[i, 4] = table.cosϕ[i] * Δt
#         A[i, 5] = table.sinϕ[i] * Δt
        
#         # Observation vector
#         # This represents the along-scan measurement
#         b[i] = table.cosϕ[i] * (αₘ[i][1] - α₀) + table.sinϕ[i] * (δₘ[i][1] - δ₀)
#     end
    
#     # Solve the weighted least squares problem
#     # (AᵀWA)x = AᵀWb
#     x = (A'*W*A) \ (A'*W*b)
    
#     five_param_sol = @SVector [
#         x[1],           # parallax (mas)
#         α₀ + x[2],     # ra (deg)
#         δ₀ + x[3],     # dec (deg)
#         x[4],          # pmra (mas/yr)
#         x[5]           # pmdec (mas/yr)
#     ]

#     # Compute log-likelihood
#     residuals = A*x - b
#     loglike = -0.5 * sum(residuals.^2 .* diag(W)) - 
#               0.5 * n_obs * log(2π*along_scan_uncertainty_mas^2)

#     return five_param_sol, loglike
# end

function _find_maxlike_gaia_model_soln(data)
    # Unpack data
    (; αₘ, δₘ, table, ref_epoch, along_scan_uncertainty_mas, catalog_rv_m_s, plx₀, ra_0, dec_0) = data
    n_obs = length(table.epoch)
    
    # Constants for unit conversion
    deg2mas = 3600.0 * 1000.0  # degrees to milliarcseconds
    
    # Initialize design matrix A and observation vector b
    A = zeros(n_obs, 5)
    b = zeros(n_obs)
    W = Diagonal(fill(1/along_scan_uncertainty_mas^2, n_obs))
    
    # Constants for RV correction
    δdist_pc_δt_sec = catalog_rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
    
    for i in 1:n_obs
        # Get scan angle components
        cos_ϕ = table.cosϕ[i]
        sin_ϕ = table.sinϕ[i]
        
        # Time difference in years from reference epoch
        Δt = (table.epoch[i] - ref_epoch)/365.25
        
        # RV correction to parallax
        Δdist_pc = δdist_pc_δt_sec * (table.epoch[i] - ref_epoch)
        plx_factor = 1.0 / (1.0 + Δdist_pc * plx₀/1000.0)
        
        # Compute offsets from reference position in mas
        Δra_mas = (αₘ[i][1] - ra_0) * deg2mas * cosd(dec_0)
        Δdec_mas = (δₘ[i][1] - dec_0) * deg2mas
        
        # Project the observed offset onto the scanning direction
        b[i] = Δra_mas * cos_ϕ + Δdec_mas * sin_ϕ
        
        # Parallax factors from Gaia's barycentric position
        # gaia_pos = @SVector [table.x[i], table.y[i], table.z[i]]  # in AU
        
        # Compute parallax effect projection
        π_ra = -table.parallaxFactorAlongScan[i] 
        π_dec = table.parallaxFactorAcrossScan[i]
        # π_factor = plx_factor * (π_ra * cos_ϕ + π_dec * sin_ϕ)
        # π_factor = plx_factor * table.parallaxFactorAlongScan[i]

        π_factor = -plx_factor * table.parallaxFactorAlongScan[i]
        
        # Design matrix - following Orbitize! approach
        A[i, 1] = π_factor                     # parallax term
        A[i, 2] = cos_ϕ                        # RA offset term
        A[i, 3] = sin_ϕ                        # Dec offset term
        # A[i, 4] = cos_ϕ * Δt                   # pmRA term
        A[i, 4] = cos_ϕ * Δt * cos(deg2rad(dec_0))  # pmRA term with cos(dec)
        A[i, 5] = sin_ϕ * Δt                   # pmDec term


        # For a single observation i
Δra_mas = (αₘ[i][1] - ra_0) * deg2mas * cos(deg2rad(dec_0))
Δdec_mas = (δₘ[i][1] - dec_0) * deg2mas
along_scan = Δra_mas * cos_ϕ + Δdec_mas * sin_ϕ
println("RA offset (mas): ", Δra_mas)
println("Dec offset (mas): ", Δdec_mas)
println("Along-scan measurement (mas): ", along_scan)
    end
    
    # Solve the weighted least squares problem
    x = (A'*W*A) \ (A'*W*b)
    
    # Convert solution back to catalog units
    five_param_sol = @SVector [
        x[1],                          # parallax (mas)
        ra_0 + x[2]/(deg2mas*cos(deg2rad(dec_0))),  # ra (deg)
        dec_0 + x[3]/deg2mas,         # dec (deg)
        x[4],                         # pmra (mas/yr)
        x[5]                          # pmdec (mas/yr)
    ]

    # Compute log-likelihood
    residuals = A*x - b
    loglike = -0.5 * sum(residuals.^2 .* diag(W)) - 
              0.5 * n_obs * log(2π*along_scan_uncertainty_mas^2)

    return five_param_sol, loglike
end

# """
# Fit the Gaia 5-parameter model to some data, and report the best fit means and covariances
# as a multivariate normal distribution.
# Uncertainties are calculated using the inverse hessian of the likelihood at the best-fit location.
# """
# function _posterior_dist_from_optimized_gaia_model(guess, data)
#     func = OptimizationFunction((args...)->-_gaia_skypath_loglikelihood(args...), AutoForwardDiff())
#     prob = OptimizationProblem(func, guess, data)
#     sol = solve(prob,
#         LBFGS(
#             m=4,
#             # linesearch=OptimizationOptimJL.Optim.LineSearches.HagerZhang(linesearchmax=50),
#             linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
#             alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
#         ),
#         # NewtonTrustRegion(), # Is a good choice too, but about 2x as slow.
#         # x_abstol=NaN,
#         # x_reltol=NaN,
#         # f_abstol=NaN,
#         # f_reltol=NaN,
#         g_abstol= 1e-5,
#         g_reltol= 1e-5,
#         # g_tol=1e-15, f_tol=NaN,  x_tol=NaN, x_reltol=NaN, 
#         allow_f_increases=true,
#         iterations=10000
#     )
#     if sol.retcode != Optimization.ReturnCode.Success
#         return nothing
#     end
#     # Compute uncertainties and covariances using inverse Hessian
#     H = ForwardDiff.hessian(u->_gaia_skypath_loglikelihood(u,data), sol.u)
#     local gaia_post_covariance
#     try
#         gaia_post_covariance = -inv(H)
#     catch 
#         return nothing
#     end
#     # We use the inverse Hessian from the gradient descent as the covariance matrix.
#     # In case of numerical error (non-hermitian matrix) we regularize it by adding 
#     # a value along the diagonals. 
#     # P = MvNormal(sol.u, Hermitian(gaia_post_covariance))
#     P = nothing
#     for e in (0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6)#, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3, 1e4)
#         try
#             P = MvNormal(sol.u, Hermitian(gaia_post_covariance))
#             return P
#         catch
#             for i in 1:length(sol.u)
#                 gaia_post_covariance[i] += e
#             end
#         end
#     end
#     return nothing
# end


# # Given catalog values, simulate
# function _simulate_gaia_5_param_model(table, plx, ra, dec, pmra, pmdec, rv_m_s, ref_epoch)
    
#     ra_apparent_deg = zeros(length(table.epoch))
#     dec_apparent_deg = zeros(length(table.epoch))

#     plx0 = plx
#     dist0 = 1000/plx0
#     δdist_pc_δt_sec = rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
#     for i in eachindex(table.epoch)
#         Δdist_pc = δdist_pc_δt_sec * (table.epoch[i] - ref_epoch)
#         dist1 = dist0 + Δdist_pc
#         plx_at_time = 1000 / dist1
#         # Calculate sky coordinate Delta at each scan epoch from the catalog position
#         # using the eath's motion ("ephemeris"), x,y,z in AU.
#         # TODO: should I rather be using the catalog Ra and Dec in here? Unclear
#         ra_apparent_Δmas = plx_at_time * (
#             table.x[i] * sind(ra) -
#             table.y[i] * cosd(ra)
#         ) + (table.epoch[i] - ref_epoch)/julian_year * pmra
#         dec_apparent_Δmas = plx_at_time * (
#             table.x[i] * cosd(ra) * sind(dec) +
#             table.y[i] * sind(ra) * sind(dec) -
#             table.z[i] * cosd(dec)
#         ) + (table.epoch[i] - ref_epoch)/julian_year * pmdec
#         ra_apparent_deg[i] = ra + ra_apparent_Δmas/(60*60*1000)/cosd(dec)
#         dec_apparent_deg[i] = dec + dec_apparent_Δmas/(60*60*1000)
#     end

#     return ra_apparent_deg, dec_apparent_deg
# end

# end
