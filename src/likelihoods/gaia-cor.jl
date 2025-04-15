

struct GaiaDifferenceLike{TCat3,TCat2,TTable} <: AbstractLikelihood
    # Source ID from each given catalog, if available
    source_id_dr3::Union{Nothing,Int}
    source_id_dr2::Union{Nothing,Int}
    # Catalog values of key parameters, if available
    dr3::TCat3
    dr2::TCat2
    table::TTable
    μ_dr2::Vector{Float64}
    Σ_dr2_dr3::Matrix{Float64}
    # GHOST predicted observations
    A_prepared_5_dr3::Matrix{Float64}
    A_prepared_5_dr2::Matrix{Float64}
end
function GaiaDifferenceLike(;
    source_id_dr2=nothing,
    source_id_dr3=nothing,
    scanlaw_table=nothing,
    dr2_error_inflation_factor=nothing,
    dr3_error_inflation_factor=nothing,
)
    if all(isnothing, (source_id_dr2,source_id_dr3))
        throw(ArgumentError("Please provide at least one of `source_id_dr1`, `source_id_dr2`, or `source_id_dr3`"))
    end

    
    # Query Gaia archive for DR3 solution
    dr3 = Octofitter._query_gaia_dr3(;gaia_id=source_id_dr3)
    dr2 = Octofitter._query_gaia_dr2(;gaia_id=source_id_dr2)
    ra_deg = dr3.ra
    dec_deg = dr3.dec


    if isnothing(dr2_error_inflation_factor)
        hgca_dr2_all = FITS(joinpath(Octofitter.datadep"HGCA_DR2", "HGCA_vDR2_corrected.fits"), "r") do fits
            Table(fits[2])
        end
        idx = findfirst(==(source_id_dr2), hgca_dr2_all.gaia_source_id)
        if isnothing(idx)
            @warn "The DR2 target was not found within the HGCA vDR2, so we have not applied an automatic error inflation to the formal uncertainties. Consider supplying a manual value on the order of `dr2_error_inflation_factor=1.5`"
            dr2_error_inflation_factor = 1
        else
            hgca_dr2 = NamedTuple(hgca_dr2_all[idx])
            dr2_error_inflation_factor = hgca_dr2.pmra_gaia_error / dr2.pmra_error
        end
    end

    if isnothing(dr3_error_inflation_factor)
        hgca_dr3_all = FITS(joinpath(joinpath(datadep"HGCA_eDR3", "HGCA_vEDR3.fits")), "r") do fits
            Table(fits[2])
        end
        idx = findfirst(==(source_id_dr3), hgca_dr3_all.gaia_source_id)
        if isnothing(idx)
            @warn "The DR3 target was not found within the HGCA vEDR3, so we have not applied an automatic error inflation to the formal uncertainties. Consider supplying a manual value on the order of `dr3_error_inflation_factor=1.3`"
            dr3_error_inflation_factor = 1
        else
            hgca_dr3 = NamedTuple(hgca_dr3_all[idx])
            dr3_error_inflation_factor = hgca_dr3.pmra_gaia_error / dr3.pmra_error
        end
    end

    # @show hgca_dr3.gaia_ra dr3.ra
    # @show (hgca_dr3.gaia_ra - hgca_dr2.gaia_ra)/(hgca_dr3.epoch_ra_gaia - hgca_dr2.epoch_ra_gaia)*60*60*1000*cosd(hgca_dr3.gaia_dec)
    # @show hgca_dr3.pmra_hg hgca_dr2.pmra_hg_error
    # @show (hgca_dr3.gaia_dec - hgca_dr2.gaia_dec)/(hgca_dr3.epoch_dec_gaia - hgca_dr2.epoch_dec_gaia)*60*60*1000

    # @show hgca_dr2.epoch_dec_gaia hgca_dr3.epoch_dec_gaia

    @info "Formal proper motion uncertainty inflations:" dr2_error_inflation_factor dr3_error_inflation_factor

    # We might want to inflate the uncertainties by the same factors as the HGCA (DR2 and eDR3 respectively)
    dr2 = (;
        dr2...,
        pmra_error=dr2.pmra_error*dr2_error_inflation_factor,
        pmdec_error=dr2.pmdec_error*dr2_error_inflation_factor
    )
    dr3 = (;
        dr3...,
        pmra_error=dr3.pmra_error*dr3_error_inflation_factor,
        pmdec_error=dr3.pmdec_error*dr3_error_inflation_factor
    )

    if isnothing(scanlaw_table)
        # @warn "No scan law table provided. We will fetch an approximate solution from the GHOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
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

    # Determine fraction of epochs in DR2 that overlap with DR3

    istart_dr2 = findfirst(>=(meta_gaia_DR2.start_mjd), vec(table.epoch))
    iend_dr2 = findlast(<=(meta_gaia_DR2.stop_mjd), vec(table.epoch))
    if isnothing(istart_dr2)
        istart_dr2 = 1
    end
    if isnothing(iend_dr2)
        iend_dr2 = length(table.epoch)
    end

    istart_dr3 = findfirst(>=(meta_gaia_DR3.start_mjd), vec(table.epoch))
    iend_dr3 = findlast(<=(meta_gaia_DR3.stop_mjd), vec(table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(table.epoch)
    end


    min_epoch = +Inf
    max_epoch = -Inf

    
    # DR2

    min_epoch = min(min_epoch,meta_gaia_DR2.start_mjd)
    max_epoch = max(max_epoch,meta_gaia_DR2.stop_mjd)


    μ_dr2 = [
        # dr2.parallax,
        dr2.ra,# deg
        dr2.dec,# deg
        dr2.pmra, 
        dr2.pmdec,
    ]
    σ_dr2 = [
        # dr2.parallax_error,
        dr2.ra_error/ 60/60/1000,
        dr2.dec_error / 60/60/1000 / cosd(dr2.dec),
        dr2.pmra_error ,
        dr2.pmdec_error,
    ]
    C_dr2 = [
        # plx                   ra                      dec                     pmra                    pmdec
        1                       dr2.ra_parallax_corr    dr2.dec_parallax_corr   dr2.parallax_pmra_corr  dr2.parallax_pmdec_corr
        dr2.ra_parallax_corr    1                       dr2.ra_dec_corr         dr2.ra_pmra_corr        dr2.ra_pmdec_corr
        dr2.dec_parallax_corr   dr2.ra_dec_corr         1                       dr2.dec_pmra_corr       dr2.dec_pmdec_corr
        dr2.parallax_pmra_corr  dr2.ra_pmra_corr        dr2.dec_pmra_corr       1                       dr2.pmra_pmdec_corr
        dr2.parallax_pmdec_corr dr2.ra_pmdec_corr       dr2.dec_pmdec_corr      dr2.pmra_pmdec_corr     1
    ]
    C_dr2 = C_dr2[2:end,2:end]
    # C_dr2 = [
    #     # pmra                    pmdec
    #     1                       dr2.pmra_pmdec_corr
    #     dr2.pmra_pmdec_corr     1
    # ]
    A_prepared_5_dr2 = prepare_A_5param(table[istart_dr2:iend_dr2], meta_gaia_DR2.ref_epoch_mjd,  meta_gaia_DR2.ref_epoch_mjd)

    
    # DR3
    min_epoch = min(min_epoch,meta_gaia_DR3.start_mjd)
    max_epoch = max(max_epoch,meta_gaia_DR3.stop_mjd)
    # Query Gaia archive for DR3 solution
    # μ_dr3 = [
    #     # dr3.parallax,
    #     # dr3.ra,# deg
    #     # dr3.dec,# deg
    #     dr3.pmra, 
    #     dr3.pmdec,
    # ]
    σ_dr3 = [
        # dr3.parallax_error,
        dr3.ra_error/ 60/60/1000,
        dr3.dec_error / 60/60/1000 / cosd(dr3.dec),
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
    C_dr3 = C_dr3[2:end,2:end]
    # C_dr3 = [
    #     # pmra                    pmdec
    #     1                       dr3.pmra_pmdec_corr
    #     dr3.pmra_pmdec_corr     1
    # ]
    A_prepared_5_dr3 = prepare_A_5param(table[istart_dr3:iend_dr3], meta_gaia_DR3.ref_epoch_mjd,  meta_gaia_DR3.ref_epoch_mjd)
    
    table = Table(table[min_epoch .<= table.epoch .<= max_epoch,:])


    # Build our overall correlation matrix -- the top left and bottom right blocks come
    # from the catalogs.
    Σ_dr2 = Diagonal(σ_dr2) * C_dr2 * Diagonal(σ_dr2)
    Σ_dr3 = Diagonal(σ_dr3) * C_dr3 * Diagonal(σ_dr3)
    ρ = sqrt(size(A_prepared_5_dr2,1)/size(A_prepared_5_dr3,1))
    K = ρ*sqrt(Σ_dr2)*sqrt(Σ_dr3)
    Σ_dr2_dr3 = [
        Σ_dr2  K 
        K      Σ_dr3
    ]



    return GaiaDifferenceLike{typeof(dr3),typeof(dr2),typeof(table),}(
        source_id_dr3,
        source_id_dr2,
        dr3,
        dr2,
        table,
        μ_dr2,
        Σ_dr2_dr3,
        A_prepared_5_dr3,
        A_prepared_5_dr2,
    )

end



function ln_like(gaialike::GaiaDifferenceLike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 
    ll, _ = simulate(gaialike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    return ll
end


function simulate(gaialike::GaiaDifferenceLike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    T = _system_number_type(θ_system)

    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol

    # Generate simulated observations from this sample draw
    # These are generated using whatever reference epoch the user has specified that corresponds
    # to their ra, dec, plx, etc variables

    ll = zero(T)

    # Now we fit a no-planet (zero mass planet) sky path model to this data.
    # These should be fit using the appropriate catalog reference epoch so 
    # that they can be compared correctly.


    # Helper functions to either get the static pmra from the orbital elements,
    # or, if using an AbsoluteVisualOrbit, get the propagated pmra at the
    # current epoch accounting for barycentric motion.
    function propagate_astrom(orbit::PlanetOrbits.AbsoluteVisualOrbit, epoch_ra, epoch_dec)
        sol_ra = orbitsolve(orbit, epoch_ra)
        cmp_ra = sol_ra.compensated
        if epoch_dec == epoch_ra
            sol_dec = sol_ra
        else
            sol_dec = orbitsolve(orbit, epoch_dec)
        end
        cmp_dec = sol_dec.compensated
        # return cmp_ra.ra2, cmp_dec.dec2, cmp_ra.pmra2, cmp_dec.pmdec2
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
    function propagate_astrom(orbit::Any, epoch_ra, epoch_dec)
        dec = θ_system.dec + θ_system.pmdec/60/60/1000/365.25*(epoch_dec-θ_system.ref_epoch)
        dec′ = θ_system.dec + θ_system.pmdec/60/60/1000/365.25*(epoch_ra-θ_system.ref_epoch)
        ra = θ_system.ra + θ_system.pmra/60/60/1000/365.25*(epoch_ra-θ_system.ref_epoch)/cosd(dec′)
        return ra, dec, θ_system.pmra, θ_system.pmdec
    end

    istart = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaialike.table.epoch))
    iend = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaialike.table.epoch))
    if isnothing(istart)
        istart = 1
    end
    if isnothing(iend)
        iend = length(gaialike.table.epoch)
    end
    Δα_mas = zeros(T, iend-istart+1)
    Δδ_mas = zeros(T, iend-istart+1)
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            gaialike.table[istart:iend], orbit,
            planet_mass_msol, 0.0,
            orbit_solutions[i_planet],
            orbit_solutions_i_epoch_start[i_planet], T
        )
    end
    # Δα, Δδ, Δμα, Δμδ = out.parameters
    out = fit_5param_prepared(gaialike.A_prepared_5_dr3, gaialike.table, Δα_mas, Δδ_mas)
    Δα, Δδ, Δμα, Δμδ = out.parameters
    
    delta_t_ra = (meta_gaia_DR3.ref_epoch_mjd - mean(@view gaialike.table.epoch[istart:iend]))/ julian_year
    delta_t_dec = (meta_gaia_DR3.ref_epoch_mjd - mean(@view gaialike.table.epoch[istart:iend]))/ julian_year
    Δα += delta_t_ra*Δμα
    Δδ += delta_t_dec*Δμδ

    α_g₀, δ_g₀, pmra_g₀, pmdec_g₀ = propagate_astrom(first(orbits), meta_gaia_DR3.ref_epoch_mjd, meta_gaia_DR3.ref_epoch_mjd)
    modelled_gaia_parameters_dr3 = [
        α_g₀ + Δα/60/60/1000,
        δ_g₀ + Δδ/60/60/1000/cosd(gaialike.dr3.dec),
        pmra_g₀+Δμα, #gaialike.dr3.pmra+Δμα,
        pmdec_g₀+Δμδ, #gaialike.dr3.pmdec+Δμδ
    ]



    istart = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaialike.table.epoch))
    iend = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaialike.table.epoch))
    if isnothing(istart)
        istart = 1
    end
    if isnothing(iend)
        iend = length(gaialike.table.epoch)
    end
    Δα_mas = zeros(T, iend-istart+1)
    Δδ_mas = zeros(T, iend-istart+1)
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            gaialike.table[istart:iend], orbit,
            planet_mass_msol, 0.0,
            orbit_solutions[i_planet],
            orbit_solutions_i_epoch_start[i_planet], T
        )
    end

    # out = fit_4param(
    #     gaialike.table[istart:iend],
    #     Δα_mas,
    #     Δδ_mas,
    #     meta_gaia_DR2.ref_epoch_mjd,
    #     meta_gaia_DR2.ref_epoch_mjd,
    # )
    # Δα, Δδ, Δμα, Δμδ = out.parameters
    out = fit_5param_prepared(gaialike.A_prepared_5_dr2, gaialike.table[istart:iend], Δα_mas, Δδ_mas)
    Δα, Δδ, Δμα, Δμδ = out.parameters

        
    # TODO: First order of business is checking the HGCA modelling to make sure we are propagating to the right epoch
    # for the perturbation
    # TODO: We need to propagate the positions we find away from the average epoch and to the requested epoch..?
    # TODO: what about pmra/pmdec? They are linear in this model, but non-linear in the other.
    # In a way, we measured the pmra at the measurement epoch, and extrapolated it linearly to the comparison epoch.
    # Maybe that's okay, maybe not? I guess it's okay since Gaia did it too! But we could get better accuracy by not
    # doing that.
    delta_t_ra = (meta_gaia_DR2.ref_epoch_mjd - mean(@view gaialike.table.epoch[istart:iend]))/ julian_year
    delta_t_dec = (meta_gaia_DR2.ref_epoch_mjd - mean(@view gaialike.table.epoch[istart:iend]))/ julian_year
    # @show delta_t_ra delta_t_dec Δμα Δμδ
    Δα += delta_t_ra*Δμα
    Δδ += delta_t_dec*Δμδ
    α_g₀, δ_g₀, pmra_g₀, pmdec_g₀ = propagate_astrom(first(orbits), meta_gaia_DR2.ref_epoch_mjd, meta_gaia_DR2.ref_epoch_mjd)
    modelled_gaia_parameters_dr2 = [
        α_g₀ + Δα/60/60/1000,
        δ_g₀ + Δδ/60/60/1000/cosd(gaialike.dr2.dec),
        pmra_g₀+Δμα, # gaialike.dr2.pmra+Δμα,
        pmdec_g₀+Δμδ, # gaialike.dr2.pmdec+Δμδ
    ] 


    μ_dr3 = @SVector [gaialike.dr3.ra, gaialike.dr3.dec, gaialike.dr3.pmra, gaialike.dr3.pmdec]


    # The Gaia DR2 reported parameter values are offset in various ways vs. DR2. 
    # Correct catalog values from DR2 for known effects:
    correction = @SVector [
        θ_system.dr2_systematic_Δra/60/60/1000,
        θ_system.dr2_systematic_Δdec/60/60/1000/cosd(gaialike.dr2.dec),
        θ_system.dr2_systematic_Δμ_ra,
        θ_system.dr2_systematic_Δμ_dec,
    ]
    μ_dr2_corrected = gaialike.μ_dr2 .+ correction
    μ_dr2_dr3 = [μ_dr2_corrected; μ_dr3]    
    
        # dist_dr2_dr3 = MvNormal(μ_dr2_dr3,Hermitian(gaialike.Σ_dr2_dr3))

        # μ_dr2_dr3_modelled = [modelled_gaia_parameters_dr2; modelled_gaia_parameters_dr3]
        # ll += logpdf(dist_dr2_dr3, μ_dr2_dr3_modelled)


          
    dist_dr2_dr3 = MvNormal(μ_dr2_dr3,Hermitian(gaialike.Σ_dr2_dr3))

    μ_dr2_dr3_modelled = [modelled_gaia_parameters_dr2; modelled_gaia_parameters_dr3]
    ll += logpdf(dist_dr2_dr3, μ_dr2_dr3_modelled)

    # @show μ_dr2_dr3
    # @show μ_dr2_dr3_modelled
    # @show θ_system.dr2_systematic_Δra θ_system.dr2_systematic_Δdec
    # @show (μ_dr2_dr3 .- μ_dr2_dr3_modelled) 
    # @show (μ_dr2_dr3 .- μ_dr2_dr3_modelled) ./ sqrt.(diag(gaialike.Σ_dr2_dr3))

    Δt = (Octofitter.meta_gaia_DR3.ref_epoch_mjd - Octofitter.meta_gaia_DR2.ref_epoch_mjd)/Octofitter.julian_year

    return convert(T,ll), (;
        modelled_gaia_parameters_dr3=modelled_gaia_parameters_dr3,
        modelled_gaia_parameters_dr2=modelled_gaia_parameters_dr2,#.-correction,
        
        pmra_dr3_model = modelled_gaia_parameters_dr3[3],
        pmdec_dr3_model = modelled_gaia_parameters_dr3[4],

        pmra_dr2_model = modelled_gaia_parameters_dr2[3],
        pmdec_dr2_model = modelled_gaia_parameters_dr2[4],
        pmra_dr32_model=((
            modelled_gaia_parameters_dr3[1]-modelled_gaia_parameters_dr2[1]
        )*60*60*1000*cosd((gaialike.dr3.dec+gaialike.dr2.dec)/2))/Δt,
        pmdec_dr32_model=((
            modelled_gaia_parameters_dr3[2]-modelled_gaia_parameters_dr2[2]
        )*60*60*1000)/Δt,
        correction,
    )
end

function build_correlation_structure(α_pos, α_pm, β_pos, δ_pos, δ_pm, γ)
    # Initialize 5x5 correlation matrix
    K = zeros(5, 5)
    
    # Parameter order: plx, ra, dec, pmra, pmdec
    
    # Set diagonal blocks
    K[1,1] = 0#α_plx  # parallax correlation
    K[2:3, 2:3] .= [α_pos δ_pos; δ_pos α_pos]  # position correlations
    K[4:5, 4:5] .= [α_pm δ_pm; δ_pm α_pm]  # proper motion correlations
    
    # Set parallax-position correlations
    K[1, 2:3] .= β_pos
    K[2:3, 1] .= β_pos
    
    # Set all other cross-correlations to γ
    # First row/column (excluding already set elements)
    K[1, 4:5] .= γ
    K[4:5, 1] .= γ
    
    # Position-PM correlations
    K[2:3, 4:5] .= γ
    K[4:5, 2:3] .= γ
    

    # TODO: shift indices down to account for removing parallax
    return K[2:end,2:end]
end


using LinearAlgebra

# """
# Compute a safe scaling factor for K to ensure positive definiteness of the block matrix
# [A    K  ]
# [K'   C  ]

# Returns:
# - f: scaling factor to apply to K
# - is_pd: whether the original matrix is already positive definite
# """
# function compute_safe_scaling(A::Matrix{Float64}, C::Matrix{Float64}, K::Matrix{Float64})
#     # Check inputs are same size
#     na, ma = size(A)
#     nc, mc = size(C)
#     nk, mk = size(K)
#     @assert na == ma "A must be square"
#     @assert nc == mc "C must be square"
#     @assert nk == na && mk == nc "K must be compatible with A and C"
    
#     # Compute K^T A^{-1} K
#     K_A_K = K' * (A \ K)
    
#     # Get eigenvalues
#     λ_c = eigvals(C)
#     λ_kak = eigvals(K_A_K)
    
#     # Check if already positive definite
#     λ_min_c = minimum(λ_c)
#     λ_max_kak = maximum(λ_kak)
    
#     is_pd = λ_min_c > λ_max_kak
    
#     # Compute safe scaling factor
#     if !is_pd && λ_max_kak > 0
#         f = sqrt(0.99 * λ_min_c / λ_max_kak)  # 0.99 for numerical safety
#     else
#         f = 1.0
#     end
    
#     return f, is_pd
# end

# """
# Apply safe scaling to ensure positive definiteness of the block matrix.
# Returns the scaled K matrix and whether scaling was needed.
# """
# function scale_K_safely(A::Matrix{Float64}, C::Matrix{Float64}, K::Matrix{Float64})
#     f, was_pd = compute_safe_scaling(A, C, K)
#     return f * K, was_pd
# end


