using LinearSolve

struct GaiaCatalogFitLikelihood{TTable,TCat,TDist,TFact} <: AbstractLikelihood
    # predicted observations from GHOST or other scanlaw
    table::TTable
    # Source ID from each given catalog, if available
    source_id::Int
    # Catalog values of key parameters, if available
    gaia_sol::TCat
    # precomputed MvNormal distribution
    dist::TDist
    A_prepared_4::TFact
    A_prepared_5::TFact
end

# TODO: add flux ratio var
function _getparams(::GaiaCatalogFitLikelihood{TTable,TCat,TDist}, θ_planet) where {TTable,TCat,TDist}
    # if fluxratio_var == :__dark
        return (;fluxratio=zero(Octofitter._system_number_type(θ_planet)))
    # end
    # fluxratio = getproperty(θ_planet, fluxratio_var)
    # return (;fluxratio)
end

function GaiaCatalogFitLikelihood(;
    gaia_id_dr2=nothing,
    gaia_id_dr3=nothing,
    scanlaw_table=nothing,
)
    # Query Gaia archive for DR3 solution
    if !isnothing(gaia_id_dr2)
        source_id = gaia_id_dr2
        gaia_sol = Octofitter._query_gaia_drs(; gaia_id=gaia_id_dr2)
        ref_epoch = meta_gaia_DR2.ref_epoch_mjd
    elseif !isnothing(gaia_id_dr3)
        source_id = gaia_id_dr3
        gaia_sol = Octofitter._query_gaia_dr3(; gaia_id=gaia_id_dr3)
        ref_epoch = meta_gaia_DR3.ref_epoch_mjd
    else
        throw(ArgumentError("Please provide at least one of `source_id_dr1`, `source_id_dr2`, or `source_id_dr3`"))
    end

    ra_deg = gaia_sol.ra
    dec_deg = gaia_sol.dec
    μ = [
        gaia_sol.parallax,
        gaia_sol.ra,# deg
        gaia_sol.dec,# deg
        gaia_sol.pmra,
        gaia_sol.pmdec,
    ]
    σ = [
        gaia_sol.parallax_error,
        gaia_sol.ra_error / 60 / 60 / 1000 / cosd(gaia_sol.dec),
        gaia_sol.dec_error / 60 / 60 / 1000,
        gaia_sol.pmra_error,
        gaia_sol.pmdec_error,
    ]
    C = [
        # plx                   ra                      dec                     pmra                    pmdec
        1 gaia_sol.ra_parallax_corr gaia_sol.dec_parallax_corr gaia_sol.parallax_pmra_corr gaia_sol.parallax_pmdec_corr
        gaia_sol.ra_parallax_corr 1 gaia_sol.ra_dec_corr gaia_sol.ra_pmra_corr gaia_sol.ra_pmdec_corr
        gaia_sol.dec_parallax_corr gaia_sol.ra_dec_corr 1 gaia_sol.dec_pmra_corr gaia_sol.dec_pmdec_corr
        gaia_sol.parallax_pmra_corr gaia_sol.ra_pmra_corr gaia_sol.dec_pmra_corr 1 gaia_sol.pmra_pmdec_corr
        gaia_sol.parallax_pmdec_corr gaia_sol.ra_pmdec_corr gaia_sol.dec_pmdec_corr gaia_sol.pmra_pmdec_corr 1
    ]
    Σ = Diagonal(σ) * C * Diagonal(σ)
    # dist = MvNormal(μ[2:end], Hermitian(Σ)[2:end,2:end])
    dist = MvNormal(μ, Hermitian(Σ))




    if isnothing(scanlaw_table)
        @info "No scan law table provided. We will fetch an approximate solution from the GHOST webservice."
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GHOST_forecast(ra_deg, dec_deg))
        forecast_table.epoch = jd2mjd.(forecast_table.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        forecast_table.scanAngle_rad = forecast_table.scanAngle_rad_
    else
        @info "Scanlaw table was provided, will not query GHOST."
        forecast_table = FlexTable(scanlaw_table)
        forecast_table.epoch = tcb_at_gaia_2mjd.(forecast_table.times)
        forecast_table.scanAngle_rad = deg2rad.(forecast_table.angles)
    end

    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π / 2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π / 2 .+ forecast_table.scanAngle_rad)

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    table = Table(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)

    # Prepare some factorized matrices for linear system solves
    A_prepared_4 = prepare_A_4param(table, ref_epoch, ref_epoch)
    A_prepared_5 = prepare_A_5param(table, ref_epoch, ref_epoch)


    return GaiaCatalogFitLikelihood(
        table,
        source_id,
        gaia_sol,
        dist,
        A_prepared_4,
        A_prepared_5,
    )
end



function ln_like(gaialike::GaiaCatalogFitLikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    ll, _ = simulate(gaialike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    return ll
end

function simulate(gaialike::GaiaCatalogFitLikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

    T = _system_number_type(θ_system)

    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass * Octofitter.mjup2msol

    sample_orbit = only(orbits)

    ll = zero(T)

    Δα_mas = zeros(T, size(gaialike.table,1))
    Δδ_mas = zeros(T, size(gaialike.table,1))

    for (orbit, θ_planet) in zip(orbits, θ_system.planets)
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        (;fluxratio) = _getparams(gaialike, θ_planet)
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            gaialike.table, orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions,
            orbit_solutions_i_epoch_start, T
        )
    end

    function propagate_astrom(orbit::PlanetOrbits.AbsoluteVisualOrbit, ref_epoch)
        sol = orbitsolve(orbit, ref_epoch)
        cmp = sol.compensated
        return cmp.ra2, cmp.dec2, cmp.pmra2, cmp.pmdec2
    end
    function propagate_astrom(orbit::Any, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end


    ref_epoch = years2mjd(gaialike.gaia_sol.ref_epoch)
    out = fit_5param(
        gaialike.table,
        Δα_mas,
        Δδ_mas,
        ref_epoch,
        ref_epoch,
    )
    Δα_g, Δδ_g, Δplx, Δpmra_g, Δpmdec_g = out.parameters
    α_g₀, δ_g₀, pmra_g₀, pmdec_g₀ = propagate_astrom(first(orbits), ref_epoch)
    modelled_gaia_parameters_dr3 = @SVector [
        gaialike.gaia_sol.parallax + Δplx,
        α_g₀ + Δα_g/60/60/1000/cosd(δ_g₀),
        δ_g₀ + Δδ_g/60/60/1000,
        pmra_g₀ + Δpmra_g,
        pmdec_g₀ + Δpmdec_g
    ]

    
    # @show resid = gaialike.dist.μ .- modelled_gaia_parameters_dr3
    ll += logpdf(gaialike.dist, modelled_gaia_parameters_dr3)

    return ll, modelled_gaia_parameters_dr3
end


function fit_5param(
    table,
    Δα_mas,
    Δδ_mas,
    reference_epoch_mjd_ra,
    reference_epoch_mjd_dec;
    include_chi2=false,
    σ_formal=0
)
    # n_obs = size(table, 1)

    # T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    # A = zeros(T, n_obs, 5)
    # b = zeros(T, n_obs)

    # for i in 1:n_obs
    #     # Position terms
    #     A[i, 1] = table.cosϕ[i]  # α
    #     A[i, 2] = table.sinϕ[i]  # δ

    #     # Parallax term
    #     A[i, 3] = -table.parallaxFactorAlongScan[i]

    #     # Proper motion terms
    #     A[i, 4] = table.cosϕ[i] * (table.epoch[i] - reference_epoch_mjd_ra)/ julian_year  # μα*
    #     A[i, 5] = table.sinϕ[i] * (table.epoch[i] - reference_epoch_mjd_dec)/ julian_year # μδ

    #     # Along-scan measurement
    #     b[i] = Δα_mas[i] * table.cosϕ[i] + Δδ_mas[i] * table.sinϕ[i]
    # end

    # # Calculate model predictions
    # # x = F \ b
    # x = A\b

    # if !include_chi2
    #     return (; parameters=x)  # degrees
    # end

    # # Analytical calculation
    # # Calculate the projection matrix P = I - A(A'A)^(-1)A'
    # # Use SVD for numerical stability

    # F = svd(A)
    # P = I - F.U * F.U'  # More stable than direct calculation

    # #####
    # # If separate  noise per epoch
    # # # Calculate Σ^(-1/2) * P * b_true
    # # # This is the "signal" contribution
    # # Σinv_sqrt = Diagonal(fill(1.0 / σ_formal, length(b)))
    # # signal_contrib = Σinv_sqrt * P * b
    # # signal_term = dot(signal_contrib, signal_contrib)

    # # # Calculate trace(Σ^(-1) * P * Σ)
    # # # This is the "noise" contribution - degrees of freedom
    # # noise_term = tr(P)  # Since Σ^(-1) * P * Σ = P for our case

    # #####
    # # If alsways same noise per epoch:
    # # Signal term simplifies to (P*b_true)^2/σ^2
    # signal_contrib = P * b
    # signal_term = dot(signal_contrib, signal_contrib) / (σ_formal^2)
    # noise_term = tr(P)

    # expected_chisq = signal_term + noise_term

    # # model_predictions = A * x
    # # residuals = b - model_predictions

    # return (;
    #     parameters=x,
    #     reference_position=(α_ref, δ_ref),  # degrees
    #     chi_squared_astro=expected_chisq,
    # )

    n_obs = size(table, 1)

    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    # Use Bumper to elide allocations
    @no_escape begin
        A = @alloc(T, n_obs, 5)
        b = @alloc(T, n_obs)
        x = @alloc(T, 5)

       
        for i in 1:n_obs
            # Position terms
            A[i, 1] = table.cosϕ[i]  # α
            A[i, 2] = table.sinϕ[i]  # δ

            #     # Parallax term
            A[i, 3] = -table.parallaxFactorAlongScan[i]

            # Proper motion terms
            A[i, 4] = table.cosϕ[i] * (table.epoch[i] - reference_epoch_mjd_ra)/ julian_year  # μα*
            A[i, 5] = table.sinϕ[i] * (table.epoch[i] - reference_epoch_mjd_dec)/ julian_year # μδ

            # Along-scan measurement
            b[i] = Δα_mas[i] * table.cosϕ[i] + Δδ_mas[i] * table.sinϕ[i]
        end

        # Calculate model predictions
        if σ_formal == 0.
            # Straight-forward solution
            # x .= A\b

            # This in-place version uses way less memory, even though it is about 10% slower.
            # Net win when using more than one thread.
            Q = qr!(A)
            ldiv!(x, Q, b)

        # Case of unceratinties provided.
        else
            # Straight-forward solution
            # d = @alloc(T, n_obs)
            # d .= 1.0 ./ (σ_formal .^ 2)
            # W = Diagonal(1.0 ./ (σ_formal .^ 2))
            # x = (A' * W * A) \ (A' * W * b) 

            A_w = @alloc(T, size(A)...)
            b_w = @alloc(T, size(b)...)
            # In-place weighted solution
            @. A_w = A ./σ_formal  # Weight design matrix in-place
            @. b_w = b ./σ_formal  # Weight observations in-place
            Q = qr!(A_w)
            ldiv!(x, Q, b_w)

            # undo for later
        end

        parameters = @SVector [x[1], x[2], x[3], x[4], x[5]]
    end

    model_predictions = A * x
    residuals = b - model_predictions
    if !include_chi2
        return (; parameters, residuals) 
    end

    # Analytical calculation
    # Calculate the projection matrix P = I - A(A'A)^(-1)A'
    # Use SVD for numerical stability

    # TODO: is SVD correct here?
    F = svd(A)
    P = I - F.U * F.U'  # More stable than direct calculation

    #####
    # If separate  noise per epoch
    # # Calculate Σ^(-1/2) * P * b_true
    # # This is the "signal" contribution
    # Σinv_sqrt = Diagonal(fill(1.0 / σ_formal, length(b)))
    # signal_contrib = Σinv_sqrt * P * b
    # signal_term = dot(signal_contrib, signal_contrib)

    # # Calculate trace(Σ^(-1) * P * Σ)
    # # This is the "noise" contribution - degrees of freedom
    # noise_term = tr(P)  # Since Σ^(-1) * P * Σ = P for our case

    #####
    # If alsways same noise per epoch:
    # Signal term simplifies to (P*b_true)^2/σ^2
    signal_contrib = P * b
    signal_term = dot(signal_contrib, signal_contrib) / (σ_formal^2)
    noise_term = tr(P)

    # @show signal_term noise_term

    expected_chisq = signal_term + noise_term


    # @show std(b)
    # @show std(residuals)
    # @show residuals


    return (;
        parameters,
        chi_squared_astro=expected_chisq,
    )
end



function prepare_A_4param(
    table,
    reference_epoch_mjd_ra,
    reference_epoch_mjd_dec,
    σ_formal=0.
)
    n_obs = size(table, 1)

    A = zeros(n_obs, 4)
    for i in 1:n_obs
        # Position terms
        A[i, 1] = table.cosϕ[i]  # α
        A[i, 2] = table.sinϕ[i]  # δ

        # Proper motion terms
        A[i, 3] = table.cosϕ[i] * (table.epoch[i] - reference_epoch_mjd_ra)/ julian_year  # μα*
        A[i, 4] = table.sinϕ[i] * (table.epoch[i] - reference_epoch_mjd_dec)/ julian_year # μδ
    end


    if σ_formal != 0.
        @. A = A .* 1 ./ σ_formal
    end
    return A
    # return qr!(A)
    # b = zeros(4)
    # prob =  LinearProblem(A, b, OperatorAssumptions(false, condition=LinearSolve.OperatorCondition.WellConditioned))
    # linsolve = init(prob)
    # return linsolve
end


function prepare_A_5param(
    table,
    reference_epoch_mjd_ra,
    reference_epoch_mjd_dec,
    σ_formal=0.
)
    n_obs = size(table, 1)

    A = zeros(n_obs, 5)
    for i in 1:n_obs
        # Position terms
        A[i, 1] = table.cosϕ[i]  # α
        A[i, 2] = table.sinϕ[i]  # δ

        # Parallax term
        A[i, 3] = -table.parallaxFactorAlongScan[i]

        # Proper motion terms
        A[i, 4] = table.cosϕ[i] * (table.epoch[i] - reference_epoch_mjd_ra)/ julian_year  # μα*
        A[i, 5] = table.sinϕ[i] * (table.epoch[i] - reference_epoch_mjd_dec)/ julian_year # μδ
    end
    if σ_formal != 0.
        @. A = A .* 1 ./ σ_formal
    end
    return A#qr!(A)
    # b = zeros(4)
    # prob =  LinearProblem(A, b, OperatorAssumptions(false, condition=LinearSolve.OperatorCondition.WellConditioned))
    # linsolve = init(prob)
    # return linsolve
end


function fit_4param_prepared(
    A_factored,
    table,
    Δα_mas,
    Δδ_mas,
    σ_formal=0.0;
    # include_chi2=false,
)
    n_obs = size(table, 1)

    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    # Use Bumper to elide allocations
    # @no_escape begin
    # begin
        b = zeros(T, n_obs)
        x = zeros(T, 4)
       
        for i in 1:n_obs
            # Along-scan measurement
            # b[i] = Δα_mas[i] * table.cosϕ[i] - Δδ_mas[i] * table.sinϕ[i]
            b[i] = Δα_mas[i] * table.cosϕ[i] + Δδ_mas[i] * table.sinϕ[i]
        end



        if σ_formal != 0.
            # weighted solution (allocates, to avoid re-writing)
            # A is already weighted at preparation time
            @. b *= 1/σ_formal  # Weight observations in-place
        end

        # differentiable_calc_x = DifferentiateWith(AutoFiniteDiff()) do b
        #     x = A \ b
        # end



        # Straight-forward solution
        x = A_factored \ b

        # differentiable_calc_x = let A = A_factored
        #     DifferentiateWith(AutoEnzyme(mode=Main.Enzyme.set_runtime_activity(Main.Enzyme.Forward), function_annotation=Main.Enzyme.Const)) do b
        #     # DifferentiateWith(AutoMooncake(config=nothing)) do b
        #     # DifferentiateWith(AutoForwardDiff()) do b
        #         x = A \ b
        #         # prob = LinearProblem(A, b,)# assumptions=OperatorAssumptions(false, condition=LinearSolve.OperatorCondition.WellConditioned))
        #         # sol = solve(prob)
        #         # return sol.u
        #     end
        # end
        # x = differentiable_calc_x(b)
        # prob = LinearProblem(A_factored, b, assumptions=OperatorAssumptions(false, condition=LinearSolve.OperatorCondition.WellConditioned))
        # sol = solve(prob)
        # x = sol.u

        # # Q = qr!(A)
        # # ldiv!(x, A_factored, b)
        # # @show typeof(A_factored)
        # # x = A_factored \ b

        # prob = LinearProblem(A_factored, b, OperatorAssumptions(false, condition=LinearSolve.OperatorCondition.WellConditioned))
        # # A_factored.b = b
        # sol = solve(prob)#, KrylovJL_LSMR())
        # x = sol.u

        parameters = @SVector [x[1], x[2], x[3], x[4]]
    # end

    # if !include_chi2
        return (; parameters) 
    # end

    # # Analytical calculation
    # # Calculate the projection matrix P = I - A(A'A)^(-1)A'
    # # Use SVD for numerical stability

    # # TODO:
    # F = svd(A_factored)
    # P = I - F.U * F.U'  # More stable than direct calculation

    # #####
    # # If separate  noise per epoch
    # # # Calculate Σ^(-1/2) * P * b_true
    # # # This is the "signal" contribution
    # # Σinv_sqrt = Diagonal(fill(1.0 / σ_formal, length(b)))
    # # signal_contrib = Σinv_sqrt * P * b
    # # signal_term = dot(signal_contrib, signal_contrib)

    # # # Calculate trace(Σ^(-1) * P * Σ)
    # # # This is the "noise" contribution - degrees of freedom
    # # noise_term = tr(P)  # Since Σ^(-1) * P * Σ = P for our case

    # #####
    # # If alsways same noise per epoch:
    # # Signal term simplifies to (P*b_true)^2/σ^2
    # signal_contrib = P * b
    # signal_term = dot(signal_contrib, signal_contrib) / (σ_formal^2)
    # noise_term = tr(P)

    # expected_chisq = signal_term + noise_term

    # # model_predictions = A * x
    # # residuals = b - model_predictions

    # return (;
    #     parameters,
    #     reference_position=(α_ref, δ_ref),  # degrees
    #     chi_squared_astro=expected_chisq,
    # )
end

using HTTP

const meta_gaia_DR1 = (;
    start_mjd=mjd("2014-07-25"),
    stop_mjd =mjd("2015-09-16"),
    ref_epoch_mjd=57023.25 # J2015.0
)

const meta_gaia_DR2 = (;
    start_mjd=mjd("2014-07-25"), # 10:30 UTC
    stop_mjd =mjd("2016-05-23"), # 11:35 UTC
    ref_epoch_mjd=57205.875 # J2015.5
)

const meta_gaia_DR3 = (;
    start_mjd=mjd("2014-07-25"), # 10:30 UTC
    stop_mjd =mjd("2017-05-28"), # 11:35 UTC
    ref_epoch_mjd=57388.5, # J2016.0
)

# OBMT (On board mission time) to MJD converter
tcb_at_gaia_2mjd(tcb_gaia) = jd2mjd(tcb_gaia+2455197.5)



"""
Given scan epochs and angles, and an AbsoluteVisual orbit describing a sky path and perturbations
from a planet, calculate a "true" model of what Gaia would have observed, including
full spherical coordiante effects, radial velocity (perspective acceleration), etc.

TODO: include effects of changing light travel time.
"""
function _simulate_skypath_observations(
    gaialike,
    sample_orbit::AbsoluteVisual,
    planet_mass_msol,
    cat_ra,
    cat_deg,
    orbit_solutions, orbit_solutions_i_epoch_start, T=Float64;
)

    # TODO: make use of potentially multi-threaded kepsolve by ensuring `epoch` column is present,
    # and using above passed-in orbit solutions.

    # Compute the Ra and Dec sky path of the star, accounting for 3D barycentric
    # motion in spherical coordinates and perturbations from planets
    α_model = zeros(T, size(gaialike.table,1))
    δ_model = zeros(T, size(gaialike.table,1))
    for i in eachindex(gaialike.table.epoch)
        sol = orbitsolve(sample_orbit, gaialike.table.epoch[i])
        cmp = sol.compensated

        # # Calculate the position of the star in cartesian coordinates
        # # Get earth position in cartesian coordaintes (in units of AU)
        # # Calculate apparent Ra and Dec.
        # x_earth_pc = gaialike.table.x[i] / PlanetOrbits.pc2au
        # y_earth_pc = gaialike.table.y[i] / PlanetOrbits.pc2au
        # z_earth_pc = gaialike.table.z[i] / PlanetOrbits.pc2au
        # x_diff_pc = cmp.x₂ - x_earth_pc
        # y_diff_pc = cmp.y₂ - y_earth_pc
        # z_diff_pc = cmp.z₂ - z_earth_pc
        # distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
        # mydtor = π / 180
        # ra_apparent_deg = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
        # arg = z_diff_pc / distance_diff
        # arg = map(arg) do arg
        #     if 1.0 < arg < 1.0 + sqrt(eps(1.0))
        #         arg = 1.0
        #     end
        #     return arg
        # end
        # dec_apparent_deg = asin(arg) / mydtor

        # @show cmp.dec2
        # @show sample_orbit.dec
        # @show sample_orbit.ref_epoch
        # @show gaialike.table.epoch[i]
        # @show sol.t
        α = cmp.ra2
        δ = cmp.dec2
        plx_at_epoch = cmp.parallax2
        α_model[i] = (α - cat_ra).*60*60*1000 *cosd(δ) + plx_at_epoch * (
            gaialike.table.x[i] * sind(α) -
            gaialike.table.y[i] * cosd(α)
        )
        δ_model[i] = (δ - cat_deg).*60*60*1000  + plx_at_epoch * (
            gaialike.table.x[i] * cosd(α) * sind(δ) +
            gaialike.table.y[i] * sind(α) * sind(δ) -
            gaialike.table.z[i] * cosd(δ)
        )

        # TODO: add perturbations from multiple planets here
        for orb in (sample_orbit,)
            sol = orbitsolve(orb, gaialike.table.epoch[i])
            # Add perturbation from planet
            α_model[i] += raoff(sol, planet_mass_msol)
            δ_model[i] += decoff(sol, planet_mass_msol)
        end
    end

    # Now we need to generate a line (of any length) that passes through these points with an
    # angle given by psi. 

    # TODO: confirm we correctly accounted for a cos(delta) when generating these line segments
    αₘ = eachrow(@. T[-1, 1]'/60/60/1000/cos(sample_orbit.dec) * gaialike.table.sinϕ)
    for i in eachindex(α_model)
        αₘ[i] .+= α_model[i]
    end

    δₘ = eachrow(@. T[1, -1]'/60/60/1000 * gaialike.table.cosϕ)
    for i in eachindex(δ_model)
        δₘ[i] .+= δ_model[i]
    end

    return (;α_model,δ_model,αₘ,δₘ)

end



function _simulate_skypath_perturbations!(
    Δα_model, Δδ_model,
    table,
    orbit::AbstractOrbit,
    planet_mass_msol,
    flux_ratio,
    orbit_solutions, orbit_solutions_i_epoch_start, T=Float64;
)
    if flux_ratio != 0
        throw(NotImplementedException())
    end
    for i in eachindex(table.epoch)
        # TODO: make use of potentially multi-threaded kepsolve by ensuring `epoch` column is present,
        # and using above passed-in orbit solutions.
        sol = orbitsolve(orbit, table.epoch[i])
        # Add perturbation from planet in mas -- simplest model, unaffected by flux of planet
        # Δα_model[i] += raoff(sol, planet_mass_msol)
        # Δδ_model[i] += decoff(sol, planet_mass_msol)

        # Add perturbation to photocentre from (possibly luminour) planet in mas
        # calculate absolute location of the primary (have this already)
        # relative location of the secondary (have this already)
        # calculate the flux-weighted offset between these, and add it to the primary
        # In the below, the following relationships hold
        # raoff(sol, planet_mass_msol) : offset of host from barycentre due to planet
        # raoff(sol) : offset of planet from host 
        ra_host_vs_bary = raoff(sol, planet_mass_msol)
        ra_planet_vs_host = raoff(sol)
        ra_planet_vs_bary = ra_host_vs_bary + ra_planet_vs_host
        ra_photocentre = ra_host_vs_bary * (1-flux_ratio) + ra_planet_vs_bary * flux_ratio
        Δα_model[i] += ra_photocentre
        dec_host_vs_bary = decoff(sol, planet_mass_msol)
        dec_planet_vs_host = decoff(sol)
        dec_planet_vs_bary = dec_host_vs_bary + dec_planet_vs_host
        dec_photocentre = dec_host_vs_bary * (1-flux_ratio) + dec_planet_vs_bary * flux_ratio
        Δδ_model[i] += dec_photocentre
    end
    return
end



function _query_gaia_dr3(;gaia_id)
    fname = "_gaia_dr3/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiaedr3.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false,
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr3")
            mkdir("_gaia_dr3")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    if length(data) < length(headers)
        error("Could not query DR3 for source")
    end
    return namedtuple(headers, data)
end

function _query_gaia_dr2(;gaia_id)
    fname = "_gaia_dr2/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr2.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr2")
            mkdir("_gaia_dr2")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    if length(data) <= 1
        error("Could not query DR2 for source")
    end
    return namedtuple(headers, data)
end

# TODO: can we add DR1 also? The challenge is that it has different source_id for the same object

function _query_gaia_dr1(;gaia_id)
    fname = "_gaia_dr1/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr1.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr1")
            mkdir("_gaia_dr1")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    if length(data) <= 1
        error("could not find source")
    end
    return namedtuple(headers, data)
end

"""
    geocentre_position_query(epoch_MJD)

Given a date+time in MJD format, return a named tuple of Earth position and velocity in AU 
on that date. The results are cached in a local `_geocentre_pos` directory for offline use.

Currently these are queried from the NASA HORIZONS online system, but this detail is subject
to change in future minor versions of the package.

The positions and velocities represent the Geocenter of the Earth relative to the solar sysytem
barycenter.
"""
function geocentre_position_query(epoch_MJD::Number)
    
    fname = "_geocentre_pos/HORIZONS-Geocenter-$(epoch_MJD)mjd.txt"
    if !isfile(fname)
        if !isdir("_geocentre_pos")
            mkdir("_geocentre_pos")
        end
        # Query the HORIZONS system for the Earth-Moon barycentre position at each specific epoch
        # of the dataset. 
        # This incurs a separate query for each epoch.
        t = DateTime(mjd2date(epoch_MJD))
        @info "Querying Earth Geocentre position from HORIZONS" epoch_MJD date=t
        HORIZONS.vec_tbl("Geocenter", t, t + Hour(1), Day(1); FILENAME=fname, CENTER="@ssb", REF_PLANE="FRAME", OUT_UNITS="AU-D", CSV_FORMAT=true, VEC_TABLE=2)
    end
        
    lines = readlines(fname)
    i = 0
    for line in lines
        i += 1
        if startswith(line, raw"$$SOE")
            break
        end
    end
    record = split(rstrip(lines[i+1], ','), ", ")
    x, y, z, vx, vy, vz = parse.(Float64, record[3:8])
    
    return (; x, y, z, vx, vy, vz)
end



"""
    forecast_table = GHOST_forecast(ra_deg,dec_deg)

Given an Ra and Dec position, retreive a forecast of Gaia observations from the GHOST tool automatically.
See tool URL here: https://gaia.esac.esa.int/gost/

Please be aware that others  might be able to discover the target coordinates you searched for
(though not who performed the search) via information leaked to the external service.
"""
function GHOST_forecast(ra_deg,dec_deg)

    fname = "GHOST-$ra_deg-$dec_deg.csv"
    if isfile(fname)
        forecast_table = CSV.read(fname, Table, normalizenames=true)
        return forecast_table
    end

    #=
    curl \
    --form "srcname=eps indi test" \
    --form "inputmode=single" \
    --form "srcra=330.84022344234336" \
    --form "srdec=-56.785978554364995" \
    --form "from=2014-07-26T00:00:00" \
    --form "to=2017-05-28T00:00:00" \
    --cookie "JSESSIONID=73EE915E5F9FF8B5376D82FC116F765D" \
    'https://gaia.esac.esa.int/gost/GostServlet' 

    curl --cookie "JSESSIONID=73EE915E5F9FF8B5376D82FC116F765D" 'https://gaia.esac.esa.int/gost/export.jsp?id=73EE915E5F9FF8B5376D82FC116F765D/1722634273&format=csv'
    =#

    # Just pick a cookie ID to use. 
    # Might be better to let the service create one for us.
    cookiejar = HTTP.CookieJar()
    @info "Contacting the GAIA scan forecast tool GHOST: https://gaia.esac.esa.int/gost/"
    resp0 = HTTP.get(
        "https://gaia.esac.esa.int/gost/",
        cookiejar=cookiejar,
    )
    if resp0.status != 200
        println(String(resp0.body))
        error("Could not contact the GAIA scan forecast tool GHOST https://gaia.esac.esa.int See above error message.")
    end
    formdata = Dict([
        "serviceCode"=>"1",
        "inputmode"=>"single",
        "srcname"=>"009",
        "srcra" => string(round(ra_deg,digits=7)),
        "srcdec" => string(round(dec_deg,digits=7)),
        "from" => "2014-07-25T10:31:26",
        "to" => "2017-05-28T00:00:00",
    ])


    @info "Retrieving forecasted GAIA scans from GHOST: https://gaia.esac.esa.int/gost/"
    resp = HTTP.post(
        "https://gaia.esac.esa.int/gost/GostServlet",
        body=HTTP.Form(formdata),
        cookiejar=cookiejar
    )
    if resp.status != 200 || contains(String(collect(resp.body)),"error")
        println(String(resp.body))
        error("Could not fetch GAIA scan forecast from GHOST. See above error message. Do you have an internet connection available?")
    end

    m = match(r"Submitted with id (\d+)", String(resp.body))
    response_id = m.captures[1]
    session_id = cookiejar.entries["gaia.esac.esa.int"]["gaia.esac.esa.int;/gost;JSESSIONID"].value
    url = "https://gaia.esac.esa.int/gost/export.jsp?id=$session_id/$response_id&format=csv"
    resp_dat = HTTP.get(
        url,
        # cookies=Dict([
        #     # "Cookie" => "$new_cookie",
        #     # "JSESSIONID" => "73EE915E5F9FF8B5376D82FC116F765D",
        #     "noagreement" => "false"
        # ]),
        cookiejar=cookiejar
    )
    @info "done"
    body = collect(resp_dat.body)
    
    # Save for offline use eg on clusters
    write(fname, body)

    io = IOBuffer(body)
    forecast_table = CSV.read(io, Table, normalizenames=true)

    return forecast_table
end


