# using LinearSolve  # Currently unused - all LinearSolve code is commented out
using SPICE
using DataDeps
using Dates

struct GaiaCatalogFitObs{TTable,TCat,TDist,TFact} <: AbstractObs
    # predicted observations from GOST or other scanlaw
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

# Backwards compatibility alias
const GaiaCatalogFitLikelihood = GaiaCatalogFitObs

export GaiaCatalogFitObs, GaiaCatalogFitLikelihood

# TODO: add flux ratio var
function _getparams(::GaiaCatalogFitObs{TTable,TCat,TDist}, θ_planet) where {TTable,TCat,TDist}
    # if fluxratio_var == :__dark
        return (;fluxratio=zero(Octofitter._system_number_type(θ_planet)))
    # end
    # fluxratio = getproperty(θ_planet, fluxratio_var)
    # return (;fluxratio)
end

function GaiaCatalogFitObs(;
    gaia_id_dr2=nothing,
    gaia_id_dr3=nothing,
    scanlaw_table=nothing,
    ref_epoch_ra=nothing,
    ref_epoch_dec=nothing
)
    # Query Gaia archive for DR3 solution
    if !isnothing(gaia_id_dr2)
        source_id = gaia_id_dr2
        gaia_sol = Octofitter._query_gaia_drs(; gaia_id=gaia_id_dr2)
        if isnothing(ref_epoch_ra)
            ref_epoch_ra = meta_gaia_DR2.ref_epoch_mjd
        end
        if isnothing(ref_epoch_dec)
            ref_epoch_dec = meta_gaia_DR2.ref_epoch_mjd
        end
    elseif !isnothing(gaia_id_dr3)
        source_id = gaia_id_dr3
        gaia_sol = Octofitter._query_gaia_dr3(; gaia_id=gaia_id_dr3)
        if isnothing(ref_epoch_ra)
            ref_epoch_ra = meta_gaia_DR3.ref_epoch_mjd
        end
        if isnothing(ref_epoch_dec)
            ref_epoch_dec = meta_gaia_DR3.ref_epoch_mjd
        end
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
        @info "No scan law table provided. We will fetch an approximate solution from the GOST webservice."
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GOST_forecast(ra_deg, dec_deg))
        forecast_table.epoch = jd2mjd.(forecast_table.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        forecast_table.scanAngle_rad = forecast_table.scanAngle_rad_
    else
        @info "Scanlaw table was provided, will not query GOST."
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

    # Now remove any known gaps -- data sourced from HTOF.py; authors G.M. Brandt et al
    gaps_dr2 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiadr2_08252020.csv"), FlexTable)
    gaps_edr23 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiaedr3_12232020.csv"), FlexTable)
    gaps = Table(
        start_mjd=obmt2mjd.(vcat(gaps_dr2.start,gaps_edr23.start)),
        stop_mjd=obmt2mjd.(vcat(gaps_dr2.end,gaps_edr23.end)),
        note=[gaps_dr2.comment; gaps_edr23.description]
    )
    table = filter(eachrow(table)) do row
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
    table = Table(map(dat->dat[], table))

    # Prepare some factorized matrices for linear system solves
    A_prepared_4 = prepare_A_4param(table, ref_epoch_ra, ref_epoch_dec)
    A_prepared_5 = prepare_A_5param(table, ref_epoch_ra, ref_epoch_dec)


    return GaiaCatalogFitObs(
        table,
        source_id,
        gaia_sol,
        dist,
        A_prepared_4,
        A_prepared_5,
    )
end



function ln_like(gaialike::GaiaCatalogFitObs, context::SystemObservationContext)
    (; θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = context
    ll, _ = simulate(gaialike, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    return ll
end

function simulate(gaialike::GaiaCatalogFitObs, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

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
        t1 = ref_epoch
        Δt = 100
        t2 = t1 + Δt
        sol′ = orbitsolve(orbit,t2)
        # This isn't right! This is double counting the proper motion which already goes into ra/dec
        # Take change in delta_time and multiply it by pmra/pmdec
        diff_lt_app_pmra = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmra2
        diff_lt_app_pmdec = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmdec2
        return cmp.ra2, cmp.dec2, cmp.pmra2+diff_lt_app_pmra, cmp.pmdec2+diff_lt_app_pmdec
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
            try
                ldiv!(x, Q, b)
            catch LAPACKException
                fill!(x, NaN)
            end

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
            try
                ldiv!(x, Q, b_w)
            catch LAPACKException
                fill!(x, NaN)
            end

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
    # P = I - A*pinv((A'*A))*A'
    # F = svd(A)
    # P = I - F.U * F.U'  # More stable than direct calculation

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
    # signal_contrib = P * b
    # signal_term = dot(signal_contrib, signal_contrib) / (σ_formal^2)
    # noise_term = tr(P)

    # expected_chisq = signal_term + noise_term


    # @show std(b)
    # @show std(residuals)
    # @show residuals

    # For uniform errors, the weighted residuals are just residuals/σ
    weighted_residuals = residuals ./ σ_formal
    
    # Chi-squared is the sum of squared weighted residuals
    chisq_residual = dot(weighted_residuals, weighted_residuals)
    expected_chisq = chisq_residual

    return (;
        parameters,
        chi_squared_astro=expected_chisq,
    )
end



function prepare_A_4param(
    table,
    reference_epoch_mjd_ra,
    reference_epoch_mjd_dec,
    # σ_formal=0.
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


    # if σ_formal != 0.
    #     @. A = A .* 1 ./ σ_formal
    # end
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
    # σ_formal=0.
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
    # if σ_formal != 0.
    #     @. A = A .* 1 ./ σ_formal
    # end

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

function fit_5param_prepared(
    A_prepared,
    table,
    Δα_mas,
    Δδ_mas,
    residuals=0.0,
    σ_formal=0.0;
    include_chi2=Val(false),
)
    n_obs = size(table, 1)

    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    # Use Bumper to elide allocations
    @no_escape begin
        b_weighted = @alloc(T, n_obs)
        A_weighted = @alloc(T, n_obs, size(A_prepared,2))
        # x = @alloc(T, 5)
       
        # for i in 1:n_obs
        #     # Along-scan measurement
        #     # b[i] = Δα_mas[i] * table.cosϕ[i] - Δδ_mas[i] * table.sinϕ[i]
        # end
        @. b_weighted = Δα_mas * table.cosϕ + Δδ_mas * table.sinϕ + residuals


        if σ_formal != 0.
            @. A_weighted = A_prepared .* 1 ./ σ_formal
            @. b_weighted *= 1/σ_formal
        else
            @. A_weighted = A_prepared
        end

        # differentiable_calc_x = DifferentiateWith(AutoFiniteDiff()) do b
        #     x = A \ b
        # end

        # Straight-forward solution
        x = A_weighted \ b_weighted

        # Q = qr!(A_weighted)
        # try
        #     ldiv!(x, Q, b_w)
        # catch LAPACKException
        #     fill!(x, NaN)
        # end

        parameters = @SVector [x[1], x[2], x[4], x[5], x[3]]

        if include_chi2 == Val(true)
            model_predictions = @alloc(T, n_obs)
            residuals = @alloc(T, n_obs)
            mul!(model_predictions, A_weighted, x)
            residuals .= b_weighted .- model_predictions
            if σ_formal == 0
                error("Asked for `include_chi2=true` but `σ_formal==0`")
            end
            # For uniform errors, the weighted residuals are just residuals/σ
        
            # Chi-squared is the sum of squared weighted residuals
            # Used for UEVA calculations
            chi_squared_astro = dot(residuals, residuals)

            # Determine parameter uncertainties

            # Compute the uncertainties based on sigma_formal
            # Compute (A^T W A)^{-1} for the covariance matrix
            # W = I/σ_formal² for uniform errors
            # AtWA = A_weighted' * A_weighted
            # cov_matrix = inv(AtWA)
            # parameter_uncertainties = sqrt.(diag(cov_matrix))

            # Chi-squared is the sum of squared weighted residuals
            chi_squared_astro = dot(residuals, residuals)

            # Determine parameter uncertainties
            n_parameters = 5
            dof = length(b_weighted) - n_parameters
            
            # Reduced chi-squared
            chi2_reduced = chi_squared_astro / dof
            
            # @show chi_squared_astro chi2_reduced inflation_factor
        end
    end

    if include_chi2 != Val(true)
        return (; parameters) 
    end
    return (;
        parameters,
        chi_squared_astro,
        chi2_reduced,
        dof
    )

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
    for i in eachindex(table.epoch)
        # TODO: make use of potentially multi-threaded kepsolve by ensuring `epoch` column is present,
        # and using above passed-in orbit solutions.
        # sol = orbitsolve(orbit, table.epoch[i])
        if orbit_solutions_i_epoch_start >= 0
            sol = orbit_solutions[orbit_solutions_i_epoch_start+i]
            @assert isapprox(table.epoch[i], PlanetOrbits.soltime(sol), rtol=1e-2)
        else
            # TODO: this is a workaround for HGCA likelihoods not having a single table
            # that can be used to precompute orbit solutions
            sol = orbitsolve(orbit, table.epoch[i])
        end

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
        ra_photocentre = (ra_host_vs_bary + ra_planet_vs_bary * flux_ratio) / (1 + flux_ratio)
        Δα_model[i] += ra_photocentre
        dec_host_vs_bary = decoff(sol, planet_mass_msol)
        dec_planet_vs_host = decoff(sol)
        dec_planet_vs_bary = dec_host_vs_bary + dec_planet_vs_host
        dec_photocentre = (dec_host_vs_bary + dec_planet_vs_bary * flux_ratio) / (1 + flux_ratio)
        Δδ_model[i] += dec_photocentre
    end
    return
end



function _query_gaia_dr3(;gaia_id)
    fname = "_gaia_dr3_final/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr3.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false,
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr3_final")
            mkdir("_gaia_dr3_final")
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

# Global variable to track if SPICE kernels are loaded
const _SPICE_KERNELS_LOADED = Ref(false)

"""
    _ensure_spice_kernels_loaded()

Ensure that SPICE kernels are loaded for Earth barycentric position calculations.
Downloads DE440 ephemeris and leap seconds kernel if not already present.
"""
function _ensure_spice_kernels_loaded()
    if _SPICE_KERNELS_LOADED[]
        return
    end
    
    # Get data directory from DataDeps
    data_dir = @datadep_str "DE440_Ephemeris"
    
    # Load leap seconds kernel
    lsk_path = joinpath(data_dir, "naif0012.tls")
    if isfile(lsk_path)
        furnsh(lsk_path)
    else
        error("Leap seconds kernel not found at $lsk_path")
    end
    
    # Load planetary ephemeris kernel
    spk_path = joinpath(data_dir, "de440.bsp")
    if isfile(spk_path)
        furnsh(spk_path)
    else
        error("DE440 ephemeris kernel not found at $spk_path")
    end
    
    _SPICE_KERNELS_LOADED[] = true
    @info "SPICE kernels loaded successfully for Earth barycentric position calculations"
end

"""
    geocentre_position_query(epoch_MJD)

Given a date+time in MJD format, return a named tuple of Earth position and velocity in AU 
on that date. Uses SPICE.jl with JPL DE440 ephemeris data for offline calculations.

The positions and velocities represent the Geocenter of the Earth relative to the solar system
barycenter in the J2000 reference frame.
"""
function geocentre_position_query(epoch_MJD::Number)
    
    # Ensure SPICE kernels are loaded
    _ensure_spice_kernels_loaded()
    
    # Convert MJD to Julian Date
    jd = epoch_MJD + 2400000.5  # MJD to JD conversion
    
    # Convert JD to DateTime
    # JD 2451545.0 = January 1, 2000, 12:00:00 TT (J2000.0 epoch)
    # Use Dates.julian2datetime for conversion
    dt = Dates.julian2datetime(jd)
    
    # Convert to ephemeris time string for SPICE
    et = utc2et(string(dt))
    
    # Get Earth's state relative to Solar System Barycenter
    # 399 = Earth geocenter, 0 = Solar System Barycenter
    state, _ = spkez(399, et, "J2000", "NONE", 0)

    # Extract position and velocity, convert to AU and AU/day
    pos_km = state[1:3]    # Position in km
    vel_km_s = state[4:6]  # Velocity in km/s
    
    # Convert to AU (1 AU = 149,597,870.7 km)
    AU_KM = Octofitter.PlanetOrbits.au2m / 1e3
    pos_au = pos_km ./ AU_KM
    vel_au_day = vel_km_s .* 86400 ./ AU_KM  # km/s to AU/day
    
    return (; x=pos_au[1], y=pos_au[2], z=pos_au[3], 
              vx=vel_au_day[1], vy=vel_au_day[2], vz=vel_au_day[3])
end



"""
    forecast_table = GOST_forecast(ra_deg,dec_deg;baseline=:dr3)

Given an Ra and Dec position, retreive a forecast of Gaia observations from the GOST tool automatically.
See tool URL here: https://gaia.esac.esa.int/gost/

Please be aware that others  might be able to discover the target coordinates you searched for
(though not who performed the search) via information leaked to the external service.

Baseline can be :dr3, :dr4, or :dr5.
"""
function GOST_forecast(ra_deg,dec_deg;baseline=:dr3)
    if baseline == :dr3
        to = "2017-06-28T00:00:00"
    elseif baseline == :dr4
        to = "2020-01-20T00:00:00"
    elseif baseline == :dr5
        to = "2025-01-15T06:16:00"
    end

    if haskey(ENV, "OCTO_GOST_CATALOG") && !isempty(ENV["OCTO_GOST_CATALOG"])
        fname = ENV["OCTO_GOST_CATALOG"]
        @info "Using provided Gaia scan forecast database $fname"
        forecast_table = CSV.read(fname, Table, normalizenames=true)
        # mask = isapprox.(forecast_table.ra_rad_, deg2rad(ra_deg)) .& isapprox.(forecast_table.dec_rad_, deg2rad(dec_deg))
        themin, idx = findmin(hypot.(
            (forecast_table.ra_rad_ .- deg2rad(ra_deg)) .*60 .*60 .*1000 .* cos(dec_deg),
            (forecast_table.dec_rad_ .- deg2rad(dec_deg)).*60 .*60 .*1000
        ))
        if themin > 500
            error("Could not find this target within the provided Gaia scan forecast database file set through OCTO_GOST_CATALOG=$fname Closest target: $themin [mas]")
        end
        ra_rad = forecast_table.ra_rad_[idx]
        dec_rad = forecast_table.dec_rad_[idx]
        mask = isapprox.(forecast_table.ra_rad_, ra_rad) .& isapprox.(forecast_table.dec_rad_, dec_rad)
        @info "Found forecasted visibility windows" windows=count(mask)
        if isempty(mask)
            error("Invalid condition: no visibility windows.")
        end
        return forecast_table[mask,:]
    end

    fname = "GOST-$ra_deg-$dec_deg-$baseline.csv"
    if isfile(fname)
        @info "Using cached Gaia scan forecast $fname"
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
    @info "Contacting the GAIA scan forecast tool GOST: https://gaia.esac.esa.int/gost/"
    resp0 = HTTP.get(
        "https://gaia.esac.esa.int/gost/",
        cookiejar=cookiejar,
    )
    if resp0.status != 200
        println(String(resp0.body))
        error("Could not contact the GAIA scan forecast tool GOST https://gaia.esac.esa.int See above error message.")
    end
    formdata = Dict([
        "serviceCode"=>"1",
        "inputmode"=>"single",
        "srcname"=>"009",
        "srcra" => string(round(ra_deg,digits=7)),
        "srcdec" => string(round(dec_deg,digits=7)),
        "from" => "2014-07-25T10:31:26",
        "to" => to,
    ])


    @info "Retrieving forecasted GAIA scans from GOST: https://gaia.esac.esa.int/gost/"
    resp = HTTP.post(
        "https://gaia.esac.esa.int/gost/GostServlet",
        body=HTTP.Form(formdata),
        cookiejar=cookiejar
    )
    if resp.status != 200 || contains(String(collect(resp.body)),"error")
        println(String(resp.body))
        error("Could not fetch GAIA scan forecast from GOST. See above error message. Do you have an internet connection available?")
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
    if bytesavailable(io) == 0
        error("Empty response from GOST service. Rate limited?")
    end
    forecast_table = CSV.read(io, Table, normalizenames=true)

    return forecast_table
end


