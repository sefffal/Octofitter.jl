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

abstract type GaiaCatalogLikelihood <: AbstractLikelihood end
struct GaiaDivergenceLikelihood_3{TCat3,TCat2,TCat1,TTable,TDist3,TDist2,TDist1} <: GaiaCatalogLikelihood
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
    scanlaw_table=nothing
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
        dist_dr3 = MvNormal(μ,Hermitian(Σ))
    end

    if isnothing(scanlaw_table)
        @warn "No scan law table provided. We will fetch an approximate solution from the GHOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GHOST_forecast(ra_deg,dec_deg))
        forecast_table.epoch = jd2mjd.(forecast_table.var" ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]")
        forecast_table.scanAngle_rad = forecast_table.var" scanAngle[rad]"
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




# function derive_IAD!(gaialike::GaiaDivergenceLikelihood_3, table=gaialike.table; guess=nothing)

#     # Simulate positions an unperturbed star would follow based on the simplified GAIA
#     # 5-parameter, rectilinear model. This is just a starting point.
#     # N_reported_variables = 0
#     if !isnothing(gaialike.dr3)
#         initial_dr = gaialike.dr3
#         intiial_ref_epoch = meta_gaia_DR3.ref_epoch_mjd
#         # N_reported_variables += 19
#     elseif !isnothing(gaialike.dr2)
#         initial_dr = gaialike.dr2
#         intiial_ref_epoch = meta_gaia_DR3.ref_epoch_mjd
#         # N_reported_variables += 19
#     else
#         initial_dr = gaialike.dr1
#         intiial_ref_epoch = meta_gaia_DR3.ref_epoch_mjd
#         # N_reported_variables += 19
#     end
#     initial_α, initial_δ = _simulate_gaia_5_param_model(
#         table,
#         initial_dr.parallax,
#         initial_dr.ra,
#         initial_dr.dec,
#         initial_dr.pmra,
#         initial_dr.pmdec,
#         getrv(initial_dr)*1e3,
#         intiial_ref_epoch
#     )


#     N = length(table.epoch)
#     if isnothing(guess)
#         initial_along_scan_residuals = zeros(N)
#         guess = [
#             initial_along_scan_residuals;
#         ]
#         if !isnothing(gaialike.dr3) 
#             push!(guess, gaialike.dr3.astrometric_excess_noise*sqrt(
#                 gaialike.dr3.visibility_periods_used/
#                     (gaialike.dr3.astrometric_n_good_obs_al + gaialike.dr3.astrometric_n_good_obs_al)
#             ));
#         end
#         if !isnothing(gaialike.dr2) 
#             push!(guess, gaialike.dr2.astrometric_excess_noise*sqrt(
#                 gaialike.dr2.visibility_periods_used/
#                     (gaialike.dr2.astrometric_n_good_obs_al + gaialike.dr2.astrometric_n_good_obs_al)
#             ));
#         end
#         if !isnothing(gaialike.dr1) 
#             push!(guess, gaialike.dr1.astrometric_excess_noise*sqrt(10));# just a guess, they don't report it
#         end
#     end

#     data = (
#         gaialike,table,initial_α,initial_δ,N,
#     )

#     initial_objective, P_dr3_initial, P_dr2_initial = _data_retrieval_objective(guess,data)

#     @show initial_objective

#     # # If the system is under-determined, we add some regularization 
#     # # using the L1 norm.
#     # if length(gaialike.table) > N_reported_variables
#     #     @info "There are more visibility windows ($(length(gaialike.table)) than total reported variables ($N_reported_variables). L1 regularization will be used."
#     #     func = OptimizationFunction(AutoFiniteDiff()) do guess, data
#     #         gaialike = data[1]
#     #         N = data[4]
#     #         objective = _data_retrieval_objective(guess, data)[1]
#     #         # L1-norm
#     #         # regul = 0.0001sum(abs, view(guess,1:N))

#     #         # Minimize the L1 distance between subsequent scans within 30 degrees of each other
#     #         regul = zero(eltype(guess))
#     #         for i in 2:N
#     #             angle_i = gaialike.table.var" scanAngle[rad]"[i]
#     #             for j in i-1:-1:1
#     #                 angle_j = gaialike.table.var" scanAngle[rad]"[j]
#     #                 angular_distance = angle_j - angle_i
#     #                 angular_distance = (angular_distance + pi) % 2pi - pi
#     #                 if abs(angular_distance) < pi/6
#     #                     regul += abs(guess[i] - guess[j])
#     #                     break
#     #                 end
#     #             end
#     #         end

#     #         return -objective + 0.001*regul
#     #     end
#     # else
#     func = OptimizationFunction((guess, data)->-_data_retrieval_objective(guess, data)[1], AutoFiniteDiff())
#     # end
#     prob = OptimizationProblem(func, guess, data)
#     sol = solve(prob,
#         # LBFGS(
#         #     m=4,
#         #     linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
#         #     alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
#         # ),
#         # NelderMead(),
#         # iterations=1_000_000,
#         # x_abstol=NaN,
#         # x_reltol=NaN,
#         # f_abstol=NaN,
#         # f_reltol=NaN,
#         # g_abstol= 5e-5,
#         # g_reltol= 5e-5,
#         # SimulatedAnnealing(),
#         ParticleSwarm(),
#         iterations=500,
#         allow_f_increases=false,
#         show_trace=true,
#         show_every=5,
#     )

#     # Store derived IAD in likelihood object
#     along_scan_residuals_mas = sol.u[1:N]
#     derived_along_scan_resid .= along_scan_residuals_mas
#     derived_initial_α .= initial_α
#     derived_initial_δ .= initial_δ

#     α✱ₐ = along_scan_residuals_mas/60/60/1000 .* gaialike.cosϕ .+ initial_α
#     δₐ = along_scan_residuals_mas/60/60/1000 .* gaialike.sinϕ .+ initial_δ

#     # Generate a line at this point, perpendicular to the scan angle.
#     # The line is defined by a pair of (αₘ, δₘ) points.
#     derived_αₘ = eachrow(@. [-1, 1]'/60/60/1000/cos(initial_δ) * gaialike.sinϕ)
#     for i in eachindex(initial_α)
#         derived_αₘ[i] .+= α✱ₐ[i]
#     end
#     derived_δₘ = eachrow(@. [1, -1]'/60/60/1000 * gaialike.cosϕ)
#     for i in eachindex(initial_δ)
#         derived_δₘ[i] .+= δₐ[i]
#     end
#     gaialike.derived_αₘ .= SVector{2,Float64}.(derived_αₘ)
#     gaialike.derived_δₘ .= SVector{2,Float64}.(derived_δₘ)
    
   
#     minimized_objective, P_dr3, P_dr2, P_dr1 = _data_retrieval_objective(sol.u,data)

#     return (;sol, initial_α, initial_δ, initial_objective, minimized_objective, P_dr3_initial, P_dr2_initial, P_dr3, P_dr2, P_dr1)
# end

# function _data_retrieval_objective(parameters,data)

#     (gaialike,table,initial_α,initial_δ,N) = data

#     T = eltype(parameters)
#     along_scan_residuals_mas = @view parameters[1:N]
#     # along_scan_uncertainty_mas = @view parameters[N+1:2N]
#     # along_scan_uncertainty_mas = fill(parameters[N+1],N)
#     along_scan_uncertainty_mas = zeros(N)
#     i_param = N+1
#     # along_scan_uncertainty_mas_dr3 = fill(parameters[N+1],N)
#     # along_scan_uncertainty_mas_dr2 = fill(parameters[N+2],N)
#     # @show along_scan_uncertainty_mas[1]

#     α✱ₐ = along_scan_residuals_mas/60/60/1000 .* gaialike.table.cosϕ .+ initial_α
#     δₐ = along_scan_residuals_mas/60/60/1000 .* gaialike.table.sinϕ .+ initial_δ

#     # Generate a line at this point, perpendicular to the scan angle.
#     # The line is defined by a pair of (αₘ, δₘ) points.
#     αₘ = eachrow(@. T[-1, 1]'/60/60/1000/cos(initial_δ) * gaialike.table.sinϕ)
#     for i in eachindex(initial_α)
#         αₘ[i] .+= α✱ₐ[i]
#     end
#     δₘ = eachrow(@. T[1, -1]'/60/60/1000 * gaialike.table.cosϕ)
#     for i in eachindex(initial_δ)
#         δₘ[i] .+= δₐ[i]
#     end

#     # We then fit a Gaia 5-parameter rectilinear models to these points
#     # Our metric to optimize against is the KL-divergence between these 
#     # model fits and the published catalog distributions (means, uncertainties,
#     # and covariances).

#     # TODO: should we weight the different releases?
#     objective = zero(T)
#     P_dr3 = P_dr2 = P_dr1 = nothing

#     ## Model 1: DR3
#     if !isnothing(gaialike.dr3)
#         istart = findfirst(>=(meta_gaia_DR3.start_mjd), gaialike.table.epoch)
#         iend = findlast(<(meta_gaia_DR3.stop_mjd), gaialike.table.epoch)
#         if isnothing(istart)
#             istart = 1
#         end
#         if isnothing(iend)
#             iend = length(gaialike.table.epoch)
#         end
#         # if iend-istart+1 != gaialike.dr3.visibility_periods_used
#             # @show iend-istart+1  gaialike.dr3.visibility_periods_used
#             # error("Number of visibility periods don't match for DR3")
#         # end
#         astrometric_excess_noise = parameters[i_param]
#         i_param += 1
#         data = (;
#             αₘ=view(αₘ, istart:iend),
#             δₘ=view(δₘ, istart:iend),
#             table=@view(gaialike.table[istart:iend]),
#             ref_epoch=meta_gaia_DR3.ref_epoch_mjd,
#             along_scan_uncertainty_mas,
#             # along_scan_uncertainty_mas=along_scan_uncertainty_mas_dr3,
#             catalog_rv_m_s=getrv(gaialike.dr3)*1e3,
#             full3D=true,
#             astrometric_excess_noise#=
#                 gaialike.dr3.astrometric_excess_noise*sqrt(
#                     gaialike.dr3.visibility_periods_used/
#                         (gaialike.dr3.astrometric_n_good_obs_al + gaialike.dr3.astrometric_n_good_obs_al)
#                     )=#
#         )
#         # parallax_guess = 10.
#         # ra_guess = mean(α)
#         # dec_guess = mean(δ)
#         # pmra_guess = (maximum(α) - minimum(α))/((gaialike.table.epoch[iend] - gaialike.table.epoch[istart])/julian_year)
#         # pmdec_guess = (maximum(δ) - minimum(δ))/((gaialike.table.epoch[iend] - gaialike.table.epoch[istart])/julian_year)
#         # @show guess = [parallax_guess,ra_guess,dec_guess,pmra_guess,pmdec_guess]
#         guess = [
#             gaialike.dr3.parallax,
#             gaialike.dr3.ra,
#             gaialike.dr3.dec,
#             gaialike.dr3.pmra,
#             gaialike.dr3.pmdec,

#         ]
        
#         # TODO: not great guesses.
#         # Compute the posterior distribution Gaia would have reported from this data
#         P_dr3 = _posterior_dist_from_optimized_gaia_model(guess, data)
#         if isnothing(P_dr3)
#             # Optimization can in some cases fail or fail to report covariances
#             # This should be improved to happen less often.
#             return -Inf,nothing
#         end
#         Q = gaialike.dist_dr3
#         kl = -Distributions.kldivergence(P_dr3,Q)
#         objective += kl
#     end

#     if !isnothing(gaialike.dr2)
#         istart = findfirst(>=(meta_gaia_DR2.start_mjd), gaialike.table.epoch)
#         iend = findlast(<=(meta_gaia_DR2.stop_mjd), gaialike.table.epoch)
#         if isnothing(istart)
#             istart = 1
#         end
#         if isnothing(iend)
#             iend = length(gaialike.table.epoch)
#         end
#         # j = 0
#         # while iend-istart+1 <= gaialike.dr2.visibility_periods_used && j <3
#         #     j += 1
#         #     iend += 1
#         # end
#         astrometric_excess_noise = parameters[i_param]
#         i_param += 1
#         data = (;
#             αₘ=view(αₘ, istart:iend),
#             δₘ=view(δₘ, istart:iend),
#             table=@view(gaialike.table[istart:iend]),
#             ref_epoch=meta_gaia_DR2.ref_epoch_mjd,
#             # along_scan_uncertainty_mas=along_scan_uncertainty_mas_dr2,
#             along_scan_uncertainty_mas,
#             catalog_rv_m_s=getrv(gaialike.dr2)*1e3,
#             full3D=true,
#             astrometric_excess_noise#=
#                 gaialike.dr2.astrometric_excess_noise*sqrt(
#                     gaialike.dr2.visibility_periods_used/
#                         (gaialike.dr2.astrometric_n_good_obs_al + gaialike.dr2.astrometric_n_good_obs_al)
#                     )=#
#         )
#         # parallax_guess = 10.
#         # ra_guess = mean(α)
#         # dec_guess = mean(δ)
#         # pmra_guess = (maximum(α) - minimum(α))/((gaialike.table.epoch[iend] - gaialike.table.epoch[istart])/julian_year)
#         # pmdec_guess = (maximum(δ) - minimum(δ))/((gaialike.table.epoch[iend] - gaialike.table.epoch[istart])/julian_year)
#         # @show guess = [parallax_guess,ra_guess,dec_guess,pmra_guess,pmdec_guess]
#         guess = [
#             gaialike.dr2.parallax,
#             gaialike.dr2.ra,
#             gaialike.dr2.dec,
#             gaialike.dr2.pmra,
#             gaialike.dr2.pmdec,

#         ]
        
#         # TODO: not great guesses.
#         # Compute the posterior distribution Gaia would have reported from this data
#         P_dr2 = _posterior_dist_from_optimized_gaia_model(guess, data)
#         if isnothing(P_dr2)
#             # Optimization can in some cases fail or fail to report covariances
#             # This should be improved to happen less often.
#             return -Inf,nothing
#         end
#         Q = gaialike.dist_dr2
#         kl = -Distributions.kldivergence(P_dr2,Q)
#         objective += kl
#     end

#     if !isnothing(gaialike.dr1)
#         istart = findfirst(>=(meta_gaia_DR1.start_mjd), gaialike.table.epoch)
#         iend = findlast(<=(meta_gaia_DR1.stop_mjd), gaialike.table.epoch)
#         if isnothing(istart)
#             istart = 1
#         end
#         if isnothing(iend)
#             iend = length(gaialike.table.epoch)
#         end
#         # j = 0
#         # while iend-istart+1 <= gaialike.dr1.visibility_periods_used && j <3
#         #     j += 1
#         #     iend += 1
#         # end
#         astrometric_excess_noise = parameters[i_param]
#         i_param += 1
#         data = (;
#             αₘ=view(αₘ, istart:iend),
#             δₘ=view(δₘ, istart:iend),
#             table=@view(gaialike.table[istart:iend]),
#             ref_epoch=meta_gaia_DR1.ref_epoch_mjd,
#             # along_scan_uncertainty_mas=along_scan_uncertainty_mas_dr1,
#             along_scan_uncertainty_mas,
#             catalog_rv_m_s=0,
#             full3D=true,
#             astrometric_excess_noise#=
#                 gaialike.dr1.astrometric_excess_noise*sqrt(
#                     gaialike.dr1.visibility_periods_used/
#                         (gaialike.dr1.astrometric_n_good_obs_al + gaialike.dr1.astrometric_n_good_obs_al)
#                     )=#
#         )
#         # parallax_guess = 10.
#         # ra_guess = mean(α)
#         # dec_guess = mean(δ)
#         # pmra_guess = (maximum(α) - minimum(α))/((gaialike.table.epoch[iend] - gaialike.table.epoch[istart])/julian_year)
#         # pmdec_guess = (maximum(δ) - minimum(δ))/((gaialike.table.epoch[iend] - gaialike.table.epoch[istart])/julian_year)
#         # @show guess = [parallax_guess,ra_guess,dec_guess,pmra_guess,pmdec_guess]
#         guess = [
#             gaialike.dr1.parallax,
#             gaialike.dr1.ra,
#             gaialike.dr1.dec,
#             gaialike.dr1.pmra,
#             gaialike.dr1.pmdec,

#         ]
        
#         # TODO: not great guesses.
#         # Compute the posterior distribution Gaia would have reported from this data
#         P_dr1 = _posterior_dist_from_optimized_gaia_model(guess, data)
#         if isnothing(P_dr1)
#             # Optimization can in some cases fail or fail to report covariances
#             # This should be improved to happen less often.
#             return -Inf,nothing
#         end
#         Q = gaialike.dist_dr1
#         kl = -Distributions.kldivergence(P_dr1,Q)
#         objective += kl
#     end

#     return (objective,P_dr3,P_dr2,P_dr1)

# end







# function ln_like(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, num_epochs::Val{L}=Val(length(gaialike.table))) where L
#     T = _system_number_type(θ_system)
#     ll = zero(T)

#     along_scan_residuals_mas = θ_system.along_scan_residuals_mas
#     if length(along_scan_residuals_mas) != length(gaialike.table.epochs)
#         error("The number of along-scan residuals must match the number of visiblity windows in the model (`length(gaialike.table.epoch)==$(length(gaialike.table.epoch))`")
#     end

#     # We have two components to our likelihood:
#     # P(orbit + motion | along-scan "measurements") and P(DR2,3 results | along-scan "measurements").

#     # The former will be _orbit_loglike(orbit params, along-scan "measurements", data=(gaia epochs, scan angles) 
#     # The latter will be _gaia_skypath_loglikelihood(along-scan "measurements", data=(gaia epochs, scan angles, P_dr2, P_dr3) 

#     # Here we provide the Gaia catalog parameters
#     args = (
#         θ_system.plx,
#         θ_system.ra,
#         θ_system.dec,
#         θ_system.pmra,
#         θ_system.pmdec,
#     )
#     # _gaia_skypath_loglikelihood(along-scan "measurements", data=(gaia epochs, scan angles, P_dr2, P_dr3) 


#     (
#         plx,
#         ra,
#         dec,
#         pmra,
#         pmdec,
#     ) = args

#     _gaia_skypath_loglikelihood(args,(;
#     αₘ,δₘ,
#     table,
#     ref_epoch,
#     along_scan_uncertainty_mas,
#     full3D,
#     catalog_rv_m_s,
#     astrometric_excess_noise
# ))

#     # We will also finally need support for vector-valued distributions in Octofitter. Sigh.
#     # Done.

    

#     return ll
# end

# function _orbit_loglike()


    # along_scan_uncertainty_deg = θ_system.along_scan_uncertainty_mas/60/60/1000

    # modelled_astrom = simulate(gaialike, θ_system, orbits)
    # for i in eachindex(modelled_astrom.resid)
    #     ll+= logpdf(Normal(0, along_scan_uncertainty_deg),modelled_astrom.resid[i])
    # end

# end


# function simulate(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, L=nothing)

#     T = _system_number_type(θ_system)
#     α✱_model_with_perturbation_out = zeros(T,length(gaialike.table.epoch))
#     δ_model_with_perturbation_out = zeros(T,length(gaialike.table.epoch))
#     resid_out = zeros(T,length(gaialike.table.epoch))

#     # All planets inthe system have orbits defined with the same ra, dec, and proper motion,
#     # since these are properties of the system.
#     orbit = first(orbits)
#     if length(orbits) > 1
#         for i in eachindex(orbits)[2:end]
#             if orbits[i].ra != orbit.ra ||
#                orbits[i].dec != orbit.dec ||
#                orbits[i].pmra != orbit.rpma ||
#                orbits[i].pmdec != orbit.ra pmdec
#                 error("Planet orbits do not have matching ra, dec, pmpra, and pmdec.")
#             end
#         end
#     end


#     for i in eachindex(gaialike.table.epoch)
#         sol_gaia_epoch = orbitsolve(orbit, gaialike.table.epoch[i])
#         cmp = sol_gaia_epoch.compensated
#         # GAIA considers a simplified rectilinear sky path model 
#         # but does include changing parallax from RV
#         # full3D
#         x_earth_pc = gaialike.table.x[i] / PlanetOrbits.pc2au
#         y_earth_pc = gaialike.table.y[i] / PlanetOrbits.pc2au
#         z_earth_pc = gaialike.table.z[i] / PlanetOrbits.pc2au
#         x_diff_pc = cmp.x₂ - x_earth_pc
#         y_diff_pc = cmp.y₂ - y_earth_pc
#         z_diff_pc = cmp.z₂ - z_earth_pc
#         distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
#         mydtor = π / 180
#         ra_apparent_deg = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
#         arg = z_diff_pc / distance_diff
#         arg = map(arg) do arg
#             if 1.0 < arg < 1.0 + sqrt(eps(1.0))
#                 arg = 1.0
#             end
#             return arg
#         end
#         dec_apparent_deg = asin(arg) / mydtor

#         α✱_perturbation = zero(T)
#         δ_perturbation = zero(T)

#         # Add perturbations from all planets
#         for planet_i in eachindex(orbits)
#             orbit = orbits[planet_i]
#             if orbit == orbit
#                 sol = sol_gaia_epoch
#             else
#                 sol = orbitsolve(orbit, gaialike.table.epoch[i])
#             end
#             # Add perturbation from planet
#             α✱_perturbation += raoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
#             δ_perturbation += decoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
#         end

#         # coordinates in degrees of RA and Dec
#         α✱_model_with_perturbation = ra_apparent_deg + α✱_perturbation/60/60/1000/cosd(gaialike.table.derived_initial_δ[i])
#         δ_model_with_perturbation = dec_apparent_deg + δ_perturbation/60/60/1000

#         # Calculate differences in milliarcseconds by mapping into a local tangent plane,
#         # with the local coordinate defined by the Hipparcos solution at this *epoch* (not 
#         # at the best solution)
#         point = @SVector [
#             α✱_model_with_perturbation,
#             δ_model_with_perturbation,
#         ]

#         # Two points defining a line along which the star's position was measured
#         line_point_1 = @SVector [
#             gaialike.table.derived_αₘ[i][1],
#             gaialike.table.derived_δₘ[i][1],
#         ]
#         line_point_2 = @SVector [
#             gaialike.table.derived_αₘ[i][2],
#             gaialike.table.derived_δₘ[i][2],
#         ]
#         # Distance from model star Delta position to this line where it was measured 
#         # by the satellite
#         resid = distance_point_to_line(point, line_point_1, line_point_2) # degrees

#         α✱_model_with_perturbation_out[i] = α✱_model_with_perturbation
#         δ_model_with_perturbation_out[i]  = δ_model_with_perturbation
#         resid_out[i] = resid
#     end

#     return (;
#         α✱_model_with_perturbation=α✱_model_with_perturbation_out,
#         δ_model_with_perturbation=δ_model_with_perturbation_out,
#         resid=resid_out
#     )
# end



















function ln_like(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, num_epochs::Val{L}=Val(length(gaialike.table))) where L
    ll, _ = simulate(gaialike, θ_system, orbits, num_epochs)
    return ll
end


function simulate(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, num_epochs::Val{L}=Val(length(gaialike.table))) where L

    T = _system_number_type(θ_system)

    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol
    # typical_along_scan_uncertainty_mas = θ_system.typical_along_scan_uncertainty_mas
    # along_scan_uncertainty_mas = sqrt(typical_along_scan_uncertainty_mas^2 + gaialike.dr3.astrometric_excess_noise)
    # along_scan_uncertainty_mas = θ_system.σ_scan

    sample_orbit = only(orbits)


    # TODO: Here is where I would, instead of simulating the skypath, just grab the x data from the sampled parameters
    # TODO: The guesses are not great at present, and I think I had a better code for this.


    # Generate simulated observations from this sample draw
    # These are generated using whatever reference epoch the user has specified that corresponds
    # to their ra, dec, plx, etc variables
    # (;α_model,δ_model,αₘ,δₘ) = _simulate_skypath_observations(gaialike, sample_orbit::AbsoluteVisual, planet_mass_msol, T)

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
        sol = orbitsolve(sample_orbit, data.ref_epoch)
        guess = [
            sol.compensated.parallax2,
            sol.compensated.ra2,
            sol.compensated.dec2,
            sample_orbit.pmra,
            sample_orbit.pmdec,
        ]
        modelled_gaia_parameters_dr3 = _find_maxlike_gaia_model_soln(guess,data)
        ll += logpdf(gaialike.dist_dr3, modelled_gaia_parameters_dr3)
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
        sol = orbitsolve(sample_orbit, data.ref_epoch)
        guess = [
            sol.compensated.parallax2,
            sol.compensated.ra2,
            sol.compensated.dec2,
            sample_orbit.pmra,
            sample_orbit.pmdec,
        ]
        modelled_gaia_parameters_dr2 = _find_maxlike_gaia_model_soln(guess,data)
        ll += logpdf(gaialike.dist_dr2, modelled_gaia_parameters_dr2)

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
Given scan epochs and angles, and an AbsoluteVisual orbit describing a sky path and perturbations
from a planet, calculate a "true" model of what Gaia would have observed, including
full spherical coordiante effects, radial velocity (perspective acceleration), etc.
TODO: include effects of changing light travel time.
"""
function _simulate_skypath_observations(gaialike, sample_orbit::AbsoluteVisual, planet_mass_msol, T=Float64)
    # Compute the Ra and Dec sky path of the star, accounting for 3D barycentric
    # motion in spherical coordinates and perturbations from planets
    α_model = zeros(T, size(gaialike.table,1))
    δ_model = zeros(T, size(gaialike.table,1))
    for i in eachindex(gaialike.table.epoch)
        sol = orbitsolve(sample_orbit, gaialike.table.epoch[i])
        cmp = sol.compensated

        # Calculate the position of the star in cartesian coordinates
        # Get earth position in cartesian coordaintes (in units of AU)
        # Calculate apparent Ra and Dec.
        x_earth_pc = gaialike.table.x[i] / PlanetOrbits.pc2au
        y_earth_pc = gaialike.table.y[i] / PlanetOrbits.pc2au
        z_earth_pc = gaialike.table.z[i] / PlanetOrbits.pc2au
        x_diff_pc = cmp.x₂ - x_earth_pc
        y_diff_pc = cmp.y₂ - y_earth_pc
        z_diff_pc = cmp.z₂ - z_earth_pc
        distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
        mydtor = π / 180
        ra_apparent_deg = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
        arg = z_diff_pc / distance_diff
        arg = map(arg) do arg
            if 1.0 < arg < 1.0 + sqrt(eps(1.0))
                arg = 1.0
            end
            return arg
        end
        dec_apparent_deg = asin(arg) / mydtor

        # TODO: add perturbations from multiple planets here
        for orb in (sample_orbit,)
            sol = orbitsolve(orb, gaialike.table.epoch[i])
            # Add perturbation from planet
            ra_apparent_deg += raoff(sol, planet_mass_msol)/60/60/1000/cos(sample_orbit.dec)
            dec_apparent_deg += decoff(sol, planet_mass_msol)/60/60/1000
        end

        α_model[i] = ra_apparent_deg
        δ_model[i] = dec_apparent_deg
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
    # astrometric_excess_noise
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



"""
Fit the Gaia 5-parameter model to some data, and report the best fit means and covariances
as a multivariate normal distribution.
Uncertainties are calculated using the inverse hessian of the likelihood at the best-fit location.
"""
function _find_maxlike_gaia_model_soln(guess, data)
    func = OptimizationFunction((args...)->-_gaia_skypath_loglikelihood(args...), AutoForwardDiff())
    prob = OptimizationProblem(func, guess, data)
    sol = solve(prob,
        LBFGS(
            m=4,
            # linesearch=OptimizationOptimJL.Optim.LineSearches.HagerZhang(linesearchmax=50),
            linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
            alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
        ),
        # NewtonTrustRegion(), # Is a good choice too, but about 2x as slow.
        # x_abstol=NaN,
        # x_reltol=NaN,
        # f_abstol=NaN,
        # f_reltol=NaN,
        g_abstol= 1e-5,
        g_reltol= 1e-5,
        # g_tol=1e-15, f_tol=NaN,  x_tol=NaN, x_reltol=NaN, 
        allow_f_increases=true,
        iterations=10000
    )
    return sol.u
end

"""
Fit the Gaia 5-parameter model to some data, and report the best fit means and covariances
as a multivariate normal distribution.
Uncertainties are calculated using the inverse hessian of the likelihood at the best-fit location.
"""
function _posterior_dist_from_optimized_gaia_model(guess, data)
    func = OptimizationFunction((args...)->-_gaia_skypath_loglikelihood(args...), AutoForwardDiff())
    prob = OptimizationProblem(func, guess, data)
    sol = solve(prob,
        LBFGS(
            m=4,
            # linesearch=OptimizationOptimJL.Optim.LineSearches.HagerZhang(linesearchmax=50),
            linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
            alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
        ),
        # NewtonTrustRegion(), # Is a good choice too, but about 2x as slow.
        # x_abstol=NaN,
        # x_reltol=NaN,
        # f_abstol=NaN,
        # f_reltol=NaN,
        g_abstol= 1e-5,
        g_reltol= 1e-5,
        # g_tol=1e-15, f_tol=NaN,  x_tol=NaN, x_reltol=NaN, 
        allow_f_increases=true,
        iterations=10000
    )
    if sol.retcode != Optimization.ReturnCode.Success
        return nothing
    end
    # Compute uncertainties and covariances using inverse Hessian
    H = ForwardDiff.hessian(u->_gaia_skypath_loglikelihood(u,data), sol.u)
    local gaia_post_covariance
    try
        gaia_post_covariance = -inv(H)
    catch 
        return nothing
    end
    # We use the inverse Hessian from the gradient descent as the covariance matrix.
    # In case of numerical error (non-hermitian matrix) we regularize it by adding 
    # a value along the diagonals. 
    # P = MvNormal(sol.u, Hermitian(gaia_post_covariance))
    P = nothing
    for e in (0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6)#, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3, 1e4)
        try
            P = MvNormal(sol.u, Hermitian(gaia_post_covariance))
            return P
        catch
            for i in 1:length(sol.u)
                gaia_post_covariance[i] += e
            end
        end
    end
    return nothing
end


# Given catalog values, simulate
function _simulate_gaia_5_param_model(table, plx, ra, dec, pmra, pmdec, rv_m_s, ref_epoch)
    
    ra_apparent_deg = zeros(length(table.epoch))
    dec_apparent_deg = zeros(length(table.epoch))

    plx0 = plx
    dist0 = 1000/plx0
    δdist_pc_δt_sec = rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
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
        ra_apparent_deg[i] = ra + ra_apparent_Δmas/(60*60*1000)/cosd(dec)
        dec_apparent_deg[i] = dec + dec_apparent_Δmas/(60*60*1000)
    end

    return ra_apparent_deg, dec_apparent_deg
end


# For testing, unless I can figure out how to pull this out.

# function simulate(gaialike::GaiaDivergenceLikelihood_3, θ_system, orbits, num_epochs::Val{L}=Val(length(gaialike.table))) where L

#     T = _system_number_type(θ_system)

#     # TODO: expand to work with multiple, named planets
#     planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol
#     # typical_along_scan_uncertainty_mas = θ_system.typical_along_scan_uncertainty_mas
#     # along_scan_uncertainty_mas = sqrt(typical_along_scan_uncertainty_mas^2 + gaialike.dr3.astrometric_excess_noise)
#     along_scan_uncertainty_mas = θ_system.along_scan_uncertainty_mas

#     sample_orbit = only(orbits)

#     # Generate simulated observations from this sample draw
#     # These are generated using whatever reference epoch the user has specified that corresponds
#     # to their ra, dec, plx, etc variables
#     (;α_model,δ_model,αₘ,δₘ) = _simulate_skypath_observations(gaialike, sample_orbit::AbsoluteVisual, planet_mass_msol, T)

#     ll = zero(T)

#     # Now we fit a no-planet (zero mass planet) sky path model to this data.
#     # These should be fit using the appropriate catalog reference epoch so 
#     # that they can be compared correctly.

#     ## Model 1: DR3
#     if !isnothing(gaialike.dr3)
#         istart = findfirst(>=(meta_gaia_DR3.start_mjd), gaialike.table.epoch)
#         iend = findlast(<=(meta_gaia_DR3.start_mjd), gaialike.table.epoch)
#         if isnothing(istart)
#             istart = 1
#         end
#         if isnothing(iend)
#             iend = length(gaialike.table.epoch)
#         end
#         data = (;
#             αₘ=αₘ[istart:iend],
#             δₘ=δₘ[istart:iend],
#             table=gaialike.table[istart:iend],
#             ref_epoch=meta_gaia_DR3.ref_epoch_mjd,
#             along_scan_uncertainty_mas,
#             catalog_rv_m_s=getrv(gaialike.dr3)*1e3,
#             full3D=false,
#         )
#         sol = orbitsolve(sample_orbit, data.ref_epoch)
#         guess = [
#             sol.compensated.parallax2,
#             sol.compensated.ra2,
#             sol.compensated.dec2,
#             sample_orbit.pmra,
#             sample_orbit.pmdec,
#         ]
#         # Compute the posterior distribution Gaia would have reported from this data
#         P = _posterior_dist_from_optimized_gaia_model(guess, data)
#         return P
#     end
# end

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
        error("Could not find source")
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
        forecast_table = CSV.read(fname, Table)
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
    forecast_table = CSV.read(io, Table)

    return forecast_table
end




