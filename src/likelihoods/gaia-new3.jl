
"""
This likelihood object models simple 5-parameter skypath / parallactic motion compared
to a set of along-scan measurements, eg. from Gaia.

Either the model parameters (plx,ra,dec,pmra,pmdec) can vary, or the along-scan measurements
can. In the former case, you infer the model from the data. In the later case, you infer the 
underlying data from reported model parameters.
""" 
struct ParallacticMotionLikelihood_v7{TTable,catalog_parameters,along_scan_residuals,σ_scan,excess_noise} <: AbstractLikelihood
    table::TTable
    ref_epoch::Float64
    catalog_rv_m_s::Float64
end
function ParallacticMotionLikelihood_v7(
    table;
    catalog_parameters,along_scan_residuals,σ_scan,excess_noise,ref_epoch,catalog_rv_m_s=0
)
    table = Table(table)
    return ParallacticMotionLikelihood_v7{typeof(table),catalog_parameters,along_scan_residuals,σ_scan,excess_noise}(table,ref_epoch,catalog_rv_m_s)
end

function ParallacticMotionLikelihood_DR2(;
    source_id_dr2,scanlaw_table=nothing,
    catalog_parameters,along_scan_residuals,σ_scan,excess_noise,ref_epoch,
    initial_α=nothing,
    initial_δ=nothing,
)
    # Query Gaia archive for DR3 solution
    dr2 = Octofitter._query_gaia_dr2(;gaia_id=source_id_dr2)
    ra_deg = dr2.ra
    dec_deg = dr2.dec
    if isnothing(scanlaw_table)
        @warn "No scan law table provided. We will fetch an approximate solution from the GHOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GHOST_forecast(ra_deg,dec_deg))
        forecast_table.epoch = jd2mjd.(forecast_table.var" ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]")
        forecast_table.scanAngle_rad = forecast_table.var" scanAngle[rad]"
    else
        @info "Scanlaw table from the `scanninglaw` python package was provided, will not use GHOST."
        scanlaw_table.epochs = tcb_at_gaia_2mjd.(scanlaw_table.times)
        scanlaw_table.epochs = tcb_at_gaia_2mjd.(scanlaw_table.times)
        forecast_table.scanAngle_rad = deg2rad.(scanlaw_table.angles)
    end


    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.scanAngle_rad)

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)

    catalog_rv_m_s = getrv(dr2)

    if !isnothing(initial_α) && !isnothing(initial_δ)
        # Set the initial skypath that the residuals will be compared against
        initial_α, initial_δ = _simulate_gaia_5_param_model(
            table,
            dr2.parallax,
            dr2.ra,
            dr2.dec,
            dr2.pmra,
            dr2.pmdec,
            catalog_rv_m_s*1e3,
            meta_gaia_DR2.ref_epoch_mjd
        )
    end
    table.initial_α = initial_α
    table.initial_δ = initial_δ

    table = table[meta_gaia_DR2.start_mjd .<= table.epoch .<= meta_gaia_DR2.stop_mjd,:]

    return ParallacticMotionLikelihood_v7(table;catalog_parameters,along_scan_residuals,σ_scan,excess_noise,ref_epoch,catalog_rv_m_s)
end


function ParallacticMotionLikelihood_DR3(;
    source_id_dr3,scanlaw_table=nothing,
    catalog_parameters,along_scan_residuals,σ_scan,excess_noise,ref_epoch,
)
    # Query Gaia archive for DR3 solution
    dr3 = Octofitter._query_gaia_dr3(;gaia_id=source_id_dr3)
    ra_deg = dr3.ra
    dec_deg = dr3.dec
    if isnothing(scanlaw_table)
        @warn "No scan law table provided. We will fetch an approximate solution from the GHOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GHOST_forecast(ra_deg,dec_deg))
        forecast_table.epoch = jd2mjd.(forecast_table.var" ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]")
        forecast_table.scanAngle_rad = forecast_table.var" scanAngle[rad]"
    else
        @info "Scanlaw table from the `scanninglaw` python package was provided, will not use GHOST."
        scanlaw_table.epochs = tcb_at_gaia_2mjd.(scanlaw_table.times)
        scanlaw_table.epochs = tcb_at_gaia_2mjd.(scanlaw_table.times)
        forecast_table.scanAngle_rad = deg2rad.(scanlaw_table.angles)
    end


    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.scanAngle_rad)

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)

    catalog_rv_m_s = getrv(dr3)

    # Set the initial skypath that the residuals will be compared against
    initial_α, initial_δ = _simulate_gaia_5_param_model(
        table,
        dr3.parallax,
        dr3.ra,
        dr3.dec,
        dr3.pmra,
        dr3.pmdec,
        catalog_rv_m_s*1e3,
        meta_gaia_DR3.ref_epoch_mjd
    )
    table.initial_α = initial_α
    table.initial_δ = initial_δ

    table = table[meta_gaia_DR3.start_mjd .<= table.epoch .<= meta_gaia_DR3.stop_mjd,:]

    return ParallacticMotionLikelihood_v7(table;catalog_parameters,along_scan_residuals,σ_scan,excess_noise,ref_epoch,catalog_rv_m_s)
end


function _getparams(skypathlike::ParallacticMotionLikelihood_v7{TTable,catalog_parameters,along_scan_residuals,σ_scan,excess_noise}, θ_system) where
    {TTable,catalog_parameters,along_scan_residuals,σ_scan,excess_noise}
    P = getproperty(θ_system, catalog_parameters)
    (
        plx,
        ra,# deg
        dec,# deg
        pmra, 
        pmdec,
    ) = P
    val_along_scan_residuals = getproperty(θ_system, along_scan_residuals)
    val_σ_scan = getproperty(θ_system, σ_scan)
    val_excess_noise = getproperty(θ_system, excess_noise)
    return (;plx, ra, dec, pmra, pmdec, along_scan_residuals=val_along_scan_residuals, σ_scan=val_σ_scan,excess_noise=val_excess_noise, skypathlike.ref_epoch)
end
function ln_like(skypathlike::ParallacticMotionLikelihood_v7, θ_system, orbits, num_epochs::Val{L}=Val(length(skypathlike.table))) where L
    T = _system_number_type(θ_system)
    ll = zero(T)
    (;
        plx, ra, dec, pmra, pmdec,
        along_scan_residuals,
        σ_scan,
        excess_noise,
        ref_epoch,
    ) = _getparams(skypathlike, θ_system)

    # @show plx ra dec pmra pmdec

    # TODO: use along_scan_residuals to define αₘ δₘ

    α✱ₐ = along_scan_residuals ./ 60 ./60 ./1000 .* skypathlike.table.cosϕ .+ skypathlike.table.initial_α
    δₐ = along_scan_residuals ./ 60 ./60 ./1000 .* skypathlike.table.sinϕ .+ skypathlike.table.initial_δ

    # Generate a line at this point, perpendicular to the scan angle.
    # The line is defined by a pair of (αₘ, δₘ) points.
    αₘ = eachrow(@. T[-1, 1]'/60/60/1000/cos(skypathlike.table.initial_δ) * skypathlike.table.sinϕ)
    for i in eachindex(skypathlike.table.initial_α)
        αₘ[i] .+= α✱ₐ[i]
    end
    δₘ = eachrow(@. T[1, -1]'/60/60/1000 * skypathlike.table.cosϕ)
    for i in eachindex(skypathlike.table.initial_δ)
        δₘ[i] .+= δₐ[i]
    end

    plx0 = plx
    dist0 = 1000/plx0
    δdist_pc_δt_sec = skypathlike.catalog_rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
    σ_tot = sqrt(σ_scan^2 + excess_noise^2)

    # ra_apparent_deg = zeros(length(skypathlike.table.epoch))
    # dec_apparent_deg = zeros(length(skypathlike.table.epoch))
    for i in eachindex(skypathlike.table.epoch)
        Δdist_pc = δdist_pc_δt_sec * (skypathlike.table.epoch[i] - ref_epoch)
        dist1 = dist0 + Δdist_pc
        plx_at_time = 1000 / dist1
        # Calculate sky coordinate Delta at each scan epoch from the catalog position
        # using the eath's motion ("ephemeris"), x,y,z in AU.
        # TODO: should I rather be using the catalog Ra and Dec in here? Unclear
        ra_apparent_Δmas = plx_at_time * (
            skypathlike.table.x[i] * sind(ra) -
            skypathlike.table.y[i] * cosd(ra)
        ) + (skypathlike.table.epoch[i] - ref_epoch)/julian_year * pmra
        dec_apparent_Δmas = plx_at_time * (
            skypathlike.table.x[i] * cosd(ra) * sind(dec) +
            skypathlike.table.y[i] * sind(ra) * sind(dec) -
            skypathlike.table.z[i] * cosd(dec)
        ) + (skypathlike.table.epoch[i] - ref_epoch)/julian_year * pmdec
        ra_apparent_deg = ra + ra_apparent_Δmas/(60*60*1000)/cosd(dec)
        dec_apparent_deg = dec + dec_apparent_Δmas/(60*60*1000)
        point = @SVector [ra_apparent_deg, dec_apparent_deg] 
        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
        line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
        # TODO: distance point to line should account for spherical coordinates!
        resid = distance_point_to_line(point, line_point_1, line_point_2) # mas
        # @show resid σ_tot/60/60/1000
        ll += logpdf(Normal(0,σ_tot/60/60/1000), resid)
    end

    # f,a,p = Main.scatterlines(ra_apparent_deg,dec_apparent_deg)
    # Main.scatterlines!(a,skypathlike.table.initial_α[:],skypathlike.table.initial_δ[:])
    # display(
    #     f
    # )

    # @show "a"

    return ll
end







struct StarAbsoluteAstrometryLikelihood_v2{TTable,along_scan_residuals,σ_scan} <: AbstractLikelihood
    table::TTable
end
function StarAbsoluteAstrometryLikelihood_v2(
    table;
    centroid_pos_al,
    σ_scan,
)
    table = Table(table)
    return StarAbsoluteAstrometryLikelihood_v2{typeof(table),centroid_pos_al,σ_scan}(table)
end

function _getparams(::StarAbsoluteAstrometryLikelihood_v2{TTable,centroid_pos_al,σ_scan}, θ_system) where
    {TTable,centroid_pos_al,σ_scan}
    val_centroid_pos_al = getproperty(θ_system, centroid_pos_al)
    val_σ_scan = getproperty(θ_system, σ_scan)
    return (;centroid_pos_al=val_centroid_pos_al, σ_scan=val_σ_scan)
end
function ln_like(
    absastromlike::StarAbsoluteAstrometryLikelihood_v2,
    θ_system,
    orbits,
    num_epochs::Val{L}=Val(length(absastromlike.table.epochs))
) where L
    T = _system_number_type(θ_system)
    ll = zero(T)
    (;
        centroid_pos_al,
        σ_scan,
    ) = _getparams(absastromlike, θ_system)
    planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol
    
    # All planets in the system have orbits defined with the same ra, dec, and proper motion,
    # since these are properties of the system.
    orbit = first(orbits)
    if length(orbits) > 1
        for i in eachindex(orbits)[2:end]
            if orbits[i].ra != orbit.ra ||
               orbits[i].dec != orbit.dec ||
               orbits[i].pmra != orbit.rpma ||
               orbits[i].pmdec != orbit.ra pmdec
                error("Planet orbits do not have matching ra, dec, pmpra, and pmdec.")
            end
        end
    end

    for i in eachindex(absastromlike.table.epoch)
        sol = orbitsolve(orbit, absastromlike.table.epoch[i])
        cmp = sol.compensated

        # Calculate the position of the star in cartesian coordinates
        # Get earth position in cartesian coordaintes (in units of AU)
        # Calculate apparent Ra and Dec.
        x_earth_pc = absastromlike.table.x[i] / PlanetOrbits.pc2au
        y_earth_pc = absastromlike.table.y[i] / PlanetOrbits.pc2au
        z_earth_pc = absastromlike.table.z[i] / PlanetOrbits.pc2au
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
        for orb in orbits
            sol = orbitsolve(orb, absastromlike.table.epoch[i])
            # Add perturbation from planet
            ra_apparent_deg += raoff(sol, planet_mass_msol)/60/60/1000/cos(orb.dec)
            dec_apparent_deg += decoff(sol, planet_mass_msol)/60/60/1000
        end

        point = @SVector [
            ra_apparent_deg ,
            dec_apparent_deg,
        ]

        α✱ₐ = along_scan_residuals[i] / 60 /60 /1000 .* absastromlike.table.cosϕ[i] + absastromlike.table.initial_α[i]
        δₐ = along_scan_residuals[i] / 60 /60 /1000 .* absastromlike.table.sinϕ[i] + absastromlike.table.initial_δ[i]

        # Generate a line at this point, perpendicular to the scan angle.
        # The line is defined by a pair of (αₘ, δₘ) points.
        αₘ = @. (-1, 1) / 60 /60/1000/cos(absastromlike.table.initial_δ[i]) * absastromlike.table.sinϕ[i] + α✱ₐ
        δₘ = @. (1, -1) / 60 /60/1000 * absastromlike.table.cosϕ[i] +  δₐ

        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [αₘ[1], δₘ[1]]
        line_point_2 = @SVector [αₘ[2], δₘ[2]]

        resid = distance_point_to_line(point, line_point_1, line_point_2) # mas
        ll += logpdf(Normal(0,σ_scan/60/60/1000), resid)
    end

    # So we can add the model ra and dec to each of our model outputs
    # And then compare to the data outputs + the 1991.25 position

    return ll
end


# function ln_like(absastromlike::StarAbsoluteAstrometryLikelihood_v2, θ_system, orbits, num_epochs::Val{L}=Val(length(absastromlike.table))) where L
#     T = _system_number_type(θ_system)
#     ll = zero(T)
#     (;
#         along_scan_residuals,
#         σ_scan,
#     ) = _getparams(absastromlike, θ_system)

#     α✱ₐ = along_scan_residuals ./ 60 ./60 ./1000 .* absastromlike.table.cosϕ .+ absastromlike.table.initial_α
#     δₐ = along_scan_residuals ./ 60 ./60 ./1000 .* absastromlike.table.sinϕ .+ absastromlike.table.initial_δ

#     # Generate a line at this point, perpendicular to the scan angle.
#     # The line is defined by a pair of (αₘ, δₘ) points.
#     αₘ = eachrow(@. T[-1, 1]'/60/60/1000/cos(absastromlike.table.initial_δ) * absastromlike.table.sinϕ)
#     for i in eachindex(absastromlike.table.initial_α)
#         αₘ[i] .+= α✱ₐ[i]
#     end
#     δₘ = eachrow(@. T[1, -1]'/60/60/1000 * absastromlike.table.cosϕ)
#     for i in eachindex(absastromlike.table.initial_δ)
#         δₘ[i] .+= δₐ[i]
#     end

#     plx0 = plx
#     dist0 = 1000/plx0
#     δdist_pc_δt_sec = absastromlike.catalog_rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
#     σ_tot = sqrt(σ_scan^2 + excess_noise^2)

#     # ra_apparent_deg = zeros(length(absastromlike.table.epoch))
#     # dec_apparent_deg = zeros(length(absastromlike.table.epoch))
#     for i in eachindex(absastromlike.table.epoch)
#         Δdist_pc = δdist_pc_δt_sec * (absastromlike.table.epoch[i] - ref_epoch)
#         dist1 = dist0 + Δdist_pc
#         plx_at_time = 1000 / dist1
#         # Calculate sky coordinate Delta at each scan epoch from the catalog position
#         # using the eath's motion ("ephemeris"), x,y,z in AU.
#         # TODO: should I rather be using the catalog Ra and Dec in here? Unclear
#         ra_apparent_Δmas = plx_at_time * (
#             absastromlike.table.x[i] * sind(ra) -
#             absastromlike.table.y[i] * cosd(ra)
#         ) + (absastromlike.table.epoch[i] - ref_epoch)/julian_year * pmra
#         dec_apparent_Δmas = plx_at_time * (
#             absastromlike.table.x[i] * cosd(ra) * sind(dec) +
#             absastromlike.table.y[i] * sind(ra) * sind(dec) -
#             absastromlike.table.z[i] * cosd(dec)
#         ) + (skypathlike.table.epoch[i] - ref_epoch)/julian_year * pmdec
#         ra_apparent_deg = ra + ra_apparent_Δmas/(60*60*1000)/cosd(dec)
#         dec_apparent_deg = dec + dec_apparent_Δmas/(60*60*1000)
#         point = @SVector [ra_apparent_deg, dec_apparent_deg] 
#         # Two points defining a line along which the star's position was measured
#         line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
#         line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
#         # TODO: distance point to line should account for spherical coordinates!
#         resid = distance_point_to_line(point, line_point_1, line_point_2) # mas
#         # @show resid σ_tot/60/60/1000
#         ll += logpdf(Normal(0,σ_tot/60/60/1000), resid)
#     end

#     # f,a,p = Main.scatterlines(ra_apparent_deg,dec_apparent_deg)
#     # Main.scatterlines!(a,skypathlike.table.initial_α[:],skypathlike.table.initial_δ[:])
#     # display(
#     #     f
#     # )

#     # @show "a"

#     return ll
# end

















# """
# Internal nested likelihood model. 
# This is our best reconstruction of the likelihood function the Gaia team would have tried
# to fit (for the 5-parameter solution). Notably it includes no effects of secondary perturbers.
# It does include the effect of perspective acceleration due to radial velocity, which should 
# be set to the catalog values (not eg the user's model inputs).

# We use this likelihood to reproduce their results.
# """
# function _gaia_skypath_loglikelihood(args,(;
#     αₘ,δₘ,
#     table,
#     ref_epoch,
#     along_scan_uncertainty_mas,
#     full3D,
#     catalog_rv_m_s,
#     astrometric_excess_noise
# ))
#     (
#         plx,
#         ra,
#         dec,
#         pmra,
#         pmdec,
#     ) = args

#     ll = zero(eltype(args))

#     if full3D
#         inner_loop_skypath = AbsoluteVisual{KepOrbit}(;
#             plx,ra,dec,rv=catalog_rv_m_s,pmra,pmdec,
#             ref_epoch=ref_epoch,
#             M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
#         )
#         for i in eachindex(table.epoch)
#             sol = orbitsolve(inner_loop_skypath, table.epoch[i])
#             cmp = sol.compensated
#             # GAIA considers a simplified rectilinear sky path model 
#             # but does include changing parallax from RV
#             # full3D
#             x_earth_pc = table.x[i] / PlanetOrbits.pc2au
#             y_earth_pc = table.y[i] / PlanetOrbits.pc2au
#             z_earth_pc = table.z[i] / PlanetOrbits.pc2au
#             x_diff_pc = cmp.x₂ - x_earth_pc
#             y_diff_pc = cmp.y₂ - y_earth_pc
#             z_diff_pc = cmp.z₂ - z_earth_pc
#             distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
#             mydtor = π / 180
#             ra_apparent_deg = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
#             arg = z_diff_pc / distance_diff
#             arg = map(arg) do arg
#                 if 1.0 < arg < 1.0 + sqrt(eps(1.0))
#                     arg = 1.0
#                 end
#                 return arg
#             end
#             dec_apparent_deg = asin(arg) / mydtor
#             # We don't add perturbations from planets here-- we are fitting Gaia's simplified
#             # unperturbed model to our "true" perturbed data, to see if it reproduces the catalog
#             # values
#             point = @SVector [ra_apparent_deg, dec_apparent_deg]
#             # Two points defining a line along which the star's position was measured
#             line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
#             line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
#             # TODO: distance point to line should account for spherical coordinates!
#             resid = distance_point_to_line(point, line_point_1, line_point_2) # degrees
#             uncertainty_mas = sqrt(along_scan_uncertainty_mas[i]^2 + astrometric_excess_noise^2)
#             ll += logpdf(Normal(0,uncertainty_mas/60/60/1000), resid)
#         end
#         return ll

#     # !full3D , ie fit a simplfied rectilinear model
#     else
#         plx0 = plx
#         dist0 = 1000/plx0
#         δdist_pc_δt_sec = catalog_rv_m_s / 1e3 / IAU_pc2km / PlanetOrbits.sec2day
#         for i in eachindex(table.epoch)
#             Δdist_pc = δdist_pc_δt_sec * (table.epoch[i] - ref_epoch)
#             dist1 = dist0 + Δdist_pc
#             plx_at_time = 1000 / dist1
#             # Calculate sky coordinate Delta at each scan epoch from the catalog position
#             # using the eath's motion ("ephemeris"), x,y,z in AU.
#             # TODO: should I rather be using the catalog Ra and Dec in here? Unclear
#             ra_apparent_Δmas = plx_at_time * (
#                 table.x[i] * sind(ra) -
#                 table.y[i] * cosd(ra)
#             ) + (table.epoch[i] - ref_epoch)/julian_year * pmra
#             dec_apparent_Δmas = plx_at_time * (
#                 table.x[i] * cosd(ra) * sind(dec) +
#                 table.y[i] * sind(ra) * sind(dec) -
#                 table.z[i] * cosd(dec)
#             ) + (table.epoch[i] - ref_epoch)/julian_year * pmdec
#             ra_apparent_deg = ra + ra_apparent_Δmas/(60*60*1000)/cosd(dec)
#             dec_apparent_deg = dec + dec_apparent_Δmas/(60*60*1000)
#             point = @SVector [ra_apparent_deg, dec_apparent_deg] 
#             # Two points defining a line along which the star's position was measured
#             line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
#             line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
#             # TODO: distance point to line should account for spherical coordinates!
#             resid = distance_point_to_line(point, line_point_1, line_point_2) # mas
#             uncertainty_mas = sqrt(along_scan_uncertainty_mas[i]^2 + astrometric_excess_noise^2)
#             ll += logpdf(Normal(0,uncertainty_mas/60/60/1000), resid)
#         end
#         return ll
#     end
# end