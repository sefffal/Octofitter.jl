#  Proper motion anomaly functions

# This code provides helpers for loading data from the HIPPARCOS GAIA Catalog of Accelerations.
"""
    gaia_plx(gaia_id=12123)

Get a distribution (truncated Normal) of parallax distance in mas of a source with 
GAIA catalog id `gaia_id`.
"""
function gaia_plx(; gaia_id, catalog=(datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits")

    # Load the Hipparcos-GAIA catalog of accelerations as a table
    hgca = FITS(catalog, "r") do fits
        Table(fits[2])
    end

    idx = findfirst(==(gaia_id), hgca.gaia_source_id)
    μ = Float64(hgca.parallax_gaia[idx,])
    σ = Float64(hgca.parallax_gaia_error[idx,])
    # Truncate the normal distribution explicitly to prevent numerical issues 
    # This formally constricts the range of the parallax parameter for some downstream methods
    return truncated(Normal(μ,σ), lower=μ-10σ, upper=μ+10σ)
end
export gaia_plx


# function ghca_pmra(;gaia_id)

struct HGCAInstantaneousLikelihood{TTable<:Table,THGCA,flux_vars} <: AbstractLikelihood
    table::TTable
    hgca::THGCA
    flux_vars::flux_vars
end
function likeobj_from_epoch_subset(obs::HGCAInstantaneousLikelihood, obs_inds)
    return HGCAInstantaneousLikelihood(obs.table[obs_inds,:,1], obs.hgca, obs.flux_vars)
end
export HGCAInstantaneousLikelihood


# We split the HGCA Likelihood into three 
struct TimeAveragedProperMotionLikelihood{TTable<:Table,THGCA} <: AbstractLikelihood
    table::TTable
end
function likeobj_from_epoch_subset(obs::TimeAveragedProperMotionLikelihood, obs_inds)
    return TimeAveragedProperMotionLikelihood(obs.table[obs_inds,:,1], obs.hgca)
end
export TimeAveragedProperMotionLikelihood



"""
    HGCAInstantaneousLikelihood(;gaia_id=1234,N_ave=1)

Model Hipparcos-Gaia Catalog of Accelerations (Brandt et al) data using a simplified instantaneous
model -- calculating the proper motion and positions at 1 to N epochs and averaging them. 
This is quicker than the full `HGCALikelihood` which involves fitting linear least-squares
solutions at each epoch. It should be entirely equivalent for systems without strong acceleration
*during* the Hipparcos or Gaia missions. 

Upon first load, you will be prompted to accept the download of the eDR3 version of the HGCA 
catalog.
"""
function HGCAInstantaneousLikelihood(;
    gaia_id,
    catalog=(datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits",
    N_ave=1,
    factor=1,
    flux_variables=(),
)

# hgca_all = Table([h])
    hgca_all = FITS(catalog, "r") do fits
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

    # Convert measurement epochs to MJD.
    # The HGCA doesn't say, but we assume these are actually Julian years and not decimal years.
    J2000_mjd = 51544.5 # year J2000 in MJD
    epoch_ra_hip_mjd = (hgca_all.epoch_ra_hip[idx] - 2000)*julian_year + J2000_mjd
    epoch_dec_hip_mjd = (hgca_all.epoch_dec_hip[idx] - 2000)*julian_year + J2000_mjd
    epoch_ra_gaia_mjd = (hgca_all.epoch_ra_gaia[idx] - 2000)*julian_year + J2000_mjd
    epoch_dec_gaia_mjd = (hgca_all.epoch_dec_gaia[idx] - 2000)*julian_year + J2000_mjd
    
    # Roughly over what time period were the observations made?
    dt_gaia = 1038 # EDR3: days between  Date("2017-05-28") - Date("2014-07-25")
    dt_hip = 4 * 365.25 # 4 years for Hipparcos

    # How many points over Δt should we average the proper motion and stellar position
    # at each epoch? This is because the PM is not an instantaneous measurement.
    if N_ave == 1 
        δt_hip = δt_gaia = 0.
    else
        # TODO: could use actual scan epochs from Hipparcos and GAIA
        δt_hip = range(-dt_hip / 2, dt_hip / 2, N_ave)
        δt_gaia = range(-dt_gaia / 2, dt_gaia / 2, N_ave)
    end

    # table cols: epoch meas inst
    rows = NamedTuple[]
    # Hipparcos
    for δt in δt_hip
        push!(rows, (;epoch=epoch_ra_hip_mjd + δt, meas=:ra, inst=:hip))
        push!(rows, (;epoch=epoch_dec_hip_mjd + δt, meas=:dec, inst=:hip))
    end
    # Gaia
    for δt in δt_gaia
        push!(rows, (;epoch=epoch_ra_gaia_mjd + δt, meas=:ra, inst=:gaia))
        push!(rows, (;epoch=epoch_dec_gaia_mjd + δt, meas=:dec, inst=:gaia))
    end

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
    ].*factor^2)
    # Hipparcos - GAIA epoch
    c = hgca.pmra_pmdec_hg[1] * hgca.pmra_hg_error[1] * hgca.pmdec_hg_error[1]
    dist_hg = MvNormal(@SArray[
        hgca.pmra_hg_error[1]^2 c
        c hgca.pmdec_hg_error[1]^2
    ].*factor^2)
    # GAIA epoch
    c = hgca.pmra_pmdec_gaia[1] * hgca.pmra_gaia_error[1] * hgca.pmdec_gaia_error[1]
    dist_gaia = MvNormal(@SArray[
        hgca.pmra_gaia_error[1]^2 c
        c hgca.pmdec_gaia_error[1]^2
    ].*factor^2)
    
    hgca = (;hgca...,dist_hip,dist_hg,dist_gaia)

    return HGCAInstantaneousLikelihood(Table(rows), hgca, flux_variables)
end
export HGCAInstantaneousLikelihood


"""
Specific HGCA proper motion modelling. Model the GAIA-Hipparcos/Δt proper motion
using 5 position measurements averaged at each of their epochs.
"""
function ln_like(hgca_like::HGCAInstantaneousLikelihood, θ_system, elements, orbit_solutions, orbit_solutions_i_epoch_start)
    ll = 0.0

    (;
        pmra_hip_model,
        pmdec_hip_model,
        pmra_gaia_model,
        pmdec_gaia_model,
        pmra_hg_model,
        pmdec_hg_model,
    ) = _simulate_hgca_instant(hgca_like, θ_system, elements, orbit_solutions, orbit_solutions_i_epoch_start)


    absolute_orbits = false
    for orbit in elements
        absolute_orbits |= orbit isa AbsoluteVisual
        # TODO: could check in a more user-friendly way
        # that we don't have a mismatch of different orbit types 
        # for different planets?
    end

    # The HGCA catalog values have an non-linearity correction added.
    # If we are doing our own rigorous propagation we don't need this
    # correction. We could subtract it from the measurements, but 
    # here we just add it to our model so that they match
    if absolute_orbits
        hg_nonlinear_dpmra = hgca_like.hgca.nonlinear_dpmra[1]
        hg_nonlinear_dpmdec = hgca_like.hgca.nonlinear_dpmdec[1]
        hip_nonlinear_dpmra = 2hgca_like.hgca.nonlinear_dpmra[1]
        hip_nonlinear_dpmdec = 2hgca_like.hgca.nonlinear_dpmdec[1]
    else
        hg_nonlinear_dpmra = 
        hg_nonlinear_dpmdec = 
        hip_nonlinear_dpmra = 
        hip_nonlinear_dpmdec = zero(hgca_like.hgca.nonlinear_dpmra[1])
    end

    # Hipparcos epoch
    resids_hip = @SArray[
        pmra_hip_model  - (hgca_like.hgca.pmra_hip - hip_nonlinear_dpmra),
        pmdec_hip_model - (hgca_like.hgca.pmdec_hip - hip_nonlinear_dpmdec)
    ]
    ll += logpdf(hgca_like.hgca.dist_hip, resids_hip)

    # Hipparcos - GAIA epoch
    resids_hg = @SArray[
        pmra_hg_model - (hgca_like.hgca.pmra_hg - hg_nonlinear_dpmra),
        pmdec_hg_model  - (hgca_like.hgca.pmdec_hg - hg_nonlinear_dpmdec)
    ]
    ll += logpdf(hgca_like.hgca.dist_hg, resids_hg)

    # GAIA epoch
    resids_gaia = @SArray[
        pmra_gaia_model - hgca_like.hgca.pmra_gaia,
        pmdec_gaia_model - hgca_like.hgca.pmdec_gaia
    ]
    @warn "skipping gaia epoch" maxlog=1
    ll += logpdf(hgca_like.hgca.dist_gaia, resids_gaia)

    # display(elements[1])
    # @show ll

    return ll
end


function _simulate_hgca_instant(pma, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)

    # This observation type just wraps one row from the HGCA (see hgca.jl)
    hgca = pma.hgca

    # Look at the position of the star around both epochs to calculate 
    # our modelled delta-position proper motion

    # If the user specified a AbsoluteVisual orbit, we will compute things a
    # little differently
    absolute_orbits = false
    # for orbit in orbits
    #     absolute_orbits |= orbit isa AbsoluteVisual
    #     # TODO: could check in a more user-friendly way
    #     # that we don't have a mismatch of different orbit types 
    #     # for different planets?
    # end

    deg2mas = 60 * 60 * 1000
    # First epoch: Hipparcos
    # Note: the catalog is stored as Float32, but we really want to do 
    # this math in Float64. Adding/subtracting small differences in RA and Dec,
    # stored in degrees, actually needs more precision than you would expect.
    # In Float32, we would easily get round off at the 0.3% relative error level.
    ra_hip_model = zero(T)
    dec_hip_model = zero(T)
    pmra_hip_model = zero(T)
    pmdec_hip_model = zero(T)
    epoch_hip_ra = 0.0
    epoch_hip_dec = 0.0
    N_ave_hip_ra = 0
    N_ave_hip_dec = 0
    # The model can support multiple planets
    for i_planet in eachindex(orbits)
        θ_planet = θ_system.planets[i_planet]
        orbit = orbits[i_planet]
        # TODO: a trait would be better here
        if !(
            (
                (orbit isa Visual || orbit isa AbsoluteVisual) && (orbit.parent isa KepOrbit)
            ) || orbit isa ThieleInnesOrbit
        )
            continue
        end
        for i_epoch in eachindex(pma.table.epoch, pma.table.inst, pma.table.meas)
            if pma.table.inst[i_epoch] != :hip
                continue
            end
            
            sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]
            if pma.table.meas[i_epoch] == :ra
                N_ave_hip_ra += 1
                epoch_hip_ra += pma.table.epoch[i_epoch]
                ra_hip_model += raoff(sol, θ_planet.mass * mjup2msol)
                if absolute_orbits
                    ra_hip_model += deg2mas*(sol.compensated.ra2)
                end
                pmra_hip_model += pmra(sol, θ_planet.mass * mjup2msol)
            elseif pma.table.meas[i_epoch] == :dec
                N_ave_hip_dec += 1
                epoch_hip_dec += pma.table.epoch[i_epoch]
                dec_hip_model += decoff(sol, θ_planet.mass * mjup2msol)
                if absolute_orbits
                    dec_hip_model += deg2mas*(sol.compensated.dec2)
                end
                pmdec_hip_model += pmdec(sol, θ_planet.mass * mjup2msol)
            end
        end
    end
    ra_hip_model /= N_ave_hip_ra
    dec_hip_model /= N_ave_hip_dec
    pmra_hip_model /= N_ave_hip_ra
    pmdec_hip_model /= N_ave_hip_dec
    epoch_hip_ra /= N_ave_hip_ra
    epoch_hip_dec /= N_ave_hip_dec
    pmra_hip_model += θ_system.pmra
    pmdec_hip_model += θ_system.pmdec

    # Last epoch: GAIA
    ra_gaia_model = zero(T)
    dec_gaia_model = zero(T)
    pmra_gaia_model = zero(T)
    pmdec_gaia_model = zero(T)
    epoch_gaia_ra = 0.0
    epoch_gaia_dec = 0.0
    N_ave_gaia_ra = 0
    N_ave_gaia_dec = 0
    # The model can support multiple planets
    for i_planet in eachindex(orbits)
        θ_planet = θ_system.planets[i_planet]
        orbit = orbits[i_planet]
        # TODO: a trait would be better here
        if !(
            (
                (orbit isa Visual || orbit isa AbsoluteVisual) && (orbit.parent isa KepOrbit)
            ) || orbit isa ThieleInnesOrbit
        )
            continue
        end
        for i_epoch in eachindex(pma.table.epoch, pma.table.inst, pma.table.meas)
            if pma.table.inst[i_epoch] != :gaia
                continue
            end
            
            sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]
            if pma.table.meas[i_epoch] == :ra
                N_ave_gaia_ra += 1
                epoch_gaia_ra += pma.table.epoch[i_epoch]
                ra_gaia_model += raoff(sol, θ_planet.mass * mjup2msol)
                if absolute_orbits
                    ra_gaia_model += deg2mas*(sol.compensated.ra2)
                end
                pmra_gaia_model += pmra(sol, θ_planet.mass * mjup2msol)
            elseif pma.table.meas[i_epoch] == :dec
                N_ave_gaia_dec += 1
                epoch_gaia_dec += pma.table.epoch[i_epoch]
                dec_gaia_model += decoff(sol, θ_planet.mass * mjup2msol)
                if absolute_orbits
                    dec_gaia_model += deg2mas*(sol.compensated.dec2)
                end
                pmdec_gaia_model += pmdec(sol, θ_planet.mass * mjup2msol)
            end
        end
    end
    ra_gaia_model /= N_ave_gaia_ra
    dec_gaia_model /= N_ave_gaia_dec
    pmra_gaia_model /= N_ave_gaia_ra
    pmdec_gaia_model /= N_ave_gaia_dec
    epoch_gaia_ra /= N_ave_gaia_ra
    epoch_gaia_dec /= N_ave_gaia_dec
    pmra_gaia_model += θ_system.pmra
    pmdec_gaia_model += θ_system.pmdec

    # Model the GAIA-Hipparcos delta-position velocity in mas/yr
    if absolute_orbits
        # Calculate the PM using the average epoch
        pmra_hg_model = (ra_gaia_model - ra_hip_model) / (epoch_gaia_ra - epoch_hip_ra)*julian_year # ((hgca.epoch_ra_gaia_mjd[1] - hgca.epoch_ra_hip_mjd[1])/julian_year)
        pmdec_hg_model = (dec_gaia_model - dec_hip_model) / (epoch_gaia_dec - epoch_hip_dec)*julian_year # ((hgca.epoch_dec_gaia_mjd[1] - hgca.epoch_dec_hip_mjd[1])/julian_year)
        # Cosine factor to go from alpha to alpha-star
        ave_dec = (dec_gaia_model + dec_hip_model)/2
        pmra_hg_model *= cosd(ave_dec/deg2mas)
    end

    # Simple linear approximation: don't deal with curvature & secular acceleration directly
    if !absolute_orbits
        pmra_hg_model = (ra_gaia_model - ra_hip_model) / (epoch_gaia_ra - epoch_hip_ra)*julian_year # ((hgca.epoch_ra_gaia_mjd[1] - hgca.epoch_ra_hip_mjd[1])/julian_year)
        pmdec_hg_model = (dec_gaia_model - dec_hip_model) / (epoch_gaia_dec - epoch_hip_dec)*julian_year # ((hgca.epoch_dec_gaia_mjd[1] - hgca.epoch_dec_hip_mjd[1])/julian_year)
        pmra_hg_model += θ_system.pmra
        pmdec_hg_model += θ_system.pmdec
    end

    result = (;
        pmra_hip_model,
        pmdec_hip_model,
        pmra_gaia_model,
        pmdec_gaia_model,
        pmra_hg_model,
        pmdec_hg_model,
    )

    return result
end




"""
Specific HGCA proper motion modelling. Model the GAIA-Hipparcos/Δt proper motion
using 25 position measurements averaged at each of their epochs.
"""
function generate_from_params(hgca_like::HGCAInstantaneousLikelihood, θ_system, orbits, solutions, sol_start_i)

    (;
        pmra_hip_model,
        pmdec_hip_model,
        pmra_gaia_model,
        pmdec_gaia_model,
        pmra_hg_model,
        pmdec_hg_model,
    ) = _simulate_hgca(hgca_like, θ_system, orbits, solutions, sol_start_i)

    # Merge the measurements together into a new observation and add noise according to the sigma
    # we were passed in from the original measurements
    return HGCAInstantaneousLikelihood(hgca_like.table, merge(hgca_like.hgca, (;
        pmra_hip=pmra_hip_model,
        pmdec_hip=pmdec_hip_model,
        pmra_gaia=pmra_gaia_model,
        pmdec_gaia=pmdec_gaia_model,
        pmra_hg=pmra_hg_model,
        pmdec_hg=pmdec_hg_model,
    )))

end

