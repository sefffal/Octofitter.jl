
# This likelihood partially wraps the Hipparcos and Gaia likelihoods
# We have our own ln_like function that doesn't exactly call into them,
# but we want to use them for loading data etc.
# Our table format is a bit different. For cross validation etc all likelihoods
# allow us to subset epochs. Here, the subsetting should happen at the granularity
# of either the entire Hipparcos or entire Gaia epochs.
# Hence, our table just includes those two values.

"""
    HGCALikelihood(;gaia_id=1234,N_ave=1)

Model Hipparcos-Gaia Catalog of Accelerations (Brandt et al) data using a full model of the Gaia and Hipparcos
measurement process and linear models.

Upon first load, you will be prompted to accept the download of the eDR3 version of the HGCA
catalog.
"""
struct HGCALikelihood{THGCA,THip,TGaia,fluxratio_var} <: AbstractLikelihood
    hgca::THGCA
    hiplike::THip
    gaialike::TGaia
    fluxratio_var::Symbol
    include_dr3_vel::Bool
    include_iad::Bool
end

function _getparams(::HGCALikelihood{THGCA,THip,TGaia,fluxratio_var}, θ_planet) where {THGCA,THip,TGaia,fluxratio_var}
    if fluxratio_var == :__dark
        return (;fluxratio=zero(Octofitter._system_number_type(θ_planet)))
    end
    fluxratio = getproperty(θ_planet, fluxratio_var)
    return (;fluxratio)
end

function HGCALikelihood(; gaia_id, fluxratio_var=nothing, hgca_catalog=(datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits", include_dr3_vel=true, include_iad=false)

    # Load the HCGA
    hgca = FITS(hgca_catalog, "r") do fits
        t = Table(fits[2])
        idx = findfirst(==(gaia_id), t.gaia_source_id)
        return NamedTuple(t[idx])
    end

    # Convert measurement epochs to MJD.
    # Careful: these are Julian years, not decimal years (T. Brant., private communications)
    J2000_mjd = 51544.5 # year J2000 in MJD
    hgca = (;
        hgca...,
        epoch_ra_hip_mjd=(hgca.epoch_ra_hip - 2000) * julian_year + J2000_mjd,
        epoch_dec_hip_mjd=(hgca.epoch_dec_hip - 2000) * julian_year + J2000_mjd,
        epoch_ra_gaia_mjd=(hgca.epoch_ra_gaia - 2000) * julian_year + J2000_mjd,
        epoch_dec_gaia_mjd=(hgca.epoch_dec_gaia - 2000) * julian_year + J2000_mjd,
    )

    # We can get the HIP id from the cross-matching done in the HGCA
    (; hip_id) = hgca

    # Load the Hipparcos IAD data for epochs and scan angles
    hip_like = HipparcosIADLikelihood(;
        hip_id,
        ref_epoch_ra=hgca.epoch_ra_hip_mjd,
        ref_epoch_dec=hgca.epoch_dec_hip_mjd
    )

    # Load the Gaia scanlaw etc
    gaia_like = GaiaCatalogFitLikelihood(;
        gaia_id_dr3=gaia_id,
        ref_epoch_ra=hgca.epoch_ra_gaia_mjd,
        ref_epoch_dec=hgca.epoch_dec_gaia_mjd
    )

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
    c = hgca.pmra_pmdec_hip[1] * hgca.pmra_hip_error[1] * hgca.pmdec_hip_error[1]
    dist_hip = MvNormal(
        @SVector([hgca.pmra_hip, hgca.pmdec_hip]),
        @SArray[
            hgca.pmra_hip_error[1]^2 c
            c hgca.pmdec_hip_error[1]^2
        ]
    )
    # Hipparcos - GAIA epoch
    c = hgca.pmra_pmdec_hg[1] * hgca.pmra_hg_error[1] * hgca.pmdec_hg_error[1]
    dist_hg = MvNormal(
        @SVector([hgca.pmra_hg, hgca.pmdec_hg]),
        @SArray [
            hgca.pmra_hg_error[1]^2 c
            c hgca.pmdec_hg_error[1]^2
        ]
    )
    # GAIA epoch
    c = hgca.pmra_pmdec_gaia[1] * hgca.pmra_gaia_error[1] * hgca.pmdec_gaia_error[1]
    dist_gaia = MvNormal(
        @SVector([hgca.pmra_gaia, hgca.pmdec_gaia]),
        @SArray [
            hgca.pmra_gaia_error[1]^2 c
            c hgca.pmdec_gaia_error[1]^2
        ]
    )

    hgca = (; hgca..., dist_hip, dist_hg, dist_gaia)

    if isnothing(fluxratio_var)
        fluxratio_var = :__dark
    end

    return HGCALikelihood{
        # typeof(table),
        typeof(hgca),
        typeof(hip_like),
        typeof(gaia_like),
        fluxratio_var,
    }(#=table,=# hgca, hip_like, gaia_like, fluxratio_var, include_dr3_vel, include_iad)

end

# function likeobj_from_epoch_subset(obs::HGCALikelihood, obs_inds)
#     # return HGCALikelihood(obs.table[obs_inds, :, 1], obs.hgca, obs.fluxratio_vars) # TODO
# end

function ln_like(hgca_like::HGCALikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    (;μ_g, μ_h, μ_hg) = simulate(hgca_like::HGCALikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)


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
        hgca_like.hgca.nonlinear_dpmra,
        hgca_like.hgca.nonlinear_dpmdec,
    ]

    # also have to remove the HGCA's nonlinear_dpmra/dec from the hipparcos epoch
    # Note: factor of two needed since dpmra is defined to the HG epoch, so H epoch
    # is twice as much. (T. Brandt, private communications).
    μ_h += @SVector [
        2hgca_like.hgca.nonlinear_dpmra,
        2hgca_like.hgca.nonlinear_dpmdec,
    ]

    if hgca_like.include_dr3_vel
        ll += logpdf(hgca_like.hgca.dist_gaia, μ_g)
    end
    ll += logpdf(hgca_like.hgca.dist_hip, μ_h)
    ll += logpdf(hgca_like.hgca.dist_hg, μ_hg)
    # @show μ_h μ_hg μ_g

    return ll
end


function simulate(hgca_like::HGCALikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)

    # * We compute the deviation caused by the planet(s) at each epoch of both likelihoods
    # * Then we perform the 5-param fit only on these deviations
    # * We add the computed deviations to the 5-param catalog values
    # * We compute the HG epoch from the delta positions
    # * We calculate the log-likelihood


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
        diff_lt_app_pmra = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmra2
        diff_lt_app_pmdec = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmdec2
        return cmp_ra.ra2, cmp_dec.dec2, cmp_ra.pmra2+diff_lt_app_pmra, cmp_dec.pmdec2+diff_lt_app_pmdec
    end
    function propagate_astrom(orbit::Any, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end

    @no_escape begin

        # First: Gaia

        Δα_mas = @alloc(T, size(hgca_like.gaialike.table,1))
        fill!(Δα_mas, 0)
        Δδ_mas = @alloc(T, size(hgca_like.gaialike.table,1))
        fill!(Δδ_mas, 0)

        for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
            planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
            (;fluxratio) = _getparams(hgca_like, θ_planet)
            _simulate_skypath_perturbations!(
                Δα_mas, Δδ_mas,
                hgca_like.gaialike.table, orbit,
                planet_mass_msol, fluxratio,
                orbit_solutions[i_planet],
                orbit_solutions_i_epoch_start[i_planet], T
            )
        end


        out = fit_5param_prepared(hgca_like.gaialike.A_prepared_5, hgca_like.gaialike.table, Δα_mas, Δδ_mas)
        # out = fit_4param_prepared(hgca_like.gaialike.A_prepared_4, hgca_like.gaialike.table, Δα_mas, Δδ_mas)
        Δα_g, Δδ_g, Δpmra_g, Δpmdec_g = out.parameters
        # Rigorously propagate the linear proper motion component in spherical coordinates
        # Account for within-gaia differential light travel time 
        α_g₀, δ_g₀, pmra_g₀, pmdec_g₀ = propagate_astrom(first(orbits), hgca_like.hgca.epoch_ra_gaia_mjd, hgca_like.hgca.epoch_dec_gaia_mjd)
        μ_g = @SVector [pmra_g₀ + Δpmra_g, pmdec_g₀ + Δpmdec_g]

        # Now: Hiparcos
        Δα_mas = @alloc(T, size(hgca_like.hiplike.table,1))
        fill!(Δα_mas, 0)
        Δδ_mas = @alloc(T, size(hgca_like.hiplike.table,1))
        fill!(Δδ_mas, 0)

        for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
            planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
            (;fluxratio) = _getparams(hgca_like, θ_planet)
            _simulate_skypath_perturbations!(
                Δα_mas, Δδ_mas,
                hgca_like.hiplike.table, orbit,
                planet_mass_msol, fluxratio,
                orbit_solutions[i_planet],
                orbit_solutions_i_epoch_start[i_planet], T
            )
        end


        if hgca_like.include_iad
            out = fit_5param_prepared(hgca_like.hiplike.A_prepared_5, hgca_like.hiplike.table, Δα_mas, Δδ_mas, hgca_like.hiplike.table.res, hgca_like.hiplike.table.sres)
        else
            out = fit_5param_prepared(hgca_like.hiplike.A_prepared_5, hgca_like.hiplike.table, Δα_mas, Δδ_mas)
        end
        Δα_h, Δδ_h, Δpmra_h, Δpmdec_h = out.parameters
        α_h₀, δ_h₀, pmra_h₀, pmdec_h₀ = propagate_astrom(first(orbits), hgca_like.hgca.epoch_ra_hip_mjd, hgca_like.hgca.epoch_dec_hip_mjd)
        μ_h = @SVector [pmra_h₀ + Δpmra_h, pmdec_h₀ + Δpmdec_h]
    end

    # Simple linear approximation: don't deal with curvature & secular acceleration directly
    if absolute_orbits
        Δα_hg_prop = (α_g₀ - α_h₀)*60*60*1000*cosd((δ_g₀ + δ_h₀)/2)
        Δδ_hg_prop = (δ_g₀ - δ_h₀)*60*60*1000
        pmra_hg_model = (Δα_g - Δα_h + Δα_hg_prop) / (
            hgca_like.hgca.epoch_ra_gaia_mjd - hgca_like.hgca.epoch_ra_hip_mjd
        )*julian_year
        pmdec_hg_model = (Δδ_g - Δδ_h + Δδ_hg_prop) / (
            hgca_like.hgca.epoch_dec_gaia_mjd - hgca_like.hgca.epoch_dec_hip_mjd
        )*julian_year
    else
        pmra_hg_model = (Δα_g - Δα_h) / (
                hgca_like.hgca.epoch_ra_gaia_mjd - hgca_like.hgca.epoch_ra_hip_mjd
        )*julian_year + θ_system.pmra
        pmdec_hg_model = (Δδ_g - Δδ_h) / (
            hgca_like.hgca.epoch_dec_gaia_mjd - hgca_like.hgca.epoch_dec_hip_mjd
        )*julian_year + θ_system.pmdec
    end

    μ_hg = @SVector [pmra_hg_model, pmdec_hg_model]

    return (;
        # Packaged up nicely
        μ_g,
        μ_h,
        μ_hg,
        # Individual
        pmra_hip_model=μ_h[1],
        pmdec_hip_model=μ_h[2],
        pmra_gaia_model=μ_g[1],
        pmdec_gaia_model=μ_g[2],
        pmra_hg_model=μ_hg[1],
        pmdec_hg_model=μ_hg[2],
    )


    # TODO -- how to remove the nonlinear_dpmra?
    # It's already baked into the pm_hg and pm_hip for example
    # If they've already added nonlinear_dpmra,
    # and here we account for it ourselves,
    # then we should subtract it from our results before the comparison

end
export HGCALikelihood
