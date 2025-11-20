
# This likelihood partially wraps the Hipparcos and Gaia likelihoods
# We have our own ln_like function that doesn't exactly call into them,
# but we want to use them for loading data etc.
# Our table format is a bit different. For cross validation etc all likelihoods
# allow us to subset epochs. Here, the subsetting should happen at the granularity
# of either the entire Hipparcos or entire Gaia epochs.
# Hence, our table just includes those two values.

"""
    HGCAObs(;
        gaia_id=1234,
        variables=@variables begin
            fluxratio ~ [Uniform(0, 1), Uniform(0, 1)]  # array for each companion
        end
    )

Model Hipparcos-Gaia Catalog of Accelerations (Brandt et al) data using a full model of the Gaia and Hipparcos
measurement process and linear models.

The `fluxratio` variable should be an array containing the flux ratio of each companion
in the same order as the planets in the system.

Upon first load, you will be prompted to accept the download of the eDR3 version of the HGCA
catalog.
"""
struct HGCAObs{TTable,THGCA,THip,TGaia} <: AbstractObs
    table::TTable
    hgca::THGCA
    priors::Priors
    derived::Derived
    hip_like::THip
    gaia_like::TGaia
    include_iad::Bool
end

# Backwards compatibility alias
const HGCALikelihood = HGCAObs

# Special method for HGCAObs which doesn't have an instrument_name field
likelihoodname(::HGCAObs) = "HGCA"

function HGCAObs(;
    variables::Tuple{Priors,Derived}=(@variables begin;end ),
    gaia_id, hgca_catalog=(datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits",
    include_iad=true)

    (priors,derived)=variables

    # Load the HCGA
    hgca = FITS(hgca_catalog, "r") do fits
        t = Table(fits[2])
        idx = findfirst(==(gaia_id), t.gaia_source_id)
        return NamedTuple(t[idx])
    end
    return HGCAObs(hgca, priors, derived; include_iad)
end
function HGCAObs(hgca::NamedTuple, priors::Priors, derived::Derived; include_iad)


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
        gaia_id_dr3=hgca.gaia_source_id,
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


    # This table serves only for plotting and keeping track of subsetting for cross-validation.
    # The actual likelihood calculations happen against the prepared linear system matrices above.
    table = Table(
        epoch=[
            isnothing(hip_like) ? NaN : years2mjd(hgca.epoch_ra_hip),
            isnothing(hip_like) ? NaN : years2mjd(hgca.epoch_dec_hip),
            isnothing(hip_like) ? NaN : years2mjd((hgca.epoch_ra_hip+hgca.epoch_ra_gaia)/2),
            isnothing(hip_like) ? NaN : years2mjd((hgca.epoch_dec_hip+hgca.epoch_dec_gaia)/2),
            years2mjd(hgca.epoch_ra_gaia),
            years2mjd(hgca.epoch_dec_gaia),
        ],
        start_epoch=[
            isnothing(hip_like) ? 0 : minimum(hip_like.table.epoch),
            isnothing(hip_like) ? 0 : minimum(hip_like.table.epoch),
            isnothing(hip_like) ? 0 : years2mjd(hgca.epoch_ra_hip),
            isnothing(hip_like) ? 0 : years2mjd(hgca.epoch_dec_hip),
            first(gaia_like.table.epoch[gaia_like.table.epoch.>=(meta_gaia_DR3.start_mjd)]),
            first(gaia_like.table.epoch[gaia_like.table.epoch.>=(meta_gaia_DR3.start_mjd)]),
        ],
        stop_epoch=[
            isnothing(hip_like) ? 0 : maximum(hip_like.table.epoch),
            isnothing(hip_like) ? 0 : maximum(hip_like.table.epoch),
            isnothing(hip_like) ? 0 : years2mjd(hgca.epoch_ra_gaia),
            isnothing(hip_like) ? 0 : years2mjd(hgca.epoch_dec_gaia),
            last(gaia_like.table.epoch[gaia_like.table.epoch.<=(meta_gaia_DR2.stop_mjd)]),
            last(gaia_like.table.epoch[gaia_like.table.epoch.<=(meta_gaia_DR2.stop_mjd)]),
        ],
        pm = [
            hgca.pmra_hip,
            hgca.pmdec_hip,
            hgca.pmra_hg,
            hgca.pmdec_hg, 
            hgca.pmra_gaia,
            hgca.pmdec_gaia,
        ],
        σ_pm = [
            hgca.pmra_hip_error,
            hgca.pmdec_hip_error,
            hgca.pmra_hg_error,
            hgca.pmdec_hg_error,
            hgca.pmra_gaia_error,
            hgca.pmdec_gaia_error,
        ],
        kind=[
            :ra_hip,
            :dec_hip,
            :ra_hg,
            :dec_hg,
            :ra_gaia,
            :dec_gaia,
        ],

    )


    return HGCAObs{
        typeof(table),
        typeof(hgca),
        typeof(hip_like),
        typeof(gaia_like),
    }(table, hgca, priors, derived, hip_like, gaia_like, include_iad)

end

function Octofitter.likeobj_from_epoch_subset(like::HGCAObs, obs_inds)
    (;
        table,
        hgca,
        priors,
        derived,
        hip_like,
        gaia_like,
        include_iad,
    ) = like

    table = table[obs_inds,:]
    return HGCAObs{
        typeof(table),
        typeof(hgca),
        typeof(hip_like),
        typeof(gaia_like),
    }(table, hgca, priors, derived, hip_like, gaia_like, include_iad)
end

function ln_like(obs::HGCAObs, ctx::SystemObservationContext)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    sim = simulate(obs, θ_system, θ_obs, orbits, orbit_solutions, -1)

    if isnothing(sim)
        return convert(T, -Inf)
    end
    (;μ_g, μ_h, μ_hg) = sim

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
        obs.hgca.nonlinear_dpmra,
        obs.hgca.nonlinear_dpmdec,
    ]

    # also have to remove the HGCA's nonlinear_dpmra/dec from the hipparcos epoch
    # Note: factor of two needed since dpmra is defined to the HG epoch, so H epoch
    # is twice as much. (T. Brandt, private communications).
    μ_h += @SVector [
        2obs.hgca.nonlinear_dpmra,
        2obs.hgca.nonlinear_dpmdec,
    ]

    if :ra_gaia ∈ obs.table.kind && :dec_gaia ∈ obs.table.kind
        ll += logpdf(obs.hgca.dist_gaia, μ_g)
    elseif :ra_gaia ∈ obs.table.kind
        μ, Σ = params(obs.hgca.dist_gaia)
        ll += logpdf(Normal(μ[1], sqrt(Σ[1,1])), μ_g[1])
    elseif :dec_gaia ∈ obs.table.kind
        μ, Σ = params(obs.hgca.dist_gaia)
        ll += logpdf(Normal(μ[2], sqrt(Σ[2,2])), μ_g[2])
    end
    if :ra_hip ∈ obs.table.kind && :dec_hip ∈ obs.table.kind
        ll += logpdf(obs.hgca.dist_hip, μ_h)
    elseif :ra_hip ∈ obs.table.kind
        μ, Σ = params(obs.hgca.dist_hip)
        ll += logpdf(Normal(μ[1], sqrt(Σ[1,1])), μ_h[1])
    elseif :dec_hip ∈ obs.table.kind
        μ, Σ = params(obs.hgca.dist_hip)
        ll += logpdf(Normal(μ[2], sqrt(Σ[2,2])), μ_h[2])
    end
    if :ra_hg ∈ obs.table.kind && :dec_hg ∈ obs.table.kind
        ll += logpdf(obs.hgca.dist_hg, μ_hg)
    elseif :ra_hg ∈ obs.table.kind
        μ, Σ = params(obs.hgca.dist_hg)
        ll += logpdf(Normal(μ[1], sqrt(Σ[1,1])), μ_hg[1])
    elseif :dec_hg ∈ obs.table.kind
        μ, Σ = params(obs.hgca.dist_hg)
        ll += logpdf(Normal(μ[2], sqrt(Σ[2,2])), μ_hg[2])
    end

    return ll
end


function simulate(hgca_like::HGCAObs, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
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


    if hasproperty(θ_system, :missed_transits)
        (;missed_transits) = θ_system
        if eltype(missed_transits) <: AbstractFloat
            missed_transits = Int.(missed_transits)
        end
        if length(unique(missed_transits)) < length(missed_transits)
            return nothing
        end
        ii = sort(setdiff(1:length(hgca_like.gaia_like.table.epoch), missed_transits))
        gaia_table = hgca_like.gaia_like.table[ii,:]
        A_prepared_5 = hgca_like.gaia_like.A_prepared_5[ii,:]
    else
        gaia_table = hgca_like.gaia_like.table
        A_prepared_5 = hgca_like.gaia_like.A_prepared_5
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

        Δα_mas = @alloc(T, size(A_prepared_5,1))
        fill!(Δα_mas, 0)
        Δδ_mas = @alloc(T, size(A_prepared_5,1))
        fill!(Δδ_mas, 0)

        for (i_planet, (orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
            planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
            fluxratio = hasproperty(θ_obs, :fluxratio) ? θ_obs.fluxratio[i_planet] : zero(T)
            _simulate_skypath_perturbations!(
                Δα_mas, Δδ_mas,
                gaia_table, orbit,
                planet_mass_msol, fluxratio,
                orbit_solutions[i_planet],
                orbit_solutions_i_epoch_start, T
            )
        end


        out = fit_5param_prepared(A_prepared_5, gaia_table, Δα_mas, Δδ_mas)
        # out = fit_4param_prepared(A_prepared_4, gaia_table, Δα_mas, Δδ_mas)
        Δα_g, Δδ_g, Δpmra_g, Δpmdec_g = out.parameters
        # Rigorously propagate the linear proper motion component in spherical coordinates
        # Account for within-gaia differential light travel time 
        α_g₀, δ_g₀, pmra_g₀, pmdec_g₀ = propagate_astrom(first(orbits), hgca_like.hgca.epoch_ra_gaia_mjd, hgca_like.hgca.epoch_dec_gaia_mjd)
        μ_g = @SVector [pmra_g₀ + Δpmra_g, pmdec_g₀ + Δpmdec_g]

        # Now: Hiparcos
        Δα_mas = @alloc(T, size(hgca_like.hip_like.table,1))
        fill!(Δα_mas, 0)
        Δδ_mas = @alloc(T, size(hgca_like.hip_like.table,1))
        fill!(Δδ_mas, 0)

        for (i_planet, (orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
            planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
            fluxratio = hasproperty(θ_obs, :fluxratio) ? θ_obs.fluxratio[i_planet] : zero(T)
            _simulate_skypath_perturbations!(
                Δα_mas, Δδ_mas,
                hgca_like.hip_like.table, orbit,
                planet_mass_msol, fluxratio,
                orbit_solutions[i_planet],
                orbit_solutions_i_epoch_start, T
            )
        end


        if hgca_like.include_iad
            out = fit_5param_prepared(hgca_like.hip_like.A_prepared_5, hgca_like.hip_like.table, Δα_mas, Δδ_mas, hgca_like.hip_like.table.res, hgca_like.hip_like.table.sres)
        else
            out = fit_5param_prepared(hgca_like.hip_like.A_prepared_5, hgca_like.hip_like.table, Δα_mas, Δδ_mas)
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

    # Adjust the reference frame such that, effectively, the pmra/pmdec system variables are referring to the primary
    # instead of the barycentre.
    # Specifically, the primary's proper motion at this epoch:
    μ_h  = μ_h  .- @SVector [Δpmra_g, Δpmdec_g,]
    μ_hg = μ_hg .- @SVector [Δpmra_g, Δpmdec_g,]
    μ_g  = μ_g  .- @SVector [Δpmra_g, Δpmdec_g,]

    # this drastically improves convergence for wide orbits / massive companions

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
export HGCAObs, HGCALikelihood



# Generate new astrometry observations
function generate_from_params(like::HGCAObs, ctx::SystemObservationContext; add_noise)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx

    sim = simulate(like, θ_system, θ_obs, orbits, orbit_solutions, -1)

    (;μ_g, μ_h, μ_hg) = sim

    # replace values in the HGCA with our new ones
    hgca = (;
        like.hgca...,
        pmra_hip = add_noise    ? μ_h[1]  + randn()*like.hgca.pmra_hip_error   : μ_h[1],
        pmdec_hip = add_noise   ? μ_h[2]  + randn()*like.hgca.pmdec_hip_error  : μ_h[2],
        pmra_hg = add_noise     ? μ_hg[1] + randn()*like.hgca.pmra_hg_error    : μ_hg[1],
        pmdec_hg = add_noise    ? μ_hg[2] + randn()*like.hgca.pmdec_hg_error   : μ_hg[2],
        pmra_gaia = add_noise   ? μ_g[1]  + randn()*like.hgca.pmra_gaia_error  : μ_g[1],
        pmdec_gaia = add_noise  ? μ_g[2]  + randn()*like.hgca.pmdec_gaia_error : μ_g[2],
    )
    new_hgca_like = HGCAObs(hgca, like.priors, like.derived; like.include_iad)
    # What do we do about the Hipparcos residuals?
    if like.include_iad
        @warn "Hipparcos residuals are not currently handled in the simulation"
        # We will need to simulate these too -- just run simulate sky path,
        # then do the 5 param fit, then store the residuals
    end
    # zero out any hipparcos residuals
    new_hgca_like.hip_like.table.res .= 0

    return new_hgca_like
end
