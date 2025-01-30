

struct AbsAstromLikelihood{TTable,TDist} <: AbstractLikelihood
    table::TTable
    dist::TDist
end
function ln_like(like::AbsAstromLikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 
    ll, _ = simulate(like, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    return ll
end


function simulate(like::AbsAstromLikelihood, θ_system, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

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

    istart = findfirst(>=(meta_gaia_DR3.start_mjd), vec(like.table.epoch))
    iend = findlast(<=(meta_gaia_DR3.stop_mjd), vec(like.table.epoch))
    if isnothing(istart)
        istart = 1
    end
    if isnothing(iend)
        iend = length(like.table.epoch)
    end
    Δα_mas = zeros(T, iend-istart+1)
    Δδ_mas = zeros(T, iend-istart+1)
    for (orbit, θ_planet) in zip(orbits, θ_system.planets)
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            like.table[istart:iend], orbit,
            planet_mass_msol, 0.0,
            orbit_solutions,
            orbit_solutions_i_epoch_start, T
        )
    end
    out = fit_4param(
        like.table[istart:iend],
        Δα_mas,
        Δδ_mas,
        meta_gaia_DR3.ref_epoch_mjd,
        meta_gaia_DR3.ref_epoch_mjd,
    )
    Δα, Δδ, Δμα, Δμδ = out.parameters
    modelled_gaia_parameters_dr3 = [
        like.dr3.ra + Δα/60/60/1000,
        like.dr3.dec + Δδ/60/60/1000/cosd(like.dr3.dec),
        θ_system.pmra+Δμα, #like.dr3.pmra+Δμα,
        θ_system.pmdec+Δμδ, #like.dr3.pmdec+Δμδ
    ] 


    istart = findfirst(>=(meta_gaia_DR2.start_mjd), vec(like.table.epoch))
    iend = findlast(<=(meta_gaia_DR2.stop_mjd), vec(like.table.epoch))
    if isnothing(istart)
        istart = 1
    end
    if isnothing(iend)
        iend = length(like.table.epoch)
    end
    Δα_mas = zeros(T, iend-istart+1)
    Δδ_mas = zeros(T, iend-istart+1)
    for (orbit, θ_planet) in zip(orbits, θ_system.planets)
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            like.table[istart:iend], orbit,
            planet_mass_msol, 0.0,
            orbit_solutions,
            orbit_solutions_i_epoch_start, T
        )
    end
    out = fit_4param(
        like.table[istart:iend],
        Δα_mas,
        Δδ_mas,
        meta_gaia_DR2.ref_epoch_mjd,
        meta_gaia_DR2.ref_epoch_mjd,
    )
    Δα, Δδ, Δμα, Δμδ = out.parameters
    modelled_gaia_parameters_dr2 = [
        like.dr2.ra + Δα/60/60/1000,
        like.dr2.dec + Δδ/60/60/1000/cosd(like.dr2.dec),
        θ_system.pmra+Δμα, # like.dr2.pmra+Δμα,
        θ_system.pmdec+Δμδ, # like.dr2.pmdec+Δμδ
    ] 


    # Build our overall correlation matrix -- the top left and bottom right blocks come
    # from the catalogs.
    # The off-diagonal blocks and fit using a semi-analytic structure.
    # (;α_pos, α_pm, β_pos, δ_pos, δ_pm, γ) = θ_system
    # K = build_correlation_structure(α_pos, α_pm, β_pos, δ_pos, δ_pm, γ)
    Σ_dr2 = Diagonal(like.σ_dr2) * like.C_dr2 * Diagonal(like.σ_dr2)
    Σ_dr3 = Diagonal(like.σ_dr3) * like.C_dr3 * Diagonal(like.σ_dr3)
    ρ = sqrt(22/33)
    K = ρ*sqrt(Σ_dr2)*sqrt(Σ_dr3)
    Σ_dr2_dr3 = [
        Σ_dr2  K 
        K               Σ_dr3
    ]

    # The Gaia DR2 reported parameter values are offset in various ways vs. DR2. 
    # Correct catalog values from DR2 for known effects:
    correction = @SVector [
        # θ_system.dr2_systematic_Δplx,
        θ_system.dr2_systematic_Δra/60/60/1000/cosd(like.dr2.dec),
        θ_system.dr2_systematic_Δdec/60/60/1000,
        θ_system.dr2_systematic_Δμ_ra,
        θ_system.dr2_systematic_Δμ_dec,
    ]
    μ_dr2_corrected = like.μ_dr2 .+ correction
    μ_dr2_dr3 = [μ_dr2_corrected; like.μ_dr3]
     
    # σ_dr2_dr3 = [like.σ_dr2; like.σ_dr3]
    # Σ_dr2_dr3 = Diagonal(σ_dr2_dr3) * C_dr2_dr3 * Diagonal(σ_dr2_dr3)

    # display(Σ_dr2_dr3)


    ##########
    # velocities and positions

    # # Not all values of correlation parameters are valid -- some may not produce
    # # a positive semi-definite covariance matrix.
    # dist_dr2_dr3 = try
    #     MvNormal(μ_dr2_dr3,Hermitian(Σ_dr2_dr3))
    # catch
    #     return -Inf, (nothing,nothing)
    # end

    # μ_dr2_dr3_modelled = [modelled_gaia_parameters_dr2; modelled_gaia_parameters_dr3]
    # ll += logpdf(dist_dr2_dr3, μ_dr2_dr3_modelled)


    ##########
    # Only 2 epoch velocities
    μ_dr2_dr3 = μ_dr2_dr3[[3:4;7:8]]
    Σ_dr2_dr3 = Σ_dr2_dr3[[3:4;7:8],[3:4;7:8]]
    dist_dr2_dr3 = try
        MvNormal(μ_dr2_dr3,Hermitian(Σ_dr2_dr3))
    catch
        return -Inf, (nothing,nothing)
    end

    μ_dr2_dr3_modelled = [modelled_gaia_parameters_dr2[3:4]; modelled_gaia_parameters_dr3[3:4]]
    ll += logpdf(dist_dr2_dr3, μ_dr2_dr3_modelled)


    # @show μ_dr2_dr3
    # resid = μ_dr2_dr3_modelled .- μ_dr2_dr3
    # @show resid
    # @show σ_dr2_dr3
    # @show μ_dr2_dr3 μ_dr2_dr3_modelled
    Δt = (Octofitter.meta_gaia_DR3.ref_epoch_mjd - Octofitter.meta_gaia_DR2.ref_epoch_mjd)/Octofitter.julian_year
    
    # @show θ_system.dr2_systematic_Δra
    # @show θ_system.dr2_systematic_Δμ_ra/Δt
    return ll, (;
        modelled_gaia_parameters_dr3=modelled_gaia_parameters_dr3,
        modelled_gaia_parameters_dr2=modelled_gaia_parameters_dr2,#.-correction,
        
        pmra_dr3_model = modelled_gaia_parameters_dr3[3],
        pmdec_dr3_model = modelled_gaia_parameters_dr3[4],

        # pmra_dr2_model = modelled_gaia_parameters_dr2[3] - θ_system.dr2_systematic_Δμ_ra,
        # pmdec_dr2_model = modelled_gaia_parameters_dr2[4] - θ_system.dr2_systematic_Δμ_dec,
        # pmra_dr32_model=((
        #     modelled_gaia_parameters_dr3[1]-modelled_gaia_parameters_dr2[1]
        # )*60*60*1000*cosd((like.dr3.ra-like.dr2.ra)/2) + θ_system.dr2_systematic_Δra)/Δt,
        # pmdec_dr32_model=((
        #     modelled_gaia_parameters_dr3[2]-modelled_gaia_parameters_dr2[2]
        # )*60*60*1000 + θ_system.dr2_systematic_Δdec)/Δt,

        pmra_dr2_model = modelled_gaia_parameters_dr2[3],
        pmdec_dr2_model = modelled_gaia_parameters_dr2[4],
        pmra_dr32_model=((
            modelled_gaia_parameters_dr3[1]-modelled_gaia_parameters_dr2[1]
        )*60*60*1000*cosd((like.dr3.ra-like.dr2.ra)/2))/Δt,
        pmdec_dr32_model=((
            modelled_gaia_parameters_dr3[2]-modelled_gaia_parameters_dr2[2]
        )*60*60*1000)/Δt,
        correction,

        # residual_pmra_dr2_model = μ_dr2_dr3[3] - μ_dr2_dr3_modelled[3],
        # residual_pmdec_dr2_model =  μ_dr2_dr3[4] - μ_dr2_dr3_modelled[4],
        # residual_pmra_dr3_model = μ_dr2_dr3[7] - μ_dr2_dr3_modelled[7],
        # residual_pmdec_dr3_model =  μ_dr2_dr3[8] - μ_dr2_dr3_modelled[8],
        # residual_pmra_dr32_model= (μ_dr2_dr3[1] - μ_dr2_dr3_modelled[1])*60*60*1000*cosd((like.dr3.ra-like.dr2.ra)/2)/Δt,
        # residual_pmdec_dr32_model= (μ_dr2_dr3[2] - μ_dr2_dr3_modelled[2])*60*60*1000/Δt,

    )
end