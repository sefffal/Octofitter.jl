

"""
General proper motion likelihood at any number of epochs.
Each epoch is averaged over 5 measurements at +-dt/2.
"""
function ln_like(pma::ProperMotionAnom, θ_system, elements)
    ll = 0.0

    # How many points over Δt should we average the proper motion at each
    # epoch? This is because the PM is not an instantaneous measurement.
    N_ave = 5
    
    for i in eachindex(pma.table.ra_epoch, pma.table.dec_epoch)
        pmra_star = 0.0
        pmdec_star = 0.0
        
        # The model can support multiple planets
        # for key in keys(θ_system.planets)
        for j in eachindex(elements)
            θ_planet = θ_system.planets[j]
            orbit = elements[j]

            if θ_planet.mass < 0
                return -Inf
            end

            # Average multiple observations over a timescale +- dt
            # to approximate what HIPPARCOS and GAIA would have measured.
            for δt = range(-pma.table.dt[i]/2, pma.table.dt[i]/2, N_ave)

                # RA and dec epochs are usually slightly different
                # Note the unit conversion here from jupiter masses to solar masses to 
                # make it the same unit as the stellar mass (element.mu)
                pmra_star += pmra(orbit, pma.table.ra_epoch[i]+δt, θ_planet.mass*mjup2msol)
                pmdec_star += pmdec(orbit, pma.table.dec_epoch[i]+δt, θ_planet.mass*mjup2msol)
            end

        end
        
        pmra_star/=N_ave
        pmdec_star/=N_ave

        residx = pmra_star + θ_system.pmra - pma.table.pmra[i]
        residy = pmdec_star + θ_system.pmdec - pma.table.pmdec[i]
        σ²x = pma.table.σ_pmra[i]^2
        σ²y = pma.table.σ_pmdec[i]^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π * σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π * σ²y))

        ll += χ²x + χ²y
    end

    return ll
end


"""
Specific HGCA proper motion modelling. Model the GAIA-Hipparcos/Δt proper motion
using 5 position measurements averaged at each of their epochs.
"""
function ln_like(pma::ProperMotionAnomHGCA, θ_system, elements)
    ll = 0.0

    # This observation type just wraps one row from the HGCA (see hgca.jl)
    hgca = pma.table
    # Roughly over what time period were the observations made?
    dt_gaia = 1038 # EDR3: days between  Date("2017-05-28") - Date("2014-07-25")
    dt_hip = 4*365
    # How many points over Δt should we average the proper motion and stellar position
    # at each epoch? This is because the PM is not an instantaneous measurement.
    N_ave = 25 #5

    # Look at the position of the star around both epochs to calculate 
    # our modelled delta-position proper motion

    # First epoch: Hipparcos
    ra_hip_model = 0.0
    dec_hip_model = 0.0
    pmra_hip_model = 0.0
    pmdec_hip_model = 0.0
    # The model can support multiple planets
    for i in eachindex(elements)
        θ_planet = θ_system.planets[i]
        orbit = elements[i]
        if θ_planet.mass < 0
            return -Inf
        end
        # Average multiple observations over a timescale +- dt/2
        # to approximate what HIPPARCOS would have measured.
        for δt = range(-dt_hip/2, dt_hip/2, N_ave)
            # RA and dec epochs are usually slightly different
            # Note the unit conversion here from jupiter masses to solar masses to 
            # make it the same unit as the stellar mass (element.mu)
            # TODO: we can't yet use the orbitsolve interface here for the pmra calls,
            # meaning we calculate the orbit 2x as much as we need.
            o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_hip[1])+δt)
            o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_hip[1])+δt)
            ra_hip_model += -raoff(o_ra) * θ_planet.mass*mjup2msol/orbit.M
            dec_hip_model += -decoff(o_dec) * θ_planet.mass*mjup2msol/orbit.M
            pmra_hip_model += pmra(o_ra, θ_planet.mass*mjup2msol)
            pmdec_hip_model += pmdec(o_dec, θ_planet.mass*mjup2msol)
        end
    end
    ra_hip_model/=N_ave
    dec_hip_model/=N_ave
    pmra_hip_model/=N_ave
    pmdec_hip_model/=N_ave

    # Last epoch: GAIA
    ra_gaia_model = 0.0
    dec_gaia_model = 0.0
    pmra_gaia_model = 0.0
    pmdec_gaia_model = 0.0
    # The model can support multiple planets
    for i in eachindex(elements)
        θ_planet = θ_system.planets[i]
        orbit = elements[i]
        if θ_planet.mass < 0
            return -Inf
        end
        # Average multiple observations over a timescale +- dt
        # to approximate what HIPPARCOS and GAIA would have measured.
        for δt = range(-dt_gaia/2, dt_gaia/2, N_ave)
            # RA and dec epochs are usually slightly different
            # Note the unit conversion here from jupiter masses to solar masses to 
            # make it the same unit as the stellar mass (element.M)
            o_ra = orbitsolve(orbit, years2mjd(hgca.epoch_ra_gaia[1])+δt)
            o_dec = orbitsolve(orbit, years2mjd(hgca.epoch_dec_gaia[1])+δt)
            ra_gaia_model += -raoff(o_ra) * θ_planet.mass*mjup2msol/orbit.M
            dec_gaia_model += -decoff(o_dec) * θ_planet.mass*mjup2msol/orbit.M
            pmra_gaia_model += pmra(o_ra, θ_planet.mass*mjup2msol)
            pmdec_gaia_model += pmdec(o_dec, θ_planet.mass*mjup2msol)
        end
    end
    ra_gaia_model/=N_ave
    dec_gaia_model/=N_ave
    pmra_gaia_model/=N_ave
    pmdec_gaia_model/=N_ave


    # Model the GAIA-Hipparcos delta-position velocity
    pmra_hg_model = (ra_gaia_model - ra_hip_model)/(years2mjd(hgca.epoch_ra_gaia[1]) - years2mjd(hgca.epoch_ra_hip[1]))
    pmdec_hg_model = (dec_gaia_model - dec_hip_model)/(years2mjd(hgca.epoch_dec_gaia[1]) - years2mjd(hgca.epoch_dec_hip[1]))

    # Compute the likelihood at all three epochs (Hipparcos, GAIA-Hip, GAIA)
    pmra_model = (pmra_hip_model, pmra_hg_model, pmra_gaia_model)
    pmdec_model = (pmdec_hip_model, pmdec_hg_model, pmdec_gaia_model)
    pmra_meas = (hgca.pmra_hip[1], hgca.pmra_hg[1], hgca.pmra_gaia[1])
    pmdec_meas = (hgca.pmdec_hip[1], hgca.pmdec_hg[1], hgca.pmdec_gaia[1])
    σ_pmra = (hgca.pmra_hip_error[1], hgca.pmra_hg_error[1], hgca.pmra_gaia_error[1])
    σ_pmdec = (hgca.pmdec_hip_error[1], hgca.pmdec_hg_error[1], hgca.pmdec_gaia_error[1])
    for i in 1:3
        residx = pmra_model[i] + θ_system.pmra - pmra_meas[i]
        residy = pmdec_model[i] + θ_system.pmdec - pmdec_meas[i]
        σ²x = σ_pmra[i]^2
        σ²y = σ_pmdec[i]^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π * σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π * σ²y))
        ll += χ²x + χ²y
    end

    return ll
end

"""
Radial velocity likelihood.
"""
function ln_like(rv::RadialVelocity, θ_system, elements)
    T = Float64
    ll = zero(T)

    for i in eachindex(rv.table.epoch)
        rv_star = zero(T)
        for j in eachindex(elements)
            θ_planet = θ_system.planets[j]
            orbit = elements[j]

            if θ_planet.mass < 0
                return -Inf 
            end

            rv_star += radvel(orbit, rv.table.epoch[i], θ_planet.mass*mjup2msol)
        end
        # Each measurement is tagged with a jitter and rv zero point variable.
        # We then query the system variables for them.
        # A previous implementation used symbols instead of indices but it was too slow.
        # barycentric_rv_inst = getproperty(θ_system, rv.table.rv0[i])::Float64
        # jitter_inst = getproperty(θ_system, rv.table.jitter[i])::Float64
        if !hasproperty(rv.table, :inst_idx)
            barycentric_rv_inst = θ_system.rv0
            jitter_inst = θ_system.jitter
        else
            inst_idx = rv.table.inst_idx[i]
            if inst_idx == 1
                barycentric_rv_inst = θ_system.rv0_1
                jitter_inst = θ_system.jitter_1
            elseif inst_idx == 2
                barycentric_rv_inst = θ_system.rv0_2
                jitter_inst = θ_system.jitter_2
            elseif inst_idx == 3
                barycentric_rv_inst = θ_system.rv0_3
                jitter_inst = θ_system.jitter_3
            elseif inst_idx == 4
                barycentric_rv_inst = θ_system.rv0_4
                jitter_inst = θ_system.jitter_4
            else
                error("More than four radial velocity instruments are not yet supported for performance reasons")
            end
        end
        resid = rv_star + barycentric_rv_inst - rv.table.rv[i]
        σ² = rv.table.σ_rv[i]^2 + jitter_inst^2
        χ² = -0.5resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end

    return ll
end

"""
Likelihood of there being planets in a sequence of images.
"""
function ln_like(images::Images, θ_system, all_elements)
    
    # Resolve the combination of system and planet parameters
    # as a VisualOrbit object. This pre-computes
    # some factors used in various calculations.
    # elements = construct_elements(θ_system, θ_planet)
    

    imgtable = images.table
    T = eltype(first(θ_system))
    ll = zero(T)
    for i in eachindex(imgtable.epoch)
        
        band = imgtable.band[i]

        # # Images are centered on the *star's* position. Large companions
        # # effectively shift the images around as they orbit the star.
        # # Account for the relative position of the star due to
        # # *all* planets
        star_δra =  0.0
        star_δdec = 0.0
        # for (θ_planet_i, orbit_i) in zip(θ_system.planets, all_elements)
        #     if hasproperty(θ_planet_i, :mass)
        #         o_i = orbitsolve(orbit_i, imgtable.epoch[i])
        #         star_δra += raoff(o_i, θ_planet_i.mass * mjup2msol)
        #         star_δdec += decoff(o_i, θ_planet_i.mass * mjup2msol)
        #     end
        # end

        # Once we have the star's reflex motion, go through and look
        # for each planet
        for (θ_planet_i, orbit_i) in zip(θ_system.planets, all_elements)
            # TODO: avoid calculating the orbit solution twice.
            # In 1.8 escape anlaysis might let us keep these in a list.
            o_i = orbitsolve(orbit_i, imgtable.epoch[i])

            # Note the x reversal between RA and image coordinates
            x = -(raoff(o_i) + star_δra)
            y = +(decoff(o_i) + star_δdec)

            # Get the photometry in this image at that location
            # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
            platescale = imgtable.platescale[i]
            f̃ₓ = imgtable.imageinterp[i](x/platescale, y/platescale)

            # Find the uncertainty in that photometry value (i.e. the contrast)
            if hasproperty(imgtable, :contrastmap)
                # If we have a 2D map
                σₓ = imgtable.contrastmapinterp[i](x/platescale, y/platescale)
            else
                # We have a 1D contrast curve
                r = √(x^2 + y^2)
                σₓ = imgtable.contrast[i](r / platescale)
            end

            # σₓ = θ_planet_i.σ

            # Verify the user has specified a prior or model for this band.
            if !hasproperty(θ_planet_i, band)
                error("No photometry variable for the band $band was specified.")
            end
            # TODO: verify this is type stable
            f_band = getproperty(θ_planet_i, band)


            # When we get a position that falls outside of our available
            # data (e.g. under the coronagraph) we cannot say anything
            # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
            # of zero.
            if !isfinite(σₓ) || !isfinite(f̃ₓ)
                continue
            end

            # Direct imaging likelihood.
            # Notes: we are assuming that the different images fed in are not correlated.
            # The general multivariate Gaussian likleihood is exp(-1/2 (x⃗-μ⃗)ᵀ𝚺⁻¹(x⃗-μ⃗)) + √((2π)ᵏ|𝚺|)
            # Because the images are uncorrelated, 𝚺 is diagonal and we can separate the equation
            # into a a product of univariate Gaussian likelihoods or sum of log-likelihoods.
            # That term for each image is given below.

            # Ruffio et al 2017, eqn (31)
            # Mawet et al 2019, eqn (8)

            σₓ² = σₓ^2
            ll += -1 / (2*σₓ²) * (f_band^2 - 2f_band * f̃ₓ)
        end
    end

    return ll
end

# Astrometry
function ln_like(astrom::Astrometry, θ_planet, orbit,)
    ll = 0.0
    for i in eachindex(astrom.table.epoch)

        star_δra =  0.
        star_δdec = 0.

        o = orbitsolve(orbit, astrom.table.epoch[i])
        # PA and Sep specified
        if hasproperty(astrom.table, :pa) && hasproperty(astrom.table, :sep)
            ρ = projectedseparation(o)
            pa = posangle(o)

            pa_diff = ( astrom.table.pa[i] - pa + π) % 2π - π;
            pa_diff = pa_diff < -π ? pa_diff + 2π : pa_diff;
            resid1 = pa_diff
            resid2 = astrom.table.sep[i] - ρ
            σ²1 = astrom.table.σ_pa[i ]^2
            σ²2 = astrom.table.σ_sep[i]^2
        # RA and DEC specified
        else
            x = raoff(o)# + star_δra
            y = decoff(o)# + star_δdec
            resid1 = astrom.table.ra[i] - x
            resid2 = astrom.table.dec[i] - y
            σ²1 = astrom.table.σ_ra[i ]^2
            σ²2 = astrom.table.σ_dec[i]^2
        end

        χ²1 = -(1/2)*resid1^2 / σ²1 - log(sqrt(2π * σ²1))
        χ²2 = -(1/2)*resid2^2 / σ²2 - log(sqrt(2π * σ²2))
        ll += χ²1 + χ²2
    end
    return ll
end

# Photometry
function ln_like(photometry::Photometry, θ_planet, _elements=nothing, _interior_planets=nothing)
    ll = 0.0
    for i in eachindex(photometry.table.band)
        band = photometry.table.band[i]
        phot_param = getproperty(θ_planet, band)
        phot_meas = photometry.table.phot[i]
        if !isfinite(phot_param)
            return -Inf
        end
        # Experimenting with fitting sigma phot
        σ_phot = photometry.table.σ_phot[i]
        # σₓ = θ_planet_i.σ

        resid = phot_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end

# # Overall log likelihood of the system given the parameters θ_system
# function ln_like(system::System, θ_system)
#     # Take some care to ensure type stability when using e.g. ForwardDiff
#     ll = zero(typeof(first(θ_system)))

#     # Fail fast if we have a negative stellar mass.
#     # Users should endeavour to use priors on e.g. stellar mass
#     # that are strictly positive, otherwise we are reduced to rejection sampling!
#     if hasproperty(θ_system, :M) && θ_system.M <= 0
#         return oftype(ll, -Inf)
#     end

#     # Go through each planet in the model and add its contribution
#     # to the ln-likelihood.
#     # out_of_bounds = Base.RefValue{Bool}(false)
#     elements = map(eachindex(system.planets)) do i
#         planet = system.planets[i]
#         θ_planet = θ_system.planets[i]

#         # Like negative stellar mass, users should use priors with supports
#         # that do not include these invalid values. But if we see them,
#         # give zero likelihood right away instead of an inscrutable error
#         # from some code expecting these invariants to hold.
#         if (hasproperty(θ_planet, :a) && θ_planet.a <= 0) ||
#             (hasproperty(θ_planet, :e) && !(0 <= θ_planet.e < 1))
#             out_of_bounds[] = true
#         end

#         # Resolve the combination of system and planet parameters
#         # as a orbit object. The type of orbitobject is stored in the 
#         # Planet type. This pre-computes some factors used in various calculations.
#         kep_elements = construct_elements(orbittype(planet), θ_system, θ_planet)

#         return kep_elements
#     end
#     # # Fail fast if out of bounds for one of the planets
#     # if out_of_bounds[]
#     #     return oftype(ll, -Inf) # Ensure type stability
#     # end

#     # Loop through the planets from the outside in. 
#     # Try to do this sorting in a non-allocating way.
#     # This way we have the option to account for each planets influence on the outer planets
#     # sma = map(elements) do elem
#     #     return elem.a
#     # end
#     # planet_sma_asc_ii = sortperm(SVector(sma))

#     # The above sorting is not currently used, so need to perform it.
#     planet_sma_asc_ii = 1:length(elements)


#     # Handle all observations attached to planets in order of semi-major axis
#     for j in eachindex(planet_sma_asc_ii)
#         i = planet_sma_asc_ii[j]
#         # Planet model
#         planet = system.planets[i]
#         # Parameters specific to this planet
#         θ_planet = θ_system.planets[i]
#         # Cached VisualOrbit with precomputed factors, etc.
#         planet_elements = elements[i]
#         # kep_elements, but for all planets interior to this one (given the current parameters)
#         # interior_planets = kep_elements[begin:min(end,i)]

#         # Loop through observations
#         for obs in planet.observations
#             ll += ln_like(obs, θ_planet, planet_elements, θ_system.planets, elements)
#         end
#     end

#     if !isfinite(ll)
#         return ll
#     end

#     # Loop through and add contribution of all observation types associated with this system as a whole
#     for obs in system.observations
#         ll += ln_like(obs, θ_system, elements)
#     end

#     return ll
# end


# # Overall log likelihood of the system given the parameters θ_system
# function ln_like(system::System, θ_system)
#     # Take some care to ensure type stability when using e.g. ForwardDiff
#     ll = zero(typeof(first(θ_system)))

#     # Should unroll this loop with RuntimeGeneratedFunction
#     # for i in eachindex(system.planets)
#     let i = 1
#         planet = system.planets[i]
#         θ_planet = θ_system.planets[i]

#         # Resolve the combination of system and planet parameters
#         # as a orbit object. The type of orbitobject is stored in the 
#         # Planet type. This pre-computes some factors used in various calculations.
#         planet_elements = construct_elements(orbittype(planet), θ_system, θ_planet)
#         # planet_elements = VisualOrbit((;
#         #     M=θ_system.M,
#         #     plx=θ_system.plx,
#         #     a=θ_planet.a,
#         #     i=θ_planet.i,
#         #     ω=θ_planet.ω,
#         #     Ω=θ_planet.Ω,
#         #     e=θ_planet.e,
#         #     τ=θ_planet.τ,
#         # ))
        
#         # Loop through observations ( could unroll this loop with RuntimeGeneratedFunction)
#         for obs in planet.observations
#             ll += ln_like(obs, θ_planet, planet_elements)
#         end
#     end
#     # Could unroll this loop

#     # # Loop through and add contribution of all observation types associated with this system as a whole
#     # for obs in system.observations
#     #     ll += ln_like(obs, θ_system)
#     # end

#     return ll
# end

function make_ln_like(system::System, θ_system)

    planet_exprs = Expr[]
    planet_keys = Symbol[]
    for i in eachindex(system.planets)
        planet = system.planets[i]
        θ_planet = θ_system.planets[i]
        key = Symbol("planet_$i")

        likelihood_exprs = map(eachindex(planet.observations)) do obs
            :(
                ll += ln_like(
                    system.planets[$(Meta.quot(i))].observations[$obs],
                    θ_system.planets.$i,
                    $(key)
                )
            )
        end

        planet_contruction = quote
            $key = $(construct_elements)($(orbittype(planet)), θ_system, θ_system.planets.$i)
            $(likelihood_exprs...)
        end
        push!(planet_exprs, planet_contruction)
        push!(planet_keys, key)
        break
    end


    sys_exprs = map(system.observations) do obs
        :(ll += ln_like($obs, θ_system, elems))
    end

    return @RuntimeGeneratedFunction(:(function (system::System, θ_system)

        ll = zero(first(θ_system))

        # Construct all orbit elements and evaluate all their individual observation likelihoods
        $(planet_exprs...)

        # Construct a tuple of existing planet orbital elements
        elems = tuple($(planet_keys...))
        
        $(sys_exprs...)

        return ll
    end))


end
