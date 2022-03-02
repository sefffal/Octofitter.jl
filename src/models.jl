

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
    ll = 0.0

    for i in eachindex(rv.table.epoch)
        rv_star = 0.0
        for j in eachindex(elements)
            θ_planet = θ_system.planets[j]
            orbit = elements[j]

            if θ_planet.mass < 0
                return -Inf 
            end

            rv_star += radvel(orbit, rv.table.epoch[i], θ_planet.mass*mjup2msol)
        end
        resid = rv_star + θ_system.rv - rv.table.rv[i]
        σ² = rv.table.σ_rv[i]^2 + θ_system.jitter^2
        χ² = -0.5resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end

    return ll
end

"""
Likelihood of there being planets in a sequence of images.
"""
function ln_like(images::Images, θ_system, θ_planet)
    
    # Resolve the combination of system and planet parameters
    # as a KeplerianElements object. This pre-computes
    # some factors used in various calculations.
    elements = construct_elements(θ_system, θ_planet)

    imgtable = images.table
    T = eltype(θ_planet)
    ll = zero(T)
    for i in eachindex(imgtable.epoch)
       
        # Calculate position at this epoch
        o = orbitsolve(elements, imgtable.epoch[i])
        # x must be negated to go from sky coordinates (ra increasing to left) to image coordinates (ra increasing to right).
        x = -raoff(o)
        y = decoff(o)

        # Get the photometry in this image at that location
        # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
        f̃ₓ = lookup_coord(imgtable.image[i], (x, y), imgtable.platescale[i])

        # Find the uncertainty in that photometry value (i.e. the contrast)
        r = √(x^2 + y^2)
        σₓ = imgtable.contrast[i](r / imgtable.platescale[i])

        # When we get a position that falls outside of our available
        # data (e.g. under the coronagraph) we cannot say anything
        # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
        # of zero.
        if !isfinite(σₓ) || !isfinite(f̃ₓ)
            continue
        end

        band = imgtable.band[i]

        # Verify the user has specified a prior or model for this band.
        if !hasproperty(θ_planet, band)
            error("No photometry prior for the band $band was specified, and neither was mass.")
        end
        # TODO: verify this is type stable
        f_band = getproperty(θ_planet, band)
        # Direct imaging likelihood.
        # Notes: we are assuming that the different images fed in are not correlated.
        # The general multivariate Gaussian likleihood is exp(-1/2 (x⃗-μ⃗)ᵀ𝚺⁻¹(x⃗-μ⃗)) + √((2π)ᵏ|𝚺|)
        # Because the images are uncorrelated, 𝚺 is diagonal and we can separate the equation
        # into a a product of univariate Gaussian likelihoods or sum of log-likelihoods.
        # That term for each image is given below.

        # Ruffio et al 2017, eqn (31)
        # Mawet et al 2019, eqn (8)

        σₓ² = σₓ^2
        ll += -1 / (2σₓ²) * (f_band^2 - 2f_band * f̃ₓ) 
    end

    return ll
end

# Astrometry
function ln_like(astrom::Astrometry, θ_planet, orbit, all_θ_planets, all_orbits,)
# function ln_like(astrom::Astrometry, θ_planet, orbit, all_orbits)
    ll = 0.0
    for i in eachindex(astrom.table.epoch)
        # o = orbitsolve(orbit, astrom.table.epoch[i])

        # Astrometry is measured relative to the star.
        # Account for the relative position of the star due to
        # this planet
        # if hasproperty(θ_planet, :mass)
        #     star_δra = -raoff(o) * θ_planet.mass*mjup2msol/orbit.M
        #     star_δdec = -decoff(o) * θ_planet.mass*mjup2msol/orbit.M
        # else
        #     star_δra =  0
        #     star_δdec = 0
        # end

        # Astrometry is measured relative to the star.
        # Account for the relative position of the star due to
        # *all* planets
        local o
        star_δra =  0
        star_δdec = 0
        for (θ_planet_i, orbit_i) in zip(all_θ_planets, all_orbits)
            o_i = orbitsolve(orbit, astrom.table.epoch[i])
            if orbit_i == orbit
                o = o_i
            end
            star_δra = -raoff(o) * θ_planet_i.mass * mjup2msol / orbit.M
            star_δdec = -decoff(o) * θ_planet_i.mass * mjup2msol / orbit.M
        end

        # # PA and Sep specified
        # if haskey(astrom.table, :pa) && haskey(astrom.table, :ρ)
        #     ρ = projectedseparation(o)
        #     pa = posangle(o)
        #     resid1 = rem2pi(astrom.table.pa[i] - pa, RoundNearest) # TODO: confirm this works at the wrap point
        #     resid2 = astrom.table.ρ[i] - ρ
        #     σ²1 = astrom.table.σ_pa[i ]^2
        #     σ²2 = astrom.table.σ_ρ[i]^2
        # # RA and DEC specified
        # else

        x = raoff(o)
        y = decoff(o)
        resid1 = astrom.table.ra[i] - star_δra - x
        resid2 = astrom.table.dec[i] - star_δdec - y
        σ²1 = astrom.table.σ_ra[i ]^2
        σ²2 = astrom.table.σ_dec[i]^2

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
        σ_phot = photometry.table.σ_phot[i]
        resid = phot_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end

# Overall log likelihood of the system given the parameters θ_system
function ln_like(system::System, θ_system)
    # Take some care to ensure type stability when using e.g. ForwardDiffJ
    ll = zero(typeof(first(θ_system)))

    # Fail fast if we have a negative stellar mass.
    # Users should endeavour to use priors on e.g. stellar mass
    # that are strictly positive, otherwise we are reduced to rejection sampling!
    if hasproperty(θ_system, :M) && θ_system.M <= 0
        return oftype(ll, -Inf)
    end

    # Go through each planet in the model and add its contribution
    # to the ln-likelihood.
    out_of_bounds = Ref(false)
    elements = map(eachindex(system.planets)) do i
        # planet = system.planets[i]
        θ_planet = θ_system.planets[i]

        # Like negative stellar mass, users should use priors with supports
        # that do not include these invalid values. But if we see them,
        # give zero likelihood right away instead of an inscrutable error
        # from some code expecting these invariants to hold.
        if (hasproperty(θ_planet, :a) && θ_planet.a <= 0) ||
            (hasproperty(θ_planet, :e) && !(0 <= θ_planet.e < 1))
            out_of_bounds[] = true
        end

        # Resolve the combination of system and planet parameters
        # as a KeplerianElements object. This pre-computes
        # some factors used in various calculations.
        kep_elements = construct_elements(θ_system, θ_planet)

        return kep_elements
    end
    # Fail fast if out of bounds for one of the planets
    if out_of_bounds[]
        return oftype(ll, -Inf) # Ensure type stability
    end

    # Loop through the planets from the outside in. 
    # Try to do this sorting in a non-allocating way.
    # This way we have the option to account for each planets influence on the outer planets
    sma = map(elements) do elem
        return elem.a
    end
    planet_sma_asc_ii = sortperm(SVector(sma))

    # planet_sma_asc_ii = sortperm(SVector(sma))

    # Handle all observations attached to planets in order of semi-major axis
    for j in eachindex(planet_sma_asc_ii)
        i = planet_sma_asc_ii[j]
        # Planet model
        planet = system.planets[i]
        # Parameters specific to this planet
        θ_planet = θ_system.planets[i]
        # Cached KeplerianElements with precomputed factors, etc.
        planet_elements = elements[i]
        # kep_elements, but for all planets interior to this one (given the current parameters)
        # interior_planets = kep_elements[begin:min(end,i)]

        # Loop through observations
        for obs in planet.observations
            ll += ln_like(obs, θ_planet, planet_elements, θ_system.planets, elements)
        end
    end

    if !isfinite(ll[])
        return ll[]
    end

    # Loop through and add contribution of all observation types associated with this system as a whole
    for obs in system.observations
        ll += ln_like(obs, θ_system, elements)
    end


    return ll
end







# This is a straight forward implementation that unfortunately is not type stable.
# This is because we are looping over a heterogeneous container
# function make_ln_prior(priors)
#     return function ln_prior(params)
#         lp = zero(first(params))
#         for i in eachindex(params)
#             pd = priors[i]
#             param = params[i]
#             lp += logpdf(pd, param)
#         end
#         return lp 
#     end
# end

function make_ln_prior(system::System)

    # This function uses meta-programming to unroll all the code at compile time.
    # This is a big performance win, since it avoids looping over all the different
    # types of distributions that might be specified as priors.
    # Otherwise we would have to loop through an abstract vector and do runtime dispatch!
    # This way all the code gets inlined into a single tight numberical function in most cases.

    i = 0
    prior_evaluations = Expr[]

    # System priors
    for prior_distribution in values(system.priors.priors)
        i += 1
        ex = :(
            lp += $logpdf($prior_distribution, arr[$i])
        )
        push!(prior_evaluations,ex)
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            i += 1
            # Work around for Beta distributions.
            # Outside of the range [0,1] logpdf returns -Inf.
            # This works fine, but AutoDiff outside this range causes a DomainError.
            if typeof(prior_distribution) <: Beta
                ex = :(
                    lp += 0 <= arr[$i] < 1 ? $logpdf($prior_distribution, arr[$i]) : -Inf
                )
            else
                ex = :(
                    lp += $logpdf($prior_distribution, arr[$i])
                )
            end
            push!(prior_evaluations,ex)
        end
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        lp = zero(first(arr))
        # Add contributions from planet priors
        @inbounds begin
           $(prior_evaluations...) 
        end
        return lp
    end))
end

# Same as above, but assumes the input to the log prior was sampled
# using transformed distributions from Bijectors.jl
# Uses logpdf_with_trans() instead of logpdf to make the necessary corrections.
function make_ln_prior_transformed(system::System)

    i = 0
    prior_evaluations = Expr[]

    # System priors
    for prior_distribution in values(system.priors.priors)
        i += 1
        ex = :(
            lp += $logpdf_with_trans($prior_distribution, arr[$i], true)
        )
        push!(prior_evaluations,ex)
    end

    # Planet priors
    for planet in system.planets
        # for prior_distribution in values(planet.priors.priors)
        for (key, prior_distribution) in zip(keys(planet.priors.priors), values(planet.priors.priors))
            i += 1
            ex = :(
                lp += $logpdf_with_trans($prior_distribution, arr[$i], true)
            )
            push!(prior_evaluations,ex)
        end
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        lp = zero(first(arr))
        # Add unrolled prior evaluations
        @inbounds begin
           $(prior_evaluations...) 
        end
        return lp
    end))
end


# # Replaces `θ = Bijectors.invlink.(priors_vec, θ_t)` with a type stable
# # unrolled version.
# function make_Bijector_invlinkvec(priors_vec)

#     i = 0
#     parameter_transformations = Expr[]

#     # System priors
#     for prior_distribution in priors_vec
#         i += 1
#         ex = :(
#             theta_out[$i] = $(Bijectors.invlink)($prior_distribution, arr[$i])
#         )
#         push!(parameter_transformations, ex)
#     end

#     # Here is the function we return.
#     # It maps an array of parameters into our nested named tuple structure
#     # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
#     # The RuntimeGeneratedFunctions package avoids these in all cases.
#     return @RuntimeGeneratedFunction(:(function (arr)
#         l = $i
#         theta_out = @MVector zeros(eltype(arr), l)
#         # theta_out = zeros(eltype(arr), l)
#         @boundscheck if length(arr) != l
#             error("Expected exactly $l elements in array (got $(length(arr)))")
#         end
#         # Add unrolled parameter transformations to fill theta_out
#         @inbounds begin
#            $(parameter_transformations...) 
#         end
#         return theta_out
#     end))
# end


# Replaces `θ = Bijectors.invlink.(priors_vec, θ_t)` with a type stable
# unrolled version.
function make_Bijector_invlinkvec(priors_vec)

    i = 0
    parameter_transformations = Expr[]

    # System priors
    for prior_distribution in priors_vec
        i += 1
        ex = :(
            $(Bijectors.invlink)($prior_distribution, arr[$i])
        )
        push!(parameter_transformations, ex)
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    # Note: eval() would normally work fine here, but sometimes we can hit "world age problemms"
    # The RuntimeGeneratedFunctions package avoids these in all cases.
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        # theta_out = zeros(eltype(arr), l)
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        # Add unrolled parameter transformations to fill theta_out
        @inbounds begin
            # theta_out = SVector{l,eltype(arr)}(
            # theta_out = MVector{l,eltype(arr)}(
            theta_out = tuple(
                $(parameter_transformations...) 
            )
        end
        return theta_out
    end))
end


function ln_post(θ, system::System)
    return ln_prior(θ, system) + ln_like(θ, system)
end
