

function ln_like_pma(θ_system, pma::ProperMotionAnom)
    ll = 0.0
    for i in eachindex(pma.ra_epoch, pma.dec_epoch)
        pm_ra_star = 0.0
        pm_dec_star = 0.0
        
        # The model can support multiple planets
        for key in keys(θ_system.planets)
            θ_planet = θ_system.planets[key]
            # TODO: we are creating these from scratch for each observation instead of sharing them
            elements = construct_elements(θ_system, θ_planet)

            # RA and dec epochs are usually slightly different
            pm_ra_star += propmotionanom(elements, pma.ra_epoch[i], elements.μ, θ_planet.mass)[1]
            pm_dec_star += propmotionanom(elements, pma.dec_epoch[i], elements.μ, θ_planet.mass)[2]
        end

        residx = pm_ra_star - pma.pm_ra[i]
        residy = pm_dec_star - pma.pm_dec[i]
        σ²x = pma.σ_pm_ra[i]^2
        σ²y = pma.σ_pm_dec[i]^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π * σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π * σ²y))

        ll += χ²x + χ²y
    end

    return ll
end

# TODO: image modelling for multi planet systems do not consider how "removing" one planet
# might increase the contrast of another.
function ln_like_images(θ_system, system)
    ll = 0.0
    for key in keys(θ_system.planets)
        θ_planet = θ_system.planets[key]
        elements = construct_elements(θ_system, θ_planet)

        if (elements.a <= 0 ||
            elements.e < 0 ||
            elements.plx < 0 ||
            elements.μ <= 0)
            ll += NaN
            continue
        end


        ll += ln_like_images_element(elements, θ_planet, system)
    end

    # # Connect the flux at each epoch to an overall flux in this band for this planet
    # # fᵢ = θ_band.epochs
    # # ll += -1/2 * sum(
    # #     (fᵢ .- θ_band.f).^2
    # # ) / (θ_band.σ_f² * mean(fᵢ)^2)

    # And connect that flux to a modelled Teff and mass
    # f_model = model_interpolator(θ_planet.Teff, θ_planet.mass)
    # ll += -1/2 * (f_model - θ_band)^2 /  (θ_planet.σ_f_model² * f_model^2)

    return ll
end

function ln_like_images_element(elements::DirectOrbits.AbstractElements, θ_planet, system)
    images = system.images
    T = eltype(θ_planet)
    ll = zero(T)
    for i in eachindex(images.epoch)
       
        # Calculate position at this epoch
        ra, dec = kep2cart(elements, images.epoch[i])
        x = -ra
        y = dec

        # Get the photometry in this image at that location
        # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
        f̃ₓ = lookup_coord(images.image[i], (x, y), images.platescale[i])

        # Find the uncertainty in that photometry value (i.e. the contrast)
        r = √(x^2 + y^2)
        σₓ = images.contrast[i](r / images.platescale[i])

        # When we get a position that falls outside of our available
        # data (e.g. under the coronagraph) we cannot say anything
        # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
        # of zero.
        if !isfinite(σₓ) || !isfinite(f̃ₓ)
            # if typeof(σₓ) <: AbstractFloat
                # println("$x $y")
                # println("$elements: $(round(x)) $(round(y)) $f̃ₓ $σₓ")
            # end
            continue
        end

        # Ruffio et al 2017, eqn 31
        # ll += -1/(2σₓ^2) * (θ_epoch_f^2 - 2θ_epoch_f*f̃ₓ)
        # ll += -1/(2σₓ^2) * (θ_band^2 - 2θ_band*f̃ₓ)

        band = images.band[i]
        # If we are doing inference directly on the photometry
        f_band::T = 0
        if hasproperty(θ_planet, band)
            f_band = getproperty(θ_planet, band)
        elseif hasproperty(θ_planet, :mass)
            mass = θ_planet.mass
            if hasproperty(system.models, band)
                f_band = system.models[band](mass)
            else
                error("Mass prior specified, but photometry prior for $band was not, and no matching atmosphere model found.")
            end
        else
            error("No photometry prior for the band $band was specified, and neither was mass.")
        end

        σₓ² = σₓ^2
        l = -1 / (2σₓ²) * (f_band^2 - 2f_band * f̃ₓ) - log(sqrt(2π * σₓ²))
        ll += l
    end

    return ll
end

# If there is no astrometry atached to the planet, it does not contribute anything to the likelihood function
function ln_like_astrom(θ, θ_planet, planet::Planet{<:Any,<:Any,Nothing})
    return 0.0
end

# Astrometry
function ln_like_astrom(θ_system, θ_planet, planet::Planet{<:Any,<:Any,<:Astrometry})
    ll = 0.0
    
    # TODO, we are creating these from scratch for each observation instead of sharing them
    elements = construct_elements(θ_system, θ_planet)
    astrom = astrometry(planet)
    for i in eachindex(astrom.epoch)

        x, y = kep2cart(elements, astrom.epoch[i])
        residx = astrom.ra[i] - x
        residy = astrom.dec[i] - y
        σ²x = astrom.σ_ra[i]^2
        σ²y = astrom.σ_dec[i]^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π * σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π * σ²y))
        ll += χ²x + χ²y
    end
    return ll
end


function ln_like(θ, system::System)
    ll = 0.0

    # Hierarchical parameters over multiple planets
    # if haskey(system.priors.priors, :σ_i²)
    #     # If the sampler wanders into negative variances, return early to prevent
    #     # taking square roots of negative values later on
    #     if θ.σ_i² < 0
    #         return -Inf
    #     end

    #     # hierarchical priors here
    #     sum_iᵢ = zero(θ.i)
    #     sum_iᵢθi² = zero(θ.i)
    #     for θ_planet in θ.planets
    #         sum_iᵢ += θ_planet.i
    #         sum_iᵢθi² += (θ_planet.i .- θ.i)^2
    #     end
    #     ll += -1/2 * sum_iᵢθi² / θ.σ_i²  - log(sqrt(2π * θ.σ_i²))
    # end
    # if haskey(system.priors.priors, :σ_Ω²)
    #     # If the sampler wanders into negative variances, return early to prevent
    #     # taking square roots of negative values later on
    #     if θ.σ_Ω² < 0
    #         return -Inf
    #     end

    #     # hierarchical priors here
    #     sum_Ωᵢ = zero(θ.Ω)
    #     sum_ΩᵢθΩ² = zero(θ.Ω)
    #     for θ_planet in θ.planets
    #         _, Ωᵢ, _  = get_ωΩτ(θ_system, θ_planet)
    #         sum_Ωᵢ += Ωᵢ
    #         sum_ΩᵢθΩ² += (Ωᵢ .- θ.Ω)^2
    #     end
    #     ll += -1/2 * sum_ΩᵢθΩ² / θ.σ_Ω²  - log(sqrt(2π * θ.σ_Ω²))
    # end

    if !isnothing(system.images)
        ll += ln_like_images(θ, system)
    end
    if !isnothing(system.propermotionanom)
        ll += ln_like_pma(θ, system.propermotionanom)
    end
    # for (θ_planet, planet) in zip(θ.planets, system.planets)
        # ll += ln_like_astrom(θ, θ_planet, planet)
        
    for key in eachindex(system.planets)
        ll += ln_like_astrom(θ, θ.planets[key], system.planets[key])
    end

    
    return ll
end








# # This is a straight forward implementation that unfortunately is not type stable.
# # This is because we are looping over a heterogeneous container
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

# This implementation is ~5x faster for the same result.
# It uses meta-programming to unroll the loop over the different
# prior types. Note that this also ensures the priors can't be modified
# after building the ln_prior function.
# function make_ln_prior(θ)

#     body = Expr[]

#     for i in eachindex(θ)
#         pd = θ[i]
#         key = keys(θ)[i]
#         ex = :(lp += $logpdf($pd, params.$key))
#         push!(body, ex)
#     end

#     ex = :(function (params)
#         lp = zero(eltype(params))
#         $(body...)
#         return lp
#     end)

#     ln_prior = @RuntimeGeneratedFunction(ex)
#     return ln_prior
# end


function make_ln_prior(system::System)

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
    return @RuntimeGeneratedFunction(:(function (arr)
        lp = zero(first(arr))
        # Add contributions from planet priors
        $(prior_evaluations...)
        return lp
    end))
end


# function ln_prior(θ, system::System)
#     lp = 0.0
#     lp += system.priors.ln_prior(θ)
#     # for (planet, θ_planet) in zip(system.planets, θ.planets)
#     # for (planet, θ_planet) in zip(Tuple(system.planets), Tuple(θ.planets))
#     for key in keys(system.planets)
#         lp += system.planets[key].priors.ln_prior(θ.planets[key])
#     end

#     return lp
# end

function ln_post(θ, system::System)
    return ln_prior(θ, system) + ln_like(θ, system)
end


"""
    make_arr2nt(system::System)

Returns a function that maps an array of parameter values (e.g. sampled from priors,
sampled from posterior, dual numbers, etc.) to a nested named tuple structure suitable
for passing into `ln_like`.
In the process, this evaluates all Deterministic variables defined in the model.

Example:
```julia
julia> arr2nt = DirectDetections.make_arr2nt(system);
julia> params = DirectDetections.sample_priors(system)
9-element Vector{Float64}:
  1.581678216418196
 29.019567489624023
  1.805551918844831
  4.527607604709967
 -2.4432825071288837
 10.165695048003027
 -2.8667829383306063
  0.0006589349005787781
 -1.2600092725864576
julia> arr2nt(params)
(μ = 1.581678216418196, plx = 29.019567489624023, planets = (B = (τ = 1.805551918844831, ω = 4.527607604709967, i = -2.4432825071288837, Ω = 10.165695048003027, loge = -2.8667829383306063, loga = 0.0006589349005787781, logm = -1.2600092725864576, e = 0.001358992505357737, a = 1.0015184052910453, mass = 0.054952914078000334),))
```
"""
function make_arr2nt(system::System)

    # Roll flat vector of variables e.g. drawn from priors or posterior
    # into a nested NamedTuple structure for e.g. ln_like
    i = 0
    body_sys_priors = Expr[]

    # Priors
    for key in keys(system.priors.priors)
        i += 1
        ex = :(
            $key = arr[$i]
        )
        push!(body_sys_priors,ex)
    end

    # Deterministic variables
    body_sys_determ = Expr[]
    if !isnothing(system.deterministic)
        for (key,func) in zip(keys(system.deterministic.variables), values(system.deterministic.variables))
            ex = :(
                $key = $func(sys)
            )
            push!(body_sys_determ,ex)
        end
    end

    # Planets: priors & deterministic variables
    body_planets = Expr[]
    for planet in system.planets
        
        # Priors
        body_planet_priors = Expr[]
        for key in keys(planet.priors.priors)
            i += 1
            ex = :(
                $key = arr[$i]
            )
            push!(body_planet_priors,ex)
        end

        if isnothing(planet.deterministic)
            ex = :(
                $(planet.name) = (;
                    $(body_planet_priors...)
                )
            )
        # Resolve deterministic vars.
        else
            body_planet_determ = Expr[]
            for (key,func) in zip(keys(planet.deterministic.variables), values(planet.deterministic.variables))
                ex = :(
                    $key = $func(sys_res, (;$(body_planet_priors...)))
                )
                push!(body_planet_determ,ex)
            end

            ex = :(
                $(planet.name) = (;
                    $(body_planet_priors...),
                    $(body_planet_determ...)
                )
            )
        end
        push!(body_planets,ex)
    end

    # Here is the function we return.
    # It maps an array of parameters into our nested named tuple structure
    func = @RuntimeGeneratedFunction(:(function (arr)
        # Expand system variables from priors
        sys = (;$(body_sys_priors...))
        # Resolve deterministic system variables
        sys_res = (;
            sys...,
            $(body_sys_determ...)
        )
        # Get resolved planets
        pln = (;$(body_planets...))
        # Merge planets into resolved system
        sys_res_pln = (;sys..., planets=pln)
        return sys_res_pln
    end))

    return func
end


