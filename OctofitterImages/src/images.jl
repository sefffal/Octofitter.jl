
const images_cols = (:image, :epoch, :platescale,)

"""
    ImageObs(
        table,
        name="images",
        variables=@variables begin
        end
    )

A block of images of a system. Pass a Table with the following columns:
$images_cols

For example:
```julia
image_dat = Table(;
    epoch = [1234.0, 1584.7],
    image = [readfits("img1.fits"), readfits("img2.fits")],
    platescale = [19.4, 19.4]
)

ImageObs(
    image_dat,
    name="SPHERE",
    variables=@variables begin
        flux ~ Normal(3.8, 0.5)        # Planet flux in image units
        platescale = 1.0               # Platescale multiplier [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
        northangle = 0.0               # North angle offset in radians [could use: northangle ~ Normal(0, deg2rad(1))]
    end
)
```
Contrast can be a function that returns the 1 sigma contrast of the image from a separation in mas to the same units as the image file.
Or, simply leave it out and it will be calculated for you.
Epoch is in MJD.
Platescale is in mas/px.
"""
struct ImageObs{TTable<:Table} <: Octofitter.AbstractObs
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
    function ImageObs(
        table;
        name::String="images",
        variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin end)
    )
        (priors, derived) = variables
        # Fallback to calculating contrast automatically
        if !in(:contrast, columnnames(table)) && !in(:contrastmap, columnnames(table))
            @info "Measuring contrast from image"
            contrast = contrast_interp.(table.image)
            table = Table(table, contrast=contrast)
        end
        if !issubset(images_cols, columnnames(table))
            error("Expected columns $images_cols")
        end
        # Create linear interpolators over the input images
        imageinterp = map(table.image) do img
            LinearInterpolation(parent.(dims(img)), img, extrapolation_bc=convert(eltype(img), NaN))
        end
        table = Table(table; imageinterp)
        if hasproperty(table, :contrastmap)
            # Create linear interpolators over the input contrastmaps
            contrastmapinterp = map(table.contrastmap) do img
                LinearInterpolation(parent.(dims(img)), img, extrapolation_bc=convert(eltype(img), NaN))
            end
            table = Table(table; contrastmapinterp)
        end
        return new{typeof(table)}(table, priors, derived, name)
    end
end
# Legacy constructor for backward compatibility
ImageObs(observations::NamedTuple...; kwargs...) = ImageObs(Table(observations...); kwargs...)

# Backwards compatibility alias
const ImageLikelihood = ImageObs

export ImageObs, ImageLikelihood


function Octofitter.likeobj_from_epoch_subset(obs::ImageObs, obs_inds)
    return ImageObs(obs.table[obs_inds,:,1]...)
end


"""
    contrast_interp(image; step=2)

Returns a linear interpolation on top of the results from `contrast`.
Extrapolated results return Inf.
"""
function contrast_interp(image::AstroImage; step=2)
    cont = contrast(image; step)
    mask = findfirst(isfinite, cont.contrast):findlast(isfinite, cont.contrast)
    return LinearInterpolation(cont.separation[mask], cont.contrast[mask], extrapolation_bc=Flat())
end


"""
    contrast(image; step=2)

Measure the contrast of an image, in the sense of high contrast imaging.
That is, divide the image into annuli moving outwards from the centre
(index 0,0 if offset image) and calculate the standard deviation in 
each.

Returns a vector of annulus locations in pixels and a vector of standard
deviations.

*NOTE* This is the 1σ contrast. Multiply by five to get the usual confidence
value.
"""
function contrast(image::AstroImage; step=2)
    dx = dims(image,X)
    dy = collect(dims(image,Y))
    dr = sqrt.(
        dx.^2 .+ (dy').^2
    )

    c_img = collect(image)
    
    bins = 0:step:maximum(dr)
    # bins = 30:step:100
    contrast = zeros(size(bins))
    mask = falses(size(image))
    mask2 = isfinite.(c_img)
    for i in eachindex(bins)
        bin = bins[i]
        mask .= (bin.-step/2) .< dr .< (bin.+step/2) 
        mask .&= mask2
        c = std(view(c_img, mask))
        contrast[i] = c
    end

    return (;separation=bins, contrast)
end


function imgsep(image::AstroImage)
    dx = dims(image,X)
    dy = collect(dims(image,Y))
    dr = sqrt.(
        dx.^2 .+ (dy').^2
    )
    return dr
end
    


"""
Likelihood of there being planets in a sequence of images.
"""
function Octofitter.ln_like(images::ImageObs, ctx::Octofitter.PlanetObservationContext)
    (; θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start) = ctx
    
    # Resolve the combination of system and planet parameters
    # as a Visual{KepOrbit} object. This pre-computes
    # some factors used in various calculations.
    this_orbit = orbits[i_planet]

    imgtable = images.table
    T = Octofitter._system_number_type(θ_planet)
    ll = zero(T)
    
    # Get calibration parameters from observation variables
    flux = θ_obs.flux
    platescale_multiplier = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : one(T)
    northangle_offset = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : zero(T)
    
    for i_epoch in eachindex(imgtable.epoch)

        # Account for any inner planets
        ra_host_perturbation = zero(T)
        dec_host_perturbation = zero(T)
        for i_other_planet in eachindex(orbits)
            orbit_other = orbits[i_other_planet]
            # Only account for inner planets
            if semimajoraxis(orbit_other) < semimajoraxis(this_orbit)
                θ_planet′ = θ_system.planets[i_other_planet]
                if !hasproperty(θ_planet′, :mass)
                    continue
                end
                mass_other = θ_planet′.mass*Octofitter.mjup2msol
                sol′ = orbit_solutions[i_other_planet][i_epoch + orbit_solutions_i_epoch_start]
                # Note about `total mass`: for this to be correct, user will have to specify
                # `M` at the planet level such that it doesn't include the outer planets.
                ra_host_perturbation += raoff(sol′, mass_other)
                dec_host_perturbation += decoff(sol′, mass_other)
            end
        end

        # Take the measurement, and *add* the Delta, to get what we compare to the model
        sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]

        # Apply north angle rotation and platescale correction
        ra_raw = raoff(sol) - ra_host_perturbation
        dec_raw = decoff(sol) - dec_host_perturbation
        
        # Apply north angle rotation
        cos_θ = cos(northangle_offset)
        sin_θ = sin(northangle_offset)
        ra_rotated = ra_raw * cos_θ - dec_raw * sin_θ
        dec_rotated = ra_raw * sin_θ + dec_raw * cos_θ

        # Note the x reversal between RA and image coordinates
        x = -ra_rotated
        y = +dec_rotated

        # Get the photometry in this image at that location
        # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
        platescale = imgtable.platescale[i_epoch] * platescale_multiplier
        f̃ₓ = imgtable.imageinterp[i_epoch](x/platescale, y/platescale)

        # Find the uncertainty in that photometry value (i.e. the contrast)
        if hasproperty(imgtable, :contrastmap)
            # If we have a 2D map
            σₓ = imgtable.contrastmapinterp[i_epoch](x/platescale, y/platescale)
        else
            # We have a 1D contrast curve
            r = √(x^2 + y^2)
            σₓ = imgtable.contrast[i_epoch](r / platescale)
        end

        # When we get a position that falls outside of our available
        # data (e.g. under the coronagraph) we cannot say anything
        # much about the planet.
        # We assume that we plateaued at the maximum flux
        if !isfinite(f̃ₓ)
            f̃ₓ = zero(typeof(f̃ₓ))
        end
        if !isfinite(σₓ) || iszero(σₓ)
            return convert(T, -Inf)
        end

        if !isfinite(flux)
            @warn "Flux variable is not finite" flux 
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
        ll_i = -1 / (2*σₓ²) * (flux^2 - 2*flux * f̃ₓ)
        ll += ll_i
    end

    return ll
end



# Generate new images  
function Octofitter.generate_from_params(like::ImageLikelihood, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

    # For image likelihood, we don't actually simulate new images in the standard sense
    # This function would need specific planet information to work properly
    # For now, return the original likelihood
    @warn "generate_from_params for ImageLikelihood not fully implemented with new API"
    return like
end
