module OctofitterImages

using Octofitter
using PlanetOrbits
using Tables, TypedTables


using ImageTransformations
using CoordinateTransformations
using Interpolations
using AstroImages


const images_cols = (:band, :image, :epoch, :platescale,)

"""
    ImageLikelihood(...)

A block of images of a system. Pass a vector of named tuples with the following fields:
$images_cols

For example:
```julia
ImageLikelihood(
    (; epoch=1234.0, band=:J, image=readfits("abc.fits"), platescale=19.4)
)
```
Contrast can be a function that returns the 1 sigma contrast of the image from a separation in mas to the same units as the image file.
Or, simply leave it out and it will be calculated for you.
Epoch is in MJD.
Band is a symbol which matches the one used in the planet's `Priors()` block.
Platescale is in mas/px.
"""
struct ImageLikelihood{TTable<:Table} <: Octofitter.AbstractLikelihood
    table::TTable
    function ImageLikelihood(observations...)
        table = Table(observations...)
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
        return new{typeof(table)}(table)
    end
end
Images(observations::NamedTuple...) = ImageLikelihood(observations)
export Images


"""
    contrast_interp(image; step=2)

Returns a linear interpolation on top of the results from `contrast`.
Extrapolated results return Inf.
"""
function contrast_interp(image::AstroImage; step=2)
    cont = contrast(image; step)
    return LinearInterpolation(cont.separation, cont.contrast, extrapolation_bc=Inf)
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





"""
Likelihood of there being planets in a sequence of images.
"""
function Octofitter.ln_like(images::ImageLikelihood, θ_planet, orbit)
    
    # Resolve the combination of system and planet parameters
    # as a VisualOrbit object. This pre-computes
    # some factors used in various calculations.
    # elements = construct_elements(θ_system, θ_planet)


    imgtable = images.table
    T = eltype(first(θ_planet))
    ll = zero(T)
    for i in eachindex(imgtable.epoch)

        soln = orbitsolve(orbit, imgtable.epoch[i])
            
        
        band = imgtable.band[i]

        star_δra =  0.0
        star_δdec = 0.0

        # Note the x reversal between RA and image coordinates
        x = -(raoff(soln) + star_δra)
        y = +(decoff(soln) + star_δdec)

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

        # Verify the user has specified a prior or model for this band.
        if !hasproperty(θ_planet, band)
            error("No photometry variable for the band $band was specified.")
        end
        # TODO: verify this is type stable
        f_band = getproperty(θ_planet, band)


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

    return ll
end



# Generate new images
function Octofitter.generate_from_params(like::ImageLikelihood, θ_system,  elements::Vector{<:VisualOrbit})

    newrows = map(like.table) do row
        (;band, image, platescale, epoch, psf) = row

        injected = copy(image)
        
        for i in eachindex(elements)
            θ_planet = θ_system.planets[i]
            elem = elements[i]
            # Generate new astrometry point
            os = orbitsolve(elem, epoch)

            # TODO: this does not consider the shift to the images due to the motion of the star
            ra = raoff(os)
            dec = decoff(os)

            phot = θ_planet[band]

            # TODO: verify signs
            dx = ra/platescale
            dy = -dec/platescale
            translation_tform = Translation(
                mean(axes(psf,1))-mean(axes(image,1))+mean(dims(image,1))+dx,
                mean(axes(psf,2))-mean(axes(image,2))+mean(dims(image,2))+dy
            )
            # TBD if we want to support rotations for handling negative sidelobes.

            psf_positioned = warp(psf, translation_tform, axes(image), fillvalue=0)
            psf_positioned[.! isfinite.(psf_positioned)] .= 0
            psf_scaled = psf_positioned .* phot ./ maximum(filter(isfinite, psf_positioned))
            injected .+= psf_scaled
        end

        return merge(row, (;image=injected))
    end

    return ImageLikelihood(newrows)
end


function __init__()



    return
end

end
