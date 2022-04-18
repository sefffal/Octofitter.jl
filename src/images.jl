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
    dx = dims(image,1)
    dy = dims(image,2)
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
