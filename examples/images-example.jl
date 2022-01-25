using DirectDetections, Distributions, Plots




##
using DirectImages
cd(@__DIR__)
images = map(centered, readfits("image-examples-1.fits",:))
pwd()
imshow2([
    images[1]
    images[2]
    images[3]
    images[4]
    images[5]
], cmap=:magma, clims=(-1.0, 6.0))
##


@named b = DirectDetections.Planet(
    Priors(
        a = Normal(16, 3),
        e = Beta(5, 25),
        τ = Normal(0.5, 1),
        ω = Normal(0.1, deg2rad(30.)),
        i = Normal(0.55, deg2rad(10.)),
        Ω = Normal(0.0, deg2rad(30.)),
        H = Normal(5.5, 1)
    ),
)

##

@named HD82134 = System(
    Priors(
        μ = Normal(2.0, 0.1),
        plx =Normal(45., 0.02),
    ),
    Images(
        (band=:H, image=images[1], platescale=10.0, epoch=1238.6),
        (band=:H, image=images[2], platescale=10.0, epoch=1584.7),
        (band=:H, image=images[3], platescale=10.0, epoch=3220.0),
        (band=:H, image=images[4], platescale=10.0, epoch=7495.9),
        (band=:H, image=images[5], platescale=10.0, epoch=7610.4),
    ),
    b,
)

chain, stats, extra = DirectDetections.hmc(
    HD82134,
    adaptation =  4_000,
    iterations =  6_000,
    tree_depth =     10,
);

##
plotmodel(chain, HD82134, lims=1000)
##
savefig("../docs/src/assets/images-model-plot.svg")

##using PairPlots
table = (;
a=         chain["b[a]"],
H=         chain["b[H]"],
e=         chain["b[e]"],
i=rad2deg.(chain["b[i]"]),
Ω=rad2deg.(chain["b[Ω]"]),
ω=rad2deg.(chain["b[ω]"]),
τ=         chain["b[τ]"],
)
labels=[
"a",
"H",
"e",
"i",
"\\Omega",
"\\omega",
"\\tau",
]
units = [
"(au)",
"(arb.)",
"",
"(\\degree)",
"(\\degree)",
"(\\degree)",
"",
]
corner(table, labels, units)#, hist2d_kwargs=(;nbins=15))
##
cd(@__DIR__)
savefig("../docs/src/assets/images-corner-plot.svg")


## This block generates the images we used for the example
using Random
using DirectOrbits
Random.seed!(1234)

truth = (;
    a = 16,
    τ = 0.75,
    ω = 0.1,
    Ω = 0.0,
    e = 0.25,
    i = 0.55,
    μ = 2.0,
    plx = 45.
)
truth_elements = KeplerianElements(truth)
times = sort(rand(Uniform(0, period(truth_elements)/2), 5,))
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

# Create synthetic images at each time with those points
using DirectImages, ImageFiltering
images_contrasts = map(eachrow(points)) do (ra,dec)
    x = -ra
    y = dec
    img = zeros(201,201)
    r = imgsep(img)
    img[r .< 2] .= NaN

    img = map(zip(img, r)) do (px,r)
        # px + 2000randn()/0.5r
        # px + 900randn()/0.5r
        px + 200randn()/0.5r
    end

    img = centered(img)

    # img = imfilter(img, Kernel.gaussian(5), NA())
    img_for_contrast = imfilter(img, Kernel.gaussian(5), "replicate")
    contrast = contrast_interp(img_for_contrast)

    img[round(Int,x/10), round(Int,y/10)] += 300#800

    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
# contrasts = [contrast for (img,contrast) in images_contrasts]
writefits("image-examples-1.fits", images...)
