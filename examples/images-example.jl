using Octofitter, Distributions, Plots


using DirectImages
cd(@__DIR__)
imagesc = map(centered, readfits("image-examples-1.fits",:))
pwd()
imshow2([
    imagesc[1]
    imagesc[2]
    imagesc[3]
    imagesc[4]
    imagesc[5]
], cmap=:magma, clims=(-1.0, 6.0))
##


@named b = Octofitter.Planet(
    Priors(
        a = Normal(16, 3),
        e = Beta(5, 25),
        τ = Normal(0.5, 1),
        ω = VonMises(0,1),
        i = VonMises(0,1),
        Ω = VonMises(0,1),
        H = Normal(5.5, 1)
    ),
)

system_images = Octofitter.Images(
    (band=:H, image=imagesc[1], platescale=10.0, epoch=times[1], contrast=contrasts[1]),
    (band=:H, image=imagesc[2], platescale=10.0, epoch=times[2], contrast=contrasts[2]),
    (band=:H, image=imagesc[3], platescale=10.0, epoch=times[3], contrast=contrasts[3]),
    (band=:H, image=imagesc[4], platescale=10.0, epoch=times[4], contrast=contrasts[4]),
    (band=:H, image=imagesc[5], platescale=10.0, epoch=times[5], contrast=contrasts[5]),
)

##

@named HD82134 = System(
    Priors(
        μ = TruncatedNormal(2.0, 0.1, 0, Inf),
        plx = TruncatedNormal(45., 0.02, 0, Inf),
    ),
    system_images,
    b,
)

chain = Octofitteradvancedhmc(
    HD82134, .60,
    MCMCThreads(),
    num_chains=4,
    adaptation =  1_000,
    iterations = 10_000,
    tree_depth =     12,
)

##
@show snr = mean(chain["X[H]"]) / std(chain["X[H]"])
##
Octofitter.MCMCChains.gelmandiag(chain)
##
plot(chain)
##
plotmodel(chain, alpha=0.05,  lims=1000)
##
savefig("../docs/src/assets/images-model-plot.svg")

##
using PairPlots
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
corner(table, labels, units)
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
    i = 0.1,
    μ = 2.0,
    plx = 45.
)
truth_elements = Visual{KepOrbit}(truth)
times = sort(rand(Uniform(0, period(truth_elements)/2), 5,))
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

# Create synthetic images at each time with those points
using DirectImageLikelihood, ImageFiltering
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
contrasts = [contrast for (img,contrast) in images_contrasts]
writefits("image-examples-1.fits", images...)

##

imshow(
    reduce(hcat, snrmap.(images, contrasts)),
    cmap=:magma,
    clims=(-1.0, 6.0)
)
