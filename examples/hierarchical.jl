import Random
using DirectDetections
using Distributions
using DirectImages
using Plots, PairPlots
theme(:dao)

using MCMCChains: Chains
using ImageFiltering

using ComponentArrays

using DirectOrbits: KeplerianElements, period, raoff, decoff

# We will keep these elements constant
static = (; # ComponentVector{SVector{2,Float64}}
    μ = 1.,
    plx = 45.,
    # τ = 100,
)

# Generate astrometry points using this template
# and we will try to fit the results
truth = (; # ComponentVector{SVector{7,Float64}}
    f = 20.,
    a = 12,
    τ = 0.25,
    ω = 0,
    Ω = 0,
    e = 0.2,
    i = 0.5,
)
truth_elements = KeplerianElements(merge(truth, static))
# truth_elements = KeplerianElements(ComponentArray(truth, static))


truth2 = (;
    f = 15.,
    a = 18,
    τ = 0.75,
    ω = 0.1,
    Ω = 0.0,
    e = 0.25,
    i = 0.55,
)
truth_elements2 = KeplerianElements(merge(truth2, static))

times = range(0, period(truth_elements)/4, length=4, )
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

times2 = range(0, period(truth_elements2)/4, length=4, )
points2 = hcat(
    raoff.(truth_elements2, times2),
    decoff.(truth_elements2, times2)
)

# Create synthetic images at each time with those points
Random.seed!(1234)
images_contrasts = map(zip(eachrow(points),eachrow(points2))) do ((ra,dec),(ra2,dec2))
    x = -ra
    y = dec
    x2 = -ra2
    y2 = dec2

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

    img[round(Int,x/10), round(Int,y/10)] += 5000
    # img[round(Int,x2/10), round(Int,y2/10)] += 5000

    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
contrasts = [contrast for (img,contrast) in images_contrasts]
display.(imshow2.(images, cmap=:turbo));

##
using ComponentArrays
const CV = ComponentVector
priors = CV(
    i = Normal(0.6, 0.3),
    Ω = Normal(0.0, 0.3),
    μ = Normal(1.0, 0.01),
    plx = Normal(45., 0.0001),
    planets = [
        CV(
            a = TruncatedNormal(16, 8, 4., Inf),
            e = TruncatedNormal(0.21, 0.2, 0.0, 0.9999),
            τ = Uniform(0,1),
            ω = Normal(0.0, 0.3),
            i = Normal(0.6, 0.3),
            Ω = Normal(0.0, 0.3),
            Teff = TruncatedNormal(1500, 800, 0, Inf),
            logg = TruncatedNormal(4, 1, 0, Inf),
            
            f = Uniform(0., 100.),
            # f_spread = TruncatedNormal(0, 1, 0., Inf),
            # TODO: rename f_i? Or do I want to allow astrometric offsets to be fit? North angle?

            # The flux in each epoch is separate, but these flux parameters are drawn
            # from a common distribution for this planet.
            epochs = [
                CV(f = Uniform(0., 100.)),
                CV(f = Uniform(0., 100.)),
                CV(f = Uniform(0., 100.)),
                CV(f = Uniform(0., 100.)),
            ]
        ),
        # CV(
        #     a = TruncatedNormal(16, 8, 4., Inf),
        #     e = TruncatedNormal(0.21, 0.2, 0.0, 0.9999),
        #     τ = Uniform(0,1),
        #     ω = Normal(0.0, 0.3),
        #     i = Normal(0.6, 0.3),
        #     Ω = Normal(0.0, 0.3),
        #     f = TruncatedNormal(28, 5, 0., Inf),
        #     f_spread = TruncatedNormal(0, 1, 0., Inf),
        #     epochs = [
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #         CV(f = TruncatedNormal(28, 5, 0., Inf)),
        #     ]
        # )
    ]

)


sum(logpdf.(priors, rand.(priors)))


##
@time chains = DirectDetections.mcmc(
    priors, images, contrasts, times;
    platescale=10.,
    burnin=10_000,
    numwalkers=5000,
    numsamples_perwalker=12_000,
    squash = true
);
nothing

##
corner(
    (;
        f=chains["planets[1].f"][:],
        a=chains["planets[1].a"][:],
        tau=chains["planets[1].τ"][:]
    ),
    hist_kwargs=(;nbins=10),
    hist2d_kwargs=(;nbins=10),
    plotscatter=false
)

##
corner(
    (;
        f=chains["planets[1].f"][:],
        f₁=chains["planets[1].epochs[1].f"][:],
        f₂=chains["planets[1].epochs[2].f"][:],
        f₃=chains["planets[1].epochs[3].f"][:],
    ),
    ["f_0", "f_1", "f_2", "f_3"],
    hist_kwargs=(;nbins=15),
    hist2d_kwargs=(;nbins=15),
    plotscatter=false
)


##