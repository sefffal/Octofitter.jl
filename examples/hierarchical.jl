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

# Generate astrometry points using this template
# and we will try to fit the results
truth = (; # ComponentVector{SVector{7,Float64}}
    a = 12,
    τ = 0.25,
    ω = 0,
    Ω = 0,
    e = 0.2,
    i = 0.5,
    μ = 1.,
    plx = 45.,
)
truth_elements = KeplerianElements(truth)
# truth_elements = KeplerianElements(ComponentArray(truth, static))

times = range(0, period(truth_elements)*3/4, length=4, )
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
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
        px + 700randn()/0.5r
    end

    img = centered(img)

    # img = imfilter(img, Kernel.gaussian(5), NA())
    img_for_contrast = imfilter(img, Kernel.gaussian(5), "replicate")
    contrast = contrast_interp(img_for_contrast)

    img[round(Int,x/10), round(Int,y/10)] += 5000 + 800randn()
    # img[round(Int,x2/10), round(Int,y2/10)] += 5000

    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
contrasts = [contrast for (img,contrast) in images_contrasts]
display.(imshow2.(images, cmap=:turbo, clims=(-35,35)));

##
# photometry = (;
#     L′ = (;
#         images=images,
#         times=times,
#         contrasts=contrast
#     ),
#     Ks = (;
#         images=images,
#         times=times,
#         contrasts=contrast
#     )
# )

photometry = (;
    L′ = [
        (phot=images[1], epoch=times[1], contrast=contrasts[1]),
        (phot=images[2], epoch=times[2], contrast=contrasts[2]),
        (phot=images[3], epoch=times[3], contrast=contrasts[3]),
        (phot=images[4], epoch=times[4], contrast=contrasts[4]),
    ],
    Ks = [
        (phot=images[1], epoch=times[1], contrast=contrasts[1]),
        (phot=images[2], epoch=times[2], contrast=contrasts[2]),
        (phot=images[3], epoch=times[3], contrast=contrasts[3]),
        (phot=images[4], epoch=times[4], contrast=contrasts[4]),
    ]
)

##
using ComponentArrays
priors = ComponentVector(
    # i = Normal(0.6, 0.3),
    # Ω = Normal(0.0, 0.3),
    μ = Normal(1.0, 0.01),
    plx = Normal(45., 0.0001),
    planets = [
        (
            a = Uniform(8, 25),
            e = TruncatedNormal(0.0, 0.4, 0.0, 0.9999),
            τ = Uniform(0,1),
            ω = Normal(0.0, 0.3),
            i = Normal(0.5, 0.3),
            Ω = Normal(0.0, 0.3),
            Teff = TruncatedNormal(1500, 800, 0, Inf),
            logg = TruncatedNormal(4, 1, 0, Inf),
            
            f = Uniform(0., 100.),
            # f_spread = TruncatedNormal(0, 1, 0., Inf),
            # TODO: rename f_i? Or do I want to allow astrometric offsets to be fit? North angle?

            # The flux in each epoch is separate, but these flux parameters are drawn
            # from a common distribution for this planet.
            epochs = [
                (;f = Uniform(0., 100.)),
                (;f = Uniform(0., 100.)),
                (;f = Uniform(0., 100.)),
                (;f = Uniform(0., 100.)),
            ]

            # photometry = (;
            #     L′ = (
            #         model = Uniform(0., 100.),
            #         model_spread = Normal(0.1, 0.1),
            #         epochs = [
            #             Uniform(0., 100.),
            #             Uniform(0., 100.),
            #             Uniform(0., 100.),
            #             Uniform(0., 100.),
            #     ],
            #     Ks = (
            #         model = Uniform(0., 100.),
            #         model_spread = Normal(0.1, 0.1),
            #         epochs = [
            #             Uniform(0., 100.),
            #             Uniform(0., 100.),
            #         ]
            #     )
            # ),
            # rv = Distribution{Univariate, Continuous}[],
            # astrometry = Distribution{Univariate, Continuous}[],
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
    burnin=4_000,
    numwalkers=5000,
    numsamples_perwalker=5_000,
    squash = true
);
nothing

##
corner(
    (;
        f=chains["planets[1].f"][:],
        a=chains["planets[1].a"][:],
        e=chains["planets[1].e"][:],
        tau=chains["planets[1].τ"][:]
    ),
    ["f_0", "a", "e", "\\tau",],
    hist_kwargs=(;nbins=15),
    hist2d_kwargs=(;nbins=15),
    plotscatter=false
)

##
corner(
    (;
        f=chains["planets[1].f"][:],
        f₁=chains["planets[1].epochs[1].f"][:],
        f₂=chains["planets[1].epochs[2].f"][:],
        f₃=chains["planets[1].epochs[3].f"][:],
        f₄=chains["planets[1].epochs[4].f"][:],
    ),
    ["f_0", "f_1", "f_2", "f_3", "f_4"],
    hist_kwargs=(;nbins=15),
    hist2d_kwargs=(;nbins=15),
    plotscatter=false
)

##
using StatsBase
bins = 23:1:40
f=fit(Histogram, chains["planets[1].f"][:], bins)
f₁=fit(Histogram, chains["planets[1].epochs[1].f"][:], bins)
f₂=fit(Histogram, chains["planets[1].epochs[2].f"][:], bins)
f₃=fit(Histogram, chains["planets[1].epochs[3].f"][:], bins)
f₄=fit(Histogram, chains["planets[1].epochs[4].f"][:], bins)
plot(fontfamily="", background=:black, grid=:none, minorgrid=:none)
for (f,l) in zip((f₁, f₂, f₃, f₄), ("f₁", "f₂", "f₃", "f₄"))
    plot!(f.edges[1][1:end-1].+step(f.edges[1])/2, f.weights, label=l)
end
plot!(f.edges[1][1:end-1].+step(f.edges[1])/2, f.weights, color=:white, lw=3, label="f")
ylabel!("Posterior density")
xlabel!("Flux - arb.")
vline!([mean(maximum.(filter.(isfinite, images)))], label="truth", color=:white)

##
sampled = map(rand(1:size(chains,1),500)) do i
    el = KeplerianElements(
        # i=chains[:i][i],
        # Ω=chains[:Ω][i],
        μ=chains[i,:μ,1],
        plx=chains[i,:plx,1],
        a=chains[i,"planets[1].a",1],
        e=chains[i,"planets[1].e",1],
        τ=chains[i,"planets[1].τ",1],
        ω=chains[i,"planets[1].ω",1],
        i=chains[i,"planets[1].i",1],
        Ω=chains[i,"planets[1].Ω",1],
    )
end
i = DirectImage(
    sum(images),
)
i.PLATESCALE = 10.
imshow(i, skyconvention=true, clims=(-10,40))
plot!(sampled, color=:white, alpha=0.005, label="")
xlims!(extrema(axes(i,1)).*i.PLATESCALE)
ylims!(extrema(axes(i,2)).*i.PLATESCALE)


##

# flux parameter of each epoch vs measured photometry, over contrast
-0.5*sum((fi - di)^2)/si^2
# Flux parameter of each epoch vs overall parameter for the planet flux
-0.5*sum((fi-f)^2)/(0.1*f)^2
# Overall planet flux parameter vs flux from a model grid, unceratinty in that model
-0.5*(f-fmodel)^2/smodel^2


##

# flux parameter of each epoch vs measured photometry, over contrast
-0.5*sum((fi - di)^2)/si^2
# Flux parameter of each epoch vs overall parameter for the planet flux
-0.5*sum((fi-f)^2)/(0.1*f)^2
# Overall planet flux parameter vs flux from a model grid, unceratinty in that model
-0.5*(f-fmodel)^2/smodel^2

fᵢ -> f -> (logg, Teff)