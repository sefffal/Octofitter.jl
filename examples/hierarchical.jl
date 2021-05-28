# using Plots, PairPlots
# theme(:dao)

##
import Random
using DirectDetections
using Distributions
using DirectImages


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

times = range(0, period(truth_elements)*4/5, length=9, )
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

# Create synthetic images at each time with those points
Random.seed!(1234)
images_contrasts = map(eachrow(points)) do (ra,dec)
    x = -ra
    y = dec

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

    img[round(Int,x/10), round(Int,y/10)] += 5000 + 300randn()
    # img[round(Int,x2/10), round(Int,y2/10)] += 5000

    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
contrasts = [contrast for (img,contrast) in images_contrasts]
nothing

##
# display.(imshow2.(images, cmap=:turbo, clims=(-35,35)));

##
input = (;
    phot = (;
        Keck_L′ = [
            (image=images[1], platescale=10.0, epoch=times[1], contrast=contrasts[1]),
            (image=images[2], platescale=10.0, epoch=times[2], contrast=contrasts[2]),
            (image=images[3], platescale=10.0, epoch=times[3], contrast=contrasts[3]),
            (image=images[5], platescale=10.0, epoch=times[5], contrast=contrasts[5]),
            (image=images[6], platescale=10.0, epoch=times[6], contrast=contrasts[6]),
            # (image=images[7], platescale=10.0, epoch=times[7], contrast=contrasts[7]),
            # (image=images[8], platescale=10.0, epoch=times[8], contrast=contrasts[8]),
            # (image=images[9], platescale=10.0, epoch=times[9], contrast=contrasts[9]),
        ],
        Keck_Ks = [ 
            (image=images[7], platescale=10.0, epoch=times[7], contrast=contrasts[7]),
            (image=images[8], platescale=10.0, epoch=times[8], contrast=contrasts[8]),
            (image=images[9], platescale=10.0, epoch=times[9], contrast=contrasts[9]),
        ]
    ),
    astrom = [
        # (; 
        #     epoch=mean(times),
        #     ra=raoff(truth_elements, mean(times)),
        #     dec=decoff(truth_elements, mean(times)),
        #     σ_ra=5.,
        #     σ_dec=5.,
        # )
    ]
)
nothing

##
priors = ComponentVector(
    # i = Normal(0.6, 0.3),
    # Ω = Normal(0.0, 0.3),
    # μ = Normal(1.0, 0.01),
    # plx = Normal(45., 0.0001),
    planets = [
        (
            a = Uniform(8, 25),
            e = TruncatedNormal(0.0, 0.4, 0.0, 0.9999),
            τ = Uniform(0,1),
            ω = Normal(0.0, 0.3),
            i = Normal(0.5, 0.3),
            Ω = Normal(0.0, 0.3),

            μ = Normal(1.0, 0.01),
            plx = Normal(45., 0.0001),

            Teff = TruncatedNormal(1200, 800, 200, 2400),
            mass = Uniform(0.55, 25),
            
            # f = Uniform(0., 100.),
            # f_spread = TruncatedNormal(0, 1, 0., Inf),
            # TODO: rename f_i? Or do I want to allow astrometric offsets to be fit? North angle?

            # The flux in each epoch is separate, but these flux parameters are drawn
            # from a common distribution for this planet.
            # epochs = [
            #     (;f = Uniform(0., 100.)),
            #     (;f = Uniform(0., 100.)),
            #     (;f = Uniform(0., 100.)),
            #     (;f = Uniform(0., 100.)),
            # ]

            phot = (;
                Keck_L′ = (
                    f = Uniform(0., 100.),
                    σ_f² = Truncated(InverseGamma(4,0.01), 0, 1),
                    σ_f_model² = Truncated(InverseGamma(4,0.01), 0, 1),
                    epochs = [
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                    ],
                ),
                Keck_Ks = (
                    f = Uniform(0., 100.),
                    σ_f² = Truncated(InverseGamma(4,0.01), 0, 1),
                    σ_f_model² = Truncated(InverseGamma(4,0.01), 0, 1),
                    epochs = [
                        Uniform(0., 100.),
                        Uniform(0., 100.),
                    ]
                )
            ),
        ),
    ]
)

# test:
θ = rand.(priors)
logpdf.(priors, θ) |> sum


##
@time chains = DirectDetections.mcmc(
    priors, input;
    # numwalkers=3000,
    # burnin=20_000,
    # numsamples_perwalker=25_000,
    
    numwalkers=300,
    burnin=20_00,
    numsamples_perwalker=25_00,
    squash = true
);
nothing

##
thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
nothing

##
reinterptted = reinterpret(reshape, eltype(first(thetase′)), thetase′);
chains = ComponentArray(collect(eachrow(reinterptted)), getaxes(thetase′[1]));
nothing


## What do we want?
chains.planets[1].a # -> matrix of chains concatenated together





##
table = (;
    chains.planets[1].a,
    chains.planets[1].mass,
    chains.planets[1].Teff,
)
corner(table,plotscatter=false)




##
corner(
    (;
        L=chains["planets[1].phot.Keck_L′.f"][:],
        K=chains["planets[1].phot.Keck_Ks.f"][:],
        # σ=chains["planets[1].phot.Keck_L′.σ_f²"][:],
        T=chains["planets[1].Teff"][:],
        m=chains["planets[1].mass"][:],
        a=chains["planets[1].a"][:],
        e=chains["planets[1].e"][:],
        tau=chains["planets[1].τ"][:]
    ),
    ["L\\prime", "Ks", "\\mathrm{T_{eff}}", "\\mathrm{mass}", "a", "i", "e", "\\tau",],
    # hist_kwargs=(;nbins=15),
    # hist2d_kwargs=(;nbins=15),
    plotscatter=false
)

##
corner(
    (;
        f =chains["planets[1].phot.Keck_L′.f"][:],
        f₁=chains["planets[1].phot.Keck_L′.epochs[1]"][:],
        f₂=chains["planets[1].phot.Keck_L′.epochs[2]"][:],
        f₃=chains["planets[1].phot.Keck_L′.epochs[3]"][:],
        f₄=chains["planets[1].phot.Keck_L′.epochs[4]"][:],
        σ_f =chains["planets[1].phot.Keck_L′.σ_f²"][:],
    ),
    ["f_0", "f_1", "f_2", "f_3", "f_4", "\\sigma_f"],
    hist_kwargs=(;nbins=15),
    hist2d_kwargs=(;nbins=15),
    plotscatter=false
)

##
histogram(chains["planets[1].phot.Keck_L′.σ_f²"][:])
##
using StatsBase
bins = 0:1:100
f₀=fit(Histogram, chains["planets[1].phot.Keck_L′.f"][:], bins)
f₁=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[1]"][:], bins)
f₂=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[2]"][:], bins)
f₃=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[3]"][:], bins)
f₄=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[4]"][:], bins)
f5=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[5]"][:], bins)
f6=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[6]"][:], bins)
f7=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[7]"][:], bins)
f8=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[8]"][:], bins)
f9=fit(Histogram, chains["planets[1].phot.Keck_L′.epochs[9]"][:], bins)
# plot(fontfamily="", background=:black, grid=:none, minorgrid=:none)
plot(fontfamily="")
for (f,l) in zip((f₁, f₂, f₃, f₄, f5, f6, f7, f8, f9), ("f₁", "f₂", "f₃", "f₄", "f₅", "f₆", "f₇", "f₈", "f₉"))
    plot!(f.edges[1][1:end-1].+step(f.edges[1])/2, f.weights, label=l)
end
plot!(f₀.edges[1][1:end-1].+step(f₀.edges[1])/2, f₀.weights, color=:black, lw=3, label="f", ls=:dash)
ylabel!("Posterior density")
xlabel!("Flux - arb.")
# vline!([mean(maximum.(filter.(isfinite, images)))], label="truth", color=:black, ls=:solid)

##
using ColorSchemes
# sampled = map(rand(1:size(chains,1),500)) do i
sampled = map(round.(Int, range(1, size(chains,1), length=500))) do i
    el = KeplerianElements(
        # μ=chains[i,:μ,1],
        # plx=chains[i,:plx,1],
        μ=chains[i,"planets[1].μ",1],
        plx=chains[i,"planets[1].plx",1],
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
    # images[1]
)
i.PLATESCALE = 10.
imshow(i, skyconvention=true, clims=(-10,40))
# plot(xflip=true, background=:black, grid=:none)
# plot!(sampled, color=:white, alpha=0.05, label="")
for (j,s) in enumerate(sampled)
    plot!(s, color=:white, alpha=0.01, label="")
end


ras = map(obs->obs.ra, input.astrom)
decs = map(obs->obs.dec, input.astrom)
σ_ras = map(obs->obs.σ_ra, input.astrom)
σ_decs = map(obs->obs.σ_dec, input.astrom)
scatter!(ras, decs, xerr=σ_ras, yerr=σ_decs, marker=(0, :square,:white), markerstrokecolor=:white, markerstrokewidth=1, label="")

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