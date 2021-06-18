##
using Revise
import Random
using DirectDetections
using Distributions
using DirectImages
using ImageFiltering
using ComponentArrays
using DirectOrbits: KeplerianElements, period, raoff, decoff, kep2cart
using Plots
##

# Generate astrometry points using this template
# and we will try to fit the results
truths = [
    KeplerianElements(;a = 12,τ = 0.25,ω = 0,Ω = 0,e = 0.1,i = 0.5,μ = 1.,plx = 45.,),
    KeplerianElements(;a = 18, τ = 0.75,ω = 0,Ω = 0,e = 0.15, i = 0.5,μ = 1.,plx = 45.,),
]
intensities = [
    1800,
    3000,
]

# truth_elements = KeplerianElements(ComponentArray(truth, static))

# times = range(0, period(truth_elements)*4/5, length=9, )
times = range(0, period(first(truths))/3, length=4,)


# Create synthetic images at each time with those points
Random.seed!(1234)
images_contrasts = map(times) do t
    img = 20000fakeimage((201,201),1.2)
    img = centered(img)
    for (inten, truth) in zip(intensities, truths)
        ra, dec = kep2cart(truth, t)
        x = -ra
        y = dec
        img[round(Int,x/10), round(Int,y/10)] += inten #2500 #+ 300randn()
    end

    img = imfilter(img, Kernel.gaussian(5), "replicate")
    contrast = contrast_interp(img)

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
contrasts = [contrast for (img,contrast) in images_contrasts]
nothing

##
imshow2(reduce(hcat, images), cmap=:turbo, clims=(-35,35));

##
input = (;
    phot = (;
        Keck_L′ = [
            (image=images[1], platescale=10.0, epoch=times[1], contrast=contrasts[1]),
            (image=images[2], platescale=10.0, epoch=times[2], contrast=contrasts[2]),
            (image=images[3], platescale=10.0, epoch=times[3], contrast=contrasts[3]),
            (image=images[4], platescale=10.0, epoch=times[4], contrast=contrasts[4]),
        ],
        # Keck_Ks = [ 
        #     (image=images[7], platescale=10.0, epoch=times[7], contrast=contrasts[7]),
        #     (image=images[8], platescale=10.0, epoch=times[8], contrast=contrasts[8]),
        #     (image=images[9], platescale=10.0, epoch=times[9], contrast=contrasts[9]),
        # ]
    ),
    # astrom = [
    #     # (; 
    #     #     epoch=mean(times),
    #     #     ra=raoff(truth_elements, mean(times)),
    #     #     dec=decoff(truth_elements, mean(times)),
    #     #     σ_ra=5.,
    #     #     σ_dec=5.,
    #     # )
    # ]
)
nothing

##
priors = ComponentVector(
    # i = Normal(0.6, 0.3),
    # Ω = Normal(0.0, 0.3),
    # μ = Normal(1.0, 0.01),
    # plx = Normal(45., 0.0001),
    planets = [
        # (
        #     a = Uniform(3, 25),
        #     e = TruncatedNormal(0.1, 0.1, 0.0, 0.4),
        #     τ = Uniform(0,1),
        #     ω = Normal(0.0, 0.1),
        #     i = Normal(0.5, 0.1),
        #     Ω = Normal(0.0, 0.1),

        #     μ = Normal(1.0, 0.01),
        #     plx = Normal(45., 0.0001),

        #     # σ_f_model² = Truncated(InverseGamma(4,0.01), 0, 1),
        #     phot = (;
        #         Keck_L′ = Uniform(0., 100.),
        #         # Keck_Ks = Uniform(0., 100.),
        #     ),
        # ),
        (
            a = Uniform(3, 25),
            e = TruncatedNormal(0.1, 0.1, 0.0, 0.4),
            τ = Uniform(0,1),
            ω = Uniform(-π,π),
            # i = Uniform(0,2π),
            i = TruncatedNormal(0, 0.2, -π,π),
            Ω = Uniform(-π,π),

            μ = Normal(1.0, 0.1),
            plx = Normal(45., 0.0001),

            # σ_f_model² = Truncated(InverseGamma(4,0.01), 0, 1),
            phot = (;
                Keck_L′ = Uniform(0., 100.),
                # Keck_Ks = Uniform(0., 100.),
            ),
        ),
    ]
)

# test:
# θ = rand.(priors)
# logpdf.(priors, θ) |> sum


##
@time chains = DirectDetections.mcmc(
    priors, input;
    numwalkers=1600,
    burnin=62_000,
    numsamples_perwalker=100_000,
    thinning=250,
    squash = false
);
nothing

##
@time chainsh, stats = DirectDetections.hmc(
    priors, input;
    burnin=2_000,
    numwalkers=4,
    numsamples_perwalker=4_000,
);
nothing

##
[count(getproperty.(s, :numerical_error)) for s in stats]
##
function orbsonimg(planet, images, N=100)
    sampled = sampleorbits(planet, N);
    i = DirectImage(sum(images))
    i.PLATESCALE = 10.
    imshow(i,skyconvention=true)
    plot!(sampled,color=:white,alpha=5/N,label="")
    xlims!(i.PLATESCALE.*extrema(axes(i,1)))
    ylims!(i.PLATESCALE.*extrema(axes(i,2)))
end

##
orbsonimg(chains.planets[1],images)
# orbsonimg(chainsh[2].planets[1],images)

##
plot()
for c in chainsh
    histogram!(c.planets[1].phot.Keck_L′, alpha=0.5)
end
current()
##

function snr(phot)
    med,up = quantile(reshape(phot,:),(0.5, 0.67))
    return med/(up-med)
end

##
for c in chainsh
    snr(c.planets[1].phot.Keck_L′) |> println
end
##
snr(chains.planets[1].phot.Keck_L′[:])
##
histogram(chains.planets[1].phot.Keck_L′[:])
##
DirectDetections.plotposterior(chains.planets[1], :a, 5000)
plot!(truths, color=:black, lw=2, ls=:dash, label="", alpha=0.5)
savefig("images/readme-orbits.png")

##
DirectDetections.plotposterior(chains.planets[1], :i, colorbartitle="inclination (rad)")

##
DirectDetections.plotposterior(chains.planets[1], (:phot, :Keck_L′), 500, colorbartitle="flux", cmap=:plasma, rev=false, clims=(0,30))
plot!(truths, color=:black, lw=2, ls=:dash, label="")

##
DirectDetections.plotposterior(chains.planets[1], :e, 500, colorbartitle="eccentricity")
plot!(truths, color=:black, lw=2, ls=:dash)

##

ra, dec = projectpositions(chains.planets[1], mean(times))
histogram2d(ra,dec,aspectratio=1, color=:plasma, background_inside=:black, framestyle=:box,xflip=true,dpi=200, colorbartitle=".\n\nposterior density")
xlims!(1.25.*(-1,+1).*max(abs.(extrema(ra))...))
ylims!(1.25.*(-1,+1).*max(abs.(extrema(dec))...))
xlabel!("Δ right ascension (as)")
ylabel!("Δ declination (as)")

##
mask = (-100 .< ra .< 150) .& (200 .< dec .< 700)
snr(chains.planets[1].phot.Keck_L′)
snr(chains.planets[1].phot.Keck_L′[mask])
histogram2d(ra[mask],dec[mask],aspectratio=1, color=:plasma, background_inside=:black, framestyle=:box, cscale=:log10kep2cart,dpi=200)
xlims!(1.25.*(-1,+1).*max(abs.(extrema(ra))...))
ylims!(1.25.*(-1,+1).*max(abs.(extrema(dec))...))
##
ras,decs = projectpositions(chains.planets[1],mean(times))
# histogram2d(ras[5end÷6:end], decs[5end÷6:end])
histogram2d(ras, decs)

##
plot(chains.planets[1].a,legend=nothing,lw=0.1)
##
mask = 7 .< chains.planets[1].a .< 15
snr(chains.planets[1].phot.Keck_L′[mask])
mask = chains.planets[1].a .< 7
snr(chains.planets[1].phot.Keck_L′[mask])



##
table = (;
    chains.planets[1].a,
    L=chains.planets[1].phot.Keck_L′,
    chains.planets[1].i,
    chains.planets[1].e,
    tau=chains.planets[1].τ,
)
corner(
    table,
    split(raw"a \mathrm{flux} i e \tau"),
    plotscatter=true,
    hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    scatter_kwargs=(;markersize=0.5),
    dpi=200,
    # lim_factor=2,
    # hist2d_kwargs=(;nbins=8),
)
nothing
##

