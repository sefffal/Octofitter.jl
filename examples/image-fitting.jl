##
using Revise
import Random
using DirectDetections
using Distributions
using DirectImages
using ImageFiltering
using DirectOrbits
using Plots
##

# Generate astrometry points using this template
# and we will try to fit the results
truths = [
    # KeplerianElements(;a = 12,τ = 0.25,ω = 0,Ω = 0,e = 0.1,i = 0.5,μ = 1.,plx = 45.,),
    KeplerianElements(;a = 12,τ = 0.65,ω = 0,Ω = 0,e = 0.1,i = 0.5,μ = 1.,plx = 45.,),
    # KeplerianElements(;a = 18, τ = 0.75,ω = 0,Ω = 0,e = 0.15, i = 0.5,μ = 1.,plx = 45.,),
]
intensities = [
    1800,
    # 3000,
    # 650,
    # 1000,
]

# truth_elements = KeplerianElements(ComponentArray(truth, static))

# times = range(0, period(truth_elements)*4/5, length=9, )
times = range(0, period(first(truths))/2, length=10,)


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
# imshow2(reduce(hcat, images), clims=(-35,35))

##

planet_priors = Priors(
    a = Uniform(7, 16),
    e = TruncatedNormal(0.1, 0.1, 0.0, 0.4),
    τ = Uniform(0,1),
    # τ = Uniform(0.4,0.8),
    ω = Uniform(-π,π),
    # ω = Uniform(-π/4,π/4),
    i = Uniform(-π,π),
    # i = Uniform(-π/4,π/4),
    # i = TruncatedNormal(0, 0.2, -π,π),
    Ω = Uniform(-π,π),
    # Ω = Uniform(-π/4,π/4),

    # mass = Uniform(0mjup2msol,70mjup2msol),

    J = Uniform(0,20)
)
planet = reparameterize(DirectDetections.Planet(planet_priors))
nothing
##


system_images = DirectDetections.Images(
    (band=:J, image=images[1], platescale=10.0, epoch=times[1]),#, contrast=contrasts[1]),
    (band=:J, image=images[2], platescale=10.0, epoch=times[2]),#, contrast=contrasts[2]),
    (band=:J, image=images[3], platescale=10.0, epoch=times[3]),#, contrast=contrasts[3]),
    (band=:J, image=images[4], platescale=10.0, epoch=times[4]),#, contrast=contrasts[4]),
)


system_priors = Priors(
    μ = TruncatedNormal(1.0, 0.01, 0, 10),
    plx = TruncatedNormal(45., 0.0001, 0, 100),
)

# system = System(system_priors, system_pma, system_images, planet_b)
system = System(system_priors, system_images, planet,)
nothing

##
@time chains_img = DirectDetections.mcmc(
    system;
    # numwalkers=800,
    # burnin=62_000,
    # numsamples_perwalker=100_000,
    # thinning=250,

    numwalkers=50,
    burnin=1000,
    numsamples_perwalker=5_000,
    thinning=1,
    
    
    squash = false
);
nothing

##
@time chainsh, stats = DirectDetections.hmc(
    system;
    burnin=10_000,
    numwalkers=10,
    numsamples_perwalker=15_000,
);
nothing

##
[count(getproperty.(s, :numerical_error)) for s in stats]
##
function orbsonimg(args...; kwargs...)
    plot()
    orbsonimg!(args...;kwargs...)
end
function orbsonimg!(planet, image, N=100; kwargs...)
    sampled = sampleorbits(planet, N);
    i = DirectImage(image)
    i.PLATESCALE = 10.
    imshow!(i;skyconvention=true,kwargs...)
    # plot!(sampled,color=:white,alpha=5/N,label="")
    plot!(sampled,color=:white,label="",)
    xlims!(i.PLATESCALE.*extrema(axes(i,1)))
    ylims!(i.PLATESCALE.*extrema(axes(i,2)))

end


##
imshow(images[3],clims=(-2,10))
# orbsonimg(chains.planets[1],images[2],clims=(-2,10))
orbsonimg(chains.planets[1],mean(images), 500,legend=nothing, clims=(-0.2,1),color=:magma,dpi=200)
plot!(truths, color=:black, lw=2, ls=:dash, label="")
# orbsonimg(chainsh[2].planets[1],images)
# writefits("tmp.fits", images[2])
##
plot()
orbsonimg!(chainsh[1].planets[1],mean(images),1,legend=nothing, clims=(-0.2,1),color=:magma,dpi=200,)
for c in chainsh
    plot!(sampleorbits(c.planets[1],125),color=:white,label="",alpha=0.1)
end
    plot!(truths, color=:black, lw=2, ls=:dash, label="",colorbar=nothing)

##


function snr(phot)
    low,med,up = quantile(reshape(phot,:),(0.33, 0.5, 0.67))
    st = std(reshape(phot,:))
    return med/(up-low)
end

##
for c in chainsh
    snr(c.planets[1].phot.Keck_L′) |> println
end
##
snr(chains.planets[1].phot.Keck_L′[:])
##
histogram(chains.planets[1].phot.Keck_L′[:])
histogram(chainsh[1].planets[1].phot.Keck_L′[:])
##
DirectDetections.plotposterior(chains_img, 1, :a, 500)
plot!(truths, color=:black, lw=2, ls=:dash, label="", alpha=0.5)
# savefig("images/readme-orbits.png")

##
DirectDetections.plotposterior(chains_img, chains_img.planets[1], :i, 500, colorbartitle="inclination (rad)",dpi=200)
# savefig("images/readme-post-i.png")

##
DirectDetections.plotposterior(chains_img, chains_img.planets[1], (:phot, :Keck_L′), 500, colorbartitle="flux", cmap=:plasma, rev=false, clims=(0,12),dpi=200)
plot!(truths, color=:black, lw=2, ls=:dash, label="")
# savefig("images/readme-post-f.png")

##
DirectDetections.plotposterior(chains.planets[1], :e, 500, colorbartitle="eccentricity",dpi=200)
plot!(truths, color=:black, lw=2, ls=:dash)
# savefig("images/readme-post-e.png")

##
DirectDetections.plotposterior(chains.planets[1],  lw=3, :e, 500, colorbartitle="eccentricity",dpi=200, cmap=:plasma, rev=false,)
plot!(truths, color=:black, lw=2, ls=:dash)
# savefig("images/readme-post-e2.png")

##

ra, dec = projectpositions(chains.planets[1], mean(times))
histogram2d(ra,dec,aspectratio=1, color=:plasma, background_inside=:black, framestyle=:box,xflip=true,dpi=200, colorbartitle=".\n\nposterior density")
xlims!(1.25.*(-1,+1).*max(abs.(extrema(ra))...))
ylims!(1.25.*(-1,+1).*max(abs.(extrema(dec))...))
xlabel!("Δ right ascension (as)")
ylabel!("Δ declination (as)")

##
mask = (-200 .< ra .< 250) .& (300 .< dec .< 800)
snr(chains.planets[1].phot.Keck_L′)
snr(chains.planets[1].phot.Keck_L′[mask])
histogram2d(ra[mask],dec[mask],aspectratio=1, color=:plasma, background_inside=:black, framestyle=:box, cscale=:log10kep2cart,dpi=200)
histogram2d(ra[mask],dec[mask],aspectratio=1, color=:plasma, background_inside=:black, framestyle=:box, cscale=:log10kep2cart,dpi=200)
# xlims!(1.25.*(-1,+1).*max(abs.(extrema(ra))...))
# ylims!(1.25.*(-1,+1).*max(abs.(extrema(dec))...))
##
ras,decs = projectpositions(chains.planets[1],mean(times))
# histogram2d(ras[5end÷6:end], decs[5end÷6:end])
histogram2d(ras, decs,aspectratio=1,xflip=true,color=:plasma)

##
mask = 7 .< chains.planets[1].a .< 15
snr(chains.planets[1].phot.Keck_L′[mask])
mask = chains.planets[1].a .< 7
snr(chains.planets[1].phot.Keck_L′[mask])



##
table = (;
    L=chains.planets[1].phot.Keck_L′,
    chains.planets[1].a,
    chains.planets[1].i,
    chains.planets[1].e,
    tau=chains.planets[1].τ,
)
corner(
    table,
    split(raw"\mathrm{flux} a i e \tau"),
    plotscatter=false,
    hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    scatter_kwargs=(;markersize=0.5),
    dpi=200,
    # lim_factor=2,
    # hist2d_kwargs=(;nbins=8),
)

##
ra, dec = projectpositions(chains.planets[1], mean(times))
pa = atan.(dec,ra)
# mask = (pa .< 0.5) .& 
#     (600 .< ra .< 800) .&
#     (-200 .< dec .< -100)
mask = vec((chains.planets[1].phot.Keck_L′ .< 6))
table = (;
    ra=ra[mask],
    dec=dec[mask],
    # pa=rad2deg.(pa[mask]),
    a=(chains.planets[1].a)[mask],
    L=(chains.planets[1].phot.Keck_L′)[mask],
    i=rad2deg.((chains.planets[1].i))[mask],
    e=(chains.planets[1].e)[mask],
)
corner(
    table,
    # split(raw"a \mathrm{flux} i e \tau"),
    plotscatter=false,
    # hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    # scatter_kwargs=(;markersize=0.5),
    # dpi=200,
    # lim_factor=2,
    hist_kwargs=(;nbins=50),
    hist2d_kwargs=(;nbins=50),
)



## Chains
table = (;
    a=(chains.planets[1].a),
    L=(chains.planets[1].phot.Keck_L′),
    i=rad2deg.((chains.planets[1].i)),
    e=(chains.planets[1].e),
)
corner(
    table,
    # split(raw"a \mathrm{flux} i e \tau"),
    plotscatter=false,
    # hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    # scatter_kwargs=(;markersize=0.5),
    # dpi=200,
    # lim_factor=2,
    hist_kwargs=(;nbins=50),
    hist2d_kwargs=(;nbins=50),
)

## Chainsh
table = (;
    a=reduce(vcat, [c.planets[1].a for c in chainsh]),
    L=reduce(vcat, [c.planets[1].phot.Keck_L′ for c in chainsh]),
    i=rad2deg.(reduce(vcat, [c.planets[1].i for c in chainsh])),
    e=reduce(vcat, [c.planets[1].e for c in chainsh]),
)
corner(
    table,
    # split(raw"a \mathrm{flux} i e \tau"),
    plotscatter=false,
    # hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    # scatter_kwargs=(;markersize=0.5),
    # dpi=200,
    # lim_factor=2,
    hist_kwargs=(;nbins=25),
    hist2d_kwargs=(;nbins=25),
)