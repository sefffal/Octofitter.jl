##
using Random: Random
using Octofitter
using Distributions
using DirectImages
using ImageFiltering
using DirectOrbits
using Plots
##

# Generate astrometry points using this template
# and we will try to fit the results
truths = [
    VisualOrbit(; a = 12, τ = 0.15, ω = 0.5, Ω = 0, e = 0.2, i = 0.15pi, M = 1.0, plx = 45.0),

]
intensities = [
    # 1800,
    # 3500*3,
    0
]

# truth_elements = VisualOrbit(ComponentArray(truth, static))

# times = range(0, period(truth_elements)*4/5, length=9, )
# times = range(0, 0.6period(first(truths)), length = 2,)
times = range(0, 0.3period(first(truths)), length = 4,)


# Create synthetic images at each time with those points
Random.seed!(1234)
images_contrasts = map(times) do t
    img = 20000fakeimage((201, 201), 1.2)
    img = centered(img)
    img_copy = copy(img)
    for (inten, truth) in zip(intensities, truths)
        o = orbitsolve(truth, t)
        x = -raoff(o)
        y = decoff(o)
        img[round(Int, x / 10), round(Int, y / 10)] += inten #2500 #+ 300randn()
    end

    img = imfilter(img, Kernel.gaussian(5), "replicate")
    img_copy = imfilter(img_copy, Kernel.gaussian(5), "replicate")

    sep = centered(DirectImages.imgsep(img))
    img[sep.<12] .= NaN
    img_copy[sep.<12] .= NaN

    contrast = contrast_interp(img_copy)
    img, contrast
end
images = [img for (img, contrast) in images_contrasts]
contrasts = [contrast for (img, contrast) in images_contrasts]


psf = let 
    image = copy(first(images))
    fill!(image, 0)
    cx,cy = round.(Int, mean.(axes(image)))
    image[cx, cy] = 1
    p = centered(imfilter(image, Kernel.gaussian(5), "replicate")[-10:10,-10:10])
    p ./= maximum(p)
    p
end

##
imshow2(reduce(vcat, images), clims = (-35, 35))
##
imshow2(psf, clims = (-1, 1))

##

@named b = Planet{VisualOrbit}(
    Variables(
        a = Normal(12, 3),
        e = Beta(10, 70),
        i = 0.15π,
        ωx = Normal(),
        ωy = Normal(),
        ω = (sys, pl) -> atan(pl.ωy, pl.ωx),
        Ωx = Normal(),
        Ωy = Normal(),
        Ω = (sys, pl) -> atan(pl.Ωy, pl.Ωx),
        τx = Normal(),
        τy = Normal(),
        τ = (sys, pl) -> atan(pl.τy, pl.τx)/2π,
        J = LogUniform(1, 100),
    )
)


system_images = Octofitter.ImageLikelihood(
    [
        (;band = :J, image = images[i], platescale = 10.0, epoch = times[i], contrast = contrasts[i], psf)
        for i in eachindex(images)
    ]...
    # (band=:J, image=images[1], platescale=10.0, epoch=times[1]),#, contrast=contrasts[1]),
    # (band=:J, image=images[2], platescale=10.0, epoch=times[2]),#, contrast=contrasts[2]),
    # (band=:J, image=images[3], platescale=10.0, epoch=times[3]),#, contrast=contrasts[3]),
    # (band=:J, image=images[4], platescale=10.0, epoch=times[4]),#, contrast=contrasts[4]),
)


system_vars = Variables(
    M = Normal(1.0, 0.01),
    plx = Normal(45.0, 0.0001),
)


@named system = System(system_vars, system_images, b)


##
θ = drawfrompriors(system)
sysg = generate(system, θ)
sysg.observations[1].table.image[1]|>imshow2
##
imshow2(snrmap((sysg.observations[1].table.image[4])))
##
out = Octofitteradvancedhmc(
    sysg, 0.65;
    # adaptation = 5_000,
    adaptation =  4_000,
    iterations = 15_000,
    # step_size=7e-4,
);

##
# out = Octofitter.temperedhmc(
#     system, 0.65;
#     # adaptation = 5_000,
#     adaptation = 500,
#     iterations = 5_000,
#     initial_samples=5_000,
#     # step_size=1e-3,
# );
# nothing

# Big differences between:
# GeneralisedNoUTurn
# SliceTS
# MultinomialTS
# StrictGeneralisedNoUTurn
##
(;chain, stats, adaptor) = out
plotmodel(chain, system)#, lims=1000)

##
ra, dec = projectpositions(chain, :b, last(times))
i = DirectImage(images[end])
i.platescale = 10
imshow(i, skyconvention=true)
ii = rand(eachindex(ra), 500)
scatter!(ra[ii], dec[ii], markerstrokewidth=0, color=:blue, ms=3, alpha=0.5)
##
function snr(phot)
    low, med, up = quantile(reshape(phot, :), (0.33, 0.5, 0.67))
    st = std(reshape(phot, :))
    return med / (up - low)
end
snr(chain["b[J]"])

##




##
ra, dec = projectpositions(chain, :b, first(times))
pa1 = atan.(dec, ra)
sep1 = sqrt.(ra .^ 2 .+ dec .^ 2)
##
ra, dec = projectpositions(chain, :b, last(times))
pa2 = atan.(dec, ra)
sep2 = sqrt.(ra .^ 2 .+ dec .^ 2)

##
table = (;
    # pa1 = rad2deg.(pa1),
    # sep1,
    # pa2 = rad2deg.(pa2),
    # sep2,
    a = vec(chain["b[a]"]),
    J = vec(chain["b[J]"]),
    e= vec(chain["b[e]"]),
    i=rad2deg.(rem.(vec(chain["b[i]"]), 2pi, RoundDown)),
    Ω=rad2deg.(rem.(vec(chain["b[Ω]"]), 1pi, RoundDown)),
    ω=rad2deg.(rem.(vec(chain["b[ω]"]), 2pi, RoundDown)),
    τ=         rem.(vec(chain["b[τ]"]),1,RoundDown),
)
labels = [
    # "\\theta_1",
    # "r_1",
    # "\\theta_2",
    # "r_2",
    "a",
    "J",
    "e",
    "i",
    "\\Omega",
    "\\omega",
    "\\tau",
]
units = [
    # "(\\degree)",
    # "(mas)",
    # "(\\degree)",
    # "(mas)",
    "(au)",
    "(arb.)",
    "",
    "(\\degree)",
    "(\\degree)",
    "(\\degree)",
    "",
]
# labels, units, 
using PairPlots
# PairPlots.corner(table, plotscatter=false, hist2d_kwargs=(;color=:white, nbins=15), hist_kwargs=(;nbins=15,framestyle=:none))
PairPlots.corner(table, labels, units)

##
using StatsBase
h = fit(Histogram, (dec,ra), (-1000:100:1000,-1000:100:1000))
contour(h.weights)#, levels=[quantile(h.weights[:], 0.39)])
##
##
[count(getproperty.(s, :numerical_error)) for s in stats]
[mean(getproperty.(s, :acceptance_rate)) for s in stats]
##
function orbsonimg(args...; kwargs...)
    plot()
    orbsonimg!(args...; kwargs...)
end
function orbsonimg!(planet, image, N = 100; kwargs...)
    sampled = sampleorbits(planet, N)
    i = DirectImage(image)
    i.PLATESCALE = 10.0
    imshow!(i; skyconvention = true, kwargs...)
    # plot!(sampled,color=:white,alpha=5/N,label="")
    plot!(sampled, color = :white, label = "",)
    xlims!(i.PLATESCALE .* extrema(axes(i, 1)))
    ylims!(i.PLATESCALE .* extrema(axes(i, 2)))

end


##
imshow(images[3], clims = (-2, 10))
# orbsonimg(chains.planets[1],images[2],clims=(-2,10))
orbsonimg(chains.planets[1], mean(images), 500, legend = nothing, clims = (-0.2, 1), color = :magma, dpi = 200)
plot!(truths, color = :black, lw = 2, ls = :dash, label = "")
# orbsonimg(chainsh[2].planets[1],images)
# writefits("tmp.fits", images[2])
##
plot()
orbsonimg!(chainsh[1].planets[1], mean(images), 1, legend = nothing, clims = (-0.2, 1), color = :magma, dpi = 200,)
for c in chainsh
    plot!(sampleorbits(c.planets[1], 125), color = :white, label = "", alpha = 0.1)
end
plot!(truths, color = :black, lw = 2, ls = :dash, label = "", colorbar = nothing)

##


function snr(phot)
    low, med, up = quantile(reshape(phot, :), (0.33, 0.5, 0.67))
    st = std(reshape(phot, :))
    return med / (up - low)
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
Octofitter.plotposterior(chains_img, 1, :a, 500)
plot!(truths, color = :black, lw = 2, ls = :dash, label = "", alpha = 0.5)
# savefig("images/readme-orbits.png")

##
Octofitter.plotposterior(chains_img, chains_img.planets[1], :i, 500, colorbartitle = "inclination (rad)", dpi = 200)
# savefig("images/readme-post-i.png")

##
Octofitter.plotposterior(chains_img, chains_img.planets[1], (:phot, :Keck_L′), 500, colorbartitle = "flux", cmap = :plasma, rev = false, clims = (0, 12), dpi = 200)
plot!(truths, color = :black, lw = 2, ls = :dash, label = "")
# savefig("images/readme-post-f.png")

##
Octofitter.plotposterior(chains.planets[1], :e, 500, colorbartitle = "eccentricity", dpi = 200)
plot!(truths, color = :black, lw = 2, ls = :dash)
# savefig("images/readme-post-e.png")

##
Octofitter.plotposterior(chains.planets[1], lw = 3, :e, 500, colorbartitle = "eccentricity", dpi = 200, cmap = :plasma, rev = false,)
plot!(truths, color = :black, lw = 2, ls = :dash)
# savefig("images/readme-post-e2.png")

##

ra, dec = projectpositions(chains.planets[1], mean(times))
histogram2d(ra, dec, aspectratio = 1, color = :plasma, background_inside = :black, framestyle = :box, xflip = true, dpi = 200, colorbartitle = ".\n\nposterior density")
xlims!(1.25 .* (-1, +1) .* max(abs.(extrema(ra))...))
ylims!(1.25 .* (-1, +1) .* max(abs.(extrema(dec))...))
xlabel!("Δ right ascension (as)")
ylabel!("Δ declination (as)")

##
mask = (-200 .< ra .< 250) .& (300 .< dec .< 800)
snr(chains.planets[1].phot.Keck_L′)
snr(chains.planets[1].phot.Keck_L′[mask])
histogram2d(ra[mask], dec[mask], aspectratio = 1, color = :plasma, background_inside = :black, framestyle = :box, cscale = :log10kep2cart, dpi = 200)
histogram2d(ra[mask], dec[mask], aspectratio = 1, color = :plasma, background_inside = :black, framestyle = :box, cscale = :log10kep2cart, dpi = 200)
# xlims!(1.25.*(-1,+1).*max(abs.(extrema(ra))...))
# ylims!(1.25.*(-1,+1).*max(abs.(extrema(dec))...))
##
ras, decs = projectpositions(chains.planets[1], mean(times))
# histogram2d(ras[5end÷6:end], decs[5end÷6:end])
histogram2d(ras, decs, aspectratio = 1, xflip = true, color = :plasma)

##
mask = 7 .< chains.planets[1].a .< 15
snr(chains.planets[1].phot.Keck_L′[mask])
mask = chains.planets[1].a .< 7
snr(chains.planets[1].phot.Keck_L′[mask])



##
table = (;
    L = chains.planets[1].phot.Keck_L′,
    chains.planets[1].a,
    chains.planets[1].i,
    chains.planets[1].e,
    tau = chains.planets[1].τ
)
corner(
    table,
    split(raw"\mathrm{flux} a i e \tau"),
    plotscatter = false,
    hist2d_kwargs = (; ylims = (0, NaN), xlims = (0, NaN)),
    scatter_kwargs = (; markersize = 0.5),
    dpi = 200,
    # lim_factor=2,
    # hist2d_kwargs=(;nbins=8),
)

##
ra, dec = projectpositions(chains.planets[1], mean(times))
pa = atan.(dec, ra)
# mask = (pa .< 0.5) .& 
#     (600 .< ra .< 800) .&
#     (-200 .< dec .< -100)
mask = vec((chains.planets[1].phot.Keck_L′ .< 6))
table = (;
    ra = ra[mask],
    dec = dec[mask],
    # pa=rad2deg.(pa[mask]),
    a = (chains.planets[1].a)[mask],
    L = (chains.planets[1].phot.Keck_L′)[mask],
    i = rad2deg.((chains.planets[1].i))[mask],
    e = (chains.planets[1].e)[mask]
)
corner(
    table,
    # split(raw"a \mathrm{flux} i e \tau"),
    plotscatter = false,
    # hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    # scatter_kwargs=(;markersize=0.5),
    # dpi=200,
    # lim_factor=2,
    hist_kwargs = (; nbins = 50),
    hist2d_kwargs = (; nbins = 50),
)



## Chains
table = (;
    a = (chains.planets[1].a),
    L = (chains.planets[1].phot.Keck_L′),
    i = rad2deg.((chains.planets[1].i)),
    e = (chains.planets[1].e)
)
corner(
    table,
    # split(raw"a \mathrm{flux} i e \tau"),
    plotscatter = false,
    # hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    # scatter_kwargs=(;markersize=0.5),
    # dpi=200,
    # lim_factor=2,
    hist_kwargs = (; nbins = 50),
    hist2d_kwargs = (; nbins = 50),
)

## Chainsh
table = (;
    a = reduce(vcat, [c.planets[1].a for c in chainsh]),
    L = reduce(vcat, [c.planets[1].phot.Keck_L′ for c in chainsh]),
    i = rad2deg.(reduce(vcat, [c.planets[1].i for c in chainsh])),
    e = reduce(vcat, [c.planets[1].e for c in chainsh])
)
corner(
    table,
    # split(raw"a \mathrm{flux} i e \tau"),
    plotscatter = false,
    # hist2d_kwargs=(;ylims=(0,NaN),xlims=(0,NaN)),
    # scatter_kwargs=(;markersize=0.5),
    # dpi=200,
    # lim_factor=2,
    hist_kwargs = (; nbins = 25),
    hist2d_kwargs = (; nbins = 25),
)
