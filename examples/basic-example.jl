

using DirectDetections


@named planet_e = DirectDetections.Planet(
    Priors(
        a = TruncatedNormal(15, 2, 0, Inf),
        e = TruncatedNormal(0.1, 0.1, 0, 0.5),
        τ = Normal(0.5,0.5),
        ω = Normal(π,π),
        i = Normal(deg2rad(30), deg2rad(4)),
        Ω = Normal(deg2rad(160), deg2rad(40)),
        Keck_L′ = TruncatedNormal(0,2e-4,0,Inf),
    ),
    Astrometry(
        (epoch=56578.,  ra=382.6sind(265.13 ), dec=382.6cosd(265.13 ), σ_ra=2., σ_dec=2.),
        (epoch=57650.,  ra=384.8sind(281.68 ), dec=384.8cosd(281.68 ), σ_ra=2., σ_dec=2.),
    )
)

system_images = DirectDetections.Images(
    (
        band=:Keck_L′,
        image=centered(readfits("abc.fits")),
        platescale=9.971,
        epoch=51231.0,
    )
)
system_priors = Priors(
    μ = Normal(1.53, 0.01),
    plx =Normal(24.2, 0.02),
)

@named hr8799 = System(
    system_priors,  
    system_images,
    planet_e,
)

##
chains, stats = DirectDetections.hmc(
    hr8799;
    numwalkers=1,
    burnin=2_000,
    numsamples_perwalker=10_00
);

plotmodel(chains[1], hr8799)
