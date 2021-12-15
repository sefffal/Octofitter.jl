using DirectDetections, Distributions
using MCMCChains

# Create a function mapping (age_Myr, mass_Mjup) -> temp_K
const cooling_tracks = sonora_cooling_interpolator()

# Create functions mapping (temp_K, mass_Mjup) -> absolute magnitude
const sonora_temp_mass_Z = sonora_photometry_interpolator(:MKO_Z)
const sonora_temp_mass_J = sonora_photometry_interpolator(:MKO_J)
const sonora_temp_mass_L = sonora_photometry_interpolator(:Keck_L′)

##
@named b = DirectDetections.Planet(
    Priors(
        mass = Uniform(0.6, 75),
        a = Uniform(1, 30),
        e = Beta(1.2, 10),
        i = Normal(0.1, 0.2),
        ω = Uniform(0, 2pi),
        Ω = Uniform(0, pi),
        τ = Uniform(0, 1),
    ),
    # Flux will be calculated from planet mass and system age
    Derived(
        Z = (sys, pl) -> sonora_temp_mass_Z(cooling_tracks(sys.age, pl.mass), pl.mass),
        J = (sys, pl) -> sonora_temp_mass_J(cooling_tracks(sys.age, pl.mass), pl.mass),
        L = (sys, pl) -> sonora_temp_mass_L(cooling_tracks(sys.age, pl.mass), pl.mass),
        # We inculde teff in the chains to ease analysis (though it can be recalculated after the fact using age and mass)
        teff = (sys, pl) -> cooling_tracks(sys.age, pl.mass),
    ),
    Photometry(
        (band = :Z, phot=15.0, σ_phot=3.),
        (band = :J, phot=13.5, σ_phot=0.5),
        (band = :L, phot=11.0, σ_phot=1.0)
    ),
    Astrometry(
        (epoch=mjd("2016-12-15"), ra=0.133*1e3, dec=-0.174*1e3, σ_ra=0.007*1e3, σ_dec=0.007*1e3),
        (epoch=mjd("2017-03-12"), ra=0.126*1e3, dec=-0.176*1e3, σ_ra=0.004*1e3, σ_dec=0.004*1e3),
        (epoch=mjd("2017-03-13"), ra=0.127*1e3, dec=-0.172*1e3, σ_ra=0.004*1e3, σ_dec=0.004*1e3),
        (epoch=mjd("2018-02-08"), ra=0.083*1e3, dec=-0.133*1e3, σ_ra=0.010*1e3, σ_dec=0.010*1e3),
        (epoch=mjd("2018-11-28"), ra=0.058*1e3, dec=-0.122*1e3, σ_ra=0.010*1e3, σ_dec=0.020*1e3),
        (epoch=mjd("2018-12-15"), ra=0.056*1e3, dec=-0.104*1e3, σ_ra=0.008*1e3, σ_dec=0.008*1e3),
    )
)

@named HD12345 = System(
    Priors(
        # age = TruncatedNormal(40, 15, 5, Inf),
        age = Uniform(5, 100),
        μ = Normal(1.61, 0.05),
        plx = gaia_plx(gaia_id=756291174721509376),
    ),
    # ProperMotionAnomHGCA(gaia_id=756291174721509376),
    ProperMotionAnomHGCA(gaia_id=2832463659640297472), #HR8799...
    b,
)

##
chains = DirectDetections.hmc(
    HD12345, 0.85,
    adaptation =   1_000,
    iterations = 20_000,
)


##
gelmandiag(chain)
##
plotmodel(chain,color=:mass,pmascatter=:mass)
##
using PairPlots
PairPlots.corner(chain)

##
table = (
    age=chain["age"],
    mass=chain["b[mass]"],
    sma=chain["b[a]"],
    Z=chain["b[Z]"],
    J=chain["b[J]"],
    L=chain["b[L]"],
)
PairPlots.corner(table)
# PairPlots.corner(table, filterscatter=false, plotcontours=false)


##
p = plot(yflip=true,yguide="abs. mag")
for planet in chain.info.model.planets
    phot = planet.photometry
    N= 1500
    ii = rand(1:size(chain,1)*size(chain,3),N)
    Zs = chain["b[Z]"][ii]
    Js = chain["b[J]"][ii]
    Ls = chain["b[L]"][ii]
    plot!(vcat(Zs', Js', Ls'), color=:black, lw=1,  alpha=0.05, label="")
    # scatter!(Ls', color=:black, lw=1,  alpha=0.05, label="")
    scatter!(vcat(Zs', Js', Ls'), color=:black, lw=1, ms=4, alpha=0.01, label="")
    scatter!(
        collect(phot.phot),
        yerr=collect(phot.σ_phot),
        xticks=(collect(eachindex(phot.phot)), string.(phot.band)),
        color=:red,
        markerstrokecolor=:red,
        label="Measured"
    )
end
p 