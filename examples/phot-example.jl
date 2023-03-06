using Octofitter, Distributions, Plots
using MCMCChains

# Create a function mapping (age_Myr, mass_Mjup) -> temp_K
const cooling_tracks = sonora_cooling_interpolator()

# Create functions mapping (temp_K, mass_Mjup) -> absolute magnitude
const sonora_temp_mass_Z = sonora_photometry_interpolator(:MKO_Z)
const sonora_temp_mass_J = sonora_photometry_interpolator(:MKO_J)
const sonora_temp_mass_L = sonora_photometry_interpolator(:Keck_L′)

##
@named b = Planet{VisualOrbit}(
    Variables(
        Z = (sys, pl) -> sonora_temp_mass_Z(cooling_tracks(sys.age, pl.mass), pl.mass),
        J = (sys, pl) -> sonora_temp_mass_J(cooling_tracks(sys.age, pl.mass), pl.mass),
        L = (sys, pl) -> sonora_temp_mass_L(cooling_tracks(sys.age, pl.mass), pl.mass),
    ),
    Priors(
        # mass = Uniform(0.53, 95)
        # mass = TruncatedNormal(7, 7, 1, 80)
        # mass = NNor(1, 50)
        mass = Uniform(2, 60)
    ),
    # PhotometryLikelihood 
    PhotometryLikelihood(
        (band = :Z, phot=16.0, σ_phot=3.),
        (band = :J, phot=17.0, σ_phot=2.),
        (band = :L, phot=12.0, σ_phot=3.0)
    )
)

@named HD12345 = System(
    Priors(
        age = TruncatedNormal(40, 15, 5, Inf)
    ),
    b,
)

##
chain = Octofitteradvancedhmc(
    HD12345, 0.85,
    adaptation =   1_000,
    iterations =  50_000,
)


##
gelmandiag(chain)
##
using PairPlots
PairPlots.corner(chain)

##
table = (
    mass=chain["b[mass]"],
    age=chain["age"],
)
PairPlots.corner(table)


##
p = plot()
for planet in chain.info.model.planets
    phot = planet.photometry
    N= 1500
    ii = rand(1:size(chain,1)*size(chain,3),N)
    Zs = chain["b[Z]"][ii]
    Js = chain["b[J]"][ii]
    Ls = chain["b[L]"][ii]
    plot!(vcat(Zs', Js', Ls'), color=:black, lw=1,  alpha=0.05, label="")
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
