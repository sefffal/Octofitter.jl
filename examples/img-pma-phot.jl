using DirectDetections, Distributions
using MCMCChains

# Create a function mapping (age_Myr, mass_Mjup) -> temp_K
const cooling_tracks = sonora_cooling_interpolator()

# Create functions mapping (temp_K, mass_Mjup) -> absolute magnitude
const sonora_temp_mass_Z = sonora_photometry_interpolator(:MKO_Z)
const sonora_temp_mass_J = sonora_photometry_interpolator(:MKO_J)
const sonora_temp_mass_L = sonora_photometry_interpolator(:Keck_L′)

##
@named f = DirectDetections.Planet(
    Priors(
        τ = Uniform(0, 1),
    ),
    # Flux will be calculated from planet mass and system age
    Derived(
        e = Returns(0),
        i = Returns(deg2rad(30)),
        ω = Returns(0),
        Ω = Returns(deg2rad(90)),
        a = (sys,pl) -> ∛((0.5sys.P1)^2),
        mass = (sys,pl) -> sys.mass,
    ),  
)

@named e = DirectDetections.Planet(
    Priors(
        τ = Uniform(0, 1),
    ),
    # Flux will be calculated from planet mass and system age
    Derived(
        e = Returns(0),
        i = Returns(deg2rad(30)),
        ω = Returns(0),
        Ω = Returns(deg2rad(90)),
        a = (sys,pl) -> ∛(sys.P1^2),
        mass = (sys,pl) -> sys.mass,
    ),  
)

@named d = DirectDetections.Planet(
    Priors(
        τ = Uniform(0, 1),
    ),
    # Flux will be calculated from planet mass and system age
    Derived(
        e = Returns(0),
        i = Returns(deg2rad(30)),
        ω = Returns(0),
        Ω = Returns(deg2rad(90)),
        a = (sys,pl) -> ∛((2sys.P1)^2),
        mass = (sys,pl) -> sys.mass,
    ),  
)

@named c = DirectDetections.Planet(
    Priors(
        τ = Uniform(0, 1),
    ),
    # Flux will be calculated from planet mass and system age
    Derived(
        e = Returns(0),
        i = Returns(deg2rad(30)),
        ω = Returns(0),
        Ω = Returns(deg2rad(90)),
        a = (sys,pl) -> ∛((4sys.P1)^2),
        mass = (sys,pl) -> sys.mass,
    ),  
)

@named b = DirectDetections.Planet(
    Priors(
        τ = Uniform(0, 1),
    ),
    # Flux will be calculated from planet mass and system age
    Derived(
        e = Returns(0),
        i = Returns(deg2rad(30)),
        ω = Returns(0),
        Ω = Returns(deg2rad(90)),
        a = (sys,pl) -> ∛((8sys.P1)^2),
        mass = (sys,pl) -> 0.5*sys.mass,
    ),  
)

@named HR8799 = System(
    Priors(
        # age = TruncatedNormal(40, 15, 5, Inf),
        # age = Uniform(5, 100),
        μ = Normal(1.5, 0.05),
        plx = gaia_plx(gaia_id=2832463659640297472),
        # mass = Uniform(0.6, 75),
        mass = TruncatedNormal(7,1,0,Inf),
        P1 = Uniform(40, 100)
    ),
    ProperMotionAnomHGCA(gaia_id=2832463659640297472),
    f, e, d, c, b
)

##
chains = DirectDetectionsadvancedhmc(
    HR8799, 0.85,
    MCMCThreads(),
    num_chains=4,
    adaptation =   1_000,
    iterations = 30_000,
)


##


# Histograms
margin_hist(chains, "e[a]", "e[τ]"; color=:grey)

##
margin_hist(chains, "d[a]", "mass"; color=:grey)

##
margin_hist(chains, "e[a]", "mass"; color=:grey)
margin_hist(chains, "d[a]", "mass"; color=:grey)

##


fig =  Figure(resolution=(700,700))
ax_hist = Axis(
    fig[1,1],
    ylabel="posterior density",
    xlabel="e[a]",
)
hideydecorations!(ax_hist)
hist_prop!(ax_hist, chains["e[a]"], color=:black, label="e[a]")
hist_prop!(ax_hist, chains["f[a]"], color=:red, label="f[a]")
Legend(fig[1,2], ax_hist)
xlims!(ax_hist, low=0)
fig
##
plot_all(chains)