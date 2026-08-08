# Fit with a Thiele-Innes Basis

This example shows how to fit relative astrometry using a Thiele-Innes orbital parameterization instead of the traditional Campbell parameterization used in other tutorials. Thiele-Innes replaces the semi-major axis `a` and the angular elements (inclination `i`, longitude of ascending node `Ω`, and argument of periastron `ω`) with the four Thiele-Innes constants (A, B, F, G). This avoids the coordinate singularities where `ω`, `Ω`, and `tp` become poorly defined as eccentricity and/or inclination approach zero.

!!! tip "When to consider Thiele-Innes"
    The Thiele-Innes parameterization may be useful when:
    - The orbit is nearly circular (e < 0.1) and/or nearly face-on (i near 0° or 180°)
    - You want to avoid angular coordinate singularities

    Both parameterizations should give consistent results for the physical orbital parameters. Choose based on your preference or specific analysis needs.

There is no separate orbit type for this: Thiele-Innes is written as what it is — four
sampled constants plus one derived line converting them into the elements PlanetOrbits
uses. [`PlanetOrbits.ThieleInnes`](https://sefffal.github.io/PlanetOrbits.jl/dev/api/)
returns a NamedTuple `(; a, i, ω, Ω)`; you still supply the shape (`e`) and phase
elements yourself.

At the end, we will convert our results back into the Campbell parameterization to compare.

```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions
using PlanetOrbits

astrom_dat = Table(;
    epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
    ra    = [-505.7637580573554, -502.570356287689, -498.2089148883798, -492.67768482682357, -485.9770335870402, -478.1095526888573, -469.0801731788123, -458.89628893460525],
    dec   = [-66.92982418533026, -37.47217527025044, -7.927548139010479, 21.63557115669823, 51.147204404903704, 80.53589069730698, 109.72870493064629, 138.65128697876773],
    σ_ra  = [10, 10, 10, 10, 10, 10, 10, 10],
    σ_dec = [10, 10, 10, 10, 10, 10, 10, 10],
    cor   = [0, 0, 0, 0, 0, 0, 0, 0]
)

A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.2, 0.1), lower=0.1)
    end
)

planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0
        e ~ Uniform(0.0, 0.5)
        # Thiele-Innes constants A, B, F, G are in milliarcseconds (not AU like semi-major axis).
        # Set the prior width to encompass the expected angular separation of your target.
        # A rough guide: if your astrometry spans ~500 mas, use Normal(0, 1000) or similar.
        A ~ Normal(0, 1000) # milliarcseconds
        B ~ Normal(0, 1000) # milliarcseconds
        F ~ Normal(0, 1000) # milliarcseconds
        G ~ Normal(0, 1000) # milliarcseconds

        # Convert to the size and orientation elements. `plx` is needed because
        # the constants are in angular units.
        ti = PlanetOrbits.ThieleInnes(; A, B, F, G, plx=system.plx)
        a = ti.a
        i = ti.i
        ω = ti.ω
        Ω = ti.Ω

        θ ~ UniformCircular()
        epoch = 50000.0   # reference epoch for θ. Choose an MJD date near your data.
    end
)

astrom_obs = RelAstromObs(
    astrom_dat;
    target = planet_b,
    ref = A,
    name = "GPI",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)

sys = System(
    name="TutoriaPrime",
    bodies=[A, planet_b],
    observations=[astrom_obs],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
```

!!! note "The node is ambiguous by ±180°"
    `(ω, Ω)` and `(ω+π, Ω+π)` give *identical* Thiele-Innes constants, so the inverse is
    genuinely two-valued: astrometry alone cannot distinguish the ascending node from the
    descending one. `PlanetOrbits.ThieleInnes` returns the branch with `Ω ∈ [0, π)`. Radial
    velocities break the tie; if you have them and they prefer the other node, use
    `ω + π` and `Ω + π`.

    If you do not know the branch in advance, the ±180° choice is *discrete* and HMC
    cannot cross it, so do not try to sample it. Fit the two branches as two models —
    identical but for the `+ π` — and compare them with
    [Bayesian evidence](@ref bayesian-evidence). Note also that once radial velocities
    are in the fit, the singularity that motivates Thiele-Innes in the first place is
    usually no longer the binding constraint, so sampling Campbell elements directly is
    often the simpler model.

Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain = initialize!(model)
octoplot(model, init_chain)
```

We now sample from the model as usual:
```@example 1
results = octofit(model)
```

!!! note
    The Thiele-Innes parameterization may reveal more complex posterior structure (e.g., multimodality) that Campbell masks through its angular parameterization. If your corner plot shows unexpected bimodality in the `A`, `B`, `F`, `G` parameters, this may reflect genuine orbital ambiguities rather than sampling issues.

We now display the results:
```@example 1
octoplot(model,results)
```

```@example 1
octocorner(model, results, small=false)
```

## Working in Campbell Elements

Because the Campbell elements are ordinary derived variables of this model, they are
already in the chain — there is no conversion step:

```@example 1
table = (;
    B_a = results[:b_a][:],
    B_e = results[:b_e][:],
    B_i = rad2deg.(results[:b_i][:]),
)
pairplot(table)
```

If you want the solved orbits themselves — for example to query them at new epochs, or to
read off elements a model did not name — rebuild the `PlanetOrbits.System` for a draw:

```@example 1
posys = construct_system(model, results, 1)
(; a = semimajoraxis(posys, 1),
   e = eccentricity(posys, 1),
   i = rad2deg(inclination(posys, 1)),
   P_days = period(posys, 1))
```

The inverse conversion is also available: `PlanetOrbits.thieleinnes(posys, 1; plx=50.0)`
returns the `(A, B, F, G)` constants of hierarchy row 1 in milliarcseconds.
