# [Priors](@id priors)

All parameters to your model must have a prior defined.
You may provide any continuous, univariate distribution from the Distributions.jl.
A few useful distributions include:

* `Normal`
* `Uniform`
* `LogNormal`
* `LogUniform`
* `truncated(Normal(μ, σ), lower=…, upper=…)`
* `VonMises`

This package also defines the [`Sine`](@ref) distribution for e.g. inclination priors and [`UniformCircular`](@ref) for periodic variables.
Internally, `UniformCircular()` creates two standard normal variables and finds the angle between them using `arctan`. This allows the sampler to smoothly cycle past the ends of the domain. You can specify a different circular domain than (0,2pi) by passing the size of the domain e.g. `x = UniformCircular(1.0)`.

!!! note
    A `UniformCircular` variable `Ω` appears in the chain as three columns: the two
    helper variables `Ωx` and `Ωy`, and the angle `Ω` derived from them. If you want to
    supply a starting guess to `initialize!` you must give `Ωx` and `Ωy`, not `Ω` —
    or replace the prior with `Ω ~ Uniform(0, 2pi)`.

The VonMises distribution is notable but not commonly used. It is the analog of a normal distribution defined on a circular domain (-π, +π). If you have a Gaussian prior on an angular parameter, a Von Mises distribution is probably more appropriate.

Behind the scenes, Octofitter remaps your parameters to unconstrained domains using the Bijectors.jl (and corrects the priors accordingly). This is essential for good sampling efficiency with HMC based samplers.

This means that e.g. if you define the eccentricity prior as `e ~ Uniform(0,0.5)`, the sampler will actually generate values across the whole real line and transform them back into the `[0,0.5]` range before evaluating the orbit.
**It is therefore essential that your priors do not include invalid domains.**

For example, setting `a ~ Normal(3,2)` will result in poor sampling efficiency as sometimes negative values for semi-major axis will be drawn (especially if you're using the parallel tempered sampler).

Instead, for parameters like semi-major axis, eccentricity, parallax, and masses, you should truncate any distributions that have negative tails.
This can easily be accomplished with `truncated(dist, lower=…, upper=…)` for any arbitrary distribution.


## Kernel Density Estimate Priors

Octofitter has support for sampling from smoothed kernel density estimate priors. These are non-parametric distributions fit to a 1D dataset consisting of random draws. This is one way to include the output of a different model as the prior to a new model. That said, it's usually best to try and incorporate the model directly into the code. There are a few examples on GitHub of this, including atmosphere model grids, cooling tracks, etc.

### Using a KDE
First, we will generate some data. In the real world, you would load this data eg. from a CSV file.
```@example 1
using Octofitter, Distributions

# create a smoothed KDE estimate of the samples from a 10+-1 gaussian
kde = Octofitter.KDEDist(randn(1000).+10)
```

Note that in Octofitter the KDE will have its support truncated to the minimum and maximum values that occur in your dataset, ie. it doesn't allow for infinite long tails.

Now add it to your model as a prior:
```@example 1
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
        a ~ kde # Sample from the KDE here
        e ~ Uniform(0.0, 0.99)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 50000.0
    end
)

sys = System(
    name="Tutoria",
    bodies=[A, planet_b],
    observations=[],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
chain = octofit(model)
```

We now examine the posterior and verify that it matches our KDE prior:
```@example 1
using Statistics
dat = chain[:b_a][:]
@show mean(dat) std(dat)
nothing # hide
```

## Priors as `@variables` lines

Any line in a `@variables` block of the form `expression ~ distribution` is a prior term.
That includes expressions of several variables, which is how you constrain a *derived*
quantity — for example the mutual inclination of two planets:

```julia
sys = System(name="coplanar", bodies=[A, b, c], variables=@variables begin
    plx ~ Uniform(1, 100)
    mut_inc = acos(cos(b.i)*cos(c.i) + sin(b.i)*sin(c.i)*cos(b.Ω - c.Ω))
    mut_inc ~ truncated(Normal(0, deg2rad(10)), lower=0)
end)
```

Lines like this are *prior-shaped terms*: they carry no data, so they are excluded from
cross-validation and pointwise likelihood accounting. Reach for an observation type
(e.g. [`PhotometryObs`](@ref)) instead when the numbers you are comparing against are
measurements.

## Dynamical priors

Several multi-body constraints are available as prior terms, which go in the **system's**
`observations=` list:

```julia
observations = [
    astrom,
    OrbitOrderPrior(b, c),                          # a_b < a_c < …
    NonCrossingPrior(bodies=(b, c)),                # apsides may not cross
    LimitClosestApproachAUPrior(0.2, 1.0; bodies=(b, c)),
    HillStabilityPrior(bodies=(b, c)),
]
```

Naming `bodies=` restricts each prior to the hierarchy rows that place those bodies. The
default — no list — covers *every* row, which is worth thinking about in a hierarchical
system: a 2+2 quadruple's rows include the wide binary orbit, and
comparing an inner pair's apsides against it is not a meaningful test. See
[Resonant Co-Planar Model](@ref fit-coplanar) for a worked example.

## Observable Based Priors

Octofitter implements observable-based priors from O'Neil 2019. Wrap an existing
observation in [`ObsPriorONeil2019`](@ref) and list **only the wrapper** in
`observations=`:

```julia
using Octofitter, Distributions

astrom_dat = Table(
    epoch=[mjd("2020-12-20")], 
    ra=[400.0], 
    σ_ra=[5.0], 
    dec=[400.0], 
    σ_dec=[5.0]
)

A = Body(name="A", variables=@variables begin
    mass = system.M_tot
end)

planet_b = Body(name="b", about=A, variables=@variables begin
    mass_jup ~ LogUniform(0.01, 100)          # [Mjup]
    mass = mass_jup * mjup                    # [M⊙] — masses are solar throughout
    # For using with ObsPriors: sample the period, derive the semi-major axis.
    # `P` itself is an orbit element, so the sampled variable needs another name.
    P_yr ~ Uniform(0.001, 1000)              # years
    a = cbrt(system.M_tot * P_yr^2)          # [AU]

    e ~ Uniform(0.0, 1.0)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()

    θ ~ UniformCircular()
    epoch = 58849.0   # reference epoch for θ. Choose an MJD date near your data.
end)

astrom_obs = RelAstromObs(astrom_dat; target=planet_b, ref=A, name="rel astrom. 1")

sys = System(
    name="System1",
    bodies=[A, planet_b],
    observations=[ObsPriorONeil2019(astrom_obs)],
    variables=@variables begin
        plx ~ Normal(21.219, 0.060)
        M_tot ~ truncated(Normal(1.1, 0.2), lower=0.1)
    end
)
```

See [Observable-Based Priors](@ref obs-priors) for a complete worked example, including the case of
wrapping a radial velocity observation, where you must name the orbit explicitly.
