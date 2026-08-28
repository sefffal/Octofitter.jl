# Fit GRAVITY-WIDE Data

## Background
Octofitter has support for directly fitting GRAVITY-WIDE closure phase data, in the OI-FITS format emitted by the pipeline.
The closure phases are mapped to a set of non-redundant kernel phases. All spectral channels are modelled separately per exposure.

!!! note
    GRAVITY modelling is supported in Octofitter via the extension package OctofitterInterferometry.
    While v9 is an unregistered prerelease it must be installed from the `v2` branch, in the same `Pkg.add` as Octofitter itself —
    see [Installation](@ref install).

The astrophysical model is a sum of point sources:

```math
V(u,v) = \frac{\sum_j f_j \, e^{-2\pi i (u\,\Delta\alpha^*_j + v\,\Delta\delta_j)}}{\sum_j f_j}
```

over the bodies named in `targets`, at whatever offsets the trajectory puts them. There is no primary — the host is one source among the others, with its own `flux_<band>` variable — and a source may orbit any body, so a moon or a wide component of a hierarchical system is expressible.

!!! note "`GRAVITYWideKPObs` is a preset, not a type"
    The kernel-phase machinery is an option on
    [Fitting Interferometric Observables](@ref)'s `InterferometryObs`, and
    `GRAVITYWideKPObs(...)` is a *function* that returns one with
    `kernel_phases=true, fiber_coupling=true` set. The equivalent spelling on the
    merged type is

    ```julia
    InterferometryObs(data; kernel_phases=true, fiber_coupling=true, targets=..., ref=..., band=...)
    ```

    `kp_correlation=false` drops the spectral correlation within a kernel phase while
    keeping the projection.

    Consequence for code that tested the type: `obs isa GRAVITYWideKPObs` no longer
    compiles. Test `obs isa InterferometryObs` instead.

Interferometer data is almost always multi-modal, so a single gradient-based chain will settle into one mode.
Multi-wavelength GRAVITY-WIDE data with multiple epochs is fairly expensive to model (can take on the order of 1ms per likelihood evaluation), so after running some tests locally, one should consider using a compute cluster.
You will probably want on the order of 30 cores and 1-5 days, depending on the scale of the problem.

## Process



```julia
using Octofitter
using OctofitterInterferometry
using Distributions
using CairoMakie
using PairPlots
```

To model orbits / brightness of a companion from GRAVITY-WIDE data, use the following observation.
For this example, we won't consider a full orbit. We will just sample from 2D separation and position angle coordinates, so we start with the bodies:

```julia
const OBS_EPOCH = 60676.00776748842

A = Body(
    name="A",
    variables=@variables begin
        mass = 1.0            # M⊙
        flux_K = 1.0          # the host is an ordinary source; pinning it to 1
    end                       # makes every other body's flux a contrast ratio
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        flux_K ~ Uniform(0, 1)   # contrast ratio against the host

        sep ~ Uniform(0, 10)     # mas
        pa ~ Uniform(0, 2pi)     # radians

        # A face-on circular orbit placed at (sep, pa) at the epoch of the
        # exposure: a fixed position, written as a degenerate orbit.
        a = sep / system.plx
        e = 0.0
        i = 0.0
        ω = 0.0
        Ω = 0.0
        θ = pa
        epoch = $OBS_EPOCH
    end
)

vis_obs = GRAVITYWideKPObs(
    (;
        filename="./GRAVI.2025-01-01T00:11:11.111_dualscivis.fits",
        epoch=OBS_EPOCH,
        wavelength_min_meters=2.025e-6,
        wavelength_max_meters=2.15e-6,
        jitter=:kp_jit,
        kp_Cy=:kp_Cy,
    ),
    # Add more exposures / epochs here if desired...
    targets = (A, b),      # every source in the visibility sum, host included
    ref     = A,           # phase centre
    band    = :K,          # which `flux_<band>` to read
    variables=@variables begin
        kp_jit ~ Uniform(0, 180)   # kernel phase jitter, degrees
        kp_Cy  ~ Uniform(-1, 1)    # spectral correlation parameter
    end
)

sys = System(
    name="sys",
    bodies=[A, b],
    observations=[vis_obs],
    variables=@variables begin
        plx = 173.5740
    end
)

model = Octofitter.LogDensityModel(sys, verbosity=4)
```

- `filename` is the path from your current working directory to the GRAVITY OI-FITS file.
- `epoch` is the average time of the exposure in MJD (not the start of the exposure!).
- `wavelength_min_meters` and `wavelength_max_meters` are wavelength cutoffs for the data if you want to restrict to only include some channels from the file (optional)
- `jitter` is a symbol giving the name of a kernel phase jitter variable to use for this exposure. This can be defined in the observation's `variables` block or in the system's.
- `kp_Cy` is a symbol giving the name of the spectral correlation variable to use for this exposure (optional). This can also live in either block.

The **fluxes are body variables**, not a `flux` vector on the observation. A second companion is a third entry in `targets` plus its own `flux_K` on its own body.

!!! note "OI-FITS files are read directly"
    OI-FITS is read through FITSIO rather than through OIFITS.jl, whose 2.0 release moved
    to AstroFITS and dropped the `read(OIDataBlock, ::FITSIO.HDU)` method Octofitter used.
    Data-block selection prefers `EXTVER = 10` — GRAVITY's science-combiner convention —
    falling back to the first block of each kind, and is exposed as an optional
    `oi_extver` field on the input row.

## Two things to be aware of

!!! note "How fibre coupling is computed"
    `GRAVITYWideKPObs` sets `fiber_coupling=true`. Each source's throughput is evaluated
    at its own offset from `fiber_pointing`, which is the quantity injection efficiency
    actually depends on.

    `fiber_pointing` defaults to `Photocentre(band)`. For most GRAVITY-WIDE
    observations the fibre is actually on the host, which is spelled

    ```julia
    GRAVITYWideKPObs(...; targets=(A, b), ref=A, band=:K, fiber_pointing=A)
    ```

    For a faint companion the old and new throughputs differ by roughly the full
    coupling loss at the companion's separation, so this is not a cosmetic change.
    Passing `fiber_pointing` without `fiber_coupling=true` is an error rather than a
    silent no-op.

!!! note "Choosing `ref`"
    Closure phases, kernel phases and squared visibilities are invariant to the phase
    centre, so `ref` is a free choice — but only modulo 360°: baseline phases are
    folded into (−180°, 180°] and the triangle sum is not. Keep `ref` near the flux
    centroid: `Barycentre` for a faint companion, or `Photocentre(:K)` in general.
    `ref=A` puts the phase centre on the host.

## Sampling

It is recommended to use Pigeons parallel tempered sampling, and to use the SliceSampler explorer. This non-default option avoids calculating the gradient of the model, which is expensive in this case.
```julia
init_chain = initialize!(model)

using Pigeons
chain, pt = octofit_pigeons(model, n_chains=8, n_chains_variational=0, n_rounds=9, explorer=SliceSampler())
```

!!! tip "Parallel tempering, and avoiding the gradient"
    Two things recommend [`octofit_pigeons`](@ref) here. Its `SliceSampler` explorer
    only ever evaluates the log density, so the expensive gradient is never formed;
    and its tempered ladder can move between the widely separated position solutions
    that a single interferometric epoch admits, which HMC cannot.

    `n_chains_variational=0` drops the variational leg: a Gaussian reference is a poor
    description of this posterior, so it buys nothing here.

    [`octofit`](@ref) will run if you prefer HMC, but it computes gradients and does not
    jump between modes. For a single-epoch GRAVITY-WIDE data set the grid search below
    remains the more trustworthy analysis either way.


```julia
fig = Figure()
ax = Axis(
    fig[1,1],
    xreversed=true,
    autolimitaspect=1
)
xlims!(ax, 10,-10)
ylims!(ax, -10,10)
x = vec(chain[:b_sep]) .* sin.(vec(chain[:b_pa]))
y = vec(chain[:b_sep]) .* cos.(vec(chain[:b_pa]))
scatter!(ax,x,y)
fig
```



Since we are only considering a single epoch, we can also go ahead and generate a detection map by performing a grid-search over positions.

The flat parameter vector `model.ℓπcallback` consumes is ordered: system priors first, then each body's priors in declaration order, then each observation's. For the model above that is `(b_flux_K, b_sep, b_pa, kp_jit, kp_Cy)` — `plx`, `A.mass` and everything derived are not free parameters. You can always check the ordering with

```julia
θ = Octofitter.sample_priors(model.system)
model.arr2nt(θ)     # shows where each flat entry landed
```

```julia
ks = 0.5:0.5#0.4:0.1:0.6 #0.01:0.1:1.0
xs =  (-10:0.25:10) .+ 1e-6
ys = (-10:0.25:10)

ks_ = reshape(ks, 1,1,:)
seps = sqrt.(xs.^2 .+ ys'.^2) .+ 0ks_
pas = atan.(xs, ys') .+ 0ks_
ks__ = ks_ .+ 0seps
LL = fill(NaN,length(xs), length(ys), length(ks))
jit = 12.6
Cy = 0.02
@time Threads.@threads for i in eachindex(LL)
    sep = seps[i]
    if sep > 10
        continue
    end
    pa = rem2pi(pas[i], RoundDown)
    K = ks__[i]
    LL[i] = model.ℓπcallback(model.link((K,sep,pa,jit,Cy)),sampled=false)
end
```

```julia
fig = Figure()
ax = Axis(
    fig[1,1],
    xreversed=true,
    autolimitaspect=1,
    backgroundcolor="#222",
    title="spec corr = 0.02"
)

N_dat = length(vis_obs.table.epoch)*length(vis_obs.table.eff_wave[1])*3 # 3 kern phases
N_param = 3
χ²_max = maximum(LL,dims=3)[:,:] ./ (N_dat + N_param)
h = heatmap!(ax,
    xs,ys, χ²_max,
    # LL[:,:,1],
    colormap=:magma,
    colorrange=(quantile(filter(isfinite,χ²_max),0.85), maximum(filter(isfinite,χ²_max)))
    # colorrange=(quantile(filter(isfinite,LL),0.85), maximum(filter(isfinite,LL)))
    # colorrange=(maximum(filter(isfinite,χ²_max))-3, maximum(filter(isfinite,χ²_max)))
)


Colorbar(fig[1,2],h,label="Log-Posterior Density")

fig
```



This single-epoch model can then be extended by replacing the fixed-position parameterization with a real Keplerian orbit — which just means declaring a different set of variables on the same body:
```julia
b_orbit = Body(
    name="b",
    about=A,
    variables=@variables begin
        flux_K ~ Uniform(0, 1)
        a ~ Uniform(0, 0.1)
        e ~ Uniform(0.0, 0.99)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 60676.0
    end
)
```
Give `A` a real mass prior at the same time (`mass ~ truncated(Normal(1.0, 0.1), lower=0.1)`), since an orbit's period depends on it.
