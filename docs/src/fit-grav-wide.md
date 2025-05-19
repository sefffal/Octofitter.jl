# Fit GRAVITY-WIDE Data

## Background
Octofitter has support for directly fitting GRAVITY-WIDE closure phase data, in the OI-FITS format emitted by the pipeline.
The closure phases are mapped to a set of non-redundant kernel phases. All spectral channels are modelled separately per exposure.

!!! note
    GRAVITY modelling is supported in Octofitter via the extension package OctofitterInterferometry. To install it, run 
    `pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterInterferometry`

The only supported astrophysical sources at this time are zero or more point sources orbitting a primary body.

Interferometer data is almost always multi-modal, requiring the use of parallel tempering.
Multi-wavelength GRAVITY-WIDE data with multiple epochs is fairly expensive to model (can take on the order of 1ms per likelihood evaluation), so one after running some tests locally, one should consider using a compute cluster.
You will probably want on the order of 30 cores and 1-5 days, depending on the scale of the problem.

## Process



```@julia
using Octofitter
using OctofitterInterferometry
using Distributions
using CairoMakie
using PairPlots
```

To model orbits / brightness of a companion from, GRAVITY-WIDE data you will want to use the following Likelihood object:

```julia
vis_like = GRAVITYWideKPLikelihood(
    (;
        filename="./GRAVI.2025-01-01T00:11:11.111_dualscivis.fits",
        epoch=60676.00776748842,
        wavelength_min_meters=2.025e-6,
        wavelength_max_meters=2.15e-6,
        spectrum_var=:K_band_spectrum,
        jitter=:kp_jit,
        kp_Cy=:kp_Cy,
    )
    # Add more exposures / epochs here if desired...
)
```

- `filename` is the path from your current working directory to the GRAVITY OI-FITS file.
- `epoch` is the average time of the exposure in MJD (not the start of the exposure!).
- `wavelength_min_meters` and `wavelength_max_meters` are wavelength cutoffs for the data if you want to restrict to only include some channels from the file (optional)
- `spectrum_var` is the variable name of the spectrum in the model. This is an array of contrast-ratios between the companion and the primary.
- `jitter` is a symbol giving the name of a  kernel phase jitter variable in the model to use for this exposure.
- `kp_Cy` is a symbol giving the name of the spectral correlation variable to use for this exposure (optional)



For this example, we won't consider a full orbit. We will just sample from 2D separation and position angle coordinates. To this end, we will use the `FixedPosition` parameterization for the planet:
```@julia

@planet b Visual{Octofitter.FixedPosition} begin
    sep ~ Uniform(0, 10) # mas
    pa ~ Uniform(0,2pi)

    # We will sample using a flat contrast ratio versus spectrum
    K ~ Uniform(0, 1) # contrast ratio variable
    K_band_spectrum = fill(b.K, $(length(vis_like.table.eff_wave[1])))

    # You could alternatively let the flux vary by wavelength:
    # K_band_spectrum = Product(fill(Uniform(0,1),length(vis_like.table.eff_wave[1])))
end

@system sys begin
    plx = 173.5740
    # Kernel phase jitters per epoch
    kp_jit ~ Uniform(0,180)# 12.8
    kp_Cy ~ Uniform(-1,1)
end vis_like b;

model = Octofitter.LogDensityModel(sys, verbosity=4)
```

It is recommended to use Pigeons parallel tempered sampling, and to use the SliceSampler explorer. This non-default option avoids calculating the gradient of the model, which is expensive in this case.
```@julia
Octofitter.default_initializer!(model,ntries=0)

using Pigeons
chain, pt = octofit_pigeons(model, n_chains=8, n_chains_variational=0, n_rounds=9, explorer=SliceSampler())
```


```@julia
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

```@julia
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
    LL[i] = model.ℓπcallback(model.link((jit,Cy,sep,pa,K)),sampled=false)
end
```

```@julia
fig = Figure()
ax = Axis(
    fig[1,1],
    xreversed=true,
    autolimitaspect=1,
    backgroundcolor="#222",
    title="spec corr = 0.02"
)

N_dat = length(vis_like.table.epoch)*length(vis_like.table.eff_wave[1])*3 # 3 kern phases
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



This single-epoch model can then be extended by replacing the `FixedPosition` parameterization with an orbit type like `KepOrbit`:
```julia
@planet b Visual{Octofitter.FixedPosition} begin
    a ~ Uniform(0, 0.1)
    e ~ Uniform(0.0, 0.99)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,60676)

    # We will sample using a flat contrast ratio versus spectrum
    K ~ Uniform(0, 1) # contrast ratio variable
    K_band_spectrum = fill(b.K, $(length(vis_like.table.eff_wave[1])))

end
```

