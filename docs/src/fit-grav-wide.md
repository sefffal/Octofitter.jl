# Fit GRAVITY-WIDE Data

## Background
Octofitter has support for directly fitting GRAVITY-WIDE closure phase data, in the OI-FITS format emitted by the pipeline.
The closure phases are actually mapped to a set of non-redundant kernel phases. All spectral channels are modelled separately per exposure.

!!! note
    GRAVITY modelling is supported in Octofitter via the extension package OctofitterInterferometry. To install it, run 
    `pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterInterferometry`

The only supported astrophysical sources at this time are zero or more point sources orbitting a primary body.

Interferometer data is almost always multi-modal, requiring the use of parallel tempering.
Multi-wavelength GRAVITY-WIDE data with multiple epochs is fairly expensive to model (can take on the order of 1ms per likelihood evaluation), so one after running some tests locally, one should consider using a compute cluster.
You will probably want on the order of 30 cores and 1-5 days, depending on the scale of the problem.

At present, these examples don't include any output plots etc because of the long runtime that would be required to generate them as part of our documentation build process.

## Process

To model orbits / brightness of a companion from, GRAVITY-WIDE data you will want to use the following Likelihood object:

```julia
vis_like = GRAVITYWideCPLikelihood(
    (;filename="./GRAVI.2025-01-01T00:11:11.111_dualscivis.fits", epoch=60676.00776748842, jitter=:kp_jit1, kp_Cy=:kp_Cy1, band=:K),
    # Add more epochs below if desired:
)
```

- `filename` is the path from your current working directory to the GRAVITY OI-FITS file.
- `epoch` is the average time of the exposure in MJD (not the start of the exposure!)
- `K` is a symbol giving the name of the contrast ratio parameter to use for this exposure
- `jitter` is a symbol giving the name of the kernel phase jitter parameter to use for this exposure
- `kp_Cy` is a symbol giving the name of the spectral correlation parameter to use for this exposure


