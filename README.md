# DirectDetections.jl


This in development package uses hierarchical Bayesian modelling to detect exoplanets and other sub-stellar companions by jointly modelling their orbit and atmosphere. Orbit modelling is handled using DirectOrbits.jl and atmospheres via grids of Sonora models.

Simply specify your priors on physical and orbital parameters, provide any direct images of the system from any bands in the Sonora grids, as well as any RV or astrometry measurements. This package will then generate a posterior distribution which can be used to assess a detection and/or constrain these parameters.
