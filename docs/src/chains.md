# [Chains](@id chains)

This page describes the format of the Monte Carlo chains created by DirectDetections.jl


The output of the samplers in DirectDetections is an MCMCChains.Chains object

A column will be present for each variable in your model, both defined in the Priors blocks or as Derived variables. 

Variables defined for the System as a whole can be accessed directly. For example:
```julia
chain["μ"]
```
This will return an array of μ values from the posterior. The format is a matrix of N-samples by N-walkers.

Variables defined for an individual Planet are grouped according to the standard MCMCChains format. For example:
```julia
chain["b[a]"]
```
This returns an array of semi-major axis values (`a`) for the planet `b` sampled from the posterior.

## Diagnostics
Printing the chains will display a number of useful summaries for each quantity, like the mean, 0.25, 0.5, and 0.75 quantiles, and convergence metrics. See MCMCChains documentation for more details.