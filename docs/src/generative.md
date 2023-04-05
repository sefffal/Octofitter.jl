# Generative Models

Given an existing model with observations and variables, you can use this package to generate a *new* model with new observations. These observations are generated from a set of parameters specified by the user or draw automatically from the priors.

This functionality is useful for a range of different Bayesian workflows, including *simulation based calibration*.

To generate a new model by drawing from priors, simply call `Octofitter.generate_from_params(system)`.


To run a single trial of simulation-based-calibration, run the following code:
```julia
settings = Dict(
    :target_accept=>0.95,
    :num_chains=>1,
    :adaptation=>1000,
    :iterations=>1000,
    :tree_depth=>14,
    :thinning => 1,
    :verbosity=>4,
)
Octofitter.sbctrial(system, settings, "sbctrial-1");
```
You should run this many (at least hundreds) of times. Each trial will 
save a chain `.fits` file and `.toml` files with the parameters drawn
and the resulting ranks of the true value in the posterior.
