# Generative Models

Given an existing model with observations and variables, you can use this package to generate a *new* model with new observations. These observations are generated from a set of parameters specified by the user or draw automatically from the priors.

This functionality is useful for a range of different Bayesian workflows, including *simulation based calibration*.

To generate a new model by drawing from priors, simply call `generate(system)`.
