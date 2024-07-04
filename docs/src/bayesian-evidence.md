# Bayesian Evidence

The Bayesian evidence, $Z$, is a normalization constant that can be used to compare models. You can estimate $Z_1/Z_2$, that is the ratio of the evidence between two models, using Octofitter's integration with the `Pigeons` sampler.

To start, define two models that describe the same data. 

As an example, we will define two models: one that allows eccentric orbits, and one which does not (consideres only circular orbits). There are innumerable other