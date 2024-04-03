# Cross-Validation


!!! note
    This tutorial is currently a stub and will be expanded in future. In the meantime, if you have questions please open an issue on GitHub.


## Calculating Pointwise Likelihoods

After you have defined a model and sampled from its posterior (eg. via `octofit`), you can see how each datapoint is influencing the posterior via the following functions.


```julia
# already have defined `model` and `chain` ...

likelihood_mat = Octofitter.pointwise_like(model, chain)
```

`likelihood_mat` is now a N_sample x N_data matrix. 

!!! note
    The columns are ordered the same as how the data are defined in the model.
    
!!! warning  
    You may see a few additional entries you didn't expect. Each `UniformCircular` and `ObsPriorAstromONeil2019` adds an additional likelihood object under the hood.


## Pareto-Smoothed Importance Sampling
After you have generated the `likelihood_mat` you can use the Julia package [ParetoSmooth.jl](https://turinglang.org/ParetoSmooth.jl/dev/) to efficiently calculate a leave-one-out cross-validataion score. This technique takes a single posterior chain and, using the pointwise likelihoods, generates `N_datapoints` new chains where each chain is adjusted as if that datapoint was held out from the model. 

In broad terms, one might say that this test verifies that no individual datapoints are overly skewing the results.

```julia
using ParetoSmooth
result = psis_loo(
    collect(likelihood_mat'),
    chain_index=ones(Int,size(chain,1))
)
```


Plot like so:
```julia
using CairoMakie

fig = Figure()

ax = Axis(
    fig[1,1],
    xlabel="data #",
    ylabel="Pareto K"
)
scatter!(ax, result.pointwise(:pareto_k))


ax = Axis(
    fig[2,1],
    xlabel="data #",
    ylabel="MCSE"
)
scatter!(ax, result.pointwise(:mcse))


ax = Axis(
    fig[3,1],
    xlabel="data #",
    ylabel="P_EFF"
)
scatter!(ax, result.pointwise(:p_eff))


fig
```