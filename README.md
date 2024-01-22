# Octofitter.jl


[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sefffal.github.io/Octofitter.jl/dev)

Octofitter is a Julia package for performing Bayesian inference 
against a wide variety of exoplanet / binary star data.
You can also use Octofitter from Python (see docs).

Tutorials and documentation are available [here](https://sefffal.github.io/Octofitter.jl/).

![](docs/src/assets/gallery.png)



### Read the paper
In addition to these documentation and tutorial pages, you can read the paper published in the [Astronomical Journal](https://dx.doi.org/10.3847/1538-3881/acf5cc) (open-access).

## Attribution
* If you use Octofitter in your work, please cite [Thompson et al](https://dx.doi.org/10.3847/1538-3881/acf5cc)
* If you use Gaia parallaxes in your work, please cite Gaia DR3 [Gaia Collaboration et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A&A...674A...1G)
* If you use Hipparcos-GAIA proper motion anomaly, please cite [Brandt 2021](https://ui.adsabs.harvard.edu/abs/2021ApJS..254...42B)
* If you use example data in one of the tutorials, please cite the sources listed [Brandt 2021](https://ui.adsabs.harvard.edu/abs/2021ApJS..254...42B)
* If you use one of the included functions for automatically retreiving data from a public dataset, eg HARPS RVBank, please cite the source as appropriate.
* If you adopt the O'Neil et al. 2019 observable based priors, please cite [O'Neil et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....158....4O).

* Please also consider citing the HMC sampler backend, [Xu et al 2020](http://proceedings.mlr.press/v118/xu20a.html)
* If you use RV phase folded plot, please consider citing Makie.jl [Danisch & Krumbiegel, (2021).](https://doi.org/10.21105/joss.03349)
* If you use TemporalGPs.jl to accelerate Gaussian processes modelling of stellar activity, please consider citing [Tebbutt et al 2021](https://proceedings.mlr.press/v161/tebbutt21a.html)
* If you use the pairplot functionality, please cite:
```
@misc{Thompson2023,
  author = {William Thompson},
  title = {{PairPlots.jl} Beautiful and flexible visualizations of high dimensional data},
  year = {2023},
  howpublished = {\url{https://sefffal.github.io/PairPlots.jl/dev}},
}
```


## Ready?


For instructions, see the documentation page:

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sefffal.github.io/Octofitter.jl/dev)
