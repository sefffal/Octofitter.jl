# Citing

To cite octofitterpy or Octiftter.jl, please use the following BibTeX entry:

```bibtex
@article{Thompson_2023,
    doi = {10.3847/1538-3881/acf5cc},
    url = {https://dx.doi.org/10.3847/1538-3881/acf5cc},
    year = {2023},
    month = {sep},
    publisher = {The American Astronomical Society},
    volume = {166},
    number = {4},
    pages = {164},
    author = {William Thompson and Jensen Lawrence and Dori Blakely and Christian Marois and Jason Wang and Mosé Giordano and Timothy Brandt and Doug Johnstone and Jean-Baptiste Ruffio and S. Mark Ammons and Katie A. Crotts and Clarissa R. Do Ó and Eileen C. Gonzales and Malena Rice},
    title = {Octofitter: Fast, Flexible, and Accurate Orbit Modeling to Detect Exoplanets},
    journal = {The Astronomical Journal},
}
```


* If you use Gaia parallaxes in your work, please cite Gaia DR3 [Gaia Collaboration et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A&A...674A...1G)
* Please cite the HMC sampler backend if you use `octofit`: [Xu et al 2020](http://proceedings.mlr.press/v118/xu20a.html)
* Please cite the [Pigeons paper](https://arxiv.org/abs/2308.09769) if you use `octofit_pigeons`.
* If you use Hipparcos-GAIA proper motion anomaly, please cite [Brandt 2021](https://ui.adsabs.harvard.edu/abs/2021ApJS..254...42B)
* If you use example data in one of the tutorials, please cite the sources listed
* If you use one of the included functions for automatically retreiving data from a public dataset, eg HARPS RVBank, please cite the source as appropriate (it will be displyed in the terminal)
* If you adopt the O'Neil et al. 2019 observable based priors, please cite [O'Neil et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....158....4O).
* If you use RV phase folded plot, please consider citing Makie.jl [Danisch & Krumbiegel, (2021).](https://doi.org/10.21105/joss.03349)
* If you use TemporalGPs.jl to accelerate Gaussian processes modelling of stellar activity, please consider citing [Tebbutt et al 2021](https://proceedings.mlr.press/v161/tebbutt21a.html)
* If you use the pairplot/cornerplot functionality, please cite:
```bibtex
@misc{Thompson2023,
  author = {William Thompson},
  title = {{PairPlots.jl} Beautiful and flexible visualizations of high dimensional data},
  year = {2023},
  howpublished = {\url{https://sefffal.github.io/PairPlots.jl/dev}},
}
```