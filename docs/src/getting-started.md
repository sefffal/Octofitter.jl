# Getting Started

The first step to using DirectDetections.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. Currently, this is 1.6.2. Click the `[help]` links next to your operating system if you require more detailed instructions.

## Installing DirectDetections
Normally, Julia packages are installed from the General registry. Since DirectDetections isn't quite ready for prime time, it requires one extra step to add an additional registry.

1. Start julia in a terminal by running `julia`
2. Type `]` to enter package-mode (see Julia documentation for more details)
3.  Type `up` to setup the General registry if this is your first time using Julia.
4. Type `registry add https://github.com/sefffal/DirectRegistry`
5. Type `add DirectDetections Distributions`

You will need the Distributions package added above so that you can specify priors for different parameters in your models.

If you would like to visualize your results, you can also install the Plots package:
4. Type `add Plots`

For loading images to sample, add the DirectImages package 
4. Type `add DirectImages`
Note: it's possible to use this package without DirectImages. This just simplifies the process of loading FITS files, and creating a centered OffsetArray with the star at index (0,0).

Finally, you may wish to `add PairPlots` to examine the correlation between variables in more detail.

This will take a little while to download all the required packages and precompile for your system.


## Fitting your first model
Start with the [Fit Astrometry](@ref fit-astrometry) tutorial which shows how to model of one planet with some astrometry points.
