## Overview

The `muppet` package executes Measurement and Uncertainty Preserving ParamETric (MUPPET) modeling. MUPPET is a modular approach to model fitting, using Bayesian approaches to modeling and estimation. Models and analyses are conceived of in terms of fragments. Each may be fit to data or analyzed, possibly conditional on antecedent fragments. In that case, the results from the antecedent fragments are carried forward to those that depend on them. 

## Installation

You can install the current public version of `muppet` from GitHub. It depends on other packages for certain fuctionality. `muppet` uses functions that are part of the `mcmcplots` package, developed by S. McKay Curtis, to plot results from the chains. However, `mcmcplots` is no longer availabe on CRAN. A verion of `mcmcplots` may be installed from GitHub, and then `muppet` can also be installed from GitHub. This may be done as follows:

```r
# install.packages("remotes")  # if you don't have remotes installed

# Install mcmcplots from maintained fork
remotes::install_github("Roy-Levy/mcmcplots")

# Then install muppet
remotes::install_github("Roy-Levy/muppet_pub")
```

## Examples

The [Examples and Tips](https://roy-levy.github.io/muppet_pub/articles/index.html) page has some examples and explanations of some of the basic functionality. 

## Contact

Have questions or suggestions for features? Found a bug? 
Please contact me at Roy.Levy@asu.edu

