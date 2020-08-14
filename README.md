mobr [![Travis build status](https://travis-ci.org/MoBiodiv/mobr.svg?branch=master)](https://travis-ci.org/MoBiodiv/mobr) [![DOI](https://zenodo.org/badge/36759899.svg)](https://zenodo.org/badge/latestdoi/36759899)
============

# Measurement of Biodiversity in R 

This repository hosts an R package that is actively being developed for 
estimating biodiversity and the components of its change. The key innovations of
this R package over other R packages that also carry out rarefaction (e.g.,
`vegan`, `iNext`) is that `mobr` is focused on 1) making empirical comparisons between 
treatments or gradients, and 2) our framework emphasizes how changes in 
biodiversity are linked to changes in community structure: the SAD, total
abundance, and spatial aggregation. 

The concepts and methods behind this R package are described in two publications.

McGlinn, D.J. X. Xiao, F. May, S. Blowes, J. Chase, N. Gotelli, T. Knight, B. McGill, and O. Purschke. 2019. MoB (Measurement of Biodiversity): a method to separate the scale-dependent effects of species abundance distribution, density, and aggregation on diversity change. Methods in Ecology and Evolution. 10:258–269. https://doi.org/10.1111/2041-210X.13102


Chase, J.M., B. McGill, D.J. McGlinn, F. May, S.A. Blowes, X. Xiao, T. Knight. 2018. Embracing scale-dependence to achieve a deeper understanding of biodiversity and its change across communities. Ecology Letters. 21: 1737–1751. https://doi.org/10.1111/ele.13151 


## Installation

The easiest option is to install the package directly from GitHub using the package `devtools`. If you do not already have `devtools` installed then need to install it.

```r
install.packages('devtools')
library(devtools)
```

Now that `devtools` is installed you can install `mobr using the following R code:

```r
install_github('MoBiodiv/mobr')
```

Note: the installation may take some time due to numerous dependencies. We are 
working on reducing the number of dependencies. 

## Examples

The package [vignette](./vignettes/mobr_intro.pdf) provides a useful walk-through
the package tools, but below is some example code that uses the two key analyses
and related graphics. 

```r
library(mobr)
data(inv_comm)
data(inv_plot_attr)
inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
inv_stats = get_mob_stats(inv_mob_in, 'group')
plot(inv_stats)
inv_deltaS = get_delta_stats(inv_mob_in, 'group', ref_level='uninvaded',
                             type='discrete', log_scale=TRUE, n_perm = 5)
plot(inv_deltaS, 'b1')
```

## Meta

* Please [report any issues or bugs](https://github.com/mobiodiv/mobr).
* License: MIT
* Get citation information for `mobr` in R doing `citation(package = 'mobr')`
* Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its terms.
