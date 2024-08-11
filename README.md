# mobr
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/mobr)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/mobr)](https://cran.r-project.org/package=mobr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4014111.svg)](https://doi.org/10.5281/zenodo.4014111)
============

# Measurement of Biodiversity in R 

This repository hosts an R package that is actively being developed for 
estimating biodiversity and the components of its change. The key innovations of
this R package over other R packages that also carry out rarefaction (e.g.,
`vegan`, `iNext`) is that `mobr` is focused on 1) making empirical comparisons between 
treatments or gradients, and 2) our framework emphasizes how changes in 
biodiversity are linked to changes in community structure: the SAD, total
abundance, and spatial aggregation. 

The concepts and methods behind this R package are described in three publications.

McGlinn, D.J., S.A. Blowes, M. Dornelas, T. Engel, I.S. Martins, H. Shimadzu,  N.J. Gotelli,  A. Magurran,  B.J. McGill,  and J.M. Chase. accepted. Disentangling non-random structure from random placement when estimating β-diversity through space or time. Ecosphere. https://doi.org/10.1101/2023.09.19.558467 


McGlinn, D.J. X. Xiao, F. May, N.J Gotelli, T. Engel, S.A Blowes, T.M. Knight, O. Purschke, J.M Chase, and B.J. McGill. 2019. MoB (Measurement of Biodiversity): a method to separate the scale-dependent effects of species abundance distribution, density, and aggregation on diversity change. Methods in Ecology and Evolution. 10:258–269. https://doi.org/10.1111/2041-210X.13102

McGlinn, D.J. T. Engel, S.A. Blowes, N.J. Gotelli, T.M. Knight, B.J. McGill, N. Sanders, and J.M. Chase. 2020. A multiscale framework for disentangling the roles of evenness, density, and aggregation on diversity gradients. Ecology. https://doi.org/10.1002/ecy.3233

Chase, J.M., B. McGill, D.J. McGlinn, F. May, S.A. Blowes, X. Xiao, T. Knight. 2018. Embracing scale-dependence to achieve a deeper understanding of biodiversity and its change across communities. Ecology Letters. 21: 1737–1751. https://doi.org/10.1111/ele.13151 

Please cite `mobr`. Run the following to get the appropriate citation for the version you're using:

```r
citation(package = "mobr")
```

## Installation

```r
install.packages('mobr')
```

Or, install the Github version

```r
install.packages('remotes')
```

Now that `remotes` is installed you can install `mobr` using the following R code:

```r
remotes::install_github('MoBiodiv/mobr')
```

## Examples

The package [vignette](https://github.com/MoBiodiv/mobr/blob/master/vignettes/mobr_intro.pdf)
provides a useful walk-through the package tools, but below is some example code
that uses the two key analyses and related graphics. 

```r
library(mobr)
library(dplyr)

data(tank_comm)
data(tank_plot_attr)
indices <- c('N', 'S', 'S_n', 'S_C', 'S_PIE')
tank_div <- tibble(tank_comm) %>% 
  group_by(group = tank_plot_attr$group) %>% 
  group_modify(~ calc_comm_div(.x, index = indices, effort = 5,
                               extrapolate = TRUE))
plot(tank_div)
tank_mob_in <- make_mob_in(tank_comm, tank_plot_attr, coord_names = c('x', 'y'))
tank_deltaS <- get_delta_stats(tank_mob_in, 'group', ref_level='low',
                             type='discrete', log_scale=TRUE, n_perm = 5)
plot(tank_deltaS, 'b1')
```

## Meta

* Please [report any issues or bugs](https://github.com/mobiodiv/mobr).
* License: MIT
* Get citation information for `mobr` in R doing `citation(package = 'mobr')`
* Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its terms.

## Thanks

* Gregor Seyer for providing a constructive review of our CRAN submission
* Kurt Hornik for helping us keep up with CRAN changes. 
