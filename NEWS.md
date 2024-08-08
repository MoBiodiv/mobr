## version 3.0.0
---

Major update
### new features
- `calc_comm_div` now replaces `get_mob_stats`. The ancillary plotting function
`plot.mob_stats` is now replaced with `plot_comm_div`. We hope to eventually
add back in bootstrapping confidence intervals for the statistics but that is
not currently supported. If you would like to know more about this design
decision see the discussion here: https://github.com/MoBiodiv/mobr/issues/255. 
- `plot_rarefaction` is a bit more versatile with many new arguments 
that provide options for smoothing or averaging rarefaction curves
when making comparisons within or between groups respectively. 
- The R package `beta_C` by Thore Engel has now been absorbed into
the `mobr` package. Use function `calc_beta_div` or `calc_comm_div`
with `beta` in the `scales` argument to compute coverage based beta
diversity. Note the index in this case is `S_C` (i.e., richness S for a
given level of coverage C)
- a vignette was added to demonstrate how to make computations of beta
diversity using `mobr` called `beta_div_demo`

## version 2.0.2
---
### change for CRAN
- remove sort on a `data.frame` from the function `plot_abu` which does not
change how this function behaves but adheres to CRAN standards.

## version 2.0.1
---

### new features
- `ref_level` can be specified in `get_mob_stats`, `plot_rarefaction` and `plot_sad`
which allows for the reference (i.e., control) group to be specified which all
others are compared to. This is similar to the implementation in the function
`get_delta_stats`

### minor bug fixes
- typo in test_rarefaction.R
- typo in get_mob_stats documentation


## version 2.0.0
---

### new features
- Gradient multi-scale mob method added to the function `get_delta_stats` as
described in McGlinn et al. ([2020](https://doi.org/10.1101/851717)
- Added a new k-nearest neighbor algorithm for the spatial sample based rarefaction curve (sSBR)
- ggplot integrated into plotting method for `get_delta_stats` -> `plot.mob_out`
- Added automated tests
- Improved documentation throughout the package

### bug fixes
- Function `sphere_dist` did not appear to be calculating great circle distances
from latitude and longitude coordinates correctly using the Haversine formula.
This was corrected with a new function this is correct according to comparisons
with `fields::rdist.earth` and `geosphere::distHaversine`.

### omitted feature
- Removed the function `overlap_effects` which provided an overlapped image of the three components SAD, N, and agg across scale. This may be added back in the future

### code refactoring
- Refactored code surrounding the `get_delta_stats` function which computes the 
differences in richness across all possible scales. 
- Refactoring was primarily accomplished by integrating in the tidyverse tools
which should make the code easier to understand and follow
- Computational time has increased as the tidyverse tools are not necessarily 
the most computationally efficient solutions
 



## version 1.0.0
---

- provided core functionality of mob tools described in McGlinn et al. 
([2019](https://doi.org/10.1111/2041-210X.13102)) and Chase et al. 
([2018](https://doi.org/10.1111/ele.13151)) 
