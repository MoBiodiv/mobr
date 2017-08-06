# Measurement of Biodiversity in R (mobr)

This reposititory hosts an R package that is being developed for 
estimating biodiversity and the components of its change.

# How to install mobr

The easiest option is to install the package directly from GitHub using the package `devtools`. If you do not already have `devtools` installed then need to install it.

```r
install.packages('devtools')
library(devtools)
```

The package also requires the `dplyr` package which has many dependencies. If you do 
not already have the `dplyr` package installed we suggest you install it and 
all of its dependencies using :

```r
install.packages(c('bindrcpp','glue','pkgconfig','tibble','plyr','dplyr'))
```

Then check that `dplyr` can be loaded with `library(dplyr)`.
Now you should be ready to go with the `mobr` install

```r
install_github('MoBiodiv/mobr')
```

# Examples
```r
library(mobr)
data(inv_comm)
data(inv_plot_attr)
inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_group='uninvaded',
                              type='discrete', log_scale=TRUE, n_perm=100)
plot(inv_mob_out, 'invaded', 'uninvaded')
```


## How to contribute to this project
1) Fork the repo to your local GitHub account

2) Clone your forked version of the repo to your machine

`git clone git@github.com:your_user_name/mobr.git`

3) Link your local repo back to the master on MoBiodiv

`git remote add upstream git@github.com:MoBiodiv/mobr.git`

4) Create a branch for your changes

`git branch new_function`

5) Checkout your branch

`git checkout new_function`

6) Make your commits on that branch and when you are done push it to your
forked copy of the repo

`git push origin new_function`

7) Submit a pull request on the GitHub website by going to your forked copy
of the repo and clicking on the pull request button 

8) After your changes are merged with master you'll want to merge that
update to master with your copies as well. 

```
git pull upstream master
git push origin master
# delete your branch as its no longer needed
git branch -d new_function
```

Before your start work on the project in the future you'll want to repeat
step 8 so that your version of the repo does not become out-of-sync
with the main repository. 
