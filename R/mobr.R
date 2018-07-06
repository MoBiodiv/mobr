#' Create the 'mob_in' object.
#' 
#' The 'mob_in' object will be passed on for analyses of biodiversity across 
#' scales.
#' 
#' @param comm community matrix with plots as rows and species columns.
#' @param plot_attr matrix which includes the environmental attributes and
#'   spatial coordinates of the plots. Environmental attributes are mandatory,
#'   while spatial coordinates are not. If spatial coordinates are provided, the
#'   column(s) has to have names "x" and/or "y".
#' @param coord_names character vector with the names of the columns of
#'   \code{plot_attr} that specify the coordinates of the samples. Defaults to
#'   'x' and 'y'. The order the names are provided matters when working with
#'   latitude-longitude coordinates (i.e., argument \code{latlong = TRUE}, and
#'   it is expected that the column specifying the x-coordinate or the longitude
#'   is provided first, y-coordinate or latitude provided second.
#' @param binary boolean, defaults to FALSE. Whether the plot by species matrix
#'   "comm" is in abundances or presence/absence.
#' @param latlong boolean, defaults to FALSE. Whether the coordinates are
#'   latitude-longitudes or not. If TRUE, distance calculations by downstream
#'   functions are based upon great circle distances rather than Euclidean
#'   distances. Note latitude-longitudes should be in decimal degree.
#'   
#' @return a "mob_in" object with four attributes. "comm" is the plot by species
#'   matrix. "env" is the environmental attribute matrix, without the spatial
#'   coordinates. "spat" contains the spatial coordinates (1-D or 2-D). "tests"
#'   specifies whether each of the three tests in the biodiversity analyses is
#'   allowed by data.
#'   
#' @author Dan McGlinn and Xiao Xiao
#' @export
#' @examples
#'  data(inv_comm)
#'  data(inv_plot_attr)
#'  inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
make_mob_in = function(comm, plot_attr, coord_names=c('x', 'y'), binary=FALSE,
                       latlong=FALSE) {
    # possibly make group_var and ref_group mandatory arguments
    out = list(tests = list(N=TRUE, SAD=TRUE, agg=TRUE))
    # carry out some basic checks
    if (nrow(comm) < 5) {
        warning("Number of plots in community is less than five therefore only individual rarefaction will be computed")
        out$tests$N = FALSE
        out$tests$agg = FALSE
    }
    if (nrow(comm) != nrow(plot_attr))
        stop("Number of plots in community does not equal number of plots in plot attribute table")
    spat_cols = sapply(coord_names, function(x) which(x == names(plot_attr)))
    if (length(spat_cols) == 1 & latlong == TRUE)
        stop("Both latitude and longitude have to be specified")
    if (any(row.names(comm) != row.names(plot_attr)))
        warning("Row names of community and plot attributes tables do not match")
    if (binary)  {
        warning("Only spatially-explict sampled based forms of rarefaction can be computed on binary data")
        out$tests$SAD = FALSE
        out$tests$N = FALSE
    } else {
        if (max(comm) == 1)
            warning("Maximum abundance is 1 which suggests data is binary, change the binary argument to TRUE")
    }
    if (any(colSums(comm) == 0)) {
        warning("Some species have zero occurrences and will be dropped from the community table")
        comm = comm[, colSums(comm) != 0]
    }
    out$comm = data.frame(comm)
    if (length(spat_cols) > 0) {
        out$env = data.frame(plot_attr[ , -spat_cols])
        colnames(out$env) = colnames(plot_attr)[-spat_cols]
        out$spat = data.frame(plot_attr[ , spat_cols])
    }
    else {
        out$tests$agg = FALSE
        out$env = data.frame(plot_attr)
        out$spat = NULL
    }
    out$latlong = latlong
    class(out) = 'mob_in'
    return(out)
}

#' Subset the rows of the mob data input object
#'
#' This function subsets the rows of comm, env, and spat attributes of the
#' mob_in object
#'
#' @param mob_in an object of class mob_in created by \code{\link{make_mob_in}}
#' @param type specifies the type of object the argument \code{subset}
#'   specifies, may be: \code{string}, \code{integer}, or \code{logical},
#'   defaults to \code{string}
#' @param drop_levels boolean if TRUE unused levels are removed from factors in
#'   mob_in$env
#' @inheritParams base::subset
#' @export
#' @examples 
#'  data(inv_comm)
#'  data(inv_plot_attr)
#'  inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#'  subset(inv_mob_in, group == 'invaded')
#'  subset(inv_mob_in, 1:4, type='integer')
#'  subset(inv_mob_in, 1:4, type='integer', drop_levels=TRUE)
#'  sub_log = c(TRUE, FALSE, TRUE, rep(FALSE, nrow(inv_mob_in$comm) - 3))
#'  subset(inv_mob_in, sub_log, type='logical')
subset.mob_in = function(mob_in, subset, type='string', drop_levels=FALSE) {
    if (missing(subset))
        r <- rep_len(TRUE, nrow(mob_in$comm))
    if (type == 'integer')
        r <- 1:nrow(mob_in$comm) %in% subset
    if (type == 'logical')
        r <- subset
    if (type == 'string'){
        e <- substitute(subset)
        r <- eval(e, mob_in$env)
        if (!is.logical(r)) 
            stop("'subset' must be logical when type = 'string'")
   } 
   mob_in$comm = base::subset(mob_in$comm, r)
   mob_in$env = base::subset(mob_in$env, r)
   if (drop_levels)
       mob_in$env = droplevels(mob_in$env)
   if (!is.null(mob_in$spat))
       mob_in$spat = mob_in$spat[r, ]
   return(mob_in)
}
  
#' Print a shortened version of the mob_in object
#' @keywords internal
#' @export
print.mob_in = function(x, nrows=6, nsp=5) {
    if (nrow(x$comm) > nrows) 
        cat(paste('Only the first', nrows, 'rows of any matrices are printed\n'))
    else 
        nrows = nrow(x$comm)
    cat('\n$tests\n')
    print(x$tests)
    if (ncol(x$comm) > nsp)
        cat(paste('\n$comm (Only first', nsp, 'species columns are printed)\n'))
    else 
        nsp = ncol(x$comm)
    print(x$comm[1:nrows, 1:nsp])
    cat('\n$env\n')
    print(head(x$env, nrows))
    cat('\n$spat\n')
    print(head(x$spat, nrows))
    cat('\n$latlong\n')
    print(x$latlong)
}

summary.mob_out = function(...) {
   #  print summary anova style table
}

#' Internal function for distance matrix assuming inputs are lat and long 
#'   on a sphere
#' 
#' @param coords a matrix with latitudes and longitudes. The longitudes should
#' be provided in the first column (they are the x-coordinate) and the latitudes
#' should be provided in the second column (they are the y-coordinate). 
#' @description  Distance matrix between points on the unit (r = 1) sphere.
#' @return a numeric value
#' @author Xiao Xiao
#' @keywords internal
sphere_dist = function(coords){
    long = coords[ , 1]
    lat = coords[ , 2]
    # Convert degrees to radians
    deg2rad = function(deg) return(deg * pi / 180)
    delta_long = as.matrix(dist(as.matrix(deg2rad(long))))
    delta_lat = as.matrix(dist(as.matrix(deg2rad(lat))))
    hav = sin(delta_lat / 2)^2 + cos(lat) %*% t(cos(lat)) * sin(delta_long / 2)^2
    dist = 2 * asin(sqrt(hav))
    return(dist)
}

#' Rarefied Species Richness
#' 
#' The expected number of species given a particular number of individuals or
#' samples under and assumption of random sampling with replacement.
#' 
#' 
#' @param x can either be a: 1) mob_in object, 2) community matrix-like
#'  object in which rows represent plots and columns represent species, or 3)
#'  a vector which contains the abundance of each species. 
#' @param method either 'indiv', 'samp', or 'spat' for individual, sample, or
#'   sample spatially explicit based rarefaction respectively. To compute the
#'   sample-based, non-spatial rarefaction curve specify 'indiv' method with the
#'   appropriate \code{dens_ratio} (see Details).
#' @param effort optional argument to specify what number of individuals or 
#'   number of samples depending on 'method' to compute rarefied richness as. If
#'   not specified all possible values from 1 to the maximum sampling effort are
#'   used
#' @param coords an optional matrix of geographic coordinates of the samples.  
#'   Only required when using the spatial rarefaction method and this information
#'   is not already supplied by \code{x}. The first column should specify 
#'   the x-coordinate (e.g., longitude) and the second coordinate should 
#'   specify the y-coordinate (e.g., latitude)
#' @param dens_ratio the ratio of individual density between a reference group
#'   and the community data (i.e., x) under consideration. This argument is
#'   used to rescale the rarefaction curve when estimating the effect of
#'   individual density on group differences in richness.
#' @param extrapolate boolean which specifies if richness should be extrapolated
#'   when effort is larger than the number of individuals using the chao1 method.
#'   Defaults to FALSE in which case it returns observed richness. Extrapolation
#'   is only implemented for individual-based rarefaction 
#'   (i.e., \code{method = 'indiv'})
#' @param return_NA boolean defaults to FALSE in which the function returns the
#'   observed S when \code{effort} is larger than the number of individuals or
#'   number of samples (depending on the method of rarefaction). If set to TRUE
#'   then NA is returned. Note that this argument is only relevant when
#'   \code{extrapolate = FALSE}.
#' @param quiet_mode boolean defaults to FALSE, if TRUE then warnings and other
#'   non-error messages are suppressed.
#' @inheritParams make_mob_in 
#'   
#' @details The analytical formulas of Cayuela et al. (2015) are used to compute
#'   the random sampling expectation for the individual and sampled based
#'   rarefaction methods. The spatially constrained rarefaction curve (Chiarucci
#'   et al. 2009) also known as the sample-based accumulation curve (Gotelli and
#'   Colwell 2001) is computed by sampling each plot in the order of their
#'   spatial proximity. If plots have the same distance from the focal plot then
#'   one is chosen randomly to be sampled first. Each plot in the dataset is
#'   treated as a starting point and then the mean of these n possible
#'   accumulation curves is computed.
#' 
#' For individual-based rarefaction if effort is greater than the number of
#' individuals and \code{extrapolate = TRUE} then the Chao1 method is used 
#' (Chao 1984, 1987). The code used to perform the extrapolation was ported
#' from \code{iNext::D0.hat} found at \url{https://github.com/JohnsonHsieh/iNEXT}. 
#' T. C. Hsieh, K. H. Ma and Anne Chao are the original authors of the
#' \code{iNEXT} package. 
#' 
#' If effort is greater than sample size and \code{extrapolate = FALSE} then the 
#' observed number of species is returned. 
#' 
#' @return A vector of rarefied species richness values
#' @author Dan McGlinn and Xiao Xiao
#' @references 
#' Cayuela, L., N.J. Gotelli, & R.K. Colwell (2015) Ecological and 
#'  biogeographic null hypotheses for comparing rarefaction curves. Ecological
#'  Monographs, 85, 437-454. Appendix A: 
#'  http://esapubs.org/archive/mono/M085/017/appendix-A.php
#'  
#' Chao, A. (1984) Nonparametric estimation of the number of classes in a
#'  population. Scandinavian Journal of Statistics, 11, 265-270.
#'  
#' Chao, A. (1987) Estimating the population size for capture-recapture data
#'  with unequal catchability. Biometrics, 43, 783-791.
#'  
#' Chiarucci, A., G. Bacaro, D. Rocchini, C. Ricotta, M. Palmer, & S. Scheiner 
#'  (2009) Spatially constrained rarefaction: incorporating the autocorrelated
#'  structure of biological communities into sample-based rarefaction. Community
#'  Ecology, 10, 209-214.
#'  
#' Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures
#' and pitfalls in the measurement and comparison of species richness. Ecology
#' Letters, 4, 379–391.
#'
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' sad = colSums(inv_comm)
#' inv_mob_in= make_mob_in(inv_comm, inv_plot_attr)
#' # rarefaction can be performed on different data inputs
#' # all three give same answer
#' # 1) the raw community site-by-species matrix
#' rarefaction(inv_comm, method='indiv', effort=1:10)
#' # 2) the SAD of the community
#' rarefaction(inv_comm, method='indiv', effort=1:10)
#' # 3) a mob_in class object
#' # rescaling of individual based rarefaction 
#' # when the density ratio is 1 the richness values are 
#' # identical to unscale rarefaction
#' rarefaction(inv_comm, method='indiv', effort=1:10, dens_ratio=1)
#' # however the curve is either shrunk when density is higher than 
#' # the reference value (i.e., dens_ratio < 1)
#' rarefaction(inv_comm, method='indiv', effort=1:10, dens_ratio=0.5)
#' # the curve is stretched when density is lower than the 
#' # reference value (i.e., dens_ratio > 1)
#' rarefaction(inv_comm, method='indiv', effort=1:10, dens_ratio=1.5)
#' # sample based rarefaction under random sampling
#' rarefaction(inv_comm, method='samp')
#' # sampled based rarefaction under spatially explicit nearest neighbor sampling
#' rarefaction(inv_comm, method='spat', coords=inv_plot_attr[ , c('x','y')],
#'             latlong=FALSE)
#' # the syntax is simplier if suppling a mob_in object
#' rarefaction(inv_mob_in, method='spat')
rarefaction = function(x, method, effort=NULL, coords=NULL, latlong=NULL, 
                       dens_ratio=1, extrapolate=FALSE, return_NA = FALSE, 
                       quiet_mode=FALSE) {
    if (!any(method %in% c('indiv', 'samp', 'spat')))
        stop('method must be "indiv", "samp", or "spat" for random individual, random sample, and spatial sample-based rarefaction, respectively')
    if (class(x) == 'mob_in') {
        x_mob_in = x
        x = x_mob_in$comm
        if (is.null(latlong))
            latlong = x_mob_in$latlong
        else if(latlong != x_mob_in$latlong)
            stop(paste('The "latlong" argument is set to', latlong, 
                       'but the value of x$latlong is', x_mob_in$latlong))
        if (is.null(coords))
            coords = x_mob_in$spat
    }
    if (method == 'samp' | method == 'spat') {
        if (is.null(dim(x)))
            stop('For random or spatially explicit sample based rarefaction "x" must be a site x species matrix as the input')
        else {
            x = (x > 0) * 1             
            # all sites are counted as samples even empty ones
            n = nrow(x) 
            if (method == 'samp')
                x = colSums(x)
        }
    }
    if (method == 'indiv') {
        if (!is.null(dim(x)))
            x = colSums(x)
        n = sum(x)
    }
    if (is.null(effort))
        if (n == 0)
            effort = 0
        else
            effort = 1:n
    if (any(effort > n)) {
        if (extrapolate & return_NA)
            stop('It does not make sense to set "extrapolate" and "return_NA" to both be TRUE, see documentation')
        if (!quiet_mode) {
            warning_mess = paste('"effort" larger than total number of',
                                 ifelse(method == 'indiv', 'individuals', 'samples'),
                                 'returning')
            if (extrapolate)
                warning(paste(warning_mess, 'extrapolated S using Chao1'))
            else if (return_NA)
                warning(paste(warning_mess, 'NA'))
            else
                warning(paste(warning_mess, 'S'))
        }
    } else if (extrapolate)
        if (!quiet_mode) 
            message('Richness was not extrapolated because effort less than or equal to the number of samples')
    if (method == 'spat') {
        explicit_loop = matrix(0, n, n)
        if (is.null(latlong))
            stop('For spatial rarefaction the argument "latlong" must be set TRUE or FALSE')
        if (latlong){
            # Compute distance on sphere if xy are longitudes and latitudes
            # Assume x is longitude and y is latitude
            pair_dist = sphere_dist(coords)
        } else {
            pair_dist = as.matrix(dist(coords))
        }
        for (i in 1:n) {
            dist_to_site = pair_dist[i, ]
            # Shuffle plots, so that tied grouping is not biased by original order.
            new_order = sample(1:n)  
            dist_new = dist_to_site[new_order]
            new_order = new_order[order(dist_new)]
            # Move focal site to the front
            new_order = c(i, new_order[new_order != i])
            comm_ordered = x[new_order, ]
            # 1 for absence, 0 for presence
            comm_bool = as.data.frame((comm_ordered == 0) * 1) 
            rich = cumprod(comm_bool)
            explicit_loop[ , i] = as.numeric(ncol(x) - rowSums(rich))
        }
        out = apply(explicit_loop, 1, mean)[effort]
    } 
    else {
        # drop species with no observations  
        x = x[x > 0] 
        S = length(x)
        if (dens_ratio == 1) {
            ldiv = lchoose(n, effort)
        } else {
            effort = effort[effort / dens_ratio <= n]
            ldiv = lgamma(n - effort / dens_ratio + 1) - lgamma(n + 1)
        }
        p = matrix(0, sum(effort <= n), S)
        out = rep(NA, length(effort))
        S_ext = NULL
        for (i in seq_along(effort)) {
            if (effort[i] <= n) {
                if (dens_ratio == 1) {
                    p[i, ] = ifelse(n - x < effort[i], 0, 
                                    exp(lchoose(n - x, effort[i]) - ldiv[i]))
                } else {
                    p[i, ] = ifelse(n - x < effort[i] / dens_ratio, 0, 
                                    exp(suppressWarnings(lgamma(n - x + 1)) -
                                        suppressWarnings(lgamma(n - x - effort[i] /
                                                                dens_ratio + 1)) +
                                        ldiv[i]))
                }
            } else if (extrapolate) {
                f1 = sum(x == 1)
                f2 = sum(x == 2)
                # estimation of unseen species via Chao1                
                f0_hat <- ifelse(f2 == 0, 
                                 (n - 1) / n * f1 * (f1 - 1) / 2, 
                                 (n - 1) / n * f1^2 / 2 / f2)
                A = n * f0_hat / (n * f0_hat + f1)
                S_ext = c(S_ext, ifelse(f1 == 0, S, 
                                        S + f0_hat * (1 - A ^ (effort[i] - n))))
            }
            else if (return_NA)
                S_ext = c(S_ext, NA)
            else 
                S_ext = c(S_ext, S)
        }
        out = rep(NA, length(effort))
        out[effort <= n] = rowSums(1 - p)
        out[effort > n] = S_ext
    }
    names(out) = effort
    return(out)
}

#' Compute permutation derived individual-based rarefaction curves
#' 
#' An internal function that can provide an independent derivation of 
#' the individual rarefaction curve for the purposes of testing the 
#' performance of the function \code{rarefaction}
#' 
#' @param abu a vector of species abundances
#' @param n_perm the number of permutations to average across, defaults to 100
#' @param n_indiv the number of individuals to evaluate the rarefaction curve
#' at. The default behavior is to evaluate it on a log2 interval from 1 to N 
#' @examples 
#' \donttest{ 
#' data(inv_comm)
#' sad = colSums(inv_comm)
#' ind_rare_perm(sad)
#' }
#' @keywords internal
ind_rare_perm = function(abu, n_perm=100, n_indiv=NULL) {
    if (!is.vector(abu)) {
        stop('abu must be a vector of abundances')
    } 
    calc_S = function(splist, n_indiv) {
        sapply(n_indiv, function(n) length(unique(splist[1:n])))
    }
    rand_splist = function(abu, S) {
        sample(unlist(mapply(rep, 1:S, abu)), replace=F)
    }
    S = length(abu)
    N = sum(abu)
    if (is.null(n_indiv))
        n_indiv = c(2^(seq(0, log2(N))), N)
    S_rand = replicate(n_perm, calc_S(rand_splist(abu, S), n_indiv))
    S_avg = apply(S_rand, 1, mean)
    S_qt = apply(S_rand, 1, quantile, c(0.025, 0.975))
    return(data.frame(n_indiv, S_avg, S_lo = S_qt[1, ], S_hi = S_qt[2, ]))
}

#' Compute average nearest neighbor distance
#' 
#' This function computes the average distance of the next
#' nearest sample for a given set of coordinates. This method
#' of sampling is used  by the function \code{rarefaction}
#' when building the spatial accumulation curves
#' 
#' @param coords a matrix with n-dimensional coordinates
#' @return a vector of average distances for each sequential number
#'   of accumulated nearest samples. 
#' @export
#' @examples 
#' # transect spatial arragnement
#' transect = 1:100
#' avg_nn_dist(transect)
#' grid = expand.grid(1:10, 1:10)
#' avg_nn_dist(grid)
#' par(mfrow=c(1,2)) 
#' plot(avg_nn_dist(transect), type='o', main='transect',
#'      xlab='# of samples', ylab='average distance')
#' # 2-D grid spatial arrangement
#' plot(avg_nn_dist(grid), type='o', main='grid',
#'      xlab='# of samples', ylab='average distance')
avg_nn_dist = function(coords) {
    pair_dist = as.matrix(dist(coords))
    sort_dist = apply(pair_dist, 1, sort)
    avg_dist = apply(sort_dist, 1, mean)
    return(avg_dist)
}


#' Auxiliary function for computing S and the effect on S of 
#' the three components of community structure: SAD, N, and aggregation
#' @param x can either be a: 1) mob_in object or 2) a vector which contains
#'  the abundance of each species (i.e., the SAD). All effects can be computed
#'  when x is a mob_in object but only the SAD effect can be computed when
#'  x is a vector of species abundances. 
#' @param tests what effects to compute defaults to 'SAD', 'N', and 'agg'
#' @param ind_dens the density of individuals to compare against for computing
#'  N effect
#' @keywords internal
get_delta_curves = function(x, tests=c('SAD', 'N', 'agg'),
                            inds=NULL, ind_dens=NULL) {
    if (is.null(inds) & any(c('SAD', 'N') %in% tests))
        stop('If SAD or N effect to be calculated inds must be specified')
    if (is.null(ind_dens) & 'N' %in% tests)
        stop('If N effect to be calculated ind_dens must be specified')
    if (any(c('N', 'agg') %in% tests) & class(x) != 'mob_in')
        stop('If N or agg effects to be computed x must be a mob_in object')
    out = list(NULL)
    if ('SAD' %in% tests) {
        S_SAD = rarefaction(x, 'indiv', inds)
        out$SAD = data.frame(test = 'SAD', sample = 'indiv',
                             effort = inds, S = S_SAD, effect = S_SAD)
    }
    if ('N' %in% tests) {
        comm_dens = sum(x$comm) / nrow(x$comm)
        dens_ratio = ind_dens / comm_dens
        S_N = rarefaction(x, 'indiv', inds, dens_ratio = dens_ratio)
        if (!('SAD' %in% tests))
            S_SAD = rarefaction(x, 'indiv', inds)
        effect = S_N - S_SAD
        out$N = data.frame(test = 'N', sample = 'indiv', 
                           effort = inds, S = S_N, effect)
    }
    if ('agg' %in% tests) {
        S_agg = rarefaction(x, 'spat')
        n_plots = nrow(x$comm)
        samp_effort = round((1:n_plots * sum(x$comm)) / n_plots)
        S_N = rarefaction(x, 'indiv', samp_effort)
        effect = S_agg - S_N
        out$agg = data.frame(test = 'agg', sample = 'plot', 
                             effort = as.numeric(names(S_agg)),
                             S = S_agg, effect)
    }
    return(flatten_dfr(tibble(out)))
}
        

#' @keywords internal
get_rand_sad = function(rad, N) {
  rand_samp = sample(1:length(rad), N, replace = T, prob = rad)
  rand_sad = table(factor(rand_samp, levels = 1:length(rad)))
  return(as.numeric(rand_sad))
}

#' Generate a null community matrix 
#' 
#' ## To DO: update this documentation
#' Two types of permutation can be carried out: 1) 'noagg': each individual of
#' each species is reassigned a plot randomly which removes any patterns due to
#' within and between plot spatial aggregation, but maintains species group
#' abundance and therefore observed group N, and 2) 'swapN': the total number of
#' individuals in a plot is shuffled and then that many individuals are drawn
#' randomly from the group specific species-abundance distribution for each plot
#' which provides a means of removing group differences in the total number of
#' individuals in a given sample. Note: The 'noagg' algorithm fixes the total number
#' of individuals for a given species within a group and the 'swapN' algorithm does
#' not. 
#' 
#' Replaces depreciated function `permute_comm`
#' 
#' @param comm community matrix with plots as rows and species columns.
#' @param tests 'SAD', 'N', or 'agg' for a null model that nullifies 
#' the given component
#' @param groups optional argument that is a vector of group ids which specify
#'   which group each site is associated with. If is NULL then all rows of the
#'   community matrix are assumed to be members of the same group
#'   
#' @return a site-by-species matrix
#' @export
#' @examples 
#' S = 3
#' N = 20
#' nplots = 4
#' comm = matrix(rpois(S*nplots, 1),
#'               ncol=S, nrow=nplots)
#' comm
#' groups = rep(1:2, each=2)
#' groups
#' get_null_comm(comm, 'noagg')
#' get_null_comm(comm, 'noagg', groups)
#' get_null_comm(comm, 'swapN')
#' get_null_comm(comm, 'swapN', groups)
get_null_comm = function(comm, tests, groups = NULL) {
    # the main component of all the null models is random sampling 
    # from a pooled or group-specific SAD
    if (!(is.matrix(comm) | is.data.frame(comm)))
        stop('comm must be a matrix or data.frame')
    if (is.null(groups))
        groups = rep(1, nrow(comm))   
  
    # compute N at each plot across groups
    N_plots = rowSums(comm)
    if (tests == "N") # shuffle these abundances in the N null model
        N_plots = sample(N_plots)
    if (tests == "SAD") {
        # compute relative abundance distribution of the species pool
        rad_pool = colSums(comm) / sum(comm)
        # randomly sample N individuals from the pool with replacement
        null_sads = map(N_plots, ~ get_rand_sad(rad_pool, .x))
        names(null_sads) = 1:length(null_sads)
    } else if (tests == "N" | tests == "agg") {
        # compute rad for each group
        rad_groups = data.frame(comm, groups) %>%
                     group_by(groups) %>%
                     summarize_all(sum) %>%
                     select(-one_of("groups")) %>%
                     t %>% as_tibble %>%
                     map(~ .x / sum(.x))
        # replicate these rads so that you have one rad for every
        # plot in the dataset
        rad_plots = rep(rad_groups, table(groups))
        names(rad_plots) = 1:length(rad_plots)
        # draw random sads 
        null_sads = map2(rad_plots, N_plots, get_rand_sad)
    }  
    # now convert sads to a new community matrix
    null_comm = null_sads %>% tibble %>% flatten_dfr %>% t  
    return(null_comm)
}


# Convert specified columns of a dataframe from factors to numeric
df_factor_to_numeric = function(dataframe, cols = NULL){
    if (is.null(cols)) cols = 1:ncol(dataframe)
    for (col in cols){
        if ('factor' %in% class(dataframe[, col]))
            dataframe[, col] = as.numeric(levels(dataframe[, col]))[dataframe[, col]]
    }
    return(dataframe)
}

# Auxiliary function for get_delta_stats()
# Overall checks for input values
get_delta_overall_checks = function(mob_in, type, group_var, env_var, 
                                    density_stat, tests){
    if (!(type %in% c('continuous', 'discrete')))
        stop('Type has to be discrete or continuous.')
    if (!(density_stat %in% c('mean', 'max', 'min')))
        stop('density_stat has to be set to min, max, or mean.')
    if (!(group_var %in% names(mob_in$env)))
        stop('group_var has to be one of the environmental variables in mob_in$env.')
    if (!(is.null(env_var)))
        if (!(env_var %in% names(mob_in$env)))
            stop('If env_var is defined, it has to be one of the environmental
                 variables in comm.')
    test_status = sapply(tests, function(x) 
        eval(parse(text = paste('mob_in$tests$', x, sep = ''))))
    approved_tests = tests[which(test_status == TRUE)]
    if (length(approved_tests) < length(tests)) {
        tests_string = paste(approved_tests, collapse=' and ')
        cat(paste('Based upon the attributes of the community object only the 
                  following tests will be performed:', tests_string))
    }
    return(approved_tests)
}

# Auxiliary function for get_delta_stats()
# Perform checks when type is "discrete"
get_delta_discrete_checks = function(ref_group, group_levels, group_data, env_var){
    if (is.null(ref_group)) {
        stop('For a discrete analysis you must specify a ref_group to compare groups to')
    } else if (!(as.character(ref_group) %in% group_levels)) {
        stop('Reference group does not exist.')
    }
    if (!is.null(env_var))
        warning('Environmental variable is not used in the discrete analysis.')
    if (!('factor' %in% class(group_data))) 
        warning('Grouping variable is not a factor. A group will be defined for each unique value.')
}

# Auxiliary function for get_delta_stats()
# Perform checks when type is "continuous"
get_delta_continuous_checks = function(corr, group_levels, env_raw){
    if (!(corr %in% c('spearman', 'pearson')))
        stop('corr has to be spearman or pearson.')
    if ('factor' %in% class(env_raw)) {
        env_vals = data.frame(groups = group_levels, 
                              values = as.numeric(env_raw)[match(group_levels, env_raw)])
        warning('Variable of interest is a factor but will be treated as a continous variable for the analysis with the above values')
        print(env_vals)
    } 
}

# Auxiliary function for get_delta_stats()
# Returns a vector of abundances where individual-based rarefaction 
# will be performed
get_inds = function(N_max, inds = NULL, log_scale = FALSE) {
    # across the groups what is the smallest total number of
    # individuals - this will be the largest N we can compute to
    if (is.null(inds)) {
        if (log_scale)
            ind_sample_size = unique(c(2^seq(0, floor(log2(N_max))), N_max))
        else 
            ind_sample_size = seq(N_max)
    }
    if (length(inds) == 1) { # if user specified an integer
        if (log_scale)  
            ind_sample_size = floor(exp(seq(inds) * log2(N_max) / inds))
        else 
            ind_sample_size = floor(seq(1, N_max, length.out = inds))
    }
    if (length(inds) > 1) { # if user specified a vector
        ind_sample_size = inds
        if (max(inds) > N_max)
            warning('Sample size is higher than abundance of at least one group!')
        ind_sample_size = unique(c(1, ind_sample_size)) # Force (1, 1) to be included
    }
    return(ind_sample_size)
}

#' Auxiliary function for get_delta_stats()
#' Returns the "assumed" density of individuals in 
#' a plot given whether min, max or mean is used
#' @keywords internal
get_ind_dens = function(comm, density_stat){
    if (density_stat == 'mean') {
        ind_dens = sum(comm) / nrow(comm)
    } else if (density_stat == 'max') {
        ind_dens = max(rowSums(comm))
    } else {
        ind_dens = min(rowSums(comm))
    }
    return(ind_dens)
}

#' Auxiliary function for effect_ functions
#' Compute an overall p-value for one factor in the discrete case
#' p-value is based on mean squared difference from zero summed across the scales
#' Method developed by Loosmore and Ford 2006 but algebraic simplifications 
#' used as developed by Baddeley et al. 2014 Ecological Archives M084-017-A1
#' @keywords internal
get_overall_p = function(effort, deltaS_emp, deltaS_null){
    delta_effort = c(effort[1], diff(effort))
    deltaS = rbind(deltaS_emp, deltaS_null)
    Hbarbar = apply(deltaS, 2, mean)                    # Baddeley Eq. A.10
    m = nrow(deltaS) - 1                                # number of permutations
    a = ((m + 1) / m)^2
    u = a * apply(deltaS, 1, function(x) 
                  sum((x - Hbarbar)^2 * delta_effort))  # Baddeley Eq. A.12-13
    overall_p = sum(u >= u[1]) / (m + 1)
    return(overall_p)
}


mod_sum = function(x) {
  # this is not the most memory efficient 
  # summary of the model as it does not use
  # cross-tabs however this will make it easy to play with
  # downstream I think when it comes to 
  # adding in null model results
  betas = coef(x)
  r2 = summary(x)$r.squared
  r2adj = summary(x)$adj.r.squared
  coef_type = c(paste0('b', 0:(length(betas) - 1)),
                'r2', 'r2adj')
  out = data.frame(coef_type, 
                   c(betas, r2, r2adj))
  names(out) = c('index', 'value')
  out
}

mod_sum_ct = function(x) {
  # this summary uses cross tabs
  betas = t(coef(x))
  colnames(betas) = paste0('b', 0:(length(betas) - 1))
  r2 = summary(x)$r.squared
  r2adj = summary(x)$adj.r.squared
  out = data.frame(betas, r2, r2adj)
  out
}




#' Conduct the MoB tests on drivers of biodiversity across scales.
#' 
#' There are three tests, on effects of 1. the shape of the SAD, 2.
#' treatment/group-level density, 3. degree of aggregation. The user can
#' specifically to conduct one or more of these tests.
#' 
#' @param mob_in an object of class mob_in created by make_mob_in()
#' @param group_var a character string specifying the environmental variable in
#'   mob_in$env used for grouping plots
#' @param env_var an optional character string specifying a environmental variable
#'   in mob_in$env which is used in correlation analysis in the continuous case. 
#'   It is not needed if "type" is discrete or group_var is used in correlation.
#' @param ref_group a character string used to define the reference group to
#'   which all other groups are compared with when "type" is discrete. It is not
#'   needed when "type" is continuous.
#' @param tests specifies which one or more of the three tests ('SAD', N', 'agg') 
#'   are to be performed. Default is to include all three tests.
#' @param type "discrete" or "continuous". If "discrete", pair-wise comparisons
#'   are conducted between all other groups and the reference group. If
#'   "continuous", a correlation analysis is conducted between the response
#'   variables and group_var or env_var (if defined).
#' @param inds effort size at which the individual-based rarefaction curves are
#'   to be evaluated, and to which the sample-based rarefaction curves are to be
#'   interpolated. It can take three types of values, a single integer, a vector
#'   of integers, and NULL. If inds = NULL (default), the curves are evaluated
#'   at every possible effort size, from 1 to the total number of individuals
#'   within the group (slow). If inds is a single integer, it is taken as the
#'   number of points at which the curves are evaluated; the positions of the
#'   points are determined by the "log_scale" argument. If inds is a vector of
#'   integers, it is taken as the exact points at which the curves are
#'   evaluated.
#' @param log_scale if "inds" is given a single integer, "log_scale" determines
#'   the position of the points. If log_scale is TRUE, the points are equally
#'   spaced on logarithmic scale. If it is FALSE (default), the points are
#'   equally spaced on arithmetic scale.
#' @param min_plots minimal number of plots for test 'agg', where plots are
#'   randomized within groups as null test. If it is given a value, all groups
#'   with fewer plots than min_plot are removed for this test. If it is NULL
#'   (default), all groups are kept. Warnings are issued if 1. there is only one
#'   group left and "type" is discrete, or 2. there are less than three groups
#'   left and "type" is continuous, or 3. reference group ("ref_group") is
#'   removed and "type" is discrete. In these three scenarios, the function will
#'   terminate. A different warning is issued if any of the remaining groups
#'   have less than five plots (which have less than 120 permutations), but the 
#'   test will be carried out.
#' @param density_stat reference density used in converting number of plots to
#'   numbers of individuals, a step in test "N". It can take one of the
#'   three values: "mean", "max", or "min". If it is "mean", the average
#'   plot-level abundance across plots (all plots when "type" is "continuous,
#'   all plots within the two groups for each pair-wise comparison when "type"
#'   is "discrete") are used. If it is "min" or "max", the minimum/maximum
#'   plot-level density is used.
#' @param corr which kind of correlation to use when "type" is "continuous". It
#'   can take two values, "spearman" or "pearson". "spearman" (default) is
#'   generally recommended because the relationship between the response and
#'   "env_var" may not be linear.
#' @param n_perm number of iterations to run for null tests, defaults to 1000.
#' @param overall_p boolean defaults to FALSE specifies if overall across scale 
#'  p-values for the null tests. This should be interpreted with caution because
#'  the overall p-values depend on scales of measurement yet do not explicitly 
#'  reflect significance at any particular scale. 
#' @return a "mob_out" object with attributes
#' @author Xiao Xiao and Dan McGlinn
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in= make_mob_in(inv_comm, inv_plot_attr)
#' inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_group='uninvaded',
#'                            type='discrete', log_scale=TRUE, n_perm=20)

get_delta_stats = function(mob_in, group_var, ref_group = NULL, 
                           tests = c('SAD', 'N', 'agg'),
                           type='discrete', inds = NULL, log_scale = FALSE,
                           min_plots = NULL, density_stat ='mean',
                           n_perm=1000, overall_p = FALSE) {
    #approved_tests = get_delta_overall_checks(mob_in, type, group_var, env_var, 
    #                                          density_stat, tests)
    S = ncol(mob_in$comm)
    groups = mob_in$env[ , group_var]
    if (type == 'discrete') {
        if (class(groups) != 'factor') {
            warning(paste("Converting", group_var, "to a factor because the argument type = 'discrete'"))
            groups = as.factor(groups)
        }
    } else if (type == 'continuous') {
        if (!is.numeric(groups)) {
            warning(paste("Converting", group_var, "to numeric because the argument type = 'continuous'"))
            groups = as.numeric(as.character(groups))
        }    
    } else
        stop('The argument "type" must be either "discrete" or "continuous"')
    #TODO if ref_group specified then the contrats of the group_var
    # need to be set so that the reference group is treated the y-intercept
    # 
    #if (is.null(env_var)){
    #    env_levels = as.numeric(names(sad_groups))
    #} else {
    #    env_levels = tapply(mob_in$env[, env_var],
    #                        list(groups), mean)
    #}
    N_max = min(tapply(rowSums(mob_in$comm), list(groups), sum))
    inds = get_inds(N_max, inds, log_scale)
    ind_dens = get_ind_dens(mob_in$comm, density_stat)

    out = list()
    out$type = type
    #out$tests = approved_tests
    out$log_scale = log_scale
    out$density_stat = list(density_stat = density_stat,
                            ind_dens = ind_dens)
    get_results = function(mob_in, groups, tests, inds, ind_dens, type) {
        
        # the approach taken here to get results for each group
        # is to first break the dataset up into a list of lists 
        # where this is one list per group - this is likely not 
        # the best pratice for memory but it makes the code much 
        # easier to follow - we may need to revisit this. 
        group_levels = unique(groups)
        group_rows = map(group_levels, ~ which(groups == .x))
        mob_in_groups = map(group_rows, ~ subset(mob_in, .x, type = 'integer'))
        names(mob_in_groups) = group_levels
        
        S_df = map_dfr(mob_in_groups, get_delta_curves, tests, inds, ind_dens,
                       .id = "group")
        if (type == 'discrete')
            S_df$group = as.factor(S_df$group)
        if (type == 'continuous')
            S_df$group = as.numeric(S_df$group)
         S_df = as_tibble(S_df)
         S_df = S_df %>% mutate_if(is.factor, as.character)

    
        # now that S and effects computed across scale compute
        # summary statistics at each scale 
    
        delta_mod = function(df) {
            lm(effect ~ group, data = df)
        }
        mod_sum = function(mod) {
            # this is not the most memory efficient 
            # summary of the model as it does not use
            # cross-tabs however this will make it easy to play with
            # downstream I think when it comes to 
            # adding in null model results
            betas = coef(mod)
            r2 = summary(mod)$r.squared
            r2adj = summary(mod)$adj.r.squared
            coef_type = c(paste0('b', 0:(length(betas) - 1)),
                         'r2', 'r2adj')
            out = data.frame(coef_type, 
                             c(betas, r2, r2adj))
            names(out) = c('index', 'value')
            out
        }    
        mod_df = S_df %>%
                 group_by(test, sample, effort) %>%
                 nest() %>%
                 mutate(fit = map(data, delta_mod)) %>%
                 mutate(sum = map(fit, mod_sum)) %>%
                 select(test, sample, effort, sum) %>% 
                 unnest(sum) %>%
                 mutate_if(is.factor, as.character)
        return(list(S_df = S_df, mod_df = mod_df))
    }
    out = append(out, 
                 get_results(mob_in, groups, tests, inds, ind_dens, type))
    run_null_models = function(mob_in, groups, tests, inds, ind_dens, type,
                               n_perm) {
        for (k in seq_along(tests)) {
            null_results = vector('list', length = n_perm)
            cat(paste('\nComputing null model for', tests[k], 'SAD effect\n'))
            pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
            for (i in 1:n_perm) {
                mob_in$comm = get_null_comm(mob_in$comm, tests[k], groups)
                null_results[[i]] = get_results(mob_in, groups, tests[k], inds,
                                                ind_dens, type)$mod_df
                setTxtProgressBar(pb, i)
            }
            close(pb)    
            # rbind across the null_results adding a permutation index
            null_results = tibble(null_results)
            null_df = flatten_dfr(null_results, .id = "perm")
            # compute quantiles
            null_qt = null_df %>%
                      group_by(test, sample, effort, index) %>%
                      summarize(low_value = quantile(value, 0.025, na.rm=T),
                                med_value = quantile(value, 0.5, na.rm=T), 
                                high_value = quantile(value, 0.975, na.rm=T))
            if (k == 1) 
              out = null_qt
            else 
              out = rbind(out, null_qt)
        }
        return(out)
    }
    null_results = run_null_models(mob_in, groups, tests, inds, ind_dens,
                                   type, n_perm)
    # merge the null_results into the model data.frame
    out$mod_df = left_join(out$mod_df, null_results, 
                           by = c("test", "sample", "effort", "index"))
    class(out) = 'mob_out'
    return(out)
}

table_effect_on_S = function(dat_sp, dat_plot, groups, ScaleBy = NA) {
  # Returns a data frame with the effects of SAD, N, and aggregation on diversity
  # across scales
  # not b/c of spatial ties the values will change every time 
  # this is calculated therefore best pratice may be
  # tst = replicate(20, table_effect_on_S(dat_sp, dat_plot, groups, ScaleBy), simplify=FALSE)
  # plyr::aaply(plyr::laply(tst, as.matrix), c(2, 3), mean)
  nplots = table(dat_plot$group)
  explicit_sample = sapply(groups, function(x) 
    rarefy_sample_explicit(dat_sp, dat_plot, x, 1:min(nplots)))
  overall = as.numeric(na.omit(explicit_sample[ , 2] - explicit_sample[ , 1]))
  deltaSsad = get_deltaSsad(dat_sp, dat_plot, groups)
  deltaSN = get_deltaSN(dat_sp, dat_plot, groups, ScaleBy) # why is this call diff
  deltaSagg = get_deltaSagg(dat_sp, dat_plot, groups)
  # Rarefy to desired abundances
  avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
  max_level = floor(log10(avg_dens * min(nplots)))
  out = as.data.frame(matrix(NA, 4, max_level))
  row.names(out) = c("overall", "SAD", "N", "aggregation")
  names(out) = as.character(10^(1:max_level))
  for (row in c(2, 3)) {
    deltaS = unlist(list(overall, deltaSsad, deltaSN, deltaSagg)[row])
    out[row, ] = sapply(10^(1:max_level), function(x) 
      ifelse(length(deltaS) >= x, deltaS[x], NA))
  }
  for (row in c(1, 4)) {
    deltaS = unlist(list(overall, deltaSsad, deltaSN, deltaSagg)[row])
    out_row = pracma::pchip(xi=(0:length(deltaS)) * avg_dens, yi=c(0, deltaS), 
                    x=10^(1:min(max_level, floor(log10(length(deltaS) * avg_dens)))))
    out[row, 1:length(out_row)] = out_row
  }
  out = cbind(out, c(overall[length(overall)], deltaSsad[length(deltaSsad)], 
                     deltaSN[length(deltaSN)], deltaSagg[length(deltaSagg)]))
  names(out)[max_level + 1] = length(deltaSsad)
  return(out)
}

pairwise_t = function(dat_sp, dat_plot, groups, lower_N = NA) {
  dat_plot_grps = dat_plot[dat_plot$group %in% groups, ]
  dat_sp = dat_sp[match(dat_plot_grps$plot, row.names(dat_sp)), ]
  S_list = rowSums(dat_sp > 0)
  N_list = rowSums(dat_sp)
  PIE_list = sapply(1:nrow(dat_sp), function(x) 
    N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, ]/N_list[x])^2)))
  if (is.na(lower_N)) {
    rarefied_S_list = apply(dat_sp, 1, function(x) 
      rarefaction(x, 'indiv', effort = 1:min(N_list)))
  } else {
    # Remove plots with abundance below lower_N in the analysis of rarefied S
    rarefied_S_list = apply(dat_sp, 1, function(x) 
      if (sum(x) < lower_N)
        rep(NA, lower_N)
      else 
        rarefaction(x, 'indiv', effort = 1:lower_N))
    if (any(is.na(rarefied_S_list))) 
      print("Warning: some plots are removed in rarefaction.")
  }
  out = as.data.frame(matrix(NA, 5, 4))
  stats_list = list(rarefied_S_list, N_list, PIE_list, S_list)
  for (i in 1:length(stats_list)) {
    stat = unlist(stats_list[i])
    stat_1 = stat[dat_plot$group == groups[1]]
    stat_2 = stat[dat_plot$group == groups[2]]
    stat_1 = stat_1[!is.na(stat_1)]
    stat_2 = stat_2[!is.na(stat_2)]
    out[ , i] = c(mean(stat_1), sd(stat_1), 
                  mean(stat_2), sd(stat_2), 
                  t.test(stat_1, stat_2)$p.val)
  }
  names(out) = c("S_rarefied", "N", "PIE", "S_raw")
  row.names(out) = c(paste(groups[1], "(mean)", sep = ""), 
                     paste(groups[1], "(sd)", sep = ""), 
                     paste(groups[2], "(mean)", sep = ""), 
                     paste(groups[2], "(sd)", sep = ""), "p_value")
  # Boxplots
  par(mfrow = c(2, 2))  # This is not ideal but I cannot get layout to work in Rstudio
  plot_names = c(paste("Rarified S at N=", 
                       ifelse(is.na(lower_N), min(N_list), lower_N), sep = ""),
                 "N", "PIE", "Raw S")
  plot_names = sapply(1:4, function(x) 
    paste(plot_names[x], " (p=", round(out[5, x], 6), ")", sep = ""))
  for (i in 1:length(stats_list)) {
    stat = unlist(stats_list[i])
    stat_1 = stat[dat_plot$group == groups[1]]
    stat_2 = stat[dat_plot$group == groups[2]]
    stat_1 = stat_1[!is.na(stat_1)]
    stat_2 = stat_2[!is.na(stat_2)]
    boxplot(stat_1, stat_2, names = c(groups[1], groups[2]), main = plot_names[i])
  }
  return(out)
}

#' Plot distributions of species abundance
#'
#' @param mob_in a 'mob_in' class object produced by 'make_mob_in'
#' @param env_var a string that specifies the column name in mob_in$env that
#'   specifies the grouping variable.
#' @param type either 'sad' or 'rad' for species abundance vs rank abundance
#'   distribution
#' @param pooled boolean specifying if abundances should be pooled at the group
#'   level or not
#' @param col optional vector of colors.
#' @inheritParams plot.mob_out
#' @inheritParams graphics::plot.default
#' @importFrom scales alpha
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' plot_abu(inv_mob_in, 'group', type='sad', pooled=FALSE, log='x')
#' plot_abu(inv_mob_in, 'group', type='rad', pooled=TRUE, log='x')
plot_abu = function(mob_in, env_var, type=c('sad', 'rad'),
                    pooled=FALSE, col=NULL, lwd=3, log='', leg_loc = 'topleft') {
    env_data = mob_in$env[ , env_var]
    grps = sort(unique(as.character(env_data)))
    if (is.null(col)) 
        col = c("#FFB3B5", "#78D3EC", "#6BDABD", "#C5C0FE",
                "#E2C288", "#F7B0E6", "#AAD28C")    
    else if (length(col) != length(grps))
      stop('Length of col vector must match the number of unique groups')
    title = ifelse(pooled, 'Group Scale', 'Sample Scale')
    if ('sad' == type) {
        plot(1, type = "n", xlab = "% abundance", ylab = "% species", 
             xlim = c(0.01, 1), ylim = c(0.01, 1), log = log, main = title)
        for (i in 1:length(grps)) {
            col_grp = col[i]
            comm_grp = mob_in$comm[env_data == grps[i], ]
            comm_grp = comm_grp[rowSums(comm_grp) > 0, ]
            if (pooled) {
                sad_grp = colSums(comm_grp)
                sad_sort = sort(sad_grp[sad_grp != 0])
                s_cul = 1:length(sad_sort) / length(sad_sort)
                n_cul = sapply(1:length(sad_sort), function(x)
                               sum(sad_sort[1:x]) / sum(sad_sort))
                lines(n_cul, s_cul, col = col_grp, lwd = lwd, type = "l")
            } else {
                for (j in 1:nrow(comm_grp)) {
                    sad_sort = sort(comm_grp[j, comm_grp[j, ] != 0])
                    s_cul = 1:length(sad_sort) / length(sad_sort)
                    n_cul = sapply(1:length(sad_sort), function(x)
                                   sum(sad_sort[1:x]) / sum(sad_sort))
                    lines(n_cul, s_cul, col = scales::alpha(col_grp, 0.5), 
                          lwd = lwd, type = "l")
                }
            }
        }
    } 
    if ('rad' == type) {
        plot(1:10, 1:10, type='n', xlab='rank', ylab='abundance',
             log=log, xlim=c(1, ncol(mob_in$comm)), 
             ylim=range(0.01, 1), cex.lab = 1.5, cex.axis = 1.5,
             main = title)
        for (i in 1:length(grps)) {
             col_grp = col[i]
             comm_grp = mob_in$comm[env_data == grps[i], ]
             comm_grp = comm_grp[rowSums(comm_grp) > 0, ]
             if (pooled) {
                sad_grp = colSums(comm_grp)
                sad_sort = sort(sad_grp[sad_grp != 0], dec=T)
                lines(sad_sort / sum(sad_sort), col = col_grp, lwd = lwd,
                      type = "l")
             } else {
                 for (j in 1:nrow(comm_grp)) {
                     sad_sort = sort(comm_grp[j, comm_grp[j, ] != 0], dec=T)
                     lines(1:length(sad_sort), sad_sort / sum(sad_sort),
                           col = scales::alpha(col_grp, 0.5),
                           lwd = lwd, type = "l")
                 }     
             }
        }
    }
    if (!is.na(leg_loc))
        legend(leg_loc, legend=grps, col = col, lwd = lwd, bty='n')
}
    
#' Plot rarefaction curves for each treatment group
#' 
#' @param pooled boolean specifying if samples should be pooled at the group
#'  level or not. Defaults to TRUE. This argument only applies when
#'  the individual based rarefaction is used (i.e., method = 'indiv')
#' @param ... other arguments to provide to \code{\link[mobr]{rarefaction}}
#' @inheritParams plot.mob_out
#' @inheritParams plot_abu
#' @inheritParams rarefaction
#' @importFrom scales alpha
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' # random individual based rarefaction curves
#' plot_rarefaction(inv_mob_in, 'group', 'indiv',
#'                  pooled=TRUE, leg_loc='bottomright')
#' plot_rarefaction(inv_mob_in, 'group', 'indiv',
#'                  pooled=FALSE, log='x')
#' # random sample based rarefaction curves 
#' plot_rarefaction(inv_mob_in, 'group', 'samp', log='xy')
#' # spatial sample based rarefaction curves 
#' plot_rarefaction(inv_mob_in, 'group', 'spat', log='xy',
#'                  coords = inv_mob_in$spat)                 
plot_rarefaction = function(mob_in, env_var, method, dens_ratio=1, pooled=T, 
                            col=NULL, lwd=3, log='', leg_loc = 'topleft',
                            ...) {
    if (pooled == FALSE & method != 'indiv')
        stop('Samples can only not be pooled at the treatment level when individual-based rarefaction is used (i.e., method="indiv")')
    env_data = mob_in$env[ , env_var]
    grps = sort(unique(as.character(env_data)))
    if (is.null(col)) 
        col = c("#FFB3B5", "#78D3EC", "#6BDABD", "#C5C0FE",
                "#E2C288", "#F7B0E6", "#AAD28C")    
    else if (length(col) != length(grps))
        stop('Length of col vector must match the number of unique groups')
    if (method == 'indiv')
        xlab = 'Number of individuals'
    else
        xlab = 'Number of samples'
    if (pooled) {
        Srare = lapply(grps, function(x) 
                       rarefaction(subset(mob_in, env_data == x, 'logical'),
                                   method, ...))
        xlim = c(1, max(unlist(sapply(Srare, function(x) as.numeric(names(x))))))
        ylim = c(1, max(unlist(Srare)))
        n = as.numeric(names(Srare[[1]]))
        plot(n, Srare[[1]], type = "n", main = "Group scale",
             xlab = xlab, ylab = "Species richness", 
             xlim = xlim, ylim = ylim, log = log)
        for (i in seq_along(grps)) {
            col_grp = col[i]
            n = as.numeric(names(Srare[[i]]))
            lines(n, Srare[[i]], col = col_grp, lwd = lwd, type = "l")
        }
    } else {
        Srare = lapply(grps, function(x)
                       apply(mob_in$comm[env_data == x, ], 1,
                             function(y)  rarefaction(y, method, ...)))
        xlim = c(1, max(unlist(lapply(Srare, function(x)
                                 lapply(x, function(y)
                                        as.numeric(names(y)))))))
        ylim = c(1, max(unlist(Srare)))
        n = as.numeric(names(Srare[[1]][[1]]))
        plot(n, Srare[[1]][[1]], type = "n", main = "Sample scale",
             xlab = xlab, ylab = "Species richness", 
             xlim = xlim, ylim = ylim, log = log)        
        for (i in seq_along(grps)) {
            col_grp = col[i]
            for (j in seq_along(Srare[[i]])) {
                 n = as.numeric(names(Srare[[i]][[j]]))
                 if (n[1] > 0)
                     lines(n, Srare[[i]][[j]], col = scales::alpha(col_grp, 0.5),
                           lwd = lwd, type = 'l')
            }
        }
    }
    if (!is.na(leg_loc))
        legend(leg_loc, legend=grps, col = col, lwd = lwd, bty='n')
}

#' Plot mob curves
#' 
#' @param mob_out a mob_out class object
#' @param trt_group a string that specifies the name of the treatment group  
#' @param ref_group a string that specifies the name of the reference group
#' @param same_scale a boolean if TRUE then the y-axis of the rarefaction and 
#'  ddelta S plots are scaled identically across the tested effects
#' @param display a string that specifies what graphics to display can be either
#'  'rarefaction', 'delta S', or 'ddelta S' defaults to all three options.
#' @param lwd the line width, a single positive number, defaults to \code{3},
#'  see \code{\link[graphics]{par}} for more information.
#' @param leg_loc the location of the legend. Defaults to 'topleft', see
#'   \code{\link[graphics]{legend}}. If set to NA then no legend is printed.
#' @param par_args optional argument that sets graphical parameters to set
#' @param ... Other plot arguments
#' 
#' @return plots the effect of the SAD, the number of individuals, and spatial
#'  aggregation on the difference in species richness
#'  
#' @author Xiao Xiao and Dan McGlinn
#' @inheritParams graphics::plot.default
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_group='uninvaded',
#'                               type='discrete', log_scale=TRUE, n_perm=2)
#' plot(inv_mob_out, 'invaded', 'uninvaded', display='rarefaction')
#' plot(inv_mob_out, 'invaded', 'uninvaded', display='delta S')
#' plot(inv_mob_out, 'invaded', 'uninvaded', display='ddelta S')
plot.mob_out = function(mob_out, trt_group, ref_group, same_scale=FALSE, 
                        log='', display=c('rarefaction', 'delta S', 'ddelta S'),
                        lwd=3, leg_loc='topleft', par_args=NULL, ...) {
    type = mob_out$type
    tests = mob_out$tests
    if (type == 'continuous')
        stop("Currently this plot only works for mob_out object with type discrete.")
    cols = list()
    cols$trt = "#FFB3B5"     # light red
    cols$ref = "#78D3EC"     # light blue
    cols$deltaS = "#C5C0FE"  # purple
    cols$ddeltaS = "#6BDABD" # green
    if (is.null(par_args)) {
        par_args = paste('mfrow = c(', length(display), ',',
                         length(tests), '), mgp = c(2.5, 1, 0)',  sep='')
        
    } 
    eval(parse(text=paste('par(', par_args, ')')))
    if (same_scale) {
        # not currently implemented for the delta S plots
        if ('rarefaction' %in% display) {
            if ('agg' %in% tests) 
                S_cols = c('impl_S', 'expl_S')
            else
                S_cols = 'impl_S'
            ylim_rare = range(list(mob_out$indiv_rare[ , -1],
                                   mob_out$sample_rare[ , S_cols]))
        }
    }
    x_axis_min = 1
    mob_out$sample_rare[, -1] = lapply(mob_out$sample_rare[, -1], function(x)
                                       as.numeric(as.character(x)))
    sample_rare_group = mob_out$sample_rare[mob_out$sample_rare == trt_group, ]
    sample_rare_ref = mob_out$sample_rare[mob_out$sample_rare == ref_group, ]
    groups = c(trt_group, ref_group)
    if ('rarefaction' %in% display) {
        if ('agg' %in% mob_out$tests) {
          if (!same_scale)
            ylim_rare = c(0, max(mob_out$sample_rare$expl_S))
          for (i in 1:length(groups)){
            group = groups[i]
            dat_group = mob_out$sample_rare[mob_out$sample_rare$group == group, ]
            if (i == 1)
              plot(dat_group$sample_plot, dat_group$expl_S, lwd = lwd,
                   type = 'l', xlab = 'Number of plots',
                   ylab = 'Richness (S)', col = cols$trt,
                   xlim = c(x_axis_min, max(dat_group$sample_plot)),
                   ylim = ylim_rare,
                   main = 'sSBR', cex.axis = 1.5, cex.lab = 1.5,
                   log=log, frame.plot=F, ...)
            else
              lines(dat_group$sample_plot, dat_group$expl_S,
                    lwd = lwd, col = cols$ref)
          }
          if (!is.na(leg_loc))
              legend(leg_loc, legend=as.character(groups),
                     col=as.character(unlist(cols)), lty=1, lwd=lwd, bty='n')
        }
        if ('N' %in% mob_out$tests) {
            if (!same_scale)
                ylim_rare = c(0, max(mob_out$sample_rare$impl_S))
            for (i in 1:length(groups)){
                group = groups[i]
                dat_group = mob_out$sample_rare[mob_out$sample_rare$group == group, ]
                if (i == 1)
                    plot(dat_group$sample_plot, dat_group$impl_S,
                         lwd = lwd, type = 'l', xlab = 'Number of plots',
                         ylab = 'Richness (S)', col = cols$trt, 
                         xlim = c(x_axis_min, max(dat_group$sample_plot)),
                         ylim = ylim_rare,
                         main = 'nsSBR', cex.axis = 1.5, cex.lab = 1.5,
                         log=log, frame.plot=F, ...)
                else
                    lines(dat_group$sample_plot, dat_group$impl_S,
                          lwd = lwd, col = cols$ref)
            }
        }
        if ('SAD' %in% mob_out$tests) {
          if (!same_scale)
            ylim_rare = range(mob_out$indiv_rare[, -1])
          plot(mob_out$indiv_rare$sample, mob_out$indiv_rare[, trt_group], 
               lwd = lwd, type = 'l', col = cols$trt, xlab = 'Number of individuals', 
               ylab = 'Richness (S)', main = 'IBR', 
               xlim = c(x_axis_min, max(mob_out$indiv_rare$sample)), ylim = ylim_rare, 
               cex.axis = 1.5, cex.lab = 1.5, log=log, frame.plot=F, ...)
          lines(mob_out$indiv_rare$sample, mob_out$indiv_rare[, ref_group], 
                lwd = lwd, col = cols$ref)
        }        
    }    
    if ('delta S' %in% display) {
      minN = min(nrow(sample_rare_group), nrow(sample_rare_ref))
      if ('agg' %in% mob_out$tests) {
        delta_Sspat = sample_rare_group$expl_S[1:minN] - 
          sample_rare_ref$expl_S[1:minN]
        plot(seq(minN), delta_Sspat, 
             ylim = c(min(delta_Sspat, 0), max(delta_Sspat, 0)),
             cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = lwd,
             col = cols$deltaS, xlab = 'Number of plots',
             ylab = expression(Delta * 'S due to SAD, N, & agg.'), frame.plot=F, ...)
        abline(h = 0, lwd = 1, lty = 2)
      }
        if ('N' %in% mob_out$tests) {
            delta_Ssample = sample_rare_group$impl_S[1:minN] - 
                            sample_rare_ref$impl_S[1:minN]
            plot(seq(minN), delta_Ssample, 
                 ylim = c(min(delta_Ssample, 0), max(delta_Ssample, 0)),
                 cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = lwd,
                 col = cols$deltaS, xlab = 'Number of plots', 
                 ylab = expression(Delta * 'S due to SAD & N'), frame.plot=F, ...)
            abline(h = 0, lwd = 1, lty = 2)
        }
      if ('SAD' %in% mob_out$tests) {
        # Create the plots for the three delta-S between groups
        deltaS_Sind = mob_out$indiv_rare[[trt_group]] - 
          mob_out$indiv_rare[[ref_group]]
        plot(mob_out$indiv_rare$sample, deltaS_Sind,
             ylim = c(min(deltaS_Sind, 0), max(deltaS_Sind, 0)),
             cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = lwd,
             col = cols$deltaS, xlab = 'Number of individuals', 
             ylab = expression(Delta * 'S due to SAD'), log=log,
             frame.plot=F, ...)
        abline(h = 0, lwd = 1, lty = 2)
        
      }       

    }
    if ('ddelta S' %in% display) {
        # Create the plots for the three ddelta S
         ylim_ddelta = range(lapply(mob_out[tests], function(x)
                              lapply(x[ , -(1:2)], function(y)
                                     as.numeric(as.character(y)))))
         if ('agg' %in% mob_out$tests) {
           mob_out$agg[, -1] = lapply(mob_out$agg[, -1], function(x)
             as.numeric(as.character(x))) 
           ddelta_Sspat = mob_out$agg[which(as.character(mob_out$agg$group) == as.character(trt_group)), ]
           if (!same_scale)
             ylim = range(ddelta_Sspat[ , -(1:2)])
           plot(ddelta_Sspat$effort_sample, ddelta_Sspat$ddeltaS_emp,
                ylim = ylim_ddelta, log=log,
                cex.axis = 1.5, cex.lab = 1.5, type = 'n', 
                xlab = 'Number of plots', 
                ylab = expression(Delta * 'S due to agg.'), frame.plot=F, ...)
           polygon(c(ddelta_Sspat$effort_sample,
                     rev(ddelta_Sspat$effort_sample)), 
                   c(ddelta_Sspat$ddeltaS_null_low,
                     rev(ddelta_Sspat$ddeltaS_null_high)),
                   col = '#C1CDCD', border = NA)
           abline(h = 0, lwd = 1, lty = 2)
           lines(ddelta_Sspat$effort_sample, ddelta_Sspat$ddeltaS_emp, 
                 lwd = lwd, col = cols$ddeltaS)
         }
         if ('N' %in% mob_out$tests) {
            mob_out$N[, -1] = lapply(mob_out$N[, -1], function(x)
                                     as.numeric(as.character(x))) 
            ddelta_Ssample = mob_out$N[which(as.character(mob_out$N$group) == as.character(trt_group)), ]
            if (!same_scale)
                ylim = range(ddelta_Ssample[ , -(1:2)])
            plot(ddelta_Ssample$effort_sample, ddelta_Ssample$ddeltaS_emp,
                 ylim = ylim_ddelta, log=log,
                 cex.axis = 1.5, cex.lab = 1.5, type = 'n', 
                 xlab = 'Number of individuals', 
                 ylab = expression(Delta * 'S due to N'), frame.plot=F, ...)
            polygon(c(ddelta_Ssample$effort_sample, 
                      rev(ddelta_Ssample$effort_sample)), 
                    c(ddelta_Ssample$ddeltaS_null_low, 
                      rev(ddelta_Ssample$ddeltaS_null_high)),
                    col = '#C1CDCD', border = NA)
            abline(h = 0, lwd = 1, lty = 2)
            lines(ddelta_Ssample$effort_sample, ddelta_Ssample$ddeltaS_emp,
                  lwd = lwd, col = cols$ddeltaS)
        }
         if ('SAD' %in% mob_out$tests) {
           mob_out$ind[, -1] = lapply(mob_out$ind[, -1], function(x)
             as.numeric(as.character(x))) 
           delta_Sind = mob_out$SAD[which(as.character(mob_out$SAD$group) == as.character(trt_group)), ]
           if (!same_scale)
             ylim = range(delta_Sind[ , -(1:2)])
           plot(delta_Sind$effort_ind, delta_Sind$deltaS_emp, 
                ylim = ylim_ddelta, log=log,
                cex.axis = 1.5, cex.lab = 1.5, type = 'n',
                xlab = 'Number of individuals', ylab = expression(Delta * 'S due to SAD'),
                frame.plot=F, ...)
           polygon(c(delta_Sind$effort_ind, rev(delta_Sind$effort_ind)), 
                   c(delta_Sind$deltaS_null_low, rev(delta_Sind$deltaS_null_high)),
                   col = '#C1CDCD', border = NA)
           abline(h = 0, lwd = 1, lty = 2)
           lines(delta_Sind$effort_ind, delta_Sind$deltaS_emp,
                 lwd = lwd, col = cols$ddeltaS)
         }         

    }
}


#' Plot summary graphics of the effect on species richness
#' 
#' All three treatment effect sizes due to the SAD, N, or aggregation are
#' graphed on a single plot. The treatment effect is defined as the treatment difference
#' in richness due to a particular component of community structure. The overlap
#' of the three effects is accomplished by rescaling numbers of individuals 
#' to the number of plots. The effect sizes may be plotted simply in their raw 
#' form or as stacked area plots. 
#' @param display the type of graphic to display options include: "raw" or
#'  "stacked" for plots of the raw effect sizes or a stacked area plot of the
#'  absolute value of the effect sizes. 
#' @param prop boolean if TRUE then proportions are used in the stacked area plot
#' @param rescale string that specifies how to rescale number of individuals to
#'  the number of samples. Defaults to 'max_effort' which rescales by 
#'  setting the maximum number of individuals considered in the individual
#'  rarefaction curve to the scale of the maximum number of samples considered
#'  in the spatial rarefaction curve. The other option is 'density_stat' which
#'  rescales individuals to number of samples using the density statistic 
#'  specified when the mob statistics where computed by \code{\link[mobr]{get_delta_stats}}. 
#' @param common_scale boolean defaults to FALSE. If TRUE then all the effects
#'  are truncated to only be across the same range of number of individuals.
#' @param xlabel_indiv boolean if TRUE then the top axis is 
#'  labeled with numbers of individuals
#' @param lty a vector of line types, see \link[graphics]{par}.
#' @param col the colors for lines. Three colors can be specified so that each
#'  line can be given its own color with the curves ordered as SAD, N, and
#'  aggregation.
#' @inheritParams plot.mob_out
#' @inheritParams graphics::plot.default
#' @importFrom pracma pchip
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_group='uninvaded',
#'                               type='discrete', log_scale=TRUE, n_perm=2)
#' overlap_effects(inv_mob_out, 'invaded', leg_loc=NA)
#' overlap_effects(inv_mob_out, 'invaded', display='stacked')
#' overlap_effects(inv_mob_out, 'invaded', display='stacked', prop=TRUE)
overlap_effects = function(mob_out, trt_group, display='raw', prop=FALSE,
                           rescale='max_effort', common_scale=FALSE, 
                           xlabel_indiv=TRUE, ylim=NULL, log='', lty=1, lwd=3,
                           leg_loc='topleft', col=c("#FFB3B5", "#78D3EC", "#C5C0FE")) {
    if (prop & display != 'stacked')
        stop("Proptional differences can only be used when considering stacked area graphs (i.e., display = 'stacked')")
    if (length(lty) == 1) 
        lty = rep(lty, 3)
    tests = mob_out$tests
    SAD = data.frame(type='SAD', 
                    mob_out$SAD[mob_out$SAD$group == trt_group, 
                                c('effort_ind', 'deltaS_emp')])
    N = data.frame(type='N', 
                   mob_out$N[mob_out$N$group == trt_group,
                             c('effort_sample', 'ddeltaS_emp')])
    agg = data.frame(type='agg',
                     mob_out$agg[mob_out$agg$group == trt_group,
                                 c('effort_sample', 'ddeltaS_emp')])
    names(SAD) = names(N) = names(agg) =  c('type', 'effort', 'effect')
    if (rescale == 'max_effort') {
        virt_effort = seq(min(SAD$effort), max(SAD$effort),
                          length.out = length(agg$effort))
        effort = agg$effort
    } else if (rescale == 'density_stat'){
        N_plots = max(agg$effort)
        N_indiv = max(SAD$effort)
        ind_dens = mob_out$density_stat$ind_dens
        # the next line assumes that the density stat is average (may need to generalize)
        effort = min(agg$effort) : round(N_indiv / N_plots)
        virt_effort = effort * ind_dens
    } else
        stop('rescale must be specified as "max_effort" or "density_stat" see documentation')
    SAD_interp = pracma::pchip(SAD$effort, SAD$effect, virt_effort)
    N_interp = pracma::pchip(N$effort, N$effect, virt_effort)    
    SAD = data.frame(type='SAD', effort, effect=SAD_interp)
    N = data.frame(type='N', effort=effort, effect=N_interp)
    agg = agg[effort, ]
    dat = rbind(SAD, N, agg)
    dat$abs_effect = abs(dat$effect)
    if (common_scale)
        dat = subset(dat , effort <= max(dat$effort[dat$type == 'SAD']))
    if (is.null(ylim)) {
        if (display == 'raw')
           ylim = range(dat$effect)
        else if (prop)
           ylim = c(0, 1)
    }
    if (display == 'raw') {
       plot(effect ~ effort, data=dat, subset= type == 'SAD',
            xlab = 'Number of Samples',
            ylab = 'Difference in Richness',
            frame.plot=F, ylim=ylim, type='l', col=col[1],
            lwd=lwd, lty=lty[1], log=log)
       lines(unique(dat$effort), dat$effect[dat$type == 'N'], col=col[2], 
             lwd=lwd, lty=lty[2])
       lines(unique(dat$effort), dat$effect[dat$type == 'agg'], col=col[3], 
             lwd=lwd, lty=lty[3])
       abline(h=0, lty=2, lwd=2)
       if (!is.na(leg_loc))
         legend(leg_loc, legend=c('SAD', 'N', 'Agg.'), col=col, lwd=lwd, 
                bty='n')       
    }
    if (display == 'stacked') {
        if (prop) {
            props = unlist(tapply(dat$abs_effect, dat$effort,
                                  function(x) x / sum(x)))
            dat = data.frame(type = rep(c('SAD', 'N', 'agg'), 
                                        times=length(effort)),
                            effort = rep(effort, each=3),
                            abs_effect = props)
            ylab = 'Fraction of abs(Diff. in Richness)'
         } else
            ylab = 'abs(Diff. in Richness)'
         plotStacked(unique(dat$effort),
                     data.frame(dat$abs_effect[dat$type == 'SAD'],
                                dat$abs_effect[dat$type == 'N'],
                                dat$abs_effect[dat$type == 'agg']),
                     xlab = 'Number of Samples', ylab = ylab ,
                     col = col, border=NA,
                     frame.plot=F, ylim=ylim, log=log)
         if (!is.na(leg_loc))
           legend(leg_loc, legend=c('SAD', 'N', 'Agg.'), fill=col,
                  border=0, box.lwd=0, box.col=0)         
    }
    ticks = axTicks(side=3)
    n_indices = round(seq(1, length(virt_effort), 
                          length.out=length(ticks)))
    if (xlabel_indiv) {
        axis(side=3, at=effort[n_indices], 
             labels=round(virt_effort[n_indices]))
        mtext(side=3, "Number of Individuals", padj=-4.5)
    }

}

#' Plot the relationship between the number of plots and the number of
#' individuals
#' 
#' The MoB methods assume a linear relationship between the number of 
#' plots and the number of individuals. This function provides a means of 
#' verifying the validity of this assumption
#' @param comm community matrix with sites as rows and species as columns
#' @param n_perm number of permutations to average across, defaults to 1000
#' @author Dan McGlinn
#' @export
#' @examples
#' data(inv_comm)
#' plot_N(inv_comm)
plot_N = function(comm, n_perm=1000) {
    N = rowSums(comm)
    ind_dens = mean(N)
    N_sum = apply(replicate(n_perm, cumsum(sample(N))), 1, mean)
    plot(N_sum, xlab='Number of plots', ylab='Number of Individuals')
    abline(a=0, b=ind_dens, col='red')
    legend('topleft', 'Expected line', lty=1, bty='n', col='red')
}

#' @title Stacked plot by Marc Taylor (@marchtaylor on gitHub)
#' @description \code{plotStacked} makes a stacked plot where each \code{y} 
#' series is plotted on top of each other using filled polygons.
#' @param x A vector of values
#' @param y A matrix of data series (columns) corresponding to x
#' @param order.method Method of ordering y plotting order. One of the 
#'   following: \code{c("as.is", "max", "first")}. \code{"as.is"} - plot in 
#'   order of y column. \code{"max"} - plot in order of when each y series 
#'   reaches maximum value. \code{"first"} - plot in order of when each y series
#'   first value > 0.
#' @param col Fill colors for polygons corresponding to y columns (will recycle).
#' @param border Border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#' @param lwd Border line width for polygons corresponding to y columns (will recycle)
#' @param xlab x-axis labels
#' @param ylab y-axis labels
#' @param ylim y-axis limits. If \code{ylim=NULL}, defaults to \code{c(0, 1.2*max(apply(y,1,sum)}.
#' @param ... Other plot arguments
#' 
#' @importFrom graphics plot polygon par
#' @importFrom grDevices rainbow 
#' @keywords internal
plotStacked <- function(
	x, y, 
	order.method="as.is",
	ylab="", xlab="", 
	border = NULL, lwd=1, 
	col=rainbow(length(y[1,])),
	ylim=NULL,
	...
){

	if(sum(y < 0) > 0) stop("y cannot contain negative numbers")

	if(is.null(border)) border <- par("fg")
	border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
	col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
	lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))

  if(is.null(ylim)) ylim=c(0, 1.2*max(apply(y,1,sum)))
  
	if(order.method == "max") {
		ord <- order(apply(y, 2, which.max))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	if(order.method == "first") {
		ord <- order(apply(y, 2, function(x) min(which(x>0))))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	top.old <- x*0
	polys <- vector(mode="list", ncol(y))
	for(i in seq(polys)){
		top.new <- top.old + y[,i]
		polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
		top.old <- top.new
	}

	if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
	plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
	for(i in seq(polys)){
		polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
	}

}

