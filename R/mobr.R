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
    out = list()
    out$tests = list(N=T, SAD=T, agg=T)
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
}

#' Print a shortened version of the mob_out object
#' @keywords internal
#' @export
print.mob_out = function(x) {
    cat('Only the first five rows of any matrices are printed\n')
    cat('\n$type\n')
    print(x$type)
    cat('\n$tests\n')
    print(x$tests)
    cat('\n$log_scale\n')
    print(x$log_scale)
    cat('\n$density_stat\n')
    print(x$density_stat)
    cat('\n$indiv_rare\n')
    print(head(x$indiv_rare))
    cat('\n$sample_rare\n')
    print(head(x$sample_rare))
    if (!is.null(x$overall_p)) {
        cat('\n$overall_p\n')
        print(x$overall_p)
    }        
    cat('\n$SAD\n')
    print(head(x$SAD))
    cat('\n$N\n')
    print(head(x$N))
    cat('\n$agg\n')
    print(head(x$agg))
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


#' Difference in S due to N
#' 
#' Internal function for computing the difference in species richness between 
#' individual-based and non-spatial sample-based rarefaction curves using an 
#' analytical approach of stretching (ref_dens > group_dens) or shrinking
#' (ref_dens < group_dens) the individual-based rarefaction curve to provide the
#' non-spatial rarefaction result.
#' 
#' @param comm community matrix with plots as rows and species columns.
#' @param ref_dens the reference density
#' @param inds the number of individuals to sample over
#' @description  Difference between the individual and non-spatial sample-based
#'   rarefaction curves for one group with the evaluation sample size (number of
#'   individuals) defined by ref_dens, evaluated at specified points (given by
#'   inds). The rescaling of sampling effort from number of samples to number
#'   of individuals is accomplished using the mean density of individuals 
#'   per sample.
#' @return a two column data.frame containing the number of individuals (inds)
#'   and the difference in species richness (deltaS)
#' @author Dan McGlinn and Xiao Xiao
#' @keywords internal
deltaS_N = function(comm, ref_dens, inds){
    nplots = nrow(comm)
    group_dens = sum(comm) / nplots
    dens_ratio = ref_dens / group_dens
    S_rescaled = rarefaction(comm, 'indiv', inds, dens_ratio=dens_ratio)
    S_raw = rarefaction(comm, 'indiv', inds)
    deltaS = S_rescaled - S_raw
    out = data.frame(inds = inds, deltaS = deltaS)
    return(out)
}

#' Permute community matrix within groups
#' 
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
#' @param comm community matrix with plots as rows and species columns.
#' @param method either 'noagg' for random shuffling of the individuals without
#'   maintaining the vector of sample total abundances or 'swapN' for random 
#'   shuffling of the individuals in which sample abundances are maintained 
#' @param groups optional argument that is a vector of group ids which specify
#'   which group each site is associated with. If is NULL then all rows of the
#'   community matrix are assumed to be members of the same group
#'   
#' @return a permuted site-by-species matrix
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
#' permute_comm(comm, 'noagg')
#' permute_comm(comm, 'noagg', groups)
#' permute_comm(comm, 'swapN')
#' permute_comm(comm, 'swapN', groups)
permute_comm = function(comm, method, groups=NULL) {
    if (!(is.matrix(comm) | is.data.frame(comm)))
        stop('comm must be a matrix or data.frame')
    if (is.null(groups))
        groups = rep(1, nrow(comm)) 
    group_levels = unique(groups)
    S = ncol(comm)
    N = rowSums(comm)
    comm_group_perm = matrix(0, ncol=S, nrow=nrow(comm))
    if (method == 'swapN')
        swapN = sample(N)
    for(i in seq_along(group_levels)) {
        row_indices = groups == group_levels[i]
        group_comm = comm[row_indices, ]
        sp_abu = colSums(group_comm)
        sp_freq = sp_abu[sp_abu > 0] / sum(sp_abu)
        plot_ids = 1:nrow(group_comm)
        if (method == 'noagg') {
            Ngroup = N[row_indices]
        } else if (method == 'swapN') {
            Ngroup = swapN[row_indices]
        }
        else 
            stop('The argument swap must be either "noagg" or "swapN"')
        sp_draws = sapply(plot_ids, function(x)
            sample(1:length(sp_freq), size=Ngroup[x], 
                   replace=T, prob=sp_freq))
        tmp_comm = t(sapply(plot_ids, function(x)
            table(c(sp_draws[[x]], 1:length(sp_freq))) - 1 ))
        # The following lines are necessary because tmp_comm may have more
        # columns than comm_group_perm
        comm_new = matrix(0, nrow = nrow(comm_group_perm), 
                          ncol = max(ncol(comm_group_perm), ncol(tmp_comm)))
        comm_new[, 1:ncol(comm_group_perm)] = comm_group_perm
        comm_new[row_indices, 1:ncol(tmp_comm)] = tmp_comm
        comm_group_perm = comm_new
    }  
    return(comm_group_perm)
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
get_delta_ind_sample = function(group_sad, inds, log_scale){
    group_minN = min(rowSums(group_sad))
    if (is.null(inds)){
        ind_sample_size = seq(group_minN)
    } else if (length(inds) > 1) {
        ind_sample_size = inds
        if (max(inds) > group_minN)
            warning('Sample size is higher than abundance of at least one group!')
    } else if (log_scale == T){
        ind_sample_size = floor(exp(seq(inds) * log(group_minN) / inds))
    } else {
        ind_sample_size = floor(seq(inds) * group_minN / inds)
    }
    ind_sample_size = unique(c(1, ind_sample_size)) # Force (1, 1) to be included
    return(ind_sample_size)
}

#' Auxiliary function for get_delta_stats()
#' Returns the "assumed" plot density given 
#' whether min, max or mean is used
#' @keywords internal
get_plot_dens = function(comm, density_stat){
    if (density_stat == 'mean') {
        plot_dens = sum(comm) / nrow(comm)
    } else if (density_stat == 'max') {
        plot_dens = max(rowSums(comm))
    } else {
        plot_dens = min(rowSums(comm))
    }
    return(plot_dens)   
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

#' Auxiliary function for get_delta_stats()
#' Obtain the swap curve and/or spatial curve for each group if asked
#' Directly add attributes to the input "out"
#' @keywords internal
get_sample_curves = function(mob_in, group_levels, group_data, approved_tests){
    if ('N' %in% approved_tests | 'agg' %in% approved_tests){
        sample_rare = data.frame(matrix(0, nrow = 0, ncol = 4), 
                                 stringsAsFactors = F)
        for (level in group_levels){
            comm_level = mob_in$comm[as.character(group_data) == level, ]
            nplots = nrow(comm_level)
            level_dens = sum(comm_level) / nplots
            samp_effort = round((1:nplots) * level_dens)
            impl_S = rarefaction(comm_level, 'indiv', samp_effort)
            sample_rare_level = data.frame(cbind(rep(level, length(impl_S)), 
                                                 seq(length(impl_S)), impl_S))
            if ('agg' %in% approved_tests){
                coords_level = mob_in$spat[as.character(group_data) == level, ]
                expl_S = rarefaction(comm_level, 'spat', coords = coords_level, 
                                     latlong = mob_in$latlong)
                sample_rare_level = cbind(sample_rare_level, expl_S)
            }
            sample_rare = rbind(sample_rare, sample_rare_level)
        }
        names(sample_rare)[1:3] = c('group', 'sample_plot', 'impl_S')
        sample_rare = df_factor_to_numeric(sample_rare, 2:ncol(sample_rare))
        if ('agg' %in% approved_tests){
            names(sample_rare)[4] = 'expl_S'
            sample_rare$deltaS_agg = sample_rare$expl_S - sample_rare$impl_S
        }
        return(sample_rare)
    }
}

#' Auxiliary function for get_delta_stats()
#' Effect of SAD when type is "continuous"
#' Directly add attributes to the input "out"
#' @keywords internal
effect_SAD_continuous = function(out, group_sad, env_levels, corr, n_perm){
    ind_sample_size = out$indiv_rare[, 1]
    ind_rare = out$indiv_rare[, -1]
    ind_cor = apply(ind_rare, 1, function(x) 
        cor(x, env_levels, method=corr))
    # Null test
    overall_sad_lumped = as.numeric(colSums(group_sad))
    sp_freq = overall_sad_lumped[overall_sad_lumped > 0] / sum(overall_sad_lumped)
    null_ind_r_mat = matrix(NA, n_perm, length(ind_sample_size))
    cat('\nComputing null model for SAD effect\n')
    pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
    for (i in 1:n_perm){
        # Note: necessary to convert sample to factor, so that zero counts are kept
        sad_perm = sapply(rowSums(group_sad), function(x)
            data.frame(table(factor(sample(1:length(sp_freq), x, replace = T,
                                           prob = sp_freq), 
                                    levels = 1:length(sp_freq))))[, 2])
        perm_ind_rare = apply(sad_perm, MARGIN = 2, function(x)
            rarefaction(x, 'indiv', ind_sample_size))
        null_ind_r_mat[i, ] = apply(perm_ind_rare, 1, function(x)
            cor(x, env_levels, method = corr))
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    ind_r_null_CI = apply(null_ind_r_mat, 2, function(x)
        quantile(x, c(0.025, 0.5, 0.975), na.rm = T)) # 95% CI
    out$continuous$SAD = data.frame(cbind(ind_sample_size, ind_cor, 
                                            t(ind_r_null_CI)))
    names(out$continuous$SAD) = c('effort_ind', 'r_emp', 'r_null_low', 
                                    'r_null_median', 'r_null_high')
    out$continuous$SAD = df_factor_to_numeric(out$continuous$SAD)
    return(out)
}

#' Auxiliary function for get_delta_stats()
#' Effect of SAD when type is "discrete"
#' @keywords internal
effect_SAD_discrete = function(out, group_sad, group_levels, ref_group, n_perm, 
                               overall_p){
    ind_sample_size = out$indiv_rare[, 1]
    ref_sad = group_sad[which(group_levels == as.character(ref_group)), ]
    out$SAD = data.frame(matrix(0, nrow = 0, ncol = 6), 
                                stringsAsFactors = F)

    cat('\nComputing null model for SAD effect\n')
    pb <- txtProgressBar(min = 0, max = n_perm * (length(group_levels) - 1), 
                         style = 3)
    k = 1
    for (level in group_levels[group_levels != ref_group]){
        deltaS = out$indiv_rare[, level] - out$indiv_rare[, as.character(ref_group)]
        level_sad = group_sad[which(group_levels == level), ]
        comp_sad_lumped = as.numeric(colSums(rbind(ref_sad, level_sad)))
        sp_freq = comp_sad_lumped[comp_sad_lumped > 0] / sum(comp_sad_lumped)
        
        null_ind_deltaS_mat = matrix(NA, n_perm, length(ind_sample_size))
        for (i in 1:n_perm){
            sad_perm = sapply(c(sum(level_sad), sum(ref_sad)), function(x)
                data.frame(table(factor(sample(1:length(sp_freq),x, replace = T,
                                               prob = sp_freq),
                                        levels = 1:length(sp_freq))))[, 2])
            perm_ind_rare = apply(sad_perm, MARGIN = 2, function(x)
                rarefaction(x, 'indiv', ind_sample_size))
            null_ind_deltaS_mat[i, ] = perm_ind_rare[, 1] - perm_ind_rare[, 2]
            setTxtProgressBar(pb, k)
            k = k + 1
        }
        ind_deltaS_null_CI = apply(null_ind_deltaS_mat, 2, function(x)
            quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
        ind_level = data.frame(cbind(rep(level,length(ind_sample_size)),
                                     ind_sample_size, deltaS, t(ind_deltaS_null_CI)))
        out$SAD = rbind(out$SAD, ind_level)
        if (overall_p){
            p_level = get_overall_p(out$indiv_rare[, 'sample'], 
                                    deltaS, null_ind_deltaS_mat)
            out$overall_p$SAD[out$overall_p$group == level] = p_level
        }
    }
    close(pb)
    out$SAD = df_factor_to_numeric(out$SAD, 2:ncol(out$SAD))
    names(out$SAD) = c('group', 'effort_ind', 'deltaS_emp',
                              'deltaS_null_low', 'deltaS_null_median',
                              'deltaS_null_high')
    return(out)
}

#' Auxiliary function for get_delta_stats()
#' Effect of N when type is "continuous"
#' @keywords internal
effect_N_continuous = function(mob_in, S, group_levels, env_levels, group_data, 
                               plot_dens, plot_abd, ind_sample_size, corr, 
                               n_perm){
    # TODO: checks?
    effect_N_by_group = data.frame(matrix(NA, ncol = length(group_levels) + 1,
                                          nrow = length(ind_sample_size)))
    effect_N_by_group[, 1] = ind_sample_size
    for (i in 1:length(group_levels)){
        level = group_levels[i]
        comm_level = mob_in$comm[which(as.character(group_data) == level), ]
        group_effect_N = deltaS_N(comm_level, plot_dens, ind_sample_size)
        effect_N_by_group[, i + 1] = group_effect_N$deltaS
    }
    effect_N_by_group = effect_N_by_group[complete.cases(effect_N_by_group), ]
    r_emp = apply(effect_N_by_group[ , -1], 1, function(x)
        cor(x, env_levels, method = corr))
    
    null_N_r_mat = matrix(NA, n_perm, length(r_emp))
    cat('\nComputing null model for N effect\n')
    pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
    for (i in 1:n_perm){
        comm_perm = permute_comm(mob_in$comm, 'swapN', group_data)
        effect_N_perm = data.frame(matrix(NA, ncol = length(group_levels),
                                          nrow = length(ind_sample_size)))
        for (j in 1:length(group_levels)){
            level_perm = group_levels[j]
            comm_level_perm = comm_perm[which(as.character(group_data) == level_perm), ]
            group_effect_N_perm = deltaS_N(comm_level_perm, plot_dens, 
                                           ind_sample_size[ind_sample_size <= sum(comm_level_perm)])
            # Ensure the column has the right length (filled with NA's if needed)
            effect_N_perm[, j] = group_effect_N_perm$deltaS[1:nrow(effect_N_perm)]
        }
        effect_N_perm = effect_N_perm[complete.cases(effect_N_perm), ]
        # If the output is not long enough, fill it with NA's
        null_N_r_mat[i, ] = apply(effect_N_perm, 1, function(x)
            cor(x, env_levels, method = corr))[1:ncol(null_N_r_mat)]
        setTxtProgressBar(pb, i)
    }
    close(pb)
    N_r_null_CI = apply(null_N_r_mat, 2, function(x) 
        quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
    out_N = data.frame(cbind(effect_N_by_group[, 1], r_emp, t(N_r_null_CI)))
    out_N = df_factor_to_numeric(out_N)
    names(out_N) = c('effort_ind', 'r_emp', 'r_null_low', 'r_null_median', 
                     'r_null_high')
    return(out_N)
}

#' Auxiliary function for get_delta_stats()
#' Effect of N when type is "discrete" 
#' @keywords internal
effect_N_discrete = function(out, mob_in, group_levels, ref_group, groups, 
                             density_stat, ind_sample_size, n_perm, overall_p){
    out_N = NULL
    if (overall_p)
        overallp_N = data.frame(matrix(0, nrow = 0, ncol = 2), 
                                stringsAsFactors = F)
    cat('\nComputing null model for N effect\n')
    pb <- txtProgressBar(min = 0, max = n_perm * (length(group_levels) - 1), 
                         style = 3)
    k = 1
    for (level in group_levels[group_levels != as.character(ref_group)]){
        row_bool = level == groups | as.character(ref_group) == groups
        comm_levels = mob_in$comm[row_bool, ]
        plot_dens_level = get_plot_dens(comm_levels, density_stat)
        plot_levels = groups[row_bool]
        N_eff = sapply(c(as.character(ref_group), level), function(x)
            deltaS_N(comm_levels[plot_levels == x, ], plot_dens_level, 
                     ind_sample_size)$deltaS)
        ddeltaS_level = N_eff[, 2] - N_eff[, 1]
        null_N_deltaS_mat = matrix(NA, n_perm, length(ddeltaS_level))
        for (i in 1:n_perm){
            # swap plot abu between group 1 and each other group
            comm_perm = permute_comm(comm_levels, 'swapN', plot_levels)  
            min_N = min(sum(comm_perm[plot_levels == as.character(ref_group), ]), 
                        sum(comm_perm[plot_levels == level, ]))
            N_eff_perm = sapply(c(as.character(ref_group), level), function(x) 
                deltaS_N(comm_perm[plot_levels == x, ], plot_dens_level, 
                         ind_sample_size[ind_sample_size <= min_N])$deltaS)
            ddeltaS_perm = N_eff_perm[ , 2] - N_eff_perm[ , 1]
            # Ensure the row has the right length (filled with NA's if needed)
            null_N_deltaS_mat[i, ] = ddeltaS_perm[1:ncol(null_N_deltaS_mat)]
            setTxtProgressBar(pb, k)
            k = k + 1
        }
        N_deltaS_null_CI = apply(null_N_deltaS_mat, 2, function(x)
            quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
        N_level = data.frame(level, ind_sample_size, ddeltaS_level,
                             t(N_deltaS_null_CI))
        out_N = rbind(out_N, N_level)
        if (overall_p){
            p_level = get_overall_p(ind_sample_size, ddeltaS_level, 
                                    null_N_deltaS_mat)
            out$overall_p$N[out$overall_p$group == level] = p_level
        }
    }
    close(pb)
    out_N = df_factor_to_numeric(out_N, 2:ncol(out_N))
    names(out_N) = c('group', 'effort_sample', 'ddeltaS_emp', 'ddeltaS_null_low', 
                     'ddeltaS_null_median', 'ddeltaS_null_high')
    out$N = out_N
    return(out)
}

#' Auxiliary function for get_delta_stats()
#' Effect of aggregation when type is "continuous"
#' @keywords internal
effect_agg_continuous = function(mob_in, sample_rare, group_plots, group_levels, 
                                 group_data, env_levels, corr, n_perm){
    min_plot_level = min(group_plots$Freq)
    r_emp = c()
    for (iplot in seq(min_plot_level)){
        deltaS_i = sample_rare$deltaS_agg[sample_rare$sample_plot == iplot]
        groups_i = as.character(sample_rare$group[sample_rare$sample_plot 
                                                  == iplot])
        env_i = env_levels[which(group_levels == groups_i)]
        r_emp = c(r_emp, cor(deltaS_i, env_i, method = corr))
    }

    null_agg_r_mat = matrix(NA, n_perm, min_plot_level)
    cat('\nComputing null model for aggregation effect\n')
    pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
    for (i in 1:n_perm){
        comm_perm = comm
        comm_perm$comm = permute_comm(mob_in$comm, 'noagg', group_data)
        sample_rare_perm = get_sample_curves(comm_perm, group_levels, group_data, 
                                             c('N', 'agg'))
        r_perm = c()
        for (iplot in seq(min_plot_level)){
            deltaS_i = sample_rare_perm$deltaS_agg[sample_rare_perm$sample_plot 
                                                   == iplot]
            groups_i = as.character(sample_rare_perm$group[sample_rare_perm$sample_plot 
                                                      == iplot])
            env_i = env_levels[which(group_levels == groups_i)]
            r_perm = c(r_perm, cor(deltaS_i, env_i, method = corr))
        }
        null_agg_r_mat[i, ] = r_perm
        setTxtProgressBar(pb, i)
    }
    close(pb)
    plot_levels = which(!is.na(r_emp))
    null_agg_r_mat = null_agg_r_mat[, plot_levels]
    agg_r_null_CI = apply(null_agg_r_mat, 2, function(x) 
        quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
    table_agg = data.frame(cbind(plot_levels, r_emp[plot_levels], 
                                 t(agg_r_null_CI)))
    names(table_agg) = c('effort_sample', 'r_emp', 'r_null_low', 
                         'r_null_median', 'r_null_high')
    return(table_agg)
}

#' Auxiliary function for get_delta_stats()
#' Effect of aggregation when type is "discrete"
#' @keywords internal
effect_agg_discrete = function(out, mob_in, ref_group, group_plots, group_data, 
                               group_levels, n_perm, overall_p){
    sample_rare = out$sample_rare
    ref_sample = sample_rare[which(sample_rare$group == 
                                       as.character(ref_group)), ]
    table_agg = data.frame(matrix(NA, nrow = 0, ncol = 6))
    cat('\nComputing null model for aggregation effect\n')
    pb <- txtProgressBar(min = 0, max = n_perm * (length(group_levels) - 1), 
                         style = 3)
    k = 1
    for (level in group_levels[group_levels != as.character(ref_group)]){
        min_plot_level = min(group_plots$Freq[which(group_plots$groups %in%
                                                        c(ref_group, level))])
        ddeltaS_level = sample_rare$deltaS_agg[sample_rare$group == level][1:min_plot_level] - 
            sample_rare$deltaS_agg[sample_rare$group == ref_group][1:min_plot_level]
        
        null_agg_deltaS_mat = matrix(NA, n_perm, min_plot_level)
        for (i in 1:n_perm){
            comm_perm = mob_in
            comm_perm$comm = permute_comm(mob_in$comm, 'noagg', group_data)
            sample_rare_perm = get_sample_curves(comm_perm, group_levels, group_data, 
                                                 c('N', 'agg'))
            ddeltaS_perm = sample_rare_perm$deltaS_agg[sample_rare_perm$group == level][1:min_plot_level] - 
                sample_rare_perm$deltaS_agg[sample_rare_perm$group == ref_group][1:min_plot_level]
            null_agg_deltaS_mat[i, ] = ddeltaS_perm
            setTxtProgressBar(pb, k)
            k = k + 1
        }
        agg_deltaS_null_CI = apply(null_agg_deltaS_mat, 2, function(x) 
            quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
        agg_level = data.frame(cbind(rep(level, min_plot_level), 1:min_plot_level, 
                                     ddeltaS_level, t(agg_deltaS_null_CI)))
        table_agg = rbind(table_agg, agg_level)
        if (overall_p){
            p_level = get_overall_p(1:min_plot_level, ddeltaS_level, 
                                    null_agg_deltaS_mat)
            out$overall_p$agg[out$overall_p$group == level] = p_level
        }
    }
    close(pb)
    table_agg = df_factor_to_numeric(table_agg, 2:ncol(table_agg))
    names(table_agg) = c('group', 'effort_sample', 'ddeltaS_emp', 'ddeltaS_null_low', 
                         'ddeltaS_null_median', 'ddeltaS_null_high')
    out$agg = table_agg
    return(out)
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

get_delta_stats = function(mob_in, group_var, env_var = NULL, ref_group = NULL, 
                           tests = c('SAD', 'N', 'agg'),
                           type='discrete', inds = NULL, log_scale = FALSE,
                           min_plots = NULL, density_stat ='mean',
                           corr='spearman', n_perm=1000, overall_p = FALSE) {
    
    approved_tests = get_delta_overall_checks(mob_in, type, group_var, env_var, 
                                              density_stat, tests)
    
    S = ncol(mob_in$comm)
    plot_abd = rowSums(mob_in$comm)
    group_data = mob_in$env[, group_var]
    groups = as.character(group_data)
    group_plots = data.frame(table(groups)) # Number of plots within each group
    
    group_sad = aggregate(mob_in$comm, by=list(group_data), sum)
    # Distinguish between group_levels, the grouping factor, and env_levels, 
    # the levels of the environmental factor for each group used for correlation
    # Make sure that the orders match!
    if (is.null(env_var)){
        env_raw = group_sad[, 1]
        if ('factor' %in% class(env_raw)){
            env_levels = as.numeric(env_raw)
        } else {
            env_levels = env_raw
        }
    } else {
        env_levels = tapply(mob_in$env[, env_var], list(group_data), mean)
    }
    group_levels = as.character(group_sad[, 1])
    group_sad = group_sad[, -1]
    ind_sample_size = get_delta_ind_sample(group_sad, inds, log_scale)
    plot_dens = get_plot_dens(mob_in$comm, density_stat)

    out = list()
    out$type = type
    out$tests = approved_tests
    out$log_scale = log_scale
    out$density_stat = list(density_stat = density_stat, plot_dens = plot_dens)
    ind_rare = data.frame(apply(group_sad, 1, function(x) 
        rarefaction(x, 'indiv', ind_sample_size)))
    out$indiv_rare = cbind(ind_sample_size, ind_rare)
    names(out$indiv_rare) = c('sample', group_levels)
    out$sample_rare = get_sample_curves(mob_in, group_levels, group_data, 
                                        approved_tests)
    
    if (type == 'continuous'){
        get_delta_continuous_checks(corr, group_levels, env_raw)
        if ('SAD' %in% approved_tests)
            out = effect_SAD_continuous(out, group_sad, env_levels, corr, n_perm)
        if ('N' %in% approved_tests)
            out = effect_N_continuous(out, mob_in, S, group_levels, env_levels, 
                                        group_data, plot_dens, plot_abd, 
                                        ind_sample_size, corr, n_perm)
        if ('agg' %in% approved_tests)
            out$agg = effect_agg_continuous(mob_in, out$sample_rare,
                                            group_plots, group_levels, 
                                            group_data, env_levels, corr, n_perm)
    } else if (type == 'discrete') {
        get_delta_discrete_checks(ref_group, group_levels, group_data, env_var)
        if (overall_p) {
            cat('Caution: Overall p-values depend on scales of measurement yet do not explicitly reflect significance at any particular scale. Be careful in interpretation.')
            out$overall_p = as.data.frame(matrix(NA, length(group_levels) - 1, 
                                                 1 + length(approved_tests)))
            names(out$overall_p) = c('group', approved_tests)
            out$overall_p$group = group_levels[group_levels != ref_group]
        }
        if ('SAD' %in% approved_tests)
            out = effect_SAD_discrete(out, group_sad, group_levels, ref_group,
                                      n_perm, overall_p)
        if ('N' %in% approved_tests)
            out = effect_N_discrete(out, mob_in, group_levels, ref_group, groups,
                                    density_stat,ind_sample_size, n_perm, overall_p)
        if ('agg' %in% approved_tests)
            out = effect_agg_discrete(out, mob_in, ref_group, group_plots, 
                                      group_data, group_levels, n_perm, overall_p)
    } else 
        stop('The argument "type" must be either "discrete" or "continuous"')
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
        plot_dens = mob_out$density_stat$plot_dens
        # the next line assumes that the density stat is average (may need to generalize)
        effort = min(agg$effort) : round(N_indiv / N_plots)
        virt_effort = effort * plot_dens
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
    plot_dens = mean(N)
    N_sum = apply(replicate(n_perm, cumsum(sample(N))), 1, mean)
    plot(N_sum, xlab='Number of plots', ylab='Number of Individuals')
    abline(a=0, b=plot_dens, col='red')
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

