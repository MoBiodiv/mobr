#' Create the 'mob_in' object.
#'
#' The 'mob_in' object will be passed on for analyses of biodiversity across
#' scales.
#'
#' @param comm community matrix in which rows are samples (e.g., plots) and
#'   columns are species.
#' @param plot_attr matrix which includes the environmental attributes and
#'   spatial coordinates of the plots. Environmental attributes are mandatory,
#'   while spatial coordinates are optional.
#' @param coord_names character vector with the names of the columns of
#'   \code{plot_attr} that specify the coordinates of the samples. Defaults to
#'   NULL (no coordinates). When providing coordinate names, the order the names
#'   are provided matters when working with latitude-longitude coordinates
#'   (i.e., argument \code{latlong = TRUE}, and it is expected that the column
#'   specifying the x-coordinate or the longitude is provided first, y-coordinate
#'    or latitude provided second. To provide coordinate names use the following
#'    syntax: \code{coord_names = c('longitude_col_name','latitude_col_name')}
#' @param binary Boolean, defaults to FALSE. Whether the plot by species matrix
#'   "comm" is in abundances or presence/absence.
#' @param latlong Boolean, defaults to FALSE. Whether the coordinates are
#'   latitude-longitudes. If TRUE, distance calculations by downstream functions
#'   are based upon great circle distances rather than Euclidean distances. Note
#'   latitude-longitudes should be in decimal degree.
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
#'  inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
make_mob_in = function(comm, plot_attr, coord_names = NULL, binary = FALSE,
                       latlong = FALSE) {
    # possibly make group_var and ref_group mandatory arguments
    out = list(tests = list(N = TRUE, SAD = TRUE, agg = TRUE))
    # carry out some basic checks
    if (nrow(comm) < 5) {
        warning("Number of plots in community is less than five therefore only individual rarefaction will be computed")
        out$tests$N = FALSE
        out$tests$agg = FALSE
    }

    if (nrow(comm) != nrow(plot_attr))
        stop("Number of plots in community does not equal number of plots in plot attribute table")

    if (is.null(coord_names) == FALSE) {
        spat_cols = sapply(coord_names, function(x) which(x == names(plot_attr)))

        if (length(spat_cols) == 1 & latlong == TRUE)
            stop("Both latitude and longitude have to be specified")
    }

    if (any(row.names(comm) != row.names(plot_attr)))
        warning("Row names of community and plot attributes tables do not match
                which may indicate different identities or orderings of samples")

    if (binary)  {
        warning("Only spatially-explicit sampled based forms of rarefaction can be computed on binary data")
        out$tests$SAD = FALSE
        out$tests$N = FALSE
    }
    else {
        if (max(comm) == 1)
            warning("Maximum abundance is 1 which suggests data is binary, change the binary argument to TRUE")
    }

    if (any(colSums(comm) == 0)) {
        warning("Some species have zero occurrences and will be dropped from the community table")
        comm = comm[ , colSums(comm) != 0]
    }

    out$comm = data.frame(comm)
    if (is.null(coord_names) == FALSE) {
        if (length(spat_cols) > 0) {
            out$env = data.frame(plot_attr[ , -spat_cols])
            colnames(out$env) = colnames(plot_attr)[-spat_cols]
            out$spat = data.frame(plot_attr[ , spat_cols])
       }
    }
    else {
        warning("Note: 'coord_names' was not supplied and therefore spatial aggregation will not be examined in downstream analyses")
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
#' @param x an object of class mob_in created by \code{\link{make_mob_in}}
#' @param subset expression indicating elements or rows to keep: missing values are taken as false.
#' @param type specifies the type of object the argument \code{subset}
#'   specifies, may be: \code{string}, \code{integer}, or \code{logical},
#'   defaults to \code{string}
#' @param drop_levels Boolean if TRUE unused levels are removed from factors in
#'   mob_in$env
#' @param ... parameters passed to other functions
#' @export
#' @examples
#'  data(inv_comm)
#'  data(inv_plot_attr)
#'  inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#'  subset(inv_mob_in, group == 'invaded')
#'  subset(inv_mob_in, 1:4, type='integer')
#'  subset(inv_mob_in, 1:4, type='integer', drop_levels=TRUE)
#'  sub_log = c(TRUE, FALSE, TRUE, rep(FALSE, nrow(inv_mob_in$comm) - 3))
#'  subset(inv_mob_in, sub_log, type='logical')
subset.mob_in = function(x, subset, type = 'string', drop_levels = FALSE, ...) {
    if (missing(subset))
        r <- rep_len(TRUE, nrow(x$comm))
    if (type == 'integer')
        r <- 1:nrow(x$comm) %in% subset
    if (type == 'logical')
        r <- subset
    if (type == 'string') {
        e <- substitute(subset)
        r <- eval(e, x$env)
        if (!is.logical(r))
            stop("'subset' must be logical when type = 'string'")
   }
   x$comm = base::subset(x$comm, r)
   x$env = base::subset(x$env, r)
   if (drop_levels)
       x$env = droplevels(x$env)
   if (!is.null(x$spat))
       x$spat = x$spat[r, ]
   return(x)
}

#' Print a shortened version of the mob_in object
#' @param x a mob_in class object
#' @param nrows the number of rows of each matrix to print
#' @param nsp the number of species columns to print
#' @param ... parameters passed to other functions
#' @importFrom utils head
#' @keywords internal
#' @rdname print.mob_in
#' @export
print.mob_in = function(x, nrows = 6, nsp = 5, ...) {
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
    print(utils::head(x$env, nrows))
    cat('\n$spat\n')
    print(head(x$spat, nrows))
    cat('\n$latlong\n')
    print(x$latlong)
}

#' Internal function for distance matrix assuming inputs are longitude and
#' latitudes on a spherical Earth.
#'
#' @param coords a matrix with longitudes and latitudes in decimal degrees. The
#'   longitudes should be provided in the first column (they are the
#'   x-coordinate) and the latitudes should be provided in the second column
#'   (they are the y-coordinate).
#'
#' @param r the radius of the Earth, defaults to 6378.137 km
#'
#' @returns distance matrix between all pairwise coordinates.
#'
#' @description This calculation uses the Haversine method of computing great
#'   circle distances in kilometers on a spherical Earth (r = 6378.137 km). This
#'   code was taken from fields::rdist.earth by Doug Nychka, John Paige, Florian
#'   Gerber.
#' @keywords internal
sphere_dist = function(coords, r = 6378.137){
    coslat1 <- cos((coords[ , 2] * pi) / 180)
    sinlat1 <- sin((coords[ , 2] * pi) / 180)
    coslon1 <- cos((coords[ , 1] * pi) / 180)
    sinlon1 <- sin((coords[ , 1] * pi) / 180)
    pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
          t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
    return(r * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
}

#' Rarefied Species Richness
#' 
#' The expected number of species given a particular number of individuals or
#' samples under random and spatially explicit nearest neighbor sampling
#' 
#' 
#' @param x can either be a: 1) mob_in object, 2) community matrix-like
#'  object in which rows represent plots and columns represent species, or 3)
#'  a vector which contains the abundance of each species. 
#' @param method a character string that specifies the method of rarefaction 
#'   curve construction it can be one of the following: 
#' \itemize{
#'     \item \code{'IBR'} ... individual-based rarefaction in which species
#'     are accumulated by randomly sampling individuals
#'     \item \code{'SBR'} ... sample-based rarefaction in which species are 
#'     accumulated by randomly sampling samples (i.e., plots). Note that within plot spatial 
#'     aggregation is maintained with this approach. Although this curve
#'     is implemented here, it is not used in the current version of the MoB framework
#'     \item \code{'nsSBR'} ... non-spatial, sampled-based rarefaction in which
#'     species are accumulated by randomly sampling samples that represent a 
#'     spatially random sample of individuals (i.e., no with-in plot spatial 
#'     aggregation). The argument \code{dens_ratio} must also be set otherwise 
#'     this sampling results in a curve identical to the IBR (see Details). 
#'     \item \code{'sSBR'} ... spatial sample-based rarefaction in which species 
#'     are accumulated by including spatially proximate samples first. 
#' }
#' @param effort optional argument to specify what number of individuals or 
#'   number of samples depending on 'method' to compute rarefied richness as. If
#'   not specified all possible values from 1 to the maximum sampling effort are
#'   used
#' @param coords an optional matrix of geographic coordinates of the samples.  
#'   Only required when using the spatial rarefaction method and this information
#'   is not already supplied by \code{x}. The first column should specify 
#'   the x-coordinate (e.g., longitude) and the second coordinate should 
#'   specify the y-coordinate (e.g., latitude)
#' @param latlong Boolean if coordinates are latitude-longitude decimal degrees
#' @param dens_ratio the ratio of individual density between a reference group
#'   and the community data (i.e., x) under consideration. This argument is
#'   used to rescale the rarefaction curve when estimating the effect of
#'   individual density on group differences in richness.
#' @param extrapolate Boolean which specifies if richness should be extrapolated
#'   when effort is larger than the number of individuals using the chao1 method.
#'   Defaults to FALSE in which case it returns observed richness. Extrapolation
#'   is only implemented for individual-based rarefaction 
#'   (i.e., \code{method = 'indiv'})
#' @param return_NA Boolean defaults to FALSE in which the function returns the
#'   observed S when \code{effort} is larger than the number of individuals or
#'   number of samples (depending on the method of rarefaction). If set to TRUE
#'   then NA is returned. Note that this argument is only relevant when
#'   \code{extrapolate = FALSE}.
#' @param quiet_mode Boolean defaults to FALSE, if TRUE then warnings and other
#'   non-error messages are suppressed.
#' @param spat_algo character string that can be either: \code{'kNN'} or \code{'kNCN'}
#' for k-nearest neighbor and k-nearest centroid neighbor sampling 
#' respectively. It defaults to k-nearest neighbor which is a 
#' more computationally efficient algorithm that closely approximates the 
#' potentially more correct k-NCN algo (see Details). 
#'   
#' @details The analytical formulas of Cayuela et al. (2015) are used to compute
#'   the random sampling expectation for the individual and sampled based
#'   rarefaction methods. The spatially constrained rarefaction curve (Chiarucci
#'   et al. 2009) also known as the sample-based accumulation curve (Gotelli and
#'   Colwell 2001) can be computed in one of two ways which is determined by the
#'   argument \code{spat_algo}. In the kNN approach each plot is accumulated by
#'   the order of their spatial proximity to the original focal cell. If plots
#'   have the same distance from the focal plot then one is chosen randomly to
#'   be sampled first. In the kNCN approach, a new centroid is computed after
#'   each plot is accumulated, then distances are recomputed from that new
#'   centroid to all other plots and the next nearest is sampled. The kNN is
#'   faster because the distance matrix only needs to be computed once, but the
#'   sampling of kNCN which simultaneously minimizes spatial distance and extent
#'   is more similar to an actual person searching a field for species. For both
#'   kNN and kNCN, each plot in the community matrix is treated as a starting
#'   point and then the mean of these n possible accumulation curves is
#'   computed.
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
#' Letters, 4, 379-391.
#'
#' @importFrom stats dist
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' sad = colSums(inv_comm)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' # rarefaction can be performed on different data inputs
#' # all three give same answer
#' # 1) the raw community site-by-species matrix
#' rarefaction(inv_comm, method='IBR', effort=1:10)
#' # 2) the SAD of the community
#' rarefaction(inv_comm, method='IBR', effort=1:10)
#' # 3) a mob_in class object
#' # rescaling of individual based rarefaction 
#' # when the density ratio is 1 the richness values are 
#' # identical to unscale rarefaction
#' rarefaction(inv_comm, method='IBR', effort=1:10, dens_ratio=1)
#' # however the curve is either shrunk when density is higher than 
#' # the reference value (i.e., dens_ratio < 1)
#' rarefaction(inv_comm, method='IBR', effort=1:10, dens_ratio=0.5)
#' # the curve is stretched when density is lower than the 
#' # reference value (i.e., dens_ratio > 1)
#' rarefaction(inv_comm, method='IBR', effort=1:10, dens_ratio=1.5)
#' # sample based rarefaction under random sampling
#' rarefaction(inv_comm, method='SBR')
#' \donttest{ 
#' # sampled based rarefaction under spatially explicit nearest neighbor sampling
#' rarefaction(inv_comm, method='sSBR', coords=inv_plot_attr[ , c('x','y')],
#'             latlong=FALSE)
#' # the syntax is simpler if supplying a mob_in object
#' rarefaction(inv_mob_in, method='sSBR', spat_algo = 'kNCN')
#' rarefaction(inv_mob_in, method='sSBR', spat_algo = 'kNN')
#' }
rarefaction = function(x, method, effort=NULL, coords=NULL, latlong=NULL, 
                       dens_ratio=1, extrapolate=FALSE, return_NA=FALSE, 
                       quiet_mode=FALSE, spat_algo=NULL) {
    
    if (method == 'indiv') {
        warning('method == "indiv" is depreciated and should be set to "IBR" for individual-based rarefaction')
        method = 'IBR'
    } else if (method == 'samp') {
        warning('method == "samp" is depreciated and should be set to "SBR" for sample-based rarefaction')
        method = 'SBR'
    } else if (method == 'spat') {
        warning('method == "spat" is depreciated and should be set to "sSBR" for spatial, sample-based rarefaction')
        method = 'sSBR'
    } else if (!any(method %in% c('IBR', 'SBR', 'nsSBR', 'sSBR')))
        stop('The argument "method" must be set to either "IBR", "SBR", "nsSBR",',
             ' or "sSBR" for random individual, random sample, non-spatial,', 
             ' sample-based (nsSBR), and spatial, sample-based rarefaction (sSBR)',
             ' respectively.')
    if (method == 'nsSBR' & dens_ratio == 1)
        warning('The nonspatial, sample-based rarefaction (nsSBR) curve only differs from the IBR when compared with a reference density by setting "dens_ratio" not equal to 1')
    if ('mob_in' %in% class(x)) {
        x_mob_in = x
        x = x_mob_in$comm
        if (is.null(latlong))
            latlong = x_mob_in$latlong
        else if (latlong != x_mob_in$latlong)
            stop(paste('The "latlong" argument is set to', latlong, 
                       'but the value of x$latlong is', x_mob_in$latlong))
        if (method == 'sSBR') {
            if (is.null(coords)) {
                if (is.null(x_mob_in$spat)) {
                    stop('Coordinate name value(s) must be supplied in the make_mob_in object in order to plot using sample spatially explicit based (spat) rarefaction')
                }
                coords = x_mob_in$spat
            }
        }
    }
    if (method == 'SBR' | method == 'sSBR') {
        if (is.null(dim(x)))
            stop('For random or spatially explicit sample based rarefaction "x" must be a site x species matrix as the input')
        else {
            x = (x > 0) * 1             
            # all sites are counted as samples even empty ones
            n = nrow(x) 
            if (method == 'SBR')
                x = colSums(x)
        }
    } else if (!is.null(spat_algo))
        warning("Setting spat_algo to a non-NULL value only has consequences when method = sSBR")
    if (method == 'IBR' | method == 'nsSBR') {
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
                                 ifelse(method == 'IBR', 'individuals', 'samples'),
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
    if (method == 'sSBR') {
        if (is.null(spat_algo)) 
            spat_algo = 'kNN'
        if (spat_algo == 'kNN') {
      
            explicit_loop = matrix(0, n, n)
            if (is.null(latlong))
                stop('For spatial rarefaction the argument "latlong" must be set TRUE or FALSE')
            if (latlong) {
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
        else if (spat_algo == "kNCN") 
            out = kNCN_average(x=x, coords=coords, latlong=latlong)[effort]
    } 
    else { 
        # drop species with no observations  
        x = x[x > 0] 
        S = length(x)
        if (dens_ratio == 1) {
            ldiv = lchoose(n, effort)
        } else {
            effort = effort / dens_ratio
            ldiv = lgamma(n - effort + 1) - lgamma(n + 1)
        }
        p = matrix(0, sum(effort <= n), S)
        S_ext = NULL
        for (i in seq_along(effort)) {
            if (effort[i] <= n) {
                if (dens_ratio == 1) {
                    p[i, ] = ifelse(n - x < effort[i], 0, 
                                    exp(lchoose(n - x, effort[i]) - ldiv[i]))
                } else {
                    p[i, ] = ifelse(n - x < effort[i], 0, 
                                    exp(suppressWarnings(lgamma(n - x + 1)) -
                                        suppressWarnings(lgamma(n - x - effort[i] + 1)) +
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
#' @keywords internal
ind_rare_perm = function(abu, n_perm = 100, n_indiv = NULL) {
    if (!is.vector(abu)) {
        stop('abu must be a vector of abundances')
    } 
    calc_S = function(splist, n_indiv) {
        sapply(n_indiv, function(n) length(unique(splist[1:n])))
    }
    rand_splist = function(abu, S) {
        sample(unlist(mapply(rep, 1:S, abu)), replace = FALSE)
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
#' when building the spatial, sample-based rarefaction curves (sSBR).
#' 
#' @param coords a matrix with n-dimensional coordinates
#' @return a vector of average distances for each sequential number
#'   of accumulated nearest samples. 
#' @export
#' @examples 
#' # transect spatial arrangement
#' transect = 1:100
#' avg_nn_dist(transect)
#' grid = expand.grid(1:10, 1:10)
#' avg_nn_dist(grid)
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow=c(1,2)) 
#' plot(avg_nn_dist(transect), type='o', main='transect',
#'      xlab='# of samples', ylab='average distance')
#' # 2-D grid spatial arrangement
#' plot(avg_nn_dist(grid), type='o', main='grid',
#'      xlab='# of samples', ylab='average distance')
#' par(oldpar)
avg_nn_dist = function(coords) {
    pair_dist = as.matrix(stats::dist(coords))
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
#' @inheritParams rarefaction
#' @importFrom tibble tibble
#' @keywords internal
get_delta_curves = function(x, tests=c('SAD', 'N', 'agg'), spat_algo=NULL,
                            inds=NULL, ind_dens=NULL, n_plots=NULL) {
    if (is.null(inds) & any(c('SAD', 'N') %in% tests))
        stop('If SAD or N effect to be calculated inds must be specified')
    if (is.null(ind_dens) & 'N' %in% tests)
        stop('If N effect to be calculated ind_dens must be specified')
    if (any(c('N', 'agg') %in% tests) & !('mob_in' %in% class(x)))
        stop('If N or agg effects to be computed x must be a mob_in object')
    out = list()
    if ('SAD' %in% tests) {
        S_SAD = rarefaction(x, 'IBR', inds)
        out$SAD = data.frame(test = 'SAD', sample = 'indiv',
                             effort = inds, S = S_SAD, effect = S_SAD,
                             stringsAsFactors = FALSE)
    }
    if ('N' %in% tests) {
        comm_dens = sum(x$comm) / nrow(x$comm)
        dens_ratio = ind_dens / comm_dens
        S_N = rarefaction(x, 'IBR', inds, dens_ratio = dens_ratio)
        if (!('SAD' %in% tests))
            S_SAD = rarefaction(x, 'IBR', inds)
        effect = S_N - S_SAD
        out$N = data.frame(test = 'N', sample = 'indiv', 
                           effort = inds, S = S_N, effect,
                           stringsAsFactors = FALSE)
    }
    if ('agg' %in% tests) {
        if (is.null(n_plots))
            n_plots = nrow(x$comm)
        S_agg = rarefaction(x, 'sSBR', 1:n_plots, spat_algo = spat_algo)
        ind_density = sum(x$comm) / nrow(x$comm)
        samp_effort = round(1:n_plots * ind_density)
        S_N = rarefaction(x, 'IBR', samp_effort)
        effect = S_agg - S_N
        out$agg = data.frame(test = 'agg', sample = 'plot', 
                             effort = as.numeric(names(S_agg)),
                             S = S_agg, effect, 
                             stringsAsFactors = FALSE)
    }
    return(flatten_dfr(tibble(out)))
}
        

#' Randomly sample of a relative abundance distribution (RAD)
#' to produce an expected species abundance distribution (SAD)
#' 
#' @param rad the relative abundance of each species
#' @param N the total number of individuals sampled
#' 
#' Randomly subsampling an RAD with replacement produces an SAD that is of a
#' similar functional form (Green and Plotkin 2007) but with overall species
#' richness equal to or less than the relative abundance distribution.
#' 
#' Literature Cited:
#' Green, J. L., and J. B. Plotkin. 2007. A statistical theory
#'  for sampling species abundances. Ecology Letters 10:1037-1045.
#' 
#' @keywords internal
get_rand_sad = function(rad, N) {
  rand_samp = sample(1:length(rad), N, replace = T, prob = rad)
  rand_sad = table(factor(rand_samp, levels = 1:length(rad)))
  return(as.numeric(rand_sad))
}

#' Generate a null community matrix 
#' 
#' Three  null models are implemented that randomize different components of
#' community structure while keeping other components constant.
#' 
#' 
#' @param comm community matrix of abundances with plots as rows and species columns.
#' @param null_model a string which specifies which null model to use options
#'   include: \code{'rand_SAD'}, \code{'rand_N'}, and \code{'rand_agg'}. See 
#'   Details for description of each null model.
#' @param groups optional argument that is a vector of group ids which specify
#'   which group each site is associated with. If is \code{NULL} then all rows
#'   of the community matrix are assumed to be members of the same group
#'   
#' @return a site-by-species matrix
#' 
#' @details 
#' This function implements three different nested null models. They are considered
#' nested because at the core of each null model is the random sampling 
#' with replacement of the relative abundance distribution (RAD) to generate 
#' a random sample of a species abundance distribution (SAD). Here we describe
#' each null model:
#' \itemize{
#'    \item \code{'rand_SAD'} ... A random SAD is generated using a sample with
#'    replacement of individuals from the species pool proportional to their
#'    observed relative abundance. This null model will produce an SAD that is
#'    of a similar functional form to the observed SAD (Green and Plotkin 2007).
#'    The total abundance of the random SAD is the same as the observed SAD but
#'    overall species richness will be equal to or less than the observed SAD.
#'    This algorithm ignores the \code{group} argument. This sampling algorithm
#'    is also used in the two other null models \code{'rand_N'} and
#'    \code{'rand_agg'}.
#'    
#'    \item \code{'rand_N'} ... The total number of individuals in a plot is
#'    shuffled across all plots (within and between groups). Then for each plot
#'    that many individuals are drawn randomly from the group specific relative
#'    abundance distribution with replacement for each plot (i.e., using the
#'    \code{'rand_SAD'} algorithm described above. This removes group
#'    differences in the total number of individuals in a given plot, but
#'    maintains group level differences in their SADs.
#'    
#'    \item \code{'rand_agg'} ... This null model nullifies the spatial
#'    structure of individuals (i.e., their aggregation), but it is constrained
#'    by the observed total number of individuals in each plot (in contrast to
#'    the \code{'rand_N'} null model), and the group specific SAD (in contrast
#'    to the \code{'rand_SAD'} null model). The other two null models also
#'    nullify spatial structure. The \code{'rand_agg'} null model is identical
#'    to the \code{'rand_N'} null model except that plot abundances are not 
#'    shuffled. 
#' }
#' 
#' Replaces depreciated function `permute_comm`
#' 
#' @references 
#' Green, J. L., and J. B. Plotkin. 2007. A statistical theory for sampling 
#' species abundances. Ecology Letters 10:1037-1045.
#' @importFrom vctrs vec_as_names
#' @import purrr
#' @import dplyr
#' @export
#' @examples 
#' S = 3
#' N = 20
#' nplots = 4
#' comm = matrix(rpois(S * nplots, 1), ncol = S, nrow = nplots)
#' comm
#' groups = rep(1:2, each=2)
#' groups
#' set.seed(1)
#' get_null_comm(comm, 'rand_SAD')
#' # null model 'rand_SAD' ignores groups argument
#' set.seed(1)
#' get_null_comm(comm, 'rand_SAD', groups)
#' set.seed(1)
#' get_null_comm(comm, 'rand_N')
#' # null model 'rand_N' does not ignore the groups argument
#' set.seed(1)
#' get_null_comm(comm, 'rand_N', groups)
#' # note that the 'rand_agg' null model is constrained by observed plot abundances
#' noagg = get_null_comm(comm, 'rand_agg', groups)
#' noagg
#' rowSums(comm)
#' rowSums(noagg)
get_null_comm = function(comm, null_model, groups = NULL) {
    # the main component of all the null models is random sampling 
    # from a pooled or group-specific SAD
    if (!(is.matrix(comm) | is.data.frame(comm)))
        stop('comm must be a matrix or data.frame')
    if (is.null(groups))
        groups = rep(1, nrow(comm))   
  
    # compute N at each plot across groups
    N_plots = rowSums(comm)
    if (null_model == "rand_SAD") {
        # NOTE: for SAD test it is not necessary to fix or not fix plot level
        # abundance because individual based rarefaction ignores that
        # detail when it is computed
        # compute relative abundance distribution of the species pool
        rad_pool = colSums(comm) / sum(comm)
        # compute average abundance per group 
        # randomly sample N individuals from the pool with replacement
        null_sads = map(N_plots, ~ get_rand_sad(rad_pool, .x))
        names(null_sads) = 1:length(null_sads)
    } else if (null_model == "rand_N" | null_model == "rand_agg") {
        if (null_model == "rand_N") # shuffle these abundances in the N null model
            N_plots = sample(N_plots)
        # compute rad for each group
        .x <- NULL   # book keeping for CRAN checks 
        rad_groups = data.frame(comm, groups) %>%
                     group_by(groups) %>%
                     summarize_all(sum) %>%
                     select(-one_of("groups")) %>% t %>%
                     as_tibble(.x, .name_repair = ~ vctrs::vec_as_names(..., quiet = TRUE)) %>%
                     map(~ .x / sum(.x))
        # replicate these rads so that you have one group specific rad for every
        # plot in the dataset
        rad_plots = rep(rad_groups, table(groups))
        names(rad_plots) = 1:length(rad_plots)
        # draw random sads from each group specific rad
        null_sads = map2(rad_plots, N_plots, get_rand_sad)
    }  
    # now convert sads to a new community matrix
    null_comm = null_sads %>% tibble %>% flatten_dfr %>% t  
    return(null_comm)
}


#' Auxiliary function for get_delta_stats()
#' Returns a vector of abundances where individual-based rarefaction 
#' will be performed
#' @keywords internal
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
        if (max(inds) > N_max) 
            warning(paste('Sample size is higher than abundance of at least one group, only n up to',
                          N_max, 'will be used'))
        ind_sample_size = inds
    }
    # ensure that no more than N_max individuals considered
    ind_sample_size = unique(c(ind_sample_size[ind_sample_size < N_max], N_max))
    # Force (1, 1) to be included
    ind_sample_size = unique(c(1, ind_sample_size))
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
get_overall_p = function(effort, perm, value){
    delta_effort = c(effort[1], diff(effort))[perm == 0]
    Hbarbar = tapply(value, effort, mean)  # Baddeley Eq. A.10
    m = max(as.numeric(perm))              # number of permutations
    a = ((m + 1) / m)^2
    u = tapply(value, perm, function(x)    # Baddeley Eq. A.12-13
               a * sum((x - Hbarbar)^2 * delta_effort)) 
    overall_p = sum(u >= u[1]) / (m + 1)
    return(overall_p)
}

#' Extract coefficients and metrics of fit from model
#' @param x a fitted model object
#' @param stats the statistics to output
#' @importFrom stats pf coef
#' @keywords internal
mod_sum = function(x, stats = c('betas', 'r', 'r2', 'r2adj', 'f', 'p')) {
    summary_lm = summary(x)
    out = list()
    if ('betas' %in% stats) 
        out$betas = coef(x)
    if ('r' %in% stats) {
        betas = coef(x)
        out$r = sqrt(summary_lm$r.squared) * 
                ifelse(betas[2] < 0, -1, 1)
    }
    if ('r2' %in% stats)
        out$r2 = summary_lm$r.squared
    if ('r2adj' %in% stats)
        out$r2adj = summary_lm$adj.r.squared
    if ('f' %in% stats)
        out$f = summary_lm$fstatistic[1]
    if ('p' %in% stats) { # interpreted as overall model p-value
        f = summary_lm$fstatistic
        out$p = unname(stats::pf(f[1], f[2], f[3], lower.tail = FALSE))
    }
    if ('betas' %in% stats)
        coef_type = c(paste0('b', 0:(length(out$betas) - 1)),
                      stats[stats != 'betas'])
    else 
        coef_type = stats
    out = data.frame(coef_type, unlist(out))
    names(out) = c('index', 'value')
    row.names(out) = NULL
    out
}

#' @import purrr
#' @import dplyr 
#' @importFrom tidyr nest unnest
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom stats lm
#' @keywords internal 
get_results = function(mob_in, env, groups, tests, inds, ind_dens, n_plots, type,
                       stats=NULL, spat_algo=NULL) {
  
    # the approach taken here to get results for each group
    # is to first break the dataset up into a list of lists 
    # where this is one list per group - this is likely not 
    # the best practice for memory but it makes the code much 
    # easier to follow - we may need to revisit this. 
    group_levels = unique(groups)
    group_rows = map(group_levels, ~ which(groups == .x))
    mob_in_groups = map(group_rows, ~ subset(mob_in, .x, type = 'integer'))
    names(mob_in_groups) = group_levels
    
    S_df = map_dfr(mob_in_groups, get_delta_curves, tests, spat_algo,
                   inds, ind_dens, n_plots, .id = "group")
    
    S_df = S_df %>% try(mutate_if(is.factor, as.character), silent = TRUE)
    
    # substitute the group variable for the env variable
    S_df = data.frame(env = env[match(S_df$group, groups)],
                      S_df)

    S_df = tibble::as_tibble(S_df, .name_repair = 'minimal')
  
    # now that S and effects computed across scale compute
    # summary statistics at each scale 
  
    delta_mod = function(df) {
        stats::lm(effect ~ env, data = df)
    }
    
    if (is.null(stats)) {
        if (type == 'discrete')
            stats = 'betas'
        else
            stats = c('betas', 'r', 'r2', 'r2adj', 'f')
    }
    mod_df = S_df %>%
             group_by(.data$test, sample, .data$effort) %>%
             nest() %>%
             mutate(fit = map(.data$data, delta_mod)) %>%
             mutate(sum = map(.data$fit, mod_sum, stats)) %>%
             select(.data$test, sample, .data$effort, sum) %>% 
             unnest(sum) %>% 
             ungroup() %>%
             try(mutate_if(is.factor, as.character), silent = TRUE)
    
    return(list(S_df = S_df, mod_df = mod_df))
}

#' @import purrr
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats quantile
#' @importFrom rlang .data
#' @keywords internal
run_null_models = function(mob_in, env, groups, tests, inds, ind_dens, n_plots, type,
                           stats, spat_algo, n_perm, overall_p) {
    if (overall_p)
        p_val = vector('list', length(tests))
    # it may be possible to run all tests at the same time
    for (k in seq_along(tests)) {
        null_results = vector('list', length = n_perm)
        cat(paste('\nComputing null model for', tests[k], 'effect\n'))
        # need to parallelize this process optionally
        # in get_mob_stats we have the following:
        #   F_rand = bind_rows(pbreplicate(n_perm,
        #             get_F_values(dat_samples, permute = TRUE),
        #              simplify = FALSE, cl = cl)) %>%
        #             ungroup()
        pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
        for (i in 1:n_perm) {
            null_mob_in = mob_in
            null_mob_in$comm = get_null_comm(mob_in$comm, paste0('rand_', tests[k]),
                                             groups)
            null_results[[i]] = get_results(null_mob_in, env, groups, tests[k], inds,
                                            ind_dens, n_plots, type, stats, spat_algo)
            setTxtProgressBar(pb, i)
        }
        close(pb)    
        # rbind across the null_results adding a permutation index
        null_results = transpose(null_results)
        null_df = map(null_results, function(x)
                      flatten_dfr(tibble(x), .id = "perm"))
        # compute quantiles
        null_qt = list()
        null_qt$S_df = null_df$S_df %>% 
                       group_by(env, .data$test, sample, .data$effort) %>%
                       summarize(low_effect = quantile(.data$effect, 0.025, na.rm = TRUE),
                                 med_effect = quantile(.data$effect, 0.5, na.rm = TRUE), 
                                 high_effect = quantile(.data$effect, 0.975, na.rm = TRUE))
          
        null_qt$mod_df = null_df$mod_df %>%
                         group_by(.data$test, sample, .data$effort, .data$index) %>%
                         summarize(low_value = quantile(.data$value, 0.025, na.rm = TRUE),
                                   med_value = quantile(.data$value, 0.5, na.rm = TRUE), 
                                   high_value = quantile(.data$value, 0.975, na.rm = TRUE))
          
        if (k == 1) # if only one test run
            out = null_qt
        else 
            out = map2(out, null_qt, rbind)
        # to compute p-value we need to also calculate the observed
        # results then the funct must be distributed across the
        # various stats and tests
        if (overall_p) {
            obs_df = get_results(mob_in, env, groups, tests[k], inds, ind_dens,
                                 n_plots, type, stats, spat_algo)
            obs_df = map(obs_df, function(x) data.frame(perm = 0, x))          
            null_df = map2(obs_df, null_df, rbind)
            p_val[[k]] = list(effect_p = null_df$S_df %>%
                                  group_by(.data$test, .data$group) %>%
                                  summarize(p = get_overall_p(.data$effort, .data$perm, .data$effect)),
                              mod_p = null_df$mod_df %>%
                                  subset(!is.na(.data$value)) %>% 
                                  group_by(.data$test, .data$index) %>% 
                                  summarize(p = get_overall_p(.data$effort, .data$perm, .data$value)))
        }
    }
    if (overall_p)
        attr(out, "p") = map(transpose(p_val), bind_rows)
    return(out)
}


#' Conduct the MoB tests on drivers of biodiversity across scales.
#' 
#' There are three tests, on effects of 1. the shape of the SAD, 2.
#' treatment/group-level density, 3. degree of aggregation. The user can
#' specifically to conduct one or more of these tests.
#' 
#' @param mob_in an object of class mob_in created by make_mob_in()
#' @param env_var a character string specifying the environmental variable in
#'   \code{mob_in$env} to be used for explaining the change in richness
#' @param group_var an optional character string 
#'   in \code{mob_in$env} which defines how samples are pooled. If not provided
#'   then each unique value of the argument \code{env_var} is used define the
#'   groups. 
#' @param ref_level a character string used to define the reference level of
#'   \code{env_var} to which all other groups are compared with. Only makes sense
#'   if \code{env_var} is a factor (i.e. \code{type == 'discrete'})
#' @param tests specifies which one or more of the three tests ('SAD', N',
#'   'agg') are to be performed. Default is to include all three tests.
#' @param spat_algo character string that can be either: \code{'kNN'} or
#'   \code{'kNCN'} for k-nearest neighbor and k-nearest centroid neighbor
#'   sampling respectively. It defaults to k-nearest neighbor which is a more
#'   computationally efficient algorithm that closely approximates the
#'   potentially more correct k-NCN algo (see Details of ?rarefaction).
#' @param type "discrete" or "continuous". If "discrete", pair-wise comparisons
#'   are conducted between all other groups and the reference group. If
#'   "continuous", a correlation analysis is conducted between the response
#'   variables and env_var.
#' @param stats a vector of character strings that specifies what statistics to
#'   summarize effect sizes with. Options include: \code{c('betas', 'r2',
#'   'r2adj', 'f', 'p')} for the beta-coefficients, r-squared, adjusted
#'   r-squared, F-statistic, and p-value respectively. The default value of
#'   \code{NULL} will result in only betas being calculated when \code{type ==
#'   'discrete'} and all possible stats being computed when \code{type ==
#'   'continuous'}. Note that for a discrete analysis all non-betas stats are
#'   meaningless because the model has zero degrees of freedom in this context.
#' @param inds effort size at which the individual-based rarefaction curves are
#'   to be evaluated, and to which the sample-based rarefaction curves are to be
#'   interpolated. It can take three types of values, a single integer, a vector
#'   of integers, and NULL. If \code{inds = NULL} (the default), the curves are
#'   evaluated at every possible effort size, from 1 to the total number of
#'   individuals within the group (slow). If inds is a single integer, it is
#'   taken as the number of points at which the curves are evaluated; the
#'   positions of the points are determined by the "log_scale" argument. If inds
#'   is a vector of integers, it is taken as the exact points at which the
#'   curves are evaluated.
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
#' @param n_perm number of iterations to run for null tests, defaults to 1000.
#' @param overall_p Boolean defaults to FALSE specifies if overall across scale 
#'  p-values for the null tests. This should be interpreted with caution because
#'  the overall p-values depend on scales of measurement yet do not explicitly 
#'  reflect significance at any particular scale. 
#' @return a "mob_out" object with attributes
#' @author Dan McGlinn and Xiao Xiao
#' @import dplyr
#' @import purrr
#' @importFrom stats sd
#' @export
#' @seealso	\code{\link{rarefaction}}
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_level='uninvaded',
#'                            type='discrete', log_scale=TRUE, n_perm=3)
#' plot(inv_mob_out)
get_delta_stats = function(mob_in, env_var, group_var=NULL, ref_level = NULL, 
                           tests = c('SAD', 'N', 'agg'), spat_algo = NULL,
                           type = c('continuous', 'discrete'),
                           stats = NULL, inds = NULL,
                           log_scale = FALSE, min_plots = NULL,
                           density_stat = c('mean', 'max', 'min'),
                           n_perm=1000, overall_p = FALSE) {
    # perform preliminary checks and variable assignments
    if (class(mob_in) != "mob_in")
        stop('mob_in must be output of function make_mob_in (i.e., of class mob_in')
    if (!(env_var %in% names(mob_in$env)))
        stop(paste(env_var, ' is not one of the columns in mob_in$env.'))
    if (!is.null(group_var))
        if  (!(group_var %in% names(mob_in$env)))
            stop(paste(group_var, ' is not one of the columns in mob_in$env.')) 
    tests = match.arg(tests, several.ok = TRUE)
    test_status = tests %in% names(unlist(mob_in$tests)) 
    approved_tests = tests[test_status]
    if (length(approved_tests) < length(tests)) {
        tests_string = paste(approved_tests, collapse = ' and ')
        warning(paste('Based upon the attributes of the community object only the following tests will be performed:',
                  tests_string))
        tests = approved_tests
    }
    type = match.arg(type)
    density_stat = match.arg(density_stat)
    
    env = mob_in$env[ , env_var]
    # if group_var is NULL then set all samples to same group (??)
    if (is.null(group_var))
        groups = env
    else {
        groups = mob_in$env[ , group_var]
        # check that for the defined groups all samples have same environmental value
        if (any(tapply(env, groups, stats::sd) > 0)) {
            # bc all env values not the same for a group then compute mean value
            message("Computed average environmental value for each group")
            env = tapply(env, groups, mean)
        }
    }    
    if (type == 'discrete') {
        if (class(env) != 'factor') {
            warning(paste("Converting", env_var, "to a factor with the default contrasts because the argument type = 'discrete'."))
            env = as.factor(env)
        }
        if (!is.null(ref_level)) { # need to ensure that contrasts on the reference level set
            env_levels = levels(env) 
            if (ref_level %in% env_levels) {
                if (env_levels[1] != ref_level)
                    env = factor(env, levels = c(ref_level, env_levels[env_levels != ref_level]))
            } else
                stop(paste(ref_level, "is not in", env_var))
        }    
    } else if (type == 'continuous') {
        if (!is.numeric(env)) {
            warning(paste("Converting", env_var, "to numeric because the argument type = 'continuous'"))
            env = as.numeric(as.character(env))
        }
        if (!is.null(ref_level))
            stop('Defining a reference level (i.e., ref_level) only makes sense when doing a discrete analysis (i.e., type = "discrete")')
    }
    #TODO It needs to be clear which beta coefficients apply to 
    # which factor level - this is likely most easily accomplished by appending
    # a variable name to the beta column or adding an additional column
    #
    #if (is.null(env_var)){
    #    env_levels = as.numeric(names(sad_groups))
    #} else {
    #    env_levels = tapply(mob_in$env[, env_var],
    #                        list(groups), mean)
    #}
    N_max = min(tapply(rowSums(mob_in$comm), groups, sum))
    inds = get_inds(N_max, inds, log_scale)
    ind_dens = get_ind_dens(mob_in$comm, density_stat)
    n_plots = min(tapply(mob_in$comm[ , 1], groups, length))

    out = list()
    out$env_var = env_var
    if (!is.null(group_var))
        out$group_var = group_var
    out$type = type
    out$tests = tests
    out$log_scale = log_scale
    out$density_stat = list(density_stat = density_stat,
                            ind_dens = ind_dens)
    out = append(out, 
                 get_results(mob_in, env, groups, tests, inds, ind_dens, n_plots,
                             type, stats, spat_algo))

    null_results = run_null_models(mob_in, env, groups, tests, inds, ind_dens,
                                   n_plots, type, stats, spat_algo,
                                   n_perm, overall_p)
    # merge the null_results into the model data.frame
    out$S_df = left_join(out$S_df, null_results$S_df, 
                         by = c("env", "test", "sample", "effort"))
    out$mod_df = left_join(out$mod_df, null_results$mod_df, 
                           by = c("test", "sample", "effort", "index"))
    if (overall_p)
        out$p = attr(null_results, "p")
    class(out) = 'mob_out'
    return(out)
}

#' Plot distributions of species abundance
#' @inheritParams get_mob_stats
#' 
#' @param mob_in a 'mob_in' class object produced by 'make_mob_in'
#' @param type either 'sad' or 'rad' for species abundance vs rank abundance
#'   distribution
#' @param pooled Boolean defaults to FALSE which specifies that abundances should
#'   not be pooled at the group level, TRUE specifies that they should be pooled 
#' @param col optional vector of colors.
#' @param lwd a number which specifies the width of the lines
#' @param log a string that specifies if any axes are to be log transformed, 
#'   options include 'x', 'y' or 'xy' in which either the x-axis, y-axis, or
#'   both axes are log transformed respectively
#' @param leg_loc a string that specifies the location of the legend, 
#'   options include: 'lowerleft', 'topleft', 'loweright','topright' 
#' @importFrom scales alpha
#' @importFrom graphics lines legend
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' plot_abu(inv_mob_in, 'group', 'uninvaded', type='sad', pooled=FALSE, log='x')
#' plot_abu(inv_mob_in, 'group', 'uninvaded', type='rad', pooled=TRUE, log='x')
plot_abu = function(mob_in, group_var, ref_level = NULL, type=c('sad', 'rad'),
                    pooled=FALSE, col=NULL, lwd=3, log='', leg_loc = 'topleft') {
    groups  = factor(mob_in$env[ , group_var])
    group_levels = levels(groups) 
    # ensure that proper contrasts in groups 
    if (!is.null(ref_level)) { 
        if (ref_level %in% group_levels) {
            if (group_levels[1] != ref_level)
                groups = factor(groups, levels = c(ref_level, group_levels[group_levels != ref_level]))
            group_levels = levels(groups)
        } else
            stop(paste(ref_level, "is not in", group_var))
    }
    
    if (is.null(col)) 
        col = c("#FFB3B5", "#78D3EC", "#6BDABD", "#C5C0FE",
                "#E2C288", "#F7B0E6", "#AAD28C")    
    else if (length(col) != length(group_levels))
      stop('Length of col vector must match the number of unique groups')
    title = ifelse(pooled, 'Group Scale', 'Sample Scale')
    if ('sad' == type) {
        plot(1, type = "n", xlab = "% abundance", ylab = "% species", 
             xlim = c(0.01, 1), ylim = c(0.01, 1), log = log, main = title)
        for (i in 1:length(group_levels)) {
            col_grp = col[i]
            comm_grp = mob_in$comm[groups == group_levels[i], ]
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
                    sad_sort = sort(as.numeric(comm_grp[j, comm_grp[j, ] != 0]))
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
        plot(1:10, 1:10, type = 'n', xlab = 'rank', ylab = 'abundance',
             log = log, xlim = c(1, ncol(mob_in$comm)), 
             ylim = range(0.01, 1), cex.lab = 1.5, cex.axis = 1.5,
             main = title)
        for (i in 1:length(group_levels)) {
             col_grp = col[i]
             comm_grp = mob_in$comm[groups == group_levels[i], ]
             comm_grp = comm_grp[rowSums(comm_grp) > 0, ]
             if (pooled) {
                sad_grp = colSums(comm_grp)
                sad_sort = sort(sad_grp[sad_grp != 0], decreasing = TRUE)
                lines(sad_sort / sum(sad_sort), col = col_grp, lwd = lwd,
                      type = "l")
             } else {
                 for (j in 1:nrow(comm_grp)) {
                     sad_sort = sort(as.numeric(comm_grp[j, comm_grp[j, ] != 0], decreasing = TRUE))
                     lines(1:length(sad_sort), sad_sort / sum(sad_sort),
                           col = scales::alpha(col_grp, 0.5),
                           lwd = lwd, type = "l")
                 }     
             }
        }
    }
    if (!is.na(leg_loc))
        legend(leg_loc, legend = group_levels, col = col, lwd = lwd, bty = 'n')
}
    
#' Plot rarefaction curves for each treatment group
#' 
#' @param pooled Boolean specifying if samples should be pooled at the group
#'  level or not. Defaults to TRUE. This argument only applies when
#'  the individual based rarefaction is used (i.e., \code{method = 'indiv'})
#' @param ... other arguments to provide to \code{\link[mobr]{rarefaction}}
#' @inheritParams get_mob_stats
#' @inheritParams plot.mob_out
#' @inheritParams plot_abu
#' @inheritParams rarefaction
#' @importFrom scales alpha
#' @importFrom graphics lines legend
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' # random individual based rarefaction curves
#' plot_rarefaction(inv_mob_in, 'group', 'uninvaded', 'IBR',
#'                  pooled=TRUE, leg_loc='bottomright')
#' plot_rarefaction(inv_mob_in, 'group', 'uninvaded', 'IBR',
#'                  pooled=FALSE, log='x')
#' # random sample based rarefaction curves 
#' plot_rarefaction(inv_mob_in, 'group', 'uninvaded', 'SBR', log='xy')
#' # spatial sample based rarefaction curves 
#' plot_rarefaction(inv_mob_in, 'group', 'uninvaded', 'sSBR', log='xy')
plot_rarefaction = function(mob_in, group_var, ref_level = NULL,
                            method, dens_ratio = 1, pooled = TRUE, 
                            spat_algo = NULL, col = NULL, lwd = 3, log = '',
                            leg_loc = 'topleft', ...) {
    if (pooled == FALSE & method != 'IBR')
        stop('Samples can only not be pooled at the treatment level when individual-based rarefaction is used (i.e., method="IBR")')
    groups  = factor(mob_in$env[ , group_var])
    group_levels = levels(groups) 
    # ensure that proper contrasts in groups 
    if (!is.null(ref_level)) { 
        if (ref_level %in% group_levels) {
            if (group_levels[1] != ref_level)
                groups = factor(groups, levels = c(ref_level, group_levels[group_levels != ref_level]))
            group_levels = levels(groups)
        } else
            stop(paste(ref_level, "is not in", group_var))
    }
  
    if (is.null(col)) 
        col = c("#FFB3B5", "#78D3EC", "#6BDABD", "#C5C0FE",
                "#E2C288", "#F7B0E6", "#AAD28C")    
    else if (length(col) != length(group_levels))
        stop('Length of col vector must match the number of unique groups')
    if (method == 'indiv')
        xlab = 'Number of individuals'
    else
        xlab = 'Number of samples'
    if (pooled) {
        Srare = lapply(group_levels, function(x) 
                       rarefaction(subset(mob_in, groups == x, 'logical'),
                                   method, spat_algo = spat_algo, ...))
        xlim = c(1, max(unlist(sapply(Srare, function(x) as.numeric(names(x))))))
        ylim = c(1, max(unlist(Srare)))
        n = as.numeric(names(Srare[[1]]))
        plot(n, Srare[[1]], type = "n", main = "Group scale",
             xlab = xlab, ylab = "Species richness", 
             xlim = xlim, ylim = ylim, log = log)
        for (i in seq_along(group_levels)) {
            col_grp = col[i]
            n = as.numeric(names(Srare[[i]]))
            lines(n, Srare[[i]], col = col_grp, lwd = lwd, type = "l")
        }
    } else {
        Srare = lapply(group_levels, function(x)
                       apply(mob_in$comm[groups == x, ], 1,
                             function(y)  rarefaction(y, method, ...)))
        xlim = c(1, max(unlist(lapply(Srare, function(x)
                                 lapply(x, function(y)
                                        as.numeric(names(y)))))))
        ylim = c(1, max(unlist(Srare)))
        n = as.numeric(names(Srare[[1]][[1]]))
        plot(n, Srare[[1]][[1]], type = "n", main = "Sample scale",
             xlab = xlab, ylab = "Species richness", 
             xlim = xlim, ylim = ylim, log = log)        
        for (i in seq_along(group_levels)) {
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
        legend(leg_loc, legend = group_levels, col = col, lwd = lwd, bty = 'n')
}

#' Plot the multiscale MoB analysis output generated by \code{get_delta_stats}.
#' 
#' @param x a mob_out class object
#' @param stat a character string that specifies what statistic should be used
#'   in the effect size plots. Options include: \code{c('b0', 'b1', 'r', 'r2',
#'   'r2adj', 'f')} for the beta-coefficients, person correlation coefficient,
#'   r-squared, adjusted  r-squared, and F-statistic respectively. If the
#'   explanatory variable is a factor then \code{'b1'} is the only reasonable
#'   option. The default is set to the regression slope \code{'b1'} because this
#'   appears to have the strongest statistical power.
#' @param log2 a character string specifying if the x- or y-axis should be
#'   rescale by log base 2. Only applies when \code{display == 'S ~ effort' | 'S
#'   ~ effort'}. Options include: \code{c('', 'x', 'y', 'xy')} for  no
#'   rescaling, x-axis, y-axis, and both x and y-axes respectively. Default is
#'   set to no rescaling.
#' @param scale_by a character string specifying if sampling effort should be
#'   rescaled. Options include: \code{NULL}, \code{'indiv'}, and \code{'plot'}
#'   for no rescaling, rescaling to number of individuals, and rescaling
#'   to number of plots respectively. The rescaling is carried out using
#'   \code{mob_out$density_stat}.
#' @param display a string that specifies what graphical panels to display. 
#'  Options include:
#'  \itemize{
#'      \item \code{S ~ expl} ... plot of S versus the explanatory variable
#'      \item \code{S ~ effort} ... plot of S versus sampling effort (i.e., a
#'      rarefaction curve) 
#'      \item \code{effect ~ expl} ... plot of agg., N, and SAD effect size
#'      versus explanatory variable
#'      \item \code{stat ~ effort} ... plot of summary statistic versus sampling
#'      effort
#' }
#' Defaults to \code{'S ~ effort'}, \code{'effect ~ expl'}, and \code{'stat ~ effort'}.
#' @param eff_sub_effort Boolean which determines if only a subset of efforts
#'   will be considered in the plot of effect size (i.e., when
#'   \code{display = 'effect ~ expl'}. Defaults to TRUE to declutter the plots.
#' @param eff_log_base a positive real number that determines the base of the 
#'   logarithm that efforts were be distributed across, the larger this number
#'   the fewer efforts will be displayed.
#' @param eff_disp_pts Boolean to display the raw effect points, defaults to TRUE
#' @param eff_disp_smooth Boolean to display the regressions used to summarize
#'  the linear effect of the explanatory variable on the effect sizes, defaults
#'  to FALSE
#' @param ... parameters passed to other functions
#' 
#' @return plots the effect of the SAD, the number of individuals, and spatial
#'  aggregation on the difference in species richness
#'  
#' @author Dan McGlinn and Xiao Xiao
#' @import ggplot2 
#' @importFrom egg ggarrange
#' @importFrom stats predict loess lm
#' @importFrom grDevices rgb
#' @importFrom rlang .data
#' @export
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_level='uninvaded',
#'                               type='discrete', log_scale=TRUE, n_perm=4)
#' plot(inv_mob_out, 'b1') 
#' \donttest{ 
#' plot(inv_mob_out, 'b1', scale_by = 'indiv')
#' }
plot.mob_out = function(x, stat = 'b1', log2 = '', scale_by = NULL, 
                        display = c('S ~ effort', 'effect ~ grad', 'stat ~ effort'),
                        eff_sub_effort = TRUE, eff_log_base = 2,
                        eff_disp_pts = TRUE, eff_disp_smooth = FALSE, ...) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (x$type == 'discrete') {
        if (stat != 'b1')
            warning('The only statistic that has a reasonable interpretation for a discrete explanatory variable is the difference in the group means from the reference group (i.e., set stat = "b1")')
    }
    if (!is.null(scale_by)) {
        if (scale_by == 'indiv') {
            x$S_df = mutate(x$S_df, 
                             effort = ifelse(sample == 'plot', 
                                             round(.data$effort * x$density_stat$ind_dens),
                                             .data$effort))
            x$mod_df = mutate(x$mod_df, 
                               effort = ifelse(sample == 'plot', 
                                               round(.data$effort * x$density_stat$ind_dens),
                                               .data$effort))
      
        }     
        if (scale_by == 'plot') {
            x$S_df = mutate(x$S_df, 
                             effort = ifelse(sample == 'indiv', 
                                             round(.data$effort / x$density_stat$ind_dens),
                                             .data$effort))
            x$mod_df = mutate(x$mod_df, 
                               effort = ifelse(sample == 'indiv', 
                                               round(.data$effort / x$density_stat$ind_dens),
                                               .data$effort))
        } 
    }
    
    p_list = vector('list', 4)

    if ('S ~ grad' %in% display) {
        facet_labs = c(`agg` = 'sSBR',
                       `N` = 'nsSBR',
                       `SAD` = 'IBR')      
        p_list[[1]] = ggplot(x$S_df, aes(.data$env, .data$S)) +
                          geom_smooth(aes(group = .data$effort, color = .data$effort),
                                      method = 'lm', se = FALSE) +
                          labs(x = x$env_var) +
                          facet_wrap(. ~ test, scales = "free",
                                     labeller = as_labeller(facet_labs))
    }
    
    if ('S ~ effort' %in% display) {
        facet_labs = c(`agg` = 'sSBR',
                       `N` = 'nsSBR',
                       `SAD` = 'IBR')
        p_list[[2]] = ggplot(x$S_df, aes(.data$effort, .data$S)) +
                          geom_line(aes(group = .data$group, color = .data$env)) +
                          facet_wrap(. ~ test, scales = "free",
                                     labeller = as_labeller(facet_labs)) +
                          labs(y = expression("richness (" *
                                                italic(S) * ")"),
                               color = x$env_var) 
    }
    
    if ('effect ~ grad' %in% display) {
        efforts = sort(unique(x$S_df$effort))
        if (is.logical(eff_sub_effort)) {
            if (eff_sub_effort) {
                # only show a subset of efforts for clarity
                effort_r = floor(log(range(efforts), eff_log_base))
                effort_2 = eff_log_base^(effort_r[1]:effort_r[2])
                effort_2 = effort_2[effort_2 > 1] # effort at 1 indiv uninformative
                eff_d = as.matrix(stats::dist(c(efforts, effort_2)))
                eff_d = eff_d[-((length(efforts) + 1):ncol(eff_d)),
                              -(1:length(efforts))]
                min_index = apply(eff_d, 2, function(x) which(x == min(x))[1])
                sub_effort = efforts[min_index]
                message(paste("Effect size shown at the following efforts:",
                            paste(sub_effort, collapse = ', ')))
            }
            else 
                sub_effort = efforts
        } else if (!is.null(eff_sub_effort))
            sub_effort = eff_sub_effort

        if (x$type == "continuous")
            x$S_df = x$S_df %>%
                  group_by(.data$test, .data$effort) %>%
                  mutate(low_effect = predict(loess(low_effect ~ .data$env), .data$env)) %>%
                  mutate(high_effect = predict(loess(high_effect ~ .data$env), .data$env)) 

        
        p_list[[3]] = ggplot(subset(x$S_df, x$S_df$effort %in% sub_effort),
                                    aes(.data$env, .data$effect)) +
                      #geom_smooth(aes(x=group, y = med_effect,
                      #                group = effort, color = effort),
                      #            method = 'lm', se = FALSE) +
                      #geom_ribbon(aes(ymin = low_effect, ymax = high_effect,
                      #                group = effort, color = effort,
                      #                fill = 'null'),
                      #            alpha = .25) +
                      geom_hline(yintercept = 0, linetype = 'dashed') + 
                      labs(x = x$env_var) +
                      facet_wrap(. ~ test, scales = "free_y") +
                      labs(y = expression('effect (' * italic(S) * ')')) +
                      scale_fill_manual(name = element_blank(),
                                        values = c(null = 'grey40')) +
                      scale_colour_gradient2(trans=scales::log2_trans(),
                                            low = rgb(248, 203, 173, maxColorValue = 255),
                                            mid = rgb(237,127, 52, maxColorValue = 255),
                                            high = rgb(165, 0 , 33, maxColorValue = 255),
                                            midpoint = 4)
                #"#74c476" 
        if (eff_disp_pts)
            p_list[[3]] = p_list[[3]] + geom_point(aes(group = .data$effort,
                                                       color = .data$effort))
        if (eff_disp_smooth)
            p_list[[3]] = p_list[[3]] + geom_smooth(aes(group = .data$effort, 
                                                        color = .data$effort),
                                                    method = lm, se = FALSE)
    }

    if ('stat ~ effort' %in% display) {
        if (stat == 'b0')
            ylab = expression('intercept (' * italic(beta)[0] * ')')
        if (stat == 'b1')
            ylab = expression('slope (' * italic(beta)[1] * ')')
        if (stat == 'r2')
            ylab = expression(italic(R^2))
        if (stat == 'r')
            ylab = expression(italic(r))
        if (stat == 'f')
            ylab = expression(italic(F))
        p_list[[4]] = ggplot(subset(x$mod_df, x$mod_df$index == stat),
                          aes(.data$effort, .data$value)) + 
                          geom_ribbon(aes(ymin = .data$low_value,
                                          ymax = .data$high_value, fill = 'null'),
                                      alpha = 0.25) +
                          geom_line(aes(group = .data$index, color = 'observed')) +
                          geom_hline(yintercept = 0, linetype = 'dashed') + 
                          facet_wrap(. ~ test, scales = "free_x") +
                          labs(y = ylab) +
                          scale_color_manual(name = element_blank(),
                                             values = c(observed = 'red')) +
                          scale_fill_manual(name = element_blank(),
                                            values = c(null = 'grey40'))
    }
    
    if (!is.null(scale_by)) { # change title of legend
        scale_by = ifelse(scale_by == 'indiv', '# of individuals', '# of plots')
        if (!is.null(p_list[[1]]))  
            p_list[[1]] = p_list[[1]] + labs(color = scale_by)
        if (!is.null(p_list[[3]])) 
            p_list[[3]] = p_list[[3]] + labs(color = scale_by)
        
        if (!is.null(p_list[[2]]))
            p_list[[2]] =  p_list[[2]] + labs(x = scale_by)
        if (!is.null(p_list[[4]]))
            p_list[[4]] =  p_list[[4]] + labs(x = scale_by)

    }
    
    if (grepl('x', log2)) {
        if (!is.null(p_list[[2]]))
            p_list[[2]] = p_list[[2]] + scale_x_continuous(trans = 'log2')
        if (!is.null(p_list[[4]]))
            p_list[[4]] = p_list[[4]] + scale_x_continuous(trans = 'log2')
    }
  
    if (grepl('y', log2)) {
        if (!is.null(p_list[[2]]))
            p_list[[2]] = p_list[[2]] + scale_y_continuous(trans = 'log2')
    }
    p_list = Filter(Negate(is.null), p_list)
    egg::ggarrange(plots = p_list)
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
#' @importFrom graphics abline legend 
#' @export
#' @examples
#' data(inv_comm)
#' plot_N(inv_comm)
plot_N = function(comm, n_perm=1000) {
    N = rowSums(comm)
    ind_dens = mean(N)
    N_sum = apply(replicate(n_perm, cumsum(sample(N))), 1, mean)
    plot(N_sum, xlab = 'Number of plots', ylab = 'Number of Individuals')
    abline(a = 0, b = ind_dens, col = 'red')
    legend('topleft', 'Expected line', lty = 1, bty = 'n', col = 'red')
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
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
	if (sum(y < 0) > 0) stop("y cannot contain negative numbers")

	if (is.null(border)) border <- par("fg")
	border <- as.vector(matrix(border, nrow = ncol(y), ncol = 1))
	col <- as.vector(matrix(col, nrow = ncol(y), ncol = 1))
	lwd <- as.vector(matrix(lwd, nrow = ncol(y), ncol = 1))

  if (is.null(ylim)) ylim = c(0, 1.2 * max(apply(y, 1, sum)))
  
	if (order.method == "max") {
		ord <- order(apply(y, 2, which.max))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	if (order.method == "first") {
		ord <- order(apply(y, 2, function(x) min(which(x > 0))))
		y <- y[ , ord]
		col <- col[ord]
		border <- border[ord]
	}

	top.old <- x*0
	polys <- vector(mode  = "list", ncol(y))
	for (i in seq(polys)) {
		top.new <- top.old + y[,i]
		polys[[i]] <- list(x = c(x, rev(x)), y = c(top.old, rev(top.new)))
		top.old <- top.new
	}

	if (is.null(ylim)) 
	  ylim <- range(sapply(polys, function(x) range(x$y, na.rm = TRUE)), na.rm = TRUE)
	plot(x, y[ , 1],  ylab = ylab, xlab = xlab, ylim = ylim, t = "n", ...)
	for (i in seq(polys)) {
		polygon(polys[[i]], border = border[i], col = col[i], lwd = lwd[i])
	}

}

