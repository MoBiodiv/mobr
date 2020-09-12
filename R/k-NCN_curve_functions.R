#' Internal function used by kNCN_average to compute the k-NCN algorithm
#' starting with a specified focal sample
#' 
#' @param x a mob_in object or a community site x species matrix
#' @param focal_sample an integer from 1 to the number of samples of x that 
#' specifies which sample to start with. Defaults to 1
#' @param n the number of samples to accumulate, defaults to NULL in which 
#' case all samples are accumulated
#' @param coords the spatial coordinates of the samples of x
#' @param latlong if latitude longitude arguments are supplied
#' 
#' @importFrom geosphere centroid midPoint
#' @importFrom stats runif
#' 
#' @keywords internal
centroid_accumulate = function(x, focal_sample = 1, n = NULL, coords = NULL, latlong = FALSE) {

    if ("mob_in" %in% class(x)) {
        comm = x$comm
        coords = x$spat
    } else {
        comm = x
    }
    comm = as.matrix(comm)
    coords = as.matrix(coords)
    rownames(comm) = NULL
    rownames(coords) = NULL
    if (is.null(n)) 
        n = nrow(comm)
    included = focal_sample
    pool = colSums(comm[included, , drop = FALSE])
    S_accumulated = rep(0, n)
    S_accumulated[1] = sum(pool > 0)
    candidates = cbind(coords, 1:n)
    for (i in 1:(n - 1)) {
        accumulated_plots = coords[included, , drop = FALSE]
        if (latlong) {
            if (i == 1)
                samp_centroid = accumulated_plots
            else if (i == 2)
                samp_centroid = geosphere::midPoint(accumulated_plots[1, ],
                                                    accumulated_plots[2, ])
            else if (latlong)
                samp_centroid = geosphere::centroid(accumulated_plots)
        } else
            samp_centroid = matrix(colMeans(accumulated_plots), ncol = 2) 
        candidates2 = candidates[-included, , drop = FALSE]
        # it appears that sphere_dist in mobr may not be correct
        # therefore it may be wise to use geosphere's function when lat long 
        # supplied but there are problems with this approch near the poles and
        # it fails going over the poles.
        if (latlong) {
            # compute great circle distance assuming sphereical earth
            spat_dists = geosphere::distHaversine(samp_centroid,
                                                  candidates2[ , -3, drop = FALSE])
        } else {
            # compute Euclidean distance
            # calculation is take 
            #spat_dists = pracma::distmat(samp_centroid, candidates2[ , -3, drop = FALSE])
            spat_dists = apply(outer(samp_centroid, t(candidates2[ , -3, drop = FALSE]), "-"),
                               c(1, 4), function(x) sqrt(sum(diag(x * x))))
        }
        random_tie_breaker = stats::runif(nrow(candidates2))
        nearest_index = order(spat_dists, random_tie_breaker)[1]
        nearest_plot = candidates2[nearest_index, 3, drop = TRUE]
        included = c(included, nearest_plot)
        pool = pool + comm[nearest_plot, ]
        S_accumulated[i + 1] = sum(pool > 0)
    }
    return(S_accumulated)
}

#' Construct spatially constrained sample-based rarefaction (sSBR) curve using
#' the k-Nearest-Centroid-neighbor (k-NCN) algorithm
#'
#' This function accumulates samples according their proximity to all previously
#' included samples (their centroid) as opposed to the proximity to the initial
#' focal sample. This ensures that included samples mutually close to each other
#' and not all over the place.
#'
#' Internally the function constructs one curve per sample whereby each sample
#' serves as the initial sample repetition times. Finally, the average curve is
#' returned.
#'
#' @param x a mob_in object or a community site x species matrix
#' @param n number of sites to include.
#' @param coords spatial coordinates of the samples. If x is a mob_in object,
#'   the function uses its 'spat' table as coordinates.
#' @param repetitions Number of times to repeat the procedure. Useful in
#'   situations where there are many ties in the distance matrix.
#' @param no_pb binary, if TRUE then a progress bar is not printed, defaults to
#'   TRUE
#' @param latlong if longitude latitudes are supplied
#'
#' @inheritParams pbapply::pbreplicate
#'
#' @return a numeric vector of estimated species richness
#'
#' @importFrom pbapply pbsapply pboptions
#' 
#' @export
#'
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' kNCN_average(inv_mob_in, n = 5)
#' \donttest{
#' # parallel evaluation using the parallel package 
#' # run in parallel
#' library(parallel)
#' cl = makeCluster(2L)
#' clusterEvalQ(cl, library(mobr))
#' clusterExport(cl, 'inv_mob_in')
#' S_kNCN = kNCN_average(inv_mob_in, cl=cl)
#'
#' stopCluster(cl)
#' }
kNCN_average = function(x, n = NULL, coords = NULL, repetitions = 1, 
                        no_pb = TRUE, latlong = FALSE, cl = NULL) {
    if ("mob_in" %in% class(x)) {
        sites = x$comm
        coords = x$spat
    } else {
        sites = x
    }
    rownames(sites) = NULL
    rownames(coords) = NULL
    if (is.null(n)) 
        n = dim(sites)[1]
    samples = rep(1:n, each = repetitions)
    op <- pbapply::pboptions()
    if (no_pb)
        pbapply::pboptions(type = 'none')
    out = pbapply::pbsapply(samples, function(sample) 
                            centroid_accumulate(x = x, focal_sample = sample, 
                                                coords = coords, latlong = latlong),
                            cl = cl)
    pbapply::pboptions(op) # reset options to default
    return(rowMeans(out))
}

#' Compare all sample-based curves (random, spatially constrained-k-NN,
#'  spatially constrained-k-NCN)
#'
#'This is just plotting all curves.
#'
#' @param x a mob_in object
#'
#' @return a plot 
#' @importFrom graphics lines legend
#' @export
#'
#' @examples
#' \donttest{
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' compare_samp_rarefaction(inv_mob_in)
#' }
compare_samp_rarefaction = function(x) {
    
    centroid_curve = rarefaction(x, method = "sSBR", spat_algo = "kNCN")
    spatcurve = rarefaction(x, method = "sSBR", spat_algo = "kNN")
    samplecurve = rarefaction(x, method = "SBR")
    plot(samplecurve, type = "l", xlab = "Samples", ylab = "Expected species richness")
    lines(centroid_curve, col = 2)
    lines(spatcurve, col = 3)
    legend("bottomright", legend = c("sample based", "k-NCN", "k-NN"), col = 1:3, 
        lty = 1)
}

