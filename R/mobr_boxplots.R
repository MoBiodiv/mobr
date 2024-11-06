#' Estimation of species richness
#' 
#' \code{calc_chao1} estimates the number of species at the asymptote
#' (\code{S_asymp}) of the species accumulation curve based on the methods
#' proposed in Chao (1984, 1987, 2005). 
#' 
#' This function is a trimmed version of \href{https://github.com/JohnsonHsieh/iNEXT}{\code{iNext::ChaoRichness}}.
#' T. C. Hsieh, K. H. Ma and Anne Chao are the original authors of the
#' \code{iNEXT} package. 
#' 
#' @param x a vector of species abundances or a site-by-species matrix
#' 
#' @returns a vector of species richness estimates
#' 
#' @examples 
#' data(inv_comm)
#' calc_chao1(inv_comm)
#' @references 
#' Chao, A. (1984) Nonparametric estimation of the number of classes in a
#' population. Scandinavian Journal of Statistics, 11, 265-270.
#' 
#' Chao, A. (1987) Estimating the population size for capture-recapture data with
#' unequal catchability. Biometrics, 43, 783-791.
#' 
#' Chao, A. (2005) Species estimation and applications. Pages 7907-7916 in
#' N. Balakrishnan, C. B. Read, and B. Vidakovic, editors. Encyclopedia of
#' statistical sciences. Second edition, volume 12. Wiley, New York, New York,
#' USA.
#' 
#' @export
calc_chao1 = function(x) {
    if (!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
        stop("invalid data structure")
    if (is.matrix(x) | is.data.frame(x)) {
        S_Chao1 = apply(x, 1, calc_chao1)
    } else {
        n = sum(x)
        D = sum(x > 0)
        f1 = sum(x == 1)
        f2 = sum(x == 2)
        if (f1 > 0 & f2 > 0)
            S_Chao1 = D + (n - 1) / n * f1^2 / (2 * f2)
        else if (f1 > 1 & f2 == 0) # bias corrected form
            S_Chao1 = D + (n - 1) / n * f1 * (f1 - 1) / (2 * (f2 + 1))
        else
            S_Chao1 = D
    }
    return(S_Chao1)
}

#' Calculate probability of interspecific encounter (PIE)
#' 
#' \code{calc_PIE} returns the probability of interspecific encounter (PIE)
#'  which is also known as Simpson's evenness index and Gini-Simpson index. 
#' 
#' By default, Hurlbert's (1971) sample-size corrected formula is used:
#' 
#' \eqn{PIE = N /(N - 1) * (1 - sum(p_i^2))}
#' 
#' where N is the total number of individuals and \eqn{p_i} is the relative
#' abundance of species i. This formulation uses sampling without replacement
#' (\code{replace = FALSE} ) For sampling with replacement (i.e., the sample-size
#' uncorrected version), set \code{replace = TRUE}.
#'
#' In earlier versions of \code{mobr}, there was an additional argument
#' (\code{ENS}) for the conversion into an effective number of species (i.e
#' S_PIE). Now, \code{calc_SPIE} has become its own function and the
#' (\code{ENS}) argument is no longer supported . Please, use \code{calc_SPIE}
#' instead.
#'
#' 
#' @inheritParams rarefaction
#' @param PIE_replace if TRUE, sampling with replacement is used. Otherwise,
#'   sampling without replacement (default).
#'
#' @returns either a single PIE value or vector of PIE values. 
#' 
#' @seealso \code{\link{calc_SPIE}}
#'
#' @author Dan McGlinn, Thore Engel
#' 
#' @references 
#' Hurlbert, S. H. (1971) The nonconcept of species diversity: a critique and
#'  alternative parameters. Ecology 52, 577-586.
#'  
#' @export
#' @examples 
#' data(inv_comm)
#' calc_PIE(inv_comm)
#' calc_PIE(inv_comm, PIE_replace = TRUE)
#' calc_PIE(c(23,21,12,5,1,2,3))
#' calc_PIE(c(23,21,12,5,1,2,3), PIE_replace = TRUE)
calc_PIE = function(x, PIE_replace = FALSE) {
    
    args = as.list(match.call())
    if (any(names(args) == "ENS")) 
        stop("The ENS argumet was removed from this function. Please, use calc_SPIE() for the ENS transformation of PIE. ")
    
    if ('mob_in' %in% class(x)) {
        x = x$comm
    }
    x = drop(as.matrix(x))
    if (any(x < 0, na.rm = TRUE)) 
        stop("input data must be non-negative")
    
    if (any(x %% 1 != 0, na.rm = TRUE))
        stop("input data must be integers")
    
    if (length(dim(x)) > 1) {
        total = apply(x, 1, sum)
        S = apply(x, 1, function(x) return(sum(x > 0)))
        p_i = sweep(x, 1, total, "/")
    } else {
        total = sum(x)
        S = sum(x > 0)
        p_i = x / total
    }
    p_i_sq = p_i * p_i
    if (length(dim(x)) > 1) {
        H = rowSums(p_i_sq, na.rm = TRUE)
    } else {
        H = sum(p_i_sq, na.rm = TRUE)
    }
    
    # calculate PIE without replacement (for total >= 2)
    if (PIE_replace) {
        PIE = 1 - H
    } else {
        PIE = total / (total - 1) * (1 - H)
    }
    # if sample had zero individuals set PIE to 0
    PIE[total == 0] = 0
    # if sample only contains 1 individual set PIE to NA
    if (!PIE_replace) 
        PIE[total == 1] = NA
    if (any(is.na(PIE))) 
        warning("NA was returned because the sample contains one or zero individuals.")
    
    return(PIE)
}

#' Calculate S_PIE
#'
#' S_PIE is the effective number of species transformation of the probability of
#' interspecific encounter (PIE) which is equal to the number of equally common
#' species that result in that value of PIE.
#'
#' By default the sample size corrected version is returned (\code{PIE_replace =
#' FALSE}), which is the asymptotic estimator for the Hill number of diversity order
#' q=2 (Chao et al, 2014). If \code{PIE_replace = TRUE} the uncorrected hill number is
#' returned. This is the same as vegan::diversity(x, index="invsimpson").
#'
#' 
#' @inheritParams calc_PIE
#'
#' @returns either a single S_PIE value or vector of S_PIE values. 
#' 
#' @seealso \code{\link{calc_PIE}}
#' 
#' @export
#' 
#' @references 
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R.
#' K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers:
#' A framework for sampling and estimation in species diversity studies.
#' Ecological Monographs 84(1), 45-67.
#'
#' @examples
#' data(inv_comm)
#' calc_SPIE(inv_comm)
#' calc_SPIE(inv_comm, PIE_replace = TRUE)
#' calc_SPIE(c(23,21,12,5,1,2,3), PIE_replace=TRUE)
calc_SPIE = function(x, PIE_replace = FALSE) {
    
    PIE = calc_PIE(x, PIE_replace = PIE_replace)
    SPIE = 1 / (1 - PIE)
    SPIE[sapply(PIE, function(x)
        isTRUE(all.equal(x, 0)))] = 0
    SPIE[sapply(PIE, function(x)
        isTRUE(all.equal(x, 1)))] = NA
    if (any(sapply(PIE, function(x)
        isTRUE(all.equal(x, 1))), na.rm = T))
        warning(
            "NA was returned because PIE = 1. This happens in samples where all species are singletons."
        )
    
    return(SPIE)
}

#' Compute various diversity indices from a vector of species abundances (i.e.,
#' one row of a community matrix)
#'
#' @param x is a vector of species abundances
#' @param C_target When computing coverage based richness (\code{S_C}) then 
#' this argument can be used to specify the coverage to be used for the richness
#' estimate. This defaults to \code{NA} in which case the target cover
#' is computed by \code{\link{calc_C_target}} (i.e., the largest allowable sample
#' size).
#' @param ... additional arguments that can be passed to the function
#'  \code{rarefaction} when computing \code{S_n}. 
#'
#' @inheritParams calc_comm_div
#' 
#'
#' @export
#' @examples  
#' data(inv_tank)
#' calc_div(tank_comm[1, ], 'S_n', effort = c(5, 10))
#' calc_div(tank_comm[1, ], 'S_C', C_target = 0.9)
calc_div = function(x, index, effort=NA, rare_thres = 0.05, PIE_replace = FALSE,
                    C_target = NULL, extrapolate = TRUE, ...) {
    if (index == 'N') out = sum(x)
    if (index == 'S') out = sum(x > 0)
    if (index == 'S_n') out = rarefaction(x, method = 'IBR', effort = effort,
                                          extrapolate = extrapolate, ...) 
    if (index == 'S_C') out = calc_S_C(x, C_target, extrapolate = extrapolate,
                                       interrupt = FALSE)
    if (index == 'PIE') out = calc_PIE(x, PIE_replace = PIE_replace)
    if (index == 'S_PIE') out = calc_SPIE(x, PIE_replace = PIE_replace)
    if (index == 'f_0') out = calc_div(x, 'S_asymp') - calc_div(x, 'S')
    if (index == 'S_asymp') {
        S_asymp = try(calc_chao1(x))
        if (methods::is(S_asymp, "try_error"))
            warning("The Chao richness estimator cannot be calculated for all samples.")
        else 
            S_asymp[!is.finite(S_asymp)] = NA
        out = S_asymp
    }    
    if (index == 'pct_rare') {
        S = calc_div(x, 'S') 
        if (S > 0) {
            N = calc_div(x, 'N')
            if (rare_thres == "N/S") {
                rare_thres = N / S
                out = 100 * (sum(x[x > 0] <= rare_thres) / S)
            } else 
                out = 100 * (sum(x[x > 0] <= (rare_thres * N)) / S)
        } else
            out = 0
    }
    return(out)
}
        

#' Calculate biodiversity statistics from sites by species table.
#'
#' @param abund_mat Abundance based site-by-species table. Species as columns
#' @param index The calculated biodiversity indices. The options are
#' \itemize{
#'    \item \code{N} ... Number of individuals (total abundance)
#'    \item \code{S} ... Number of species
#'    \item \code{S_n} ... Rarefied or extrapolated number of species for n individuals
#'    \item \code{S_C} ... Estimate species richness of a given level of coverage by \code{C_target_gamma}
#'    \item \code{S_asymp} ... Estimated asymptotic species richness
#'    \item \code{f_0} ... Estimated number of undetected species
#'    \item \code{pct_rare} ... The percent of rare species as defined by \code{rare_thres}
#'    \item \code{PIE} ... Hurlbert's PIE (Probability of Interspecific Encounter)
#'    \item \code{S_PIE} ... Effective number of species based on PIE
#'
#' }
#'   See \emph{Details} for additional information on the biodiversity
#'   statistics.
#'
#' @param effort The standardized number of individuals used for the calculation
#'   of rarefied species richness. This must be a single integer. 
#'
#' @param extrapolate Boolean which specifies if richness should be extrapolated
#'   when effort is larger than the number of individuals using the chao1
#'   method (Chao 1984, 1987). Defaults to TRUE.
#'
#' @param return_NA Boolean in which the rarefaction function returns the
#'   observed S when \code{effort} is larger than the number of individuals. If
#'   set to TRUE then NA is returned. Note that this argument is only relevant
#'   when \code{extrapolate = FALSE}.
#'
#' @param rare_thres The threshold that determines how the metric
#'   \code{pct_rare} is computed. It can range from (0, 1] and defaults to 0.05
#'   which specifies that any species with less than or equal to 5% of the total
#'   abundance in a sample is considered rare. It can also be specified as "N/S"
#'   which results in using average abundance as the threshold which McGill
#'   (2011) found to have the best small sample behavior.
#'
#' @param scales The scales to compute the diversity indices for:
#' \itemize{
#'     \item \code{alpha} ... for each row of the site x species community matrix
#'     \item \code{gamma} ... for the entire site x species community matrix
#'     \item \code{beta} ... the ratio of diversity at the \code{gamma} and
#'                            \code{alpha} scales.
#' } Defaults to all three scales: \code{c('alpha', 'gamma', 'beta')}
#'
#' @param avg_alpha Boolean if TRUE then the alpha values are averaged. Defaults
#' to FALSE.
#' 
#' @param PIE_replace Used for \code{PIE} and \code{SPIE}.  If TRUE, sampling with
#'   replacement is used. Otherwise, sampling without replacement (default).
#'
#' @param C_target_gamma When computing coverage based richness (\code{S_C})
#'   then this argument can be used to specify the coverage to be used for the
#'   gamma scale richness estimate. This defaults to \code{NA} in which case the
#'   target cover is computed by \code{\link{calc_C_target}} (i.e., the largest
#'   allowable sample size).
#'
#' @param ... additional arguments that can be passed to \code{\link{calc_div}}
#'
#' @details
#'
#' \strong{BIODIVERSITY INDICES}
#'
#' \strong{N: total community abundance} is the total number of individuals
#' observed across all species in the sample
#'
#' \strong{S: species richness} is the observed number of species that occurs at
#' least once in a sample
#'
#' \strong{S_n: Rarefied species richness} is the expected number of species,
#' given a defined number of sampled individuals (n) (Gotelli & Colwell 2001).
#' Rarefied richness at the alpha-scale is calculated for the values provided in
#' \code{effort_samples} as long as these values are not smaller than the
#' user-defined minimum value \code{effort_min}. In this case the minimum value
#' is used and samples with less individuals are discarded. When no values for
#' \code{effort_samples} are provided the observed minimum number of individuals
#' of the samples is used, which is the standard in rarefaction analysis
#' (Gotelli & Colwell 2001). Because the number of individuals is expected to
#' scale linearly with sample area or effort, at the gamma-scale the number of
#' individuals for rarefaction is calculated as the minimum number of samples
#' within groups multiplied by \code{effort_samples}. For example, when there
#' are 10
#' samples within each group, \code{effort_groups} equals \code{10 *
#' effort_samples}. If n is larger than the number of individuals in sample and
#' \code{extrapolate = TRUE} then the Chao1 (Chao 1984, 1987) method is used to
#' extrapolate the rarefaction curve.
#'
#' \strong{pct_rare: Percent of rare species} Is the ratio of the number of rare
#' species to the number of observed species x 100 (McGill 2011). Species are
#' considered rare in a particular sample if they have fewer individuals than
#' \code{rare_thres * N} where \code{rare_thres} can be set by the user and
#' \code{N} is the total number of individuals in the sample. The default value
#' of \code{rare_thres} of 0.05 is arbitrary and was chosen because McGill
#' (2011) found this metric of rarity performed well and was generally less
#' correlated with other common metrics of biodiversity. Essentially this metric
#' attempt to estimate what proportion of the species in the same occur in the
#' tail of the species abundance distribution and is therefore sensitive to
#' presence of rare species.
#'
#' \strong{S_asymp: Asymptotic species richness} is the expected number of
#' species given complete sampling and here it is calculated using the Chao1
#' estimator (Chao 1984, Chao 1987) see \code{\link{calc_chao1}}. Note: this
#' metric is typically highly correlated with S (McGill 2011).
#'
#' \strong{f_0: Undetected species richness} is the number of undetected species
#' or the number of species observed 0 times which is an indicator of the degree
#' of rarity in the community. If there is a greater rarity then f_0 is expected
#' to increase. This metric is calculated as \code{S_asymp - S}. This metric is
#' less correlated with S than the raw \code{S_asymp} metric.
#'
#' \strong{PIE: Probability of intraspecific encounter} represents the
#' probability that two randomly drawn individuals belong to the same species.
#' Here we use the definition of Hurlbert (1971), which considers sampling
#' without replacement. PIE is closely related to the well-known Simpson
#' diversity index, but the latter assumes sampling with replacement.
#'
#' \strong{S_PIE: Effective number of species for PIE} represents the effective
#' number of species derived from the PIE. It is calculated using the asymptotic
#' estimator for Hill numbers of diversity order 2 (Chao et al. 2014). S_PIE
#' represents the species richness of a hypothetical community with
#' equally-abundant species and infinitely many individuals corresponding to the
#' same value of PIE as the real community. An intuitive interpretation of S_PIE
#' is that it corresponds to the number of dominant (highly abundant) species in
#' the species pool.
#'
#' For species richness \code{S}, rarefied richness \code{S_n}, undetected
#' richness \code{f_0}, and the Effective Number of Species \code{S_PIE} we also
#' calculate beta-diversity using multiplicative partitioning (Whittaker 1972,
#' Jost 2007). That means for these indices we estimate beta-diversity as the
#' ratio of gamma-diversity (total diversity across all plots) divided by
#' alpha-diversity (i.e., average plot diversity).
#'
#' @returns A \code{data.frame} with four columns:
#' \itemize{
#'    \item \code{scale} ... Group label for sites
#'    \item \code{index} ... Name of the biodiversity index
#'    \item \code{sample_size} ... The number of samples used to compute the
#'     statistic, helpful for interpreting beta and gamma metrics.
#'    \item \code{effort} ... Sampling effort for rarefied richness
#'    (NA for the other indices)
#'    \item \code{gamma_coverage} ... The coverage value for that particular
#'    effort value on the gamma scale rarefaction curve. Will be \code{NA} unless
#'    coverage based richness (\code{S_C}) and/or beta diversity is computed.
#'    \item \code{value} ... Value of the biodiversity index
#' }
#'
#' @author Felix May and Dan McGlinn
#'
#' @references
#'
#' Chao, A. 1984. Nonparametric Estimation of the Number of Classes in a
#' Population. Scandinavian Journal of Statistics 11:265–270.
#'
#' Chao, A. 1987. Estimating the population size for capture-recapture data with
#' unequal catchability. Biometrics, 43, 783-791.
#'
#' Gotelli, N. J., and R. K. Colwell. 2001. Quantifying biodiversity: procedures
#' and pitfalls in the measurement and comparison of species richness. Ecology
#' Letters 4:379–391.
#'
#' Hurlbert, S. H. 1971. The nonconcept of species diversity: a critique and
#' alternative parameters. Ecology 52:577–586.
#'
#' Jost, L. 2007. Partitioning diversity into independent alpha and beta
#' components. Ecology 88:2427–2439.
#'
#' McGill, B. J. 2011. Species abundance distributions. Pages 105-122 Biological
#' Diversity: Frontiers in Measurement and Assessment, eds. A.E. Magurran and
#' B.J. McGill.
#'
#' Whittaker, R. H. 1972. Evolution and measurement of species diversity. Taxon
#' 21:213–251.
#'
#'
#' @examples
#' data(tank_comm)
#' div_metrics <- calc_comm_div(tank_comm, 'S')
#' div_metrics
#' div_metrics <- calc_comm_div(tank_comm, 'S', avg_alpha = TRUE)
#' div_metrics
#' div_metrics <- calc_comm_div(tank_comm, 'S_n', effort = 10, avg_alpha = TRUE)
#' div_metrics
#' div_metrics <- calc_comm_div(tank_comm, 'S_C', C_target_gamma = 0.75, avg_alpha = TRUE)
#' div_metrics
#' @export
calc_comm_div = function(abund_mat, index, effort = NA, 
                         extrapolate = TRUE,
                         return_NA = FALSE, rare_thres = 0.05,
                         scales = c('alpha', 'gamma', 'beta'),
                         avg_alpha = FALSE, 
                         PIE_replace = FALSE, C_target_gamma = NA, ...) {
    
    # store each calculated index into its own data.frame in a list
    out = vector('list', length = length(index))
    names(out) = index
    if (any(index == 'S_n') && any(is.na(effort))) 
        stop('effort value is needed to compute S_n')
    if (any(index == 'S_C') & is.na(C_target_gamma))
        C_target_gamma <- calc_C_target(abund_mat)
    # compute indices ---------------------------------------------------------
    for (i in seq_along(index)) {
        if (any(c('gamma', 'beta') %in% scales))
            gamma = calc_div(colSums(abund_mat), index[i], effort, rare_thres,
                         extrapolate = extrapolate, return_NA = return_NA, 
                         quiet = TRUE, PIE_replace = PIE_replace, C_target = C_target_gamma, ...)
        if (any(c('alpha', 'beta') %in% scales)) {
            if (index[i] == 'S_C' & 'beta' %in% scales) {
                effort_eff = attributes(gamma)$N
                index_eff = 'S_n'
            }
            else {
                index_eff = index[i]
                effort_eff = effort
            }    
            alpha = apply(abund_mat, 1, calc_div, index_eff, effort_eff, rare_thres,
                          extrapolate = extrapolate, return_NA = return_NA, 
                          quiet = TRUE, PIE_replace = PIE_replace,
                          C_target = C_target_gamma, ...)
        }
        if ('beta' %in% scales) {
            # compute beta
            if (index[i] == 'S_n' & length(effort) > 1) {
                beta = gamma / rowMeans(alpha, na.rm = TRUE)
            } else {
                beta = gamma / mean(alpha, na.rm = TRUE)
            }
        }  
        if (index[i] == 'S_n' | index[i] == 'S_C')
            effort_out <- effort_eff 
        else
            effort_out <- NA
        gamma_coverage <- ifelse(index[i] == 'S_C', C_target_gamma, NA)
        # compute number of finite samples used for calculation
        sample_size <- nrow(abund_mat)
        if ('alpha' %in% scales) {
            if (avg_alpha) 
                alpha <- mean(alpha, na.rm = TRUE)
            out[[i]]$alpha = data.frame(scale = 'alpha', index = index[i],
                                        sample_size = ifelse(avg_alpha, sample_size, 1),
                                        effort = effort_out,
                                        gamma_coverage = gamma_coverage,
                                        value = alpha)
        }    
        if ('gamma' %in% scales) 
            out[[i]]$gamma = data.frame(scale = 'gamma', index = index[i], 
                                        sample_size = 1, effort = effort_out,
                                        gamma_coverage = gamma_coverage,
                                        value = gamma)
        if ('beta' %in% scales & index[i] != 'N') 
            out[[i]]$beta = data.frame(scale = 'beta',
                                       index = paste('beta', index[i], sep = '_'),
                                       sample_size = 1, effort = effort_out,
                                       gamma_coverage = gamma_coverage,
                                       value = beta)
        
        if (index[i] == 'S_n' & length(effort) > 1) {
            out[[i]] = lapply(out[[i]], dplyr::arrange, effort)
        }
        out[[i]] = dplyr::bind_rows(out[[i]])
        row.names(out[[i]]) = 1:nrow(out[[i]])
    }
    out = dplyr::bind_rows(out)
    return(out)
}

#' Calculate beta diversity from sites by species table.
#' 
#' This function computes multiplicative beta diversity by computing the ratio 
#' of gamma diversity to average alpha diversity. 
#' It is a wrapper for the function \code{calc_comm_div} when that function's
#' \code{scales} is set to \code{'betas'}. 
#' 
#' @inherit calc_comm_div
#' @param ... other arguments to pass to \code{calc_comm_div}
#' @seealso \code{\link{calc_comm_div}}
#' @examples 
#' data(inv_comm)
#' beta_metrics = calc_beta_div(inv_comm, 'S_n', effort = c(5, 10))
#' beta_metrics
#' @export
calc_beta_div = function(abund_mat, index, effort = NA, C_target_gamma = NA, 
                         extrapolate = TRUE, ...) {
    out <- calc_comm_div(abund_mat, index, effort,scales = 'beta', 
                         C_target_gamma = C_target_gamma, ...)
    return(out)
}

#' Generate distribution of sampled community matrices from an original matrix. 
#' 
#' The sampled matrices are either bootstrap or leave-one-out (the default)
#' samples. The bootstrap samples represent a random set of the rows sampled
#' with replacement from the original matrix. The leave-one-out samples are
#' generated by removing each sample one at a time from the original community
#' matrix.
#' 
#' These sampled community matrices can become the input to \code{calc_comm_div_ci}
#' for computing confidence intervals of diversity metrics. 
#' 
#' Note, it is unclear which sampling algorithm (bootstrap or loo) is least
#' biased and most efficient for the purpose of generating confidence intervals.
#' 
#' @inherit calc_comm_div 
#' @param algo can be either 'boot' or 'loo' for bootstrap or leave-one-out 
#'   methods respectively. Default value is 'loo'. 
#' @param n_boot how to many boot strapped samples to create, defaults to 1000.
#'
#' @returns a list of community matrices which are sampled from the original 
#'   input matrix. 
#' 
#' @seealso \code{\link{calc_comm_div_ci}}
#' @export
#' @examples
#' data(tank_comm)
#' # 2 leave-one-out samples
#' lapply(get_samples(tank_comm)[1:2], head)
#' # 2 bootstrap samples
#' lapply(get_samples(tank_comm, algo = 'boot', n_boot = 2), head)
#' 
get_samples <- function(abund_mat, algo = 'loo', n_boot = 1000) {
  nsamp <- nrow(abund_mat)
  if (algo == 'boot') {
    # if bootstrap samples wanted
    sample_out <- replicate(n_boot, abund_mat[sample(nsamp, replace = TRUE), ], 
                            simplify = FALSE)
  } else if (algo == 'loo') {
    # if leave-one-out samples wanted
    sample_out <- lapply(1:nsamp, function(i) abund_mat[-i, ])
  }
  return(sample_out)
}

#' Compute non-parametric confidence intervals across diversity indices. 
#' 
#' This function take a list of community matrices and returns the central tendency
#' and the confidence interval for each diversity index of interest across that
#' list of communities.  
#' 
#' The measure of central tendency can be the median (default) or mean, and the
#' range of the confidence interval can be specified but defaults to a 95% 
#' confidence interval.  
#' 
#' @param samples a list of community matrices (i.e., the output of \code{get_samples})
#' @param cent_stat a string that is either 'mean' or 'median' which specifies the measure
#'   of central tendency. Defaults to 'median'.
#' @param ci a numeric vector of two numbers specifying the lower and upper
#'   quantiles of the distribution to return. The default is to report
#'   the 0.025 and 0.975 quantiles in other words a 95 percent interval.
#' @inheritParams calc_comm_div
#'
#' @returns a data.frame that has the lower, middle, and upper quantiles of the 
#'  sampled distribution of each diversity index. 
#' 
#' @seealso \code{\link{get_samples}} for generating samples, and \code{\link{calc_comm_div}}
#'   for the calculation of diversity indices. 
#' 
#' @importFrom rlang .data
#' @importFrom stats median
#' 
#' @export
#' @examples
#' data(tank_comm)
#' samples <- get_samples(tank_comm, algo = 'loo')
#' calc_comm_div_ci(samples, index = 'S_PIE')
#' samples <- get_samples(tank_comm, algo = 'boot', n_boot = 20)
#' calc_comm_div_ci(samples, index = 'S_PIE')
#' # compute ci for average diversity (rather than median)
#' calc_comm_div_ci(samples, index = 'S_PIE', cent_stat = 'avg')
calc_comm_div_ci <- function(samples, cent_stat = 'median', ci = c(0.025, 0.975),
                             index, effort = NA, extrapolate = TRUE,
                             return_NA = FALSE, rare_thres = 0.05,
                             scales = c('alpha', 'gamma', 'beta'),
                             PIE_replace = FALSE, C_target_gamma = NA, ...) {
    if (!is.numeric(ci) & length(ci) != 2)
        stop('The ci argument must be a numeric vector of length 2 specifing the lower and upper quantiles of the confidence interval to return.')
    if (any(index %in% 'S_C') & is.na(C_target_gamma)) {
        C_target_gamma <- min(sapply(samples, calc_C_target))
    }
    # compute the diversity indices on a random sample of 
    sample_div <- pbapply::pblapply(samples, function(x)
        calc_comm_div(x, index, effort, extrapolate, 
                      return_NA, rare_thres, scales, avg_alpha = TRUE,
                      PIE_replace, C_target_gamma, ...))
    # bind across the bootstrap replicates
    sample_div <- dplyr::bind_rows(sample_div, .id = 'id')
    # compute quantiles across bootstraps. 
    sample_qts <- sample_div |> 
        dplyr::group_by(scale, index) |>        
        dplyr::summarize(
            sample_size = mean(.data$sample_size),
            effort = mean(effort),
            gamma_coverage = mean(.data$gamma_coverage),
            lo_value = quantile(.data$value, ci[1], na.rm = TRUE),
            hi_value = quantile(.data$value, ci[2], na.rm = TRUE), 
            value = ifelse(cent_stat == 'avg', mean(.data$value, na.rm = TRUE),
                           median(.data$value, na.rm = TRUE)),
            .groups = 'keep')

    return(data.frame(sample_qts))
}



#' Carry out biodiversity metric comparisons between groups. 
#' 
#' This function can compute a range of biodiversity metrics, their uncertainty, 
#' and carries out a permutation test to examine if groups differ in a particular
#' metric more than would be expected due to random chance. 
#' 
#' This function is partially a wrapper for the functions
#' \code{\link{calc_comm_div}} or \code{\link{calc_comm_div_ci}} that makes
#' group comparisons easier to implement.
#' 
#' @inheritParams get_delta_stats
#' @inheritParams calc_comm_div
#' @inheritParams pbapply::pbreplicate
#' 
#' @param group_var String that specifies which variable in \code{mob_in$env} the
#'   data should be grouped by.
#' 
#' @param effort_samples An integer that specifies the standardized number of
#'   individuals used for the calculation of rarefied species richness at the
#'   alpha-scale. It must be a single integer. The default value of \code{NULL}
#'   will result in the using the minimum number of individuals found across the
#'   samples is used, when this is not smaller than \code{effort_min}.
#'   
#' @param effort_min The minimum number of individuals considered for the 
#'   calculation of rarefied richness (Default value of 5). Samples with less
#'   individuals then \code{effort_min} are excluded from the analysis with a
#'   warning. Accordingly, when \code{effort_samples} is set by the user it has
#'   to be higher than \code{effort_min}.
#'
#' @param ci boolean, if TRUE then confidence intervals are calculated. Defaults
#'  to TRUE. 
#' 
#' @param ci_cent_stat a string that is either 'mean' or 'median' which
#'   specifies the measure of central tendency. Defaults to 'median'.
#'
#' @param ci_algo can be either 'boot' or 'loo' for bootstrap or leave-one-out
#'   methods respectively. Default value is 'boot'.
#'  
#' @param quiet_mode boolean, defaults to TRUE in which case warnings from 
#'  the diversity calculations are not returned. These are typically generated
#'  when calculating the metric SPIE which is undefined when all species are
#'  singletons. 
#'
#' @details
#' 
#' See \code{\link{calc_comm_div}} for more details on the biodiversity indices. 
#' 
#' \strong{Group comparison metric and test}
#' 
#' For each metric group comparison the function computes \code{D_bar}: the
#' average absolute difference between the groups. At the alpha scale the
#' indices are averaged first before computing \code{D_bar}.
#' 
#' Permutation tests are carried out for testing differences of the biodiversity
#' statistics among the groups (Legendre & Legendre 1998). This is accomplished
#' by using \code{D_bar} as the test statistic and random shuffling the group 
#' label across the samples. The p-value indicates how many of the permutations 
#' result in a \code{D_bar} as large as the observed \code{D_bar} value.
#' 
#' @returns A list of class \code{mob_stats} that contains two objects: 
#' 1) \code{comm_div} a data.frame of each diversity metrics for each group
#' at each scale specified, and
#' 2) \code{gorup_tests} a data.frame of the average difference between groups
#' in their diversity metric (\code{D_bar}) with an associated p-value derived 
#' from the permutation test. 
#'
#' @inherit calc_comm_div references
#' 
#' 
#' @import dplyr
#' @importFrom pbapply pbreplicate
#' @importFrom rlang .data
#' @importFrom stats median
#' 
#' @export
#' @examples
#' # tank community analysis
#' data(tank_comm)
#' data(tank_plot_attr)
#' tank_mob <- make_mob_in(tank_comm, tank_plot_attr)
#' tank_stats <- get_mob_stats(tank_mob, 'group', index = c('S', 'S_PIE', 'S_C'),
#'                             n_perm = 5, ci_n_boot = 5)
#' tank_stats    
get_mob_stats <- function(mob_in, group_var, 
                          index = c("N", "S", "S_n", "S_PIE"),
                          effort_samples = NULL, effort_min = 5,
                          extrapolate = TRUE, return_NA = FALSE, 
                          rare_thres = 0.05, 
                          scales = c('alpha', 'gamma', 'beta'),
                          PIE_replace = FALSE, C_target_gamma = NA,
                          n_perm = 199, cl = NULL, 
                          ci = TRUE, ci_cent_stat = 'median', ci_algo = 'boot',
                          ci_n_boot = 1000, quiet_mode = TRUE, ...) {
  EPS <- sqrt(.Machine$double.eps)
  if (n_perm < 1) 
    stop('Set n_perm to a value greater than 1') 
  
  INDICES <- c("N", "S", "S_n", "S_C", "S_asymp", "f_0",
              "pct_rare", "PIE", "S_PIE")
  index <- match.arg(index, INDICES, several.ok = TRUE)
  
  
  groups <- factor(mob_in$env[ , group_var])
  group_levels <- levels(groups) 

  # Get rarefaction level
  samples_N <- rowSums(mob_in$comm) 
  samples_per_group <- table(groups)
  
  if (any(is.null(effort_samples)) | !is.numeric(effort_samples)) {
    N_min_sample <- min(samples_N)
    effort_samples <- N_min_sample
  } else {
    effort_samples <- floor(effort_samples)
  }
  
  if (any(effort_samples < effort_min) & 'S_n' %in% index) {
    warning(paste("The number of individuals for rarefaction analysis is too low and is set to the minimum of", effort_min,"individuals."))
    
    effort_samples[effort_samples < effort_min] <- effort_min
    effort_samples <- unique(effort_samples)
    
    print(paste("Number of individuals for rarefaction:",
                paste(effort_samples, collapse = ", ")))
  }
  

  if(ci) {
    cat('\nComputing confidence intervals\n')
    sample_list <- data.frame(groups, mob_in$comm) |> 
      group_by(groups) |> 
      group_map(~ get_samples(.x, algo = ci_algo, 
                              n_boot = ci_n_boot))
    names(sample_list) <- group_levels
    if (is.na(C_target_gamma) & "S_C" %in% index) {
      # compute C_target_gamma across groups
      C_target_gamma <- min(sapply(sample_list, function(x)
                            sapply(x, calc_C_target)))
    }
    if (quiet_mode)
      dat_div <- suppressWarnings(lapply(sample_list, calc_comm_div_ci, ci_cent_stat,
                        index = index, effort = effort_samples,
                        extrapolate = extrapolate, return_NA = return_NA,
                        rare_thres = rare_thres, scales = scales, 
                        PIE_replace = FALSE, C_target_gamma = C_target_gamma, ...))
    else
      dat_div <- lapply(sample_list, calc_comm_div_ci, ci_cent_stat,
                        index = index, effort = effort_samples,
                        extrapolate = extrapolate, return_NA = return_NA,
                        rare_thres = rare_thres, scales = scales, 
                        PIE_replace = FALSE, C_target_gamma = C_target_gamma, ...)
    dat_div <- bind_rows(dat_div, .id = group_var)
  } else {
    if (is.na(C_target_gamma) & "S_C" %in% index) {
      # compute C_target_gamma across two groups
      grpC <- data.frame(groups, mob_in$comm) |>  
                group_by(groups) |> 
                group_map(~ calc_C_target(.x))
      C_target_gamma <- min(unlist(grpC))
    }
    if (quiet_mode)
      dat_div <- suppressWarnings(data.frame(groups, mob_in$comm) |>  
        group_by(groups) |> 
        group_map(~ calc_comm_div(.x, index = index, effort = effort_samples,
                                  extrapolate = extrapolate, return_NA = return_NA,
                                  rare_thres = rare_thres, scales = scales, avg_alpha = TRUE, 
                                  PIE_replace = FALSE, C_target_gamma = C_target_gamma, ...)))
    else
      dat_div <- data.frame(groups, mob_in$comm) |>  
        group_by(groups) |> 
        group_map(~ calc_comm_div(.x, index = index, effort = effort_samples,
                                  extrapolate = extrapolate, return_NA = return_NA,
                                  rare_thres = rare_thres, scales = scales, avg_alpha = TRUE, 
                                  PIE_replace = FALSE, C_target_gamma = C_target_gamma, ...))
    names(dat_div) <- group_levels
    dat_div <- bind_rows(dat_div, .id = group_var)
  }
  # D_bar is the average difference in a diversity metric between
  # one or more groups. At the alpha scale the diversity metric is 
  # averaged within groups and then the average of the differences 
  # between groups is computed. 
  D_bar <- dat_div |>
    summarise(D_bar = mean(stats::dist(.data$value)), 
              .by = c(scale, index))
  D_obs <- D_bar
  # Significance tests -------------------------------------------------------
  cat('\nComputing permutation tests\n')
  if (quiet_mode)
    D_rand <- suppressWarnings(bind_rows(pbreplicate(n_perm,
                                  data.frame(groups = sample(groups), mob_in$comm) |> 
                                    group_by(groups) |> 
                                    group_map(~ calc_comm_div(.x, index = index,
                                      effort = effort_samples, extrapolate = extrapolate,
                                      return_NA = return_NA, rare_thres = rare_thres,
                                      scales = scales, avg_alpha = TRUE, 
                                      PIE_replace = FALSE,
                                      C_target_gamma = C_target_gamma, ...)) |>
                                    bind_rows(.id = 'id') |>
                                    summarise(D_bar = mean(stats::dist(.data$value)),
                                              .by = c(scale, index)),
                                  simplify = FALSE, cl = cl)))
  else
    D_rand <- bind_rows(pbreplicate(n_perm,
                                    data.frame(groups = sample(groups), mob_in$comm) |> 
                                      group_by(groups) |> 
                                      group_map(~ calc_comm_div(.x, index = index,
                                                                effort = effort_samples, extrapolate = extrapolate,
                                                                return_NA = return_NA, rare_thres = rare_thres,
                                                                scales = scales, avg_alpha = TRUE, 
                                                                PIE_replace = FALSE,
                                                                C_target_gamma = C_target_gamma, ...)) |>
                                      bind_rows(.id = 'id') |>
                                      summarise(D_bar = mean(stats::dist(.data$value)),
                                                .by = c(scale, index)),
                                    simplify = FALSE, cl = cl)) 
  D_tmp <- D_obs |> mutate(D_bar_obs = .data$D_bar, D_bar = NULL)
  D_rand <- left_join(D_rand, D_tmp, by = c("scale", "index"))
  
  perm_tests <- D_rand |> 
    group_by(.data$scale, .data$index) |>
    summarise(p_val = (sum(.data$D_bar >= .data$D_bar_obs - EPS) + 1) /
                (n_perm + 1)) |>
    ungroup()
  perm_tests <- left_join(D_bar, perm_tests, by = c("scale", "index"))
    
  # order output data frames by indices
  dat_div$index <- factor(dat_div$index,
                             levels = c("N",
                                        "S", "beta_S",
                                        "S_n", "beta_S_n",
                                        "S_C", "beta_S_C",
                                        "S_asymp", "beta_S_asympS",
                                        "f_0", "beta_f_0",
                                        "pct_rare", "beta_pct_rare",
                                        "PIE", "beta_PIE",
                                        "S_PIE", "beta_S_PIE"))
  dat_div <- dat_div[order(dat_div$index, dat_div$effort, dat_div[[group_var]]), ]
  dat_div$index <- factor(dat_div$index)
  perm_tests$index <- factor(perm_tests$index,
                          levels = c("N",
                                     "S", "beta_S",
                                     "S_n", "beta_S_n",
                                     "S_C", "beta_S_C",
                                     "S_asymp", "beta_S_asympS",
                                     "f_0", "beta_f_0",
                                     "pct_rare", "beta_pct_rare",
                                     "PIE", "beta_PIE",
                                     "S_PIE", "beta_S_PIE"))
  perm_tests <- perm_tests[order(perm_tests$index), ]
  perm_tests$index <- factor(perm_tests$index)

  # create output
  out <- list(comm_div = dat_div,
              group_tests = perm_tests)
  if ("pct_rare" %in% index)
    out$rare_thres <- rare_thres
  class(out) <- 'mob_stats'
  return(out)
}


#' Plot statistics for a MoB analysis
#' 
#' Plots a \code{mob_stats} object which is produced by the 
#' function \code{\link{get_mob_stats}}. The p-value for each statistic
#' is displayed in the plot title if applicable.
#' 
#' @param mob_stats a \code{mob_stats} object produced by the function
#'   \code{\link{get_mob_stats}}
#'   
#' @param group_var a character string which specifies which variable provides
#'   the groups that the diversity indices are compared across. 
#' 
#' @param group_order Optional vector of group levels the order of 
#'   which defines the order in which the groups are plotted from left-to-right. 
#'   Note that the default in R is to order groups alphabetically.
#' 
#' @param index optional string of indices that are to be plotted. If
#'   left to the default NULL value then all indices are plotted. 
#' 
#' @importFrom rlang .data
#' 
#' @export
#' @examples
#' data(tank_comm)
#' data(tank_plot_attr)
#' tank_mob <- make_mob_in(tank_comm, tank_plot_attr)
#' # for quick results we specify a low number of permutations 
#' # and bootstrap samples these should be increased in real 
#' # analyses. 
#' tank_stats <- get_mob_stats(tank_mob, 'group',
#'                             index = c('N', 'S', 'S_n', 'S_PIE', 'S_C'),
#'                             n_perm = 19, ci_n_boot = 20)
#' p <- plot(tank_stats, 'group')
#' # plot of S
#' p$S
#' # change the order of the groups on the plot
#' p <- plot(tank_stats, 'group', group_order = c('low', 'high'))
#' p$S
plot.mob_stats <- function(mob_stats, group_var, group_order = NULL, index = NULL) {
  all_stats <- dplyr::left_join(mob_stats$comm_div, mob_stats$group_tests,
                                by=c('scale', 'index'))
  if (group_var != "group") {
    # rename the group_var column to group
    # to avoid ggplot issues when plotting N 
    name_i <- which(names(all_stats) == group_var)
    names(all_stats)[name_i] <- "group"
  }
  if (!is.null(index)) {
    indices_to_plot <- index
    if ('beta' %in% all_stats$scale)
      indices_to_plot <- c(index, paste('beta_', index, sep=''))
    all_stats <- subset(all_stats, index %in% indices_to_plot)
  }
  # if users would like groups to be plotted in a different
  # order then accommodate that. 
  if (is.null(group_order)) { 
    all_stats$group <- factor(all_stats$group)
  } else {
    all_stats$group <- factor(all_stats$group, levels = group_order)
  }
  indices <- sub("beta_", "", all_stats$index)
  uni_index <- unique(indices)
  p <- vector('list', length = length(uni_index))
  names(p) <- uni_index
  
  for (i in seq_along(uni_index)) {
    tmp <- subset(all_stats, indices == uni_index[i])
    # setup y-axis label
    y_label <- switch(uni_index[i],
                      "N" = expression("Abundance (" * italic(N) * ")"),
                      "S" = expression('Richness (' * italic(S) * ')'),
                      "S_C" = expression('Richness (' * italic(S)[C] *')'),
                      "S_n" = expression('Richness (' * italic(S)[n]*')'), 
                      "S_asymp" = expression('Asympotic richness (' *
                                               italic(S[asymp]) * ')'),
                      "pct_rare" = paste("% of species in lower",
                                         ifelse(is.character(mob_stats$rare_thres),
                                                mob_stats$rare_thres,
                                                mob_stats$rare_thres * 100),
                                         "% of abundance"),
                      "f_0" = expression('Undetected richness (' * italic(f)[0] * ')'),
                      "S_PIE" = expression('ENS of PIE (' * italic(S)[PIE] * ')'))
    labs <- with(tmp, paste(scale, "\nDbar = ", round(D_bar, 1),
                            ", p = ", round(p_val, 3), sep=""))
    names(labs) <- tmp$scale
    p[[i]] <- ggplot(tmp, aes(x = group, y = value)) +
      geom_point() + 
      geom_errorbar(aes(ymin = lo_value, ymax = hi_value),
                    width = 0.2) +
      ylab(y_label) + 
      facet_wrap(~ scale, scales = 'free', 
                 labeller = as_labeller(labs)) + 
      xlab(group_var) + 
      theme_classic()
    
  }
  return(p)
}

#' Panel function for alpha-scale results
#' @importFrom graphics boxplot mtext
#' @keywords internal
#' @noRd
samples_panel1 = function(sample_dat, col, ylab = "",
                          main = expression(alpha * "-scale"), 
                          cex.axis=1.2, ...) {
  #label = substitute(paste(italic(bar(D)), ' = ', D_bar, ', ', italic(p), 
  #                         ' = ', p_val),
  #                   list(D_bar = round(samples_tests$D_bar, 2),
  #                        p_val = round(samples_tests$p_val, 3)))
  boxplot(value ~ group, data = sample_dat, main = main,
          ylab =  ylab, col = col, cex.axis = cex.axis, cex.main = 1.5,
          frame.plot = TRUE, ...)
  #groups = levels(sample_dat$group)
  #axis(side=1, at=1:length(groups), labels=groups, tick=FALSE,
  #      cex.axis=cex.axis)
  #mtext(label, side = 3, line = 0)  
}

#' Panel function for gamma-scale results
#' @importFrom graphics boxplot points mtext 
#' @keywords internal
#' @noRd
groups_panel1 = function(group_dat, col, ylab = "",
                         main = expression(gamma * "-scale"),
                         cex.axis=1.2, ...) {
  #label = substitute(paste(italic(bar(D)), ' = ', D_bar, ', ', italic(p), 
  #                         ' = ', p_val), 
  #                   list(D_bar = round(tests$D_bar, 2),
  #                        p_val = round(tests$p_val, 3)))
  boxplot(value ~ group, data = group_dat, main = main,
          ylab = ylab, boxwex = 0, 
          ylim = c(0, 1.1 * max(group_dat$value, na.rm = TRUE)),
          col = col, cex.axis = cex.axis, cex.main = 1.5, frame.plot = TRUE,
          ...)
  groups = levels(group_dat$group)
  points(value ~ group, data = group_dat, pch = 8, cex = 1.5, lwd = 2,
         col = col, ...)
  #mtext(label, side = 3, line = 0)
}

#' Panel function for gamma-scale results with confidence intervals
#' @importFrom plotrix plotCI
#' @keywords internal
#' @noRd
groups_panel2 = function(group_dat, col, ylab = "", 
                         main = expression(gamma * "-scale"),
                         cex.axis=1.2, ...) {
  boxplot(median ~ group, data = group_dat, main = main,
          ylab = ylab, boxwex = 0, ylim = c(0, 1.1 * max(group_dat$upper)),
          col = col, cex.axis = cex.axis, cex.main = 1.5, frame.plot = TRUE, 
          ...)
  groups = levels(group_dat$group)
  plotCI(1:nrow(group_dat), group_dat$median, li = group_dat$lower,
         ui = group_dat$upper, add = TRUE, pch = 19, cex = 1.5,
         sfrac = 0.02, col = col, ...)
}

#' Plot alpha-, beta-, and gamma-scale biodiversity statistics for a MoB analysis
#' 
#' Plots the community diversity metrics from produced by the function 
#' \code{calc_comm_div}. The p-value for each statistic
#' is displayed in the plot title if applicable.
#' 
#' The user may specify which results to plot or simply to plot 
#' all the results. 
#' 
#' @param comm_div a table that is output by \code{calc_comm_div} that has the
#' sample (alpha) and group (gamma) level statistics
#' 
#' @param index The biodiversity statistics that should be plotted.
#' See \code{\link{calc_comm_div}} for information on the indices. By default there
#' is one figure for each index, with panels for alpha- and gamma-scale results
#' as well as for beta-diversity when applicable. 
#' 
#' @param multi_panel A logical variable. If \code{multi_panel = TRUE} then a 
#' multipanel plot is produced, which shows observed, rarefied, and asymptotic 
#' species richness and S_PIE at the alpha- and gamma-scale.
#' This set of variables conveys a comprehensive picture of the underlying 
#' biodiversity changes. 
#' 
#' @param col a vector of colors for the groups, set to NA if no color is
#' preferred
#' 
#' @param cex.axis The magnification to be used for axis annotation relative to
#' the current setting of cex. Defaults to 1.2. 
#' 
#' @param ... additional arguments to provide to \code{boxplot}, \code{points},
#'   and confidence interval functions
#' 
#' @author Felix May, Xiao Xiao, and Dan McGlinn 
#' 
#' @importFrom rlang .data 
#' 
#' @export
#' 
#' @examples 
#' library(dplyr)
#' data(tank_comm)
#' data(tank_plot_attr)
#' indices <- c('N', 'S', 'S_C', 'S_n', 'S_PIE')
#' # the grouping variable must be called group in the 
#' # tank_div data.frame
#' tank_div <- tibble(tank_comm) |> 
#'   group_by(group = tank_plot_attr$group) |> 
#'   group_modify(~ calc_comm_div(.x, index = indices, effort = 5,
#'                                extrapolate = TRUE))
#' # plot the community metrics                                 
#' plot_comm_div(tank_div, index = "S")
#' plot_comm_div(tank_div, index = "S_n")
#' # or plot all of the indices at once with
#' plot_comm_div(tank_div)
plot_comm_div = function(comm_div, index = NULL, multi_panel = FALSE, 
                          col = c("#FFB3B5", "#78D3EC", "#6BDABD", "#C5C0FE",
                                  "#E2C288", "#F7B0E6", "#AAD28C"), 
                          cex.axis=1.2, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  # default colors derived with colorspace::rainbow_hcl(5, c=60, l=80)
  if (any(is.na(col)) & length(col) == 1) 
    col_groups = 1
  else
    col_groups = col
  if (is.null(index))
    index = as.character(unique(comm_div$index))
  INDICES = c("N", "S", "S_C", "S_n", "S_asymp", "f_0", 
              "pct_rare", "PIE", "S_PIE")
  if (multi_panel) 
    index = c("S","S_n","pct_rare","S_PIE")
  index = match.arg(index, INDICES, several.ok = TRUE)
  
  var_names = unique(comm_div$index)
  var_names2 = var_names[var_names != "beta_S" & var_names != "beta_S_PIE"]
  
  index_match = intersect(index, var_names)
  if (length(index_match) == 0)
    stop(paste("The indices", paste(index, collapse = ", "), 
               "are missing in the input. Please choose other indices or re-run get_comm_div with the indices of interest."))
  
  index_missing = setdiff(index, var_names2)
  if (length(index_missing) > 0)
    warning(paste("The indices", paste(index, collapse = ", "), 
                  "are missing in the input and cannot be plotted."))
  
  if ("S_n" %in% index_match) {
    S_n_samples = filter(comm_div, index == "S_n" & scale == "alpha")
    S_n_groups = filter(comm_div, index == "S_n" & scale == "gamma")
    S_n_len = length(unique(S_n_samples$effort))
  } else {
    S_n_len = 0
  }
  
  if (multi_panel) {
    n_rows = 3 + S_n_len
    op = par(mfrow = c(n_rows,3), las = 1, cex.lab = 1.4, 
             oma = c(0, 1, 0, 0), xpd = NA)
  } 
  
  for (var in index_match) {
    
    if (var == "N") {
      
      if (!multi_panel)
        op = par(mfrow = c(1, 3), cex.lab = 1.6,
                 oma = c(0, 2, 0, 0), mar = c(4, 3, 5, 1), xpd = NA)
      
      #op = par(mfrow = c(1, 2), las = 1, cex.lab = 1.3,
      #         oma = c(0, 2, 0, 0), mar = c(4, 3, 5, 1),
      #         xpd = NA)
      
      
      y_label = switch(var,
                       "N" = expression("Abundance (" * italic(N) * ")"),
                       "PvarIE" = "PIE")
      
      dat_samples = filter(comm_div, index == var & scale == 'alpha')
      #dat_tests = filter(comm_div$samples_tests, index == var)
      samples_panel1(dat_samples, ylab = y_label,
                     main = expression(alpha * "-scale"), col = col,
                     cex.axis = cex.axis, ...)
      
      # insert blank space b/c non beta plot for these stats
      plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '')
      
      dat_groups = filter(comm_div, index == var & scale == 'gamma')
      
      #if (is.null(comm_div$groups_tests)) {
        # then plot median and quantiles
      #  groups_panel2(dat_groups, col = col_groups, cex.axis = cex.axis, ...) 
      #} else {  
      #  tests = filter(comm_div$groups_tests, index == var)
        groups_panel1(dat_groups, col = col_groups, cex.axis = cex.axis, ...) 
      #}
    }
    
    if (var %in% c("S", "S_C", "S_asymp", "pct_rare", "f0", "PIE", "S_PIE")) {
      
      if (!multi_panel)
        op = par(mfrow = c(1, 3), cex.lab = 1.6,
                 oma = c(0, 2, 0, 0), mar = c(4, 3, 5, 1), xpd = NA)
      
      if (multi_panel) {
        if (var == "f_0") 
          par(fig = c(0, 0.33, 1 / n_rows, 2 / n_rows), new = TRUE)
        if (var == "S_PIE") 
          par(fig = c(0, 0.33, 0         , 1 / n_rows), new = TRUE)
      }
      
      y_label = switch(var,
                       "S" = expression('Richness (' * italic(S) * ')'),
                       "S_C" = expression('Richness (' * italic(S)[C] *')'),
                       "S_asymp" = expression('Asympotic richness (' *
                                                italic(S[asymp]) * ')'),
                       "pct_rare" = paste("% of species in lower",
                                          ifelse(is.character(comm_div$rare_thres),
                                                 comm_div$rare_thres,
                                                 comm_div$rare_thres * 100),
                                          "% of abundance"),
                       "f_0" = expression('Undetected richness (' * italic(f)[0] * ')'),
                       "S_PIE" = expression('ENS of PIE (' * italic(S)[PIE] * ')'))
      
      dat_samples = filter(comm_div, index == var & scale == 'alpha')
      #dat_tests = filter(comm_div$samples_tests, index == var)
      samples_panel1(dat_samples, ylab =  y_label,
                     main = expression(alpha * "-scale"), col = col, 
                     cex.axis = cex.axis, ...)
      
      if (multi_panel) {
        if (var == "pct_rare") 
          par(fig = c(0.33, 0.67, 1/n_rows, 2/n_rows), new = TRUE)
        if (var == "S_PIE") 
          par(fig = c(0.33, 0.67, 0       , 1/n_rows), new = TRUE)
      }
      
      beta_var = paste("beta", var, sep = "_")
      dat_samples = filter(comm_div, index == beta_var & scale == 'beta')
      #dat_tests = filter(comm_div$samples_tests, index == beta_var)
      samples_panel1(dat_samples, ylab =  "",
                     main = expression(beta * "-diversity (=" *
                                         gamma / alpha * ")"),
                     col = col, cex.axis = cex.axis, ...)
      
      if (multi_panel) {
        if (var == "pct_rare")
          par(fig = c(0.67, 1.0, 1 / n_rows, 2 / n_rows), new = TRUE)
        if (var == "S_PIE")
          par(fig = c(0.67, 1.0, 0         , 1 / n_rows), new = TRUE)
      }
      
      dat_groups = filter(comm_div, index == var & scale == 'gamma')
      
      #if (is.null(comm_div$groups_tests)) {
      #  groups_panel2(dat_groups, col = col_groups, ...) 
      #} else {
      #  tests = filter(comm_div$groups_tests, index == var)
        groups_panel1(dat_groups, col = col_groups, cex.axis = cex.axis,
                      ...) 
      #}
    }    
    
    if (var == "S_n") {
      
      if (!multi_panel) {
        op = par(mfrow = c(S_n_len, 3), las = 1, cex.lab = 1.6,
                 oma = c(0, 2, 0, 0), mar = c(4, 3, 5, 1), xpd = NA)
      }
      
      y_label = expression('Rarefied richness (' * italic(S[n]) * ')')
      
      effort_samples = unique(S_n_samples$effort)
      effort_groups = unique(S_n_groups$effort)
      
      for (i in seq_along(effort_samples)) {
        
        if (multi_panel)
          par(fig = c(0, 0.33, (1 + i) / n_rows, (2 + i) / n_rows),
              new = TRUE)
        
        dat_samples = filter(S_n_samples, .data$effort == effort_samples[i])
        #dat_tests = filter(comm_div$samples_tests, 
        #                   index == var & .data$effort == effort_samples[i])
        
        fig_title = substitute(paste(alpha, "-scale, n = ", n),
                               list(n = effort_samples[i]))
        
        samples_panel1(dat_samples, 
                       ylab = y_label,
                       main = '', col = col, cex.axis = cex.axis,
                       ...)
        par(new = TRUE)
        plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '',
             main = fig_title, cex.main = 1.5)
        
        if (multi_panel)
          par(fig = c(0.33, 0.67, (1 + i) / n_rows, (2 + i) / n_rows),
              new = TRUE)
        
        beta_var = paste("beta", var, sep = "_")
        dat_samples = filter(comm_div, index == beta_var & scale == 'beta')
        #dat_tests = filter(comm_div$samples_tests,
        #                   index == beta_var & .data$effort == effort_samples[i])
        n = effort_samples[i]
        fig_title = substitute(paste(beta, "-diversity (=", gamma / alpha,
                                     "), n = ", n), 
                               list(n = effort_samples[i]))
        samples_panel1(dat_samples, main = '',
                       ylab = "", col = col, cex.axis = cex.axis, ...)
        
        par(new = TRUE)
        plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '',
             main = fig_title, cex.main = 1.5)
        
        if (multi_panel)
          par(fig = c(0.67, 1.0, (1 + i) / n_rows, (2 + i) / n_rows),
              new = TRUE)
        
        dat_groups = filter(S_n_groups, .data$effort == effort_groups[i])
        fig_title = substitute(paste(gamma, "-scale, n = ", n),
                               list(n = effort_groups[i]))            
        #if (is.null(comm_div$groups_test)) {
        #  groups_panel2(dat_groups, main = '', col = col_groups,
        #                ...) 
        #  par(new = TRUE)
        #  plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '',
        #       main = fig_title, cex.main = 1.5)
        #} else {
        #  tests = filter(comm_div$groups_tests, 
        #                 index == var & .data$effort == effort_groups[i])
          groups_panel1(dat_groups, ylab = "", main = '',
                        col = col_groups, cex.axis = cex.axis, ...)
          par(new = TRUE)
          plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '',
               main = fig_title, cex.main = 1.5)
        #}
      }
      y_coords = (S_n_len:0) / S_n_len
    }
  }
  par(op)
}
