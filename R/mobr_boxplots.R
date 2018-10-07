#' Estimation of species richness
#' 
#' \code{calc_chao1} estimates the number of species at the asymptote
#' (\code{S_asymp}) of the species accumulation curve based on the methods
#' proposed in Chao (1984, 1987, 2005). 
#' 
#' This function is a trimmed version of \href{https://github.com/JohnsonHsieh/iNEXT}{\code{iNext::ChaoRichess}}.
#' T. C. Hsieh, K. H. Ma and Anne Chao are the original authors of the
#' \code{iNEXT} package. 
#' 
#' @param x a vector of species abundances or a site-by-species matrix
#' 
#' @return a vector of species richness estimates
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
#' Chao, A. (2005) Species estimation and applications. Pages 7907–7916 in
#' N. Balakrishnan, C. B. Read, and B. Vidakovic, editors. Encyclopedia of
#' statistical sciences. Second edition, volume 12. Wiley, New York, New York,
#' USA.
#' 
#' @export
calc_chao1 = function(x) {
    if (!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
        stop("invalid data structure")
    if (is.matrix(x) | is.data.frame(x)) {
        S_Chao1= apply(x, 1, calc_chao1)
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
#'  \code{calc_PIE} returns the probability of interspecific  encounter (PIE)
#'  which is also known as Simpson's evenness index and Gini-Simpson index. For \code{ENS=TRUE},
#'  PIE will be converted to an asymptotic effective number of species (S_PIE).
#' 
#' The formula of Hurlbert (1971) is used to calculate PIE:
#' 
#' \eqn{PIE = N /(N - 1) * (1 - p_i^2)}
#' 
#' where N is the total number of individuals and \eqn{p_i} is the relative abundance
#' of species i. This formulation uses sampling without replacement and it is
#' sometimes referred to as the bias corrected formulation of PIE.
#' 
#' For \code{ENS = TRUE}, S_PIE will be returned which represents the species richness of
#' a hypothetical community with equally-abundant species and infinitely many individuals
#' corresponding to the observed value of PIE. It is computed as
#' \eqn{S_PIE = 1 /(1 - PIE)}, which is equal to the
#' asymptotic estimator for Hill numbers of diversity order 2 provided by Chao et al (2014).
#' Note that S_PIE is undefined for communities with exactly one individual per species.
#'  
#' The code in this function borrows heavily from the function vegan::diversity()
#' but computes a different quantity. The function vegan::diversity() computes
#' PIE when sampling with replacement is assumed. The difference between the two 
#' formulations will decrease as N becomes large. Jari Oksanen and Bob O'Hara are
#' the original authors of the function vegan::diversity().
#' 
#' @inheritParams rarefaction
#' @param ENS boolean that determines if the effective number of species should
#' be returned or the raw PIE value. Defaults to FALSE
#'
#' @author Dan McGlinn, Thore Engel
#' 
#' @references 
#' Hurlbert, S. H. (1971) The nonconcept of species diversity: a critique and
#'  alternative parameters. Ecology 52, 577–586.
#'  
#' Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K., & Ellison, A. M. (2014).
#'  Rarefaction and extrapolation with Hill numbers: A framework for sampling and estimation in species diversity studies.
#'  Ecological Monographs 84(1), 45–67.
#'
#' @export
#' @examples 
#' data(inv_comm)
#' calc_PIE(inv_comm)
#' calc_PIE(inv_comm, ENS=TRUE)
calc_PIE = function(x, ENS=FALSE) {
    if (class(x) == 'mob_in') {
        x = x_mob_in$comm
    }
    x = drop(as.matrix(x))
    if (any(x < 0, na.rm = TRUE)) 
        stop("input data must be non-negative")
    if (length(dim(x)) > 1) {
        total = apply(x, 1, sum)
        S = apply(x, 1, function(x) return(sum(x > 0)))
        x = sweep(x, 1, total, "/")
    } else {
        total = sum(x)
        S = sum(x > 0)
        x = x / total
    }
    x = x * x
    if (length(dim(x)) > 1) {
        H = rowSums(x, na.rm = TRUE)
    } else {
        H = sum(x, na.rm = TRUE)
        }
    # calculate PIE without replacement (for total >= 2)
    H = ifelse(total < 2, NA, (total / (total - 1) * (1 - H)))
    if (ENS) {
        # convert to effective number of species (except for PIE == 1)
        H = ifelse(H==1| S == total, NA, (1/ (1-H)))
    }     
    return(H)
}

# generate a single bootstrap sample of gamma-scale biodiversity indices
boot_sample_groups = function(abund_mat, index, effort, extrapolate, return_NA,
                              rare_thres) {
    # sample rows and calculate abundance vector
    sample_dat = by(abund_mat, INDICES = abund_mat$group_id, FUN = sample_frac,
                    replace = T)
    class(sample_dat) = "list"
    sample_dat = bind_rows(sample_dat)
   
    # abundance distribution pooled in groups
    abund_group = aggregate(sample_dat[ , -1], by = list(sample_dat[ , 1]),
                            FUN = "sum")
   
    dat_groups = calc_biodiv(abund_mat = abund_group[ , -1],
                             groups = abund_group[ , 1],
                             index = index, effort = effort,
                             extrapolate = extrapolate,
                             return_NA = return_NA, 
                             rare_thres = rare_thres)
    return(dat_groups)
}


#' Calculate biodiversity statistics from sites by species table.
#' 
#' @param abund_mat Sites by species table with species abundances
#' in the respective cells
#' 
#' @param groups Vector with group labels for the sites. The length
#' of the vector has to correspond to the number of rows of the
#' sites by species table.
#' 
#' @param index The calculated biodiversity indices. The options are
#' \itemize{
#'    \item \code{N} ... Number of individuals (total abundance)
#'    \item \code{S} ... Number of species
#'    \item \code{S_n} ... Rarefied or extrapolated number of species for n individuals
#'    \item \code{S_asymp} ... Estimated asymptotic species richness
#'    \item \code{f_0} ... Estimated number of undetected species 
#'    \item \code{pct_rare} ... The percent of species with abundances below \code{rare_thres}
#'    \item \code{PIE} ... Hurlbert's PIE (Probability of Interspecific Encounter)
#'    \item \code{S_PIE} ... Effective number of species based on PIE
#' }
#' 
#' See the documentation of \code{\link{get_mob_stats}} for further details on the
#' biodiversity indices.
#' 
#' @param effort The standardized number of individuals used for the 
#'   calculation of rarefied species richness. This can a be
#'   single value or an integer vector. 
#'   
#' @param extrapolate boolean which specifies if richness should be extrapolated
#'   when effort is larger than the number of individuals using the chao1
#'   method.
#'   
#' @param rare_thres The threshold that determines how pct_rare is computed.
#'   It can range from (0, 1] and defaults to 0.05 which specifies that any 
#'   species with less than or equal to 5% of the total abundance in a sample is
#'   considered rare. It can also be specified as "N/S" which results in using
#'   average abundance as the threshold which McGill (2011) found to have the 
#'   best small sample behavior. 
#'   
#' @details This function is primarily intended as auxiliary function used in
#' \code{\link{get_mob_stats}}, but can be also used directly for data exploration.
#' 
#' 
#' @return A dataframe with four columns:
#' \itemize{
#'    \item \code{group} ... Group label for sites
#'    \item \code{index} ... Name of the biodiversity index
#'    \item \code{effort} ... Sampling effort for rarefeid richness 
#'    (NA for the other indices)
#'    \item \code{value} ... Value of the biodiversity index
#' }
#'   
#' @author Felix May and Dan McGlinn
#' 
#' @references 
#' 
#' McGill, B. J. 2011. Species abundance distributions. Pages 105–122 Biological
#' Diversity: Frontiers in Measurement and Assessment, eds. A.E. Magurran and
#' B.J. McGill.
#'
#' @export
calc_biodiv = function(abund_mat, groups, index, effort, extrapolate, return_NA,
                       rare_thres) {
    out = expand.grid(group = groups,
                      index = index[index != "S_n"],
                      effort = NA, value = NA)
    N = rowSums(abund_mat) 
   
    # Number of individuals -----------------------------------------------------
    if (any(index == "N")) {
        out$value[out$index == "N"] = N
    } 
   
    # Number of species ---------------------------------------------------------
    if (any(index == "S")) {
        out$value[out$index == "S"] = rowSums(abund_mat > 0)
    }  
   
    # Rarefied richness ---------------------------------------------------------
    if (any(index == "S_n")) {
      
        dat_S_n = expand.grid(group = groups,
                                 index = 'S_n',
                                 effort = effort, value = NA)
        out = rbind(out, dat_S_n)
        S_n  = apply(abund_mat, 1, rarefaction, method = 'indiv', effort = effort,
                     extrapolate = extrapolate, return_NA = return_NA,
                     quiet_mode = TRUE)
        out$value[out$index == "S_n"] = as.numeric(t(S_n))
    } 
    
    # Percent rare ------------------------------------------------------------
    if (any(index == "pct_rare")) {
        pct_rare = function(x, rare_thres) {
            S = sum(x > 0)
            if (S > 0) {
                N = sum(x)
                if (rare_thres == "N/S") {
                    rare_thres = N / S
                    100 * (sum(x[x > 0] <= rare_thres) / S)
                } else 
                    100 * (sum(x[x > 0] <= (rare_thres * N)) / S)
            } else
              0
        }
        out$value[out$index == "pct_rare"] = apply(abund_mat, 1, pct_rare, rare_thres)
    }
 
    # Probability of Interspecific Encounter (PIE)-------------------------------
    if (any(index == "PIE")) {
        # Hurlbert's PIE can only be calculated for two or more individuals
        plots_n01 = N <= 1 
        if (sum(plots_n01) == 1) {
            warning(paste("There is", sum(plots_n01), "plot with less than two individuals.
                          This is removed for the calculation of PIE."))
         
        }
        if (sum(plots_n01) > 1) {
            warning(paste("There are", sum(plots_n01), "plots with less than two individuals.
                          These are removed for the calculation of PIE."))
         
        }
        out$value[out$index == "PIE"] = calc_PIE(abund_mat)
    }
   
    # Effective number of species based on PIE ----------------------------------
    if (any(index == "S_PIE")) {
        plots_n01 = N <= 1
        S_PIE = calc_PIE(abund_mat, ENS=TRUE)
        S_PIE[plots_n01] = NA
        out$value[out$index == "S_PIE"] = S_PIE
    }
   
    # Asymptotic estimates species richness -------------------------------------
    if (any(index == "S_asymp" | index == "f_0")) {
        S_asymp = try(calc_chao1(abund_mat))
        if (class(S_asymp) == "try_error") 
            warning("The Chao richness estimator cannot be calculated for all groups.")
        else 
            S_asymp[!is.finite(S_asymp)] = NA
        S = apply(abund_mat, 1, function(x) sum(x > 0))
        if (any(index == "S_asymp"))
            out$value[out$index == "S_asymp"] = S_asymp
        if (any(index == "f_0"))
            out$value[out$index == "f_0"] = S_asymp - S
    }    
    return(out)
}

#Get F statistics from diversity indices and grouping vector
get_F_values = function(div_dat, permute = F) {
    if (permute)
        div_dat = div_dat %>%
                  group_by(index, effort) %>%
                  mutate(group = sample(group))
    
    models = div_dat %>%
             group_by(index, effort) %>%
             do(mod = try(lm(value ~ group, data = .), silent=TRUE))
   
    models = models %>% 
             mutate(F_val = ifelse(class(mod) == 'lm', anova(mod)$F[1], NA)) %>%
             ungroup()
   
    
    models$F_val[!is.finite(models$F_val)] = NA
   
    return(models)
}

# Get gamma-scale differences 
get_group_delta = function(abund_mat, group_id, index, effort, extrapolate,
                           return_NA, rare_thres, permute = F) {
    if (permute)
        group_id = sample(group_id)
   
    abund_group = aggregate(abund_mat, by = list(group_id), FUN = "sum")
   
    dat_groups = calc_biodiv(abund_mat = abund_group[ , -1],
                             groups = abund_group[ , 1],
                             index = index, effort = effort, 
                             extrapolate, return_NA, rare_thres = rare_thres)
    delta_groups = dat_groups %>%
                   group_by(index, effort) %>%
                   summarise(delta = mean(dist(value))) %>%
                   ungroup()
   
    return(delta_groups)
}

#' Calculate sample based and group based biodiversity statistics.
#' @inheritParams get_delta_stats
#' 
#' @param group_var String that specifies which field in \code{mob_in$env} the
#'   data should be grouped by
#' 
#' @param index The calculated biodiversity indices. The options are
#' \itemize{
#'    \item \code{N} ... Number of individuals (total abundance)
#'    \item \code{S} ... Number of species
#'    \item \code{S_n} ... Rarefied or extrapolated number of species for n individuals
#'    \item \code{S_asymp} ... Estimated asymptotic species richness
#'    \item \code{f_0} ... Estimated number of undetected species 
#'    \item \code{pct_rare} ... The percent of rare species as defined by \code{rare_thres}
#'    \item \code{PIE} ... Hurlbert's PIE (Probability of Interspecific Encounter)
#'    \item \code{S_PIE} ... Effective number of species based on PIE
#'    
#' }
#'   If index is not specified then N, S, S_n, pct_rare, and S_PIE are computed
#'   by default. See \emph{Details} for additional information on the
#'   biodiversity statistics.
#' 
#' @param effort_samples The standardized number of individuals used for the 
#'   calculation of rarefied species richness at the alpha-scale. This can a be
#'   single value or an integer vector. As default the minimum number of
#'   individuals found across the samples is used, when this is not smaller than
#'   \code{effort_min}.
#'   
#' @param effort_min The minimum number of individuals considered for the 
#'   calculation of rarefied richness (Default value of 5). Samples with less
#'   individuals then \code{effort_min} are excluded from the analysis with a
#'   warning. Accordingly, when \code{effort_samples} is set by the user it has
#'   to be higher than \code{effort_min}.
#'  
#' @param extrapolate extrapolate	boolean which specifies if richness should be
#'   extrapolated when \code{effort_samples} is larger than the number of
#'   individuals using the chao1 method. Defaults to TRUE. 
#'   
#' @param return_NA boolean defaults to FALSE in which the rarefaction function
#'   returns the observed S when \code{effort} is larger than the number of
#'   individuals . If set to TRUE then NA is returned. Note that this argument
#'   is only relevant when \code{extrapolate = FALSE}.
#'     
#' @param rare_thres The threshold that determines how pct_rare is computed.
#'   It can range from (0, 1] and defaults to 0.05 which specifies that any 
#'   species with less than or equal to 5% of the total abundance in a sample is
#'   considered rare. It can also be specified as "N/S" which results in using
#'   average abundance as the threshold which McGill (2011) found to have the 
#'   best small sample behavior. 
#'   
#' @param n_perm The number of permutations to use for testing for treatment
#'   effects.
#'   
#' @param boot_groups Use bootstrap resampling within groups to derive
#'   gamma-scale confidence intervals for all biodiversity indices. Default is
#'   \code{FALSE}. See \emph{Details} for information on the bootstrap approach.
#'   
#' @param conf_level Confidence level used for the calculation of gamma-scale 
#'   bootstrapped confidence intervals. Only used when \code{boot_groups =
#'   TRUE}.
#'   
#' @inheritParams pbapply::pbreplicate
#'
#' @details 
#' 
#' \strong{BIODIVERSITY INDICES}
#' 
#' \strong{S_n: Rarefied species richness} is the expected number of species, given a
#' defined number of sampled individuals (n) (Gotelli & Colwell 2001). Rarefied
#' richness at the alpha-scale is calculated for the values provided in 
#' \code{effort_samples} as long as these values are not smaller than the 
#' user-defined minimum value \code{effort_min}. In this case the minimum value 
#' is used and samples with less individuals are discarded. When no values for
#' \code{effort_samples} are provided the observed minimum number of individuals
#' of the samples is used, which is the standard in rarefaction analysis
#' (Gotelli & Colwell 2001). Because the number of individuals is expected to
#' scale linearly with sample area or effort, at the gamma-scale the number of
#' individuals for rarefaction is calculated as the minimum number of samples
#' within groups multiplied by \code{effort_samples}. For example, when there are 10
#' samples within each group, \code{effort_groups} equals \code{10 *
#' effort_samples}. If n is larger than the number of individuals in sample and
#' \code{extrapolate = TRUE} then the Chao1 (Chao 1984, Chao 1987) method is
#' used to extrapolate the rarefaction curve.
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
#' estimator (Chao 1984, Chao 1987) see \code{\link{calc_chao1}}. Note: this metric
#' is typically highly correlated with S (McGill 2011).
#'  
#' \strong{f_0: Undetected species richness} is the number of undetected species
#' or the number of species observed 0 times which is an indicator of the degree
#' of rarity in the community. If there is a greater rarity then f_0 is expected
#' to increase. This metric is calculated as \code{S_asymp - S}. This metric is less 
#' correlated with S than the raw \code{S_asymp} metric. 
#' 
#' \strong{PIE: Probability of intraspecific encounter} represents the probability that two randomly drawn individuals 
#' belong to the same species. Here we use the definition of Hurlbert (1971),
#' which considers sampling without replacement. PIE is closely related to the
#' well-known Simpson diversity index, but the latter assumes sampling with
#' replacement.
#' 
#' \strong{S_PIE: Effective number of species for PIE} represents the effective number of species derived from the
#' PIE. It is calculated using the asymptotic estimator for Hill numbers of diversity order 2 (Chao et al, 2014).
#' S_PIE represents the species richness of a hypothetical community with equally-abundant species
#' and infinitely many individuals corresponding to the same value of PIE as the real community.
#' An intuitive interpretation of S_PIE is that it corresponds to the number of
#' dominant (highly abundant) species in the species pool.
#' 
#' For species richness \code{S}, rarefied richness \code{S_n}, undetected
#' richness \code{f_0}, and the Effective Number of Species \code{S_PIE} we also
#' calculate beta-diversity using multiplicative partitioning (Whittaker 1972,
#' Jost 2007). That means for these indices we estimate beta-diversity as the
#' ratio of gamma-diversity (total diversity across all plots) divided by
#' alpha-diversity (i.e., average plot diversity).
#' 
#' \strong{PERMUTATION TESTS AND BOOTSTRAP}
#' 
#' For both the alpha and gamma scale analyses we summarize effect size in each
#' biodiversity index by computing \code{D_bar}: the average absolute difference
#' between the groups. At the alpha scale the indices are averaged first before
#' computing \code{D_bar}.
#'
#' We used permutation tests for testing differences of the biodiversity
#' statistics among the groups (Legendre & Legendre 1998). At the alpha-scale,
#' one-way ANOVA (i.e. F-test) is implemented by shuffling treatment group
#' labels across samples. The test statistic for this test is the F-statistic
#' which is a pivotal statistic (Legendre & Legendre 1998). At the gamma-scale
#' we carried out the permutation test by shuffling the treatment group labels
#' and using \code{D_bar} as the test statistic. We could not use the
#' F-statistic as the test statistic at the gamma scale because at this scale
#' there are no replicates and therefore the F-statistic is undefined.
#' 
#' A bootstrap approach can be used to also test differences at the gamma-scale.
#' When \code{boot_groups = TRUE} instead of the gamma-scale permutation test,
#' there will be resampling of samples within groups to derive gamma-scale
#' confidence intervals for all biodiversity indices. The function output
#' includes lower and upper confidence bounds and the median of the bootstrap
#' samples. Please note that for the richness indices sampling with replacement
#' corresponds to rarefaction to ca. 2/3 of the individuals, because the same
#' samples occur several times in the resampled data sets.
#' 
#' 
#' @return A list of class \code{mob_stats} that contains alpha-scale and 
#'   gamma-scale biodiversity statistics, as well as the p-values for
#'   permutation tests at both scales.
#'   
#'   When \code{boot_groups = TRUE} there are no p-values at the gamma-scale.
#'   Instead there is lower bound, median, and upper bound for each biodiversity
#'   index derived from the bootstrap within groups.
#
#' @author Felix May and Dan McGlinn
#' 
#' @references 
#' 
#' Chiu, C.-H., Wang, Y.-T., Walther, B.A. & Chao, A. (2014) An improved
#' nonparametric lower bound of species richness via a modified good-turing
#' frequency formula. Biometrics, 70, 671-682.
#' 
#' Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures
#' and pitfalls in the measurement and comparison of species richness. Ecology
#' letters, 4, 379-391.
#' 
#' Hurlbert, S.H. (1971) The Nonconcept of Species Diversity: A Critique and
#' Alternative Parameters. Ecology, 52, 577-586.
#' 
#' Jost, L. (2006) Entropy and diversity. Oikos, 113, 363-375.
#' 
#' Jost, L. (2007) Partitioning Diversity into Independent Alpha and Beta
#' Components. Ecology, 88, 2427–2439.
#' 
#' Legendre, P. & Legendre, L.F.J. (1998) Numerical Ecology, Volume 24, 2nd
#' Edition Elsevier, Amsterdam; Boston.
#' 
#' McGill, B.J. (2011) Species abundance distributions. 105-122 in Biological 
#' Diversity: Frontiers in Measurement and Assessment. eds. A.E. Magurran
#' B.J. McGill.
#' 
#' Whittaker, R.H. (1972) Evolution and Measurement of Species Diversity.
#' Taxon, 21, 213-251.
#' 
#' @export
#' 
#' @examples 
#' # a binary grouping variable (uninvaded or invaded)
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, group_var = "group",
#'                           n_perm = 19, effort_samples = c(5,10))
#' plot(inv_stats)
#' 
#' \donttest{
#' # Not run: 
#' # parallel evaluation using the parallel package 
#' # run in parallel
#' library(parallel)
#' cl = makeCluster(2L)
#' clusterEvalQ(cl, library(mobr))
#' clusterExport(cl, 'inv_mob_in')
#' inv_mob_stats = get_mob_stats(inv_mob_in, 'group', n_perm=999, cl=cl)
#'
#' stopCluster(cl)
#' }
get_mob_stats = function(mob_in, group_var, 
                         index = c("N", "S", "S_n", "S_PIE"),
                         effort_samples = NULL, effort_min = 5,
                         extrapolate = TRUE, return_NA = FALSE, 
                         rare_thres = 0.05, n_perm = 199, 
                         boot_groups = F, conf_level = 0.95, cl=NULL, 
                         ...) {
    EPS <- sqrt(.Machine$double.eps)
    if (n_perm < 1) 
        stop('Set nperm to a value greater than 1') 
     
    INDICES = c("N", "S", "S_n", "S_asymp", "f_0",
                "pct_rare", "PIE", "S_PIE")
    index = match.arg(index, INDICES, several.ok = TRUE)
    
    group_id  = factor(mob_in$env[, group_var])
    
    # Get rarefaction level
    samples_N = rowSums(mob_in$comm) 
    samples_per_group = table(group_id)
   
    if (any(is.null(effort_samples)) | !is.numeric(effort_samples)) {
        N_min_sample = min(samples_N)
        effort_samples = N_min_sample
    } else {
        effort_samples = floor(effort_samples)
    }
   
    if (any(effort_samples < effort_min) & 'S_n' %in% index) {
        warning(paste("The number of individuals for rarefaction analysis is too low and is set to the minimum of", effort_min,"individuals."))
      
        effort_samples[effort_samples < effort_min] = effort_min
        effort_samples = unique(effort_samples)
      
        print(paste("Number of individuals for rarefaction:",
                    paste(effort_samples, collapse = ", ")))
    }
  
    # Group-level indices
    effort_groups = c(effort_samples, effort_samples * min(samples_per_group))
   
    # Abundance distribution pooled in groups
    abund_group = aggregate(mob_in$comm, by = list(group_id), FUN = "sum")
   
    dat_groups = calc_biodiv(abund_mat = abund_group[ , -1],
                             groups = abund_group[ , 1],
                             index = index, 
                             effort = effort_groups, 
                             extrapolate = extrapolate,
                             return_NA = return_NA,
                             rare_thres = rare_thres)
   
    dat_samples = calc_biodiv(abund_mat = mob_in$comm,
                              groups = group_id,
                              index = index,
                              effort = effort_samples, 
                              extrapolate = extrapolate,
                              return_NA = return_NA,
                              rare_thres = rare_thres)
   
    # beta-diversity
   
    # Number of species ---------------------------------------------------------
    if (any(index == "S")) {
        gamma = with(dat_groups, value[index == "S"])
        alpha = with(dat_samples,  value[index == "S"])
      
        beta_S = gamma[group_id] / alpha
        beta_S[!is.finite(beta_S)] = NA
      
        dat_beta_S = data.frame(group = group_id,
                               index = "beta_S",
                               effort = NA,
                               value = beta_S)
        dat_samples = rbind(dat_samples, dat_beta_S)
    }  
   
    # Rarefied richness ---------------------------------------------------------
    if ("S_n" %in% index) {  
        for (i in seq_along(effort_samples)) {
             gamma = with(dat_groups, 
                          value[index == "S_n" & effort == effort_groups[i]])
             alpha = with(dat_samples,
                          value[index == "S_n" & effort == effort_samples[i]])
             beta_S_n = gamma[group_id] / alpha
             beta_S_n[!is.finite(beta_S_n)] = NA
         
             dat_beta_S_n = data.frame(group = group_id,
                                          index = "beta_S_n",
                                          effort = effort_samples[i],
                                          value = beta_S_n)
             dat_samples = rbind(dat_samples, dat_beta_S_n)
        }
        # clean up dat_groups by removing S_n computed using alpha n-value
        dat_groups = subset(dat_groups, !(effort %in% effort_samples))
    } # end rarefied richness

    # Number of rare species ---------------------------------------------------------
    if (any(index == "pct_rare")) {
      gamma = with(dat_groups, value[index == "pct_rare"])
      alpha = with(dat_samples,  value[index == "pct_rare"])
      
      beta_pct_rare = gamma[group_id] / alpha
      beta_pct_rare[!is.finite(beta_pct_rare)] = NA
      
      dat_beta_pct_rare = data.frame(group = group_id,
                                     index = "beta_pct_rare",
                                     effort = NA,
                                     value = beta_pct_rare)
      dat_samples = rbind(dat_samples, dat_beta_pct_rare)
    }  
    
    
    # Asymptotic estimates species richness -------------------------------------
    if ("S_asymp" %in% index) {
        gamma = with(dat_groups, value[index == "S_asymp"])
        alpha = with(dat_samples,  value[index == "S_asymp"])
        beta_S_asymp = gamma[group_id] / alpha
        beta_S_asymp[!is.finite(beta_S_asymp)] = NA
      
        dat_beta_S_asymp = data.frame(group = group_id,
                                      index = "beta_S_asymp",
                                      effort = NA,
                                      value = beta_S_asymp)
        dat_samples = rbind(dat_samples, dat_beta_S_asymp)
    }
    
    
    # Estimated # of missing species richness -------------------------------------
    if ("f_0" %in% index) {
        gamma = with(dat_groups, value[index == "f_0"])
        alpha = with(dat_samples,  value[index == "f_0"])
        beta_f_0 = gamma[group_id] / alpha
        beta_f_0[!is.finite(beta_f_0)] = NA
      
        dat_beta_f_0 = data.frame(group = group_id,
                                      index = "beta_f_0",
                                      effort = NA,
                                      value = beta_f_0)
        dat_samples = rbind(dat_samples, dat_beta_f_0)
    }
   
    # Effective number of species based on PIE ----------------------------------
    if ("S_PIE" %in% index) {
        gamma = with(dat_groups, value[index == "S_PIE"])
        alpha = with(dat_samples,  value[index == "S_PIE"])
      
        beta_S_PIE = gamma[group_id] / alpha
        beta_S_PIE[!is.finite(beta_S_PIE)] = NA
      
        dat_beta_S_PIE = data.frame(group = group_id,
                                    index = "beta_S_PIE",
                                    effort = NA,
                                    value = beta_S_PIE)
        dat_samples = rbind(dat_samples, dat_beta_S_PIE)
    }
   
    # Significance tests -------------------------------------------------------
   
    # alpha-scale
    alpha_avg = dat_samples %>% 
                group_by(group, index, effort) %>%
                summarise(alpha_avg = mean(value, na.rm=T)) %>%
                ungroup()
    D_bar = alpha_avg %>%
                group_by(index, effort) %>%
                summarise(D_bar = mean(dist(alpha_avg), na.rm=T)) %>%
                ungroup()
    F_obs = get_F_values(dat_samples, permute = F)
    cat('\nComputing null model at alpha-scale\n')
    F_rand = bind_rows(pbapply::pbreplicate(n_perm, 
                                  get_F_values(dat_samples, permute = T),
                              simplify = F, cl = cl)) %>%
             ungroup()
    F_obs = F_obs %>% mutate(F_val_obs = F_val, F_val = NULL)
    F_rand = left_join(F_rand, F_obs, by = c("index", "effort"))

    samples_tests = F_rand %>% 
                   group_by(index, effort) %>%
                   summarise(F_stat = first(F_val_obs),
                             p_val = (sum(F_val >= F_val_obs - EPS) + 1) /
                                     (n_perm + 1)) %>%
                              
                   ungroup()
    samples_tests = left_join(D_bar, samples_tests, by = c("index", "effort"))
   
    # gamma-scale
    if (boot_groups) {
        # bootstrap sampling within groups
      
        abund_dat = cbind(group_id, mob_in$comm)
        cat('\nComputing bootstrap confidence intervals at the gamma-scale\n')
        boot_repl_groups = pbapply::pbreplicate(n_perm,
                               boot_sample_groups(abund_dat,
                                                  index = index,
                                                   effort = effort_groups,
                                                   extrapolate = extrapolate,
                                                   return_NA = return_NA, 
                                                  rare_thres = rare_thres),
                            simplify = F, cl = cl)
        boot_repl_groups = bind_rows(boot_repl_groups)
      
        alpha = 1 - conf_level
        probs = c(alpha / 2, 0.5, 1 - alpha / 2)
     
        dat_groups = boot_repl_groups %>% 
                     group_by(group, index, effort) %>%
                     do(setNames(data.frame(t(quantile(.$value, probs, na.rm = T))),
                                 c("lower","median","upper")))      
    } else {
        delta_obs = get_group_delta(mob_in$comm, group_id, index,
                                    effort_groups, extrapolate,
                                    return_NA, rare_thres, permute=F)
        cat('\nComputing null model at gamma-scale\n')
        delta_rand = bind_rows(pbapply::pbreplicate(n_perm, 
                                   get_group_delta(mob_in$comm, group_id,
                                       index, effort_groups, extrapolate,
                                       return_NA, rare_thres, permute = T),
                               simplify = F, cl = cl))
        delta_obs = delta_obs %>% mutate(d_obs = delta, delta = NULL)
        delta_rand = suppressMessages(left_join(delta_rand, delta_obs))
      
        groups_tests = delta_rand %>% 
                       group_by(index, effort) %>%
                       summarise(D_bar = first(d_obs),
                                 p_val = (sum(delta >= d_obs - EPS) + 1) /
                                         (n_perm + 1)) %>%
                       ungroup()
    }
    # order output data frames by indices
    dat_samples$index = factor(dat_samples$index,
                               levels = c("N",
                                          "S", "beta_S",
                                          "S_n", "beta_S_n",
                                          "S_asymp", "beta_S_asympS",
                                          "f_0", "beta_f_0",
                                          "pct_rare", "beta_pct_rare",
                                          "PIE",
                                          "S_PIE", "beta_S_PIE"))
    dat_samples = dat_samples[order(dat_samples$index, dat_samples$effort,
                                    dat_samples$group), ]
   
    dat_groups$index = factor(dat_groups$index, levels = index)
    dat_groups = dat_groups[order(dat_groups$index, dat_groups$effort,
                                  dat_groups$group), ]
   
    #remove unused factor levels
    dat_samples$index = factor(dat_samples$index)
    if (boot_groups) {
        out = list(samples_stats = dat_samples,
                   groups_stats  = dat_groups,
                   samples_tests  = samples_tests)
    } else {      
        groups_tests$index = factor(groups_tests$index, levels = index)
        groups_tests = groups_tests[order(groups_tests$index), ]
       
        out = list(samples_stats = dat_samples,
                   groups_stats  = dat_groups,
                   samples_tests  = samples_tests,
                   groups_tests   = groups_tests)
    }
    if ("pct_rare" %in% index)
        out$rare_thres = rare_thres
    class(out) = 'mob_stats'
    return(out)
}

# Panel function for alpha-scale results
samples_panel1 = function(sample_dat, samples_tests, col, ylab = "",
                           main = expression(alpha * "-scale"), 
                          cex.axis=1.2, ...) {
   label = substitute(paste(italic(bar(D)), ' = ', D_bar, ', ', italic(p), 
                            ' = ', p_val),
                      list(D_bar = round(samples_tests$D_bar, 2),
                           p_val = round(samples_tests$p_val, 3)))
   boxplot(value ~ group, data = sample_dat, main = main,
           ylab =  ylab, col = col, cex.axis=cex.axis, cex.main = 1.5,
           frame.plot=T, ...)
   groups = levels(sample_dat$group)
   #axis(side=1, at=1:length(groups), labels=groups, tick=FALSE,
   #      cex.axis=cex.axis)
   mtext(label, side = 3, line = 0)  
}

# Panel function for gamma-scale results
groups_panel1 = function(group_dat, tests, col, ylab = "",
                         main = expression(gamma * "-scale"),
                         cex.axis=1.2, ...) {
    label = substitute(paste(italic(bar(D)), ' = ', D_bar, ', ', italic(p), 
                             ' = ', p_val), 
                       list(D_bar = round(tests$D_bar, 2),
                            p_val = round(tests$p_val, 3)))
    boxplot(value ~ group, data = group_dat, main = main,
            ylab = ylab, boxwex = 0, 
            ylim = c(0, 1.1 * max(group_dat$value, na.rm = T)),
            col = col, cex.axis=cex.axis, cex.main = 1.5, frame.plot=T,
            ...)
    groups = levels(group_dat$group)
    points(value ~ group, data = group_dat, pch = 8, cex = 1.5, lwd = 2,
           col = col, ...)
    mtext(label, side = 3, line = 0)
}

# Panel function for gamma-scale results with confidence intervals
groups_panel2 = function(group_dat, col, ylab = "", 
                         main = expression(gamma * "-scale"),
                         cex.axis=1.2, ...) {
    boxplot(median ~ group, data = group_dat, main = main,
            ylab = ylab, boxwex = 0, ylim = c(0, 1.1*max(group_dat$upper)),
            col = col, cex.axis=cex.axis, cex.main=1.5, frame.plot=T, 
            ...)
    groups = levels(group_dat$group)
    plotrix::plotCI(1:nrow(group_dat), group_dat$median, li = group_dat$lower,
                    ui = group_dat$upper, add = T, pch = 19, cex = 1.5,
                    sfrac = 0.02, col = col, ...)
}

#' Plot alpha- and gamma-scale biodiversity statistics for a MoB analysis
#' 
#' Plots a \code{mob_stats} object which is produced by the 
#' function \code{get_mob_stats}. The p-value for each statistic
#' is displayed in the plot title if applicable.
#' 
#' The user may specify which results to plot or simply to plot 
#' all the results. 
#' 
#' @param mob_stats a \code{mob_stats} object that has the samples and 
#' treatment level statistics
#' 
#' @param index The biodiversity statistics that should be plotted.
#' See \code{\link{get_mob_stats}} for information on the indices. By default there
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
#' @export
#' 
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' # without bootstrap CI for gamma-scale
#' inv_stats = get_mob_stats(inv_mob_in, group_var = "group", n_perm = 20)
#' plot(inv_stats) 
#' 
#' plot(inv_stats, multi_panel = TRUE)
#' # with bootstrap CI for gamma-scale
#' inv_stats_boot = get_mob_stats(inv_mob_in, group_var = "group", n_perm = 20,
#'                                boot_groups=TRUE)
#' plot(inv_stats_boot)
plot.mob_stats = function(mob_stats, index = NULL, multi_panel = FALSE, 
                          col = c("#FFB3B5", "#78D3EC", "#6BDABD", "#C5C0FE",
                                  "#E2C288", "#F7B0E6", "#AAD28C"), 
                          cex.axis=1.2, ...) {
    # default colors derived with colorspace::rainbow_hcl(5, c=60, l=80)
    if (any(is.na(col)) & length(col) == 1) 
        col_groups = 1
    else
        col_groups = col
    if (is.null(index))
        index = as.character(unique(mob_stats$samples_stats$index))
    INDICES = c("N", "S", "S_n", "S_asymp", "f_0", 
                "pct_rare", "PIE", "S_PIE")
    if (multi_panel) 
        index = c("S","S_n","pct_rare","S_PIE")
    index = match.arg(index, INDICES, several.ok = TRUE)
   
    var_names = levels(mob_stats$samples_stats$index)
    var_names2 = var_names[var_names != "beta_S" & var_names != "beta_S_PIE"]
   
    index_match = intersect(index, var_names)
    if (length(index_match) == 0)
        stop(paste("The indices", paste(index, collapse = ", "), 
                   "are missing in the input. Please choose other indices or re-run get_mob_stats with the indices of interest."))
   
    index_missing = setdiff(index, var_names2)
    if (length(index_missing) > 0)
        warning(paste("The indices", paste(index, collapse = ", "), 
                      "are missing in the input and cannot be plotted."))
   
    if ("S_n" %in% index_match) {
        S_n_samples = filter(mob_stats$samples_stats, index == "S_n")
        S_n_groups = filter(mob_stats$groups_stats, index == "S_n")
        S_n_len = max(length(unique(S_n_samples$effort)),
                        length(unique(S_n_groups$effort)))
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
                             "PIE" = "PIE")
         
            dat_samples = filter(mob_stats$samples_stats, index == var)
            dat_tests = filter(mob_stats$samples_tests, index == var)
            samples_panel1(dat_samples, dat_tests, ylab = y_label,
                           main = expression(alpha * "-scale"), col = col,
                           cex.axis=cex.axis, ...)
            
            # insert blank space b/c non beta plot for these stats
            plot(1:10, 1:10, type='n', axes=F, xlab='', ylab='')
         
            dat_groups = filter(mob_stats$groups_stats, index == var)
         
            if (is.null(mob_stats$groups_tests)) {
                groups_panel2(dat_groups, col = col_groups, cex.axis=cex.axis, ...) 
            } else {  
                tests = filter(mob_stats$groups_tests, index == var)
                groups_panel1(dat_groups, tests, col = col_groups, cex.axis=cex.axis, ...) 
            }
        }
      
        if (var %in% c("S", "S_asymp", "pct_rare", "f0", "PIE", "S_PIE")) {
         
            if (!multi_panel)
                op = par(mfrow = c(1, 3), cex.lab = 1.6,
                         oma = c(0, 2, 0, 0), mar = c(4, 3, 5, 1), xpd = NA)
         
            if (multi_panel) {
                if (var == "f_0") 
                    par(fig = c(0, 0.33, 1 / n_rows, 2 / n_rows), new = T)
                if (var == "S_PIE") 
                    par(fig = c(0, 0.33, 0         , 1 / n_rows), new = T)
            }
            
            y_label = switch(var,
                             "S" = expression('Richness (' * italic(S) * ')'),
                             "S_asymp" = expression('Asympotic richness (' *
                                                      italic(S[asymp]) * ')'),
                             "pct_rare" = paste("% of species in lower",
                                               ifelse(is.character(mob_stats$rare_thres),
                                                      mob_stats$rare_thres,
                                                      mob_stats$rare_thres * 100),
                                               "% of abundance"),
                             "f_0" = expression('Undetected richness (' * italic(f)[0] * ')'),
                             "S_PIE" = expression('ENS of PIE (' * italic(S)[PIE] * ')'))
         
            dat_samples = filter(mob_stats$samples_stats, index == var)
            dat_tests = filter(mob_stats$samples_tests, index == var)
            samples_panel1(dat_samples, dat_tests, p_val = dat_tests$p_val, ylab =  y_label,
                           main = expression(alpha * "-scale"), col = col, 
                           cex.axis=cex.axis, ...)
         
           if (multi_panel) {
               if (var == "pct_rare") 
                   par(fig = c(0.33, 0.67, 1/n_rows, 2/n_rows), new = T)
               if (var == "S_PIE") 
                   par(fig = c(0.33, 0.67, 0       , 1/n_rows), new = T)
           }
         
           beta_var = paste("beta", var, sep = "_")
           dat_samples = filter(mob_stats$samples_stats, index == beta_var)
           dat_tests = filter(mob_stats$samples_tests, index == beta_var)
           samples_panel1(dat_samples, dat_tests, ylab =  "",
                          main = expression(beta * "-diversity (=" *
                                            gamma / alpha * ")"),
                          col = col, cex.axis=cex.axis, ...)
         
           if (multi_panel) {
               if (var == "pct_rare")
                   par(fig = c(0.67, 1.0, 1 / n_rows, 2 / n_rows), new = T)
               if (var == "S_PIE")
                   par(fig = c(0.67, 1.0, 0         , 1 / n_rows), new = T)
           }
         
           dat_groups = filter(mob_stats$groups_stats, index == var)
         
           if (is.null(mob_stats$groups_tests)) {
               groups_panel2(dat_groups, col = col_groups, ...) 
           } else {
               tests = filter(mob_stats$groups_tests, index == var)
               groups_panel1(dat_groups, tests, col = col_groups, cex.axis=cex.axis,
                             ...) 
           }
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
                        new = T)
            
                dat_samples = filter(S_n_samples, effort == effort_samples[i])
                dat_tests = filter(mob_stats$samples_tests, 
                                   index == var & effort == effort_samples[i])
            
                fig_title = substitute(paste(alpha, "-scale, n = ", n),
                                       list(n = effort_samples[i]))
            
                samples_panel1(dat_samples, dat_tests, 
                               ylab = y_label,
                               main = '', col = col, cex.axis=cex.axis,
                               ...)
                par(new=TRUE)
                plot(1:10, 1:10, type='n', axes=F, xlab='', ylab='',
                     main=fig_title, cex.main=1.5)
            
                if (multi_panel)
                    par(fig = c(0.33, 0.67, (1 + i) / n_rows, (2 + i) / n_rows),
                        new = T)
            
                beta_var = paste("beta", var, sep = "_")
                dat_samples = filter(mob_stats$samples_stats, index == beta_var)
                dat_tests = filter(mob_stats$samples_tests,
                                   index == beta_var & effort == effort_samples[i])
                n = effort_samples[i]
                fig_title = substitute(paste(beta, "-diversity (=", gamma / alpha,
                                             "), n = ", n), 
                                       list(n = effort_samples[i]))
                samples_panel1(dat_samples, dat_tests, main = '',
                               ylab = "", col = col, cex.axis=cex.axis, ...)

                par(new=TRUE)
                plot(1:10, 1:10, type='n', axes=F, xlab='', ylab='',
                     main=fig_title, cex.main=1.5)
                                
                if (multi_panel)
                    par(fig = c(0.67, 1.0, (1 + i) / n_rows, (2 + i) / n_rows),
                        new = T)
            
                dat_groups = filter(S_n_groups, effort == effort_groups[i])
                fig_title = substitute(paste(gamma, "-scale, n = ", n),
                                       list(n = effort_groups[i]))            
                if (is.null(mob_stats$groups_test)) {
                    groups_panel2(dat_groups, main = '', col = col_groups,
                                  ...) 
                    par(new=TRUE)
                    plot(1:10, 1:10, type='n', axes=F, xlab='', ylab='',
                         main=fig_title, cex.main=1.5)
                } else {
                    tests = filter(mob_stats$groups_tests, 
                                   index == var & effort == effort_groups[i])
                    groups_panel1(dat_groups, tests, ylab = "", main = '',
                                  col = col_groups, cex.axis=cex.axis, ...)
                    par(new=TRUE)
                    plot(1:10, 1:10, type='n', axes=F, xlab='', ylab='',
                         main=fig_title, cex.main=1.5)
                }
            }
            y_coords = (S_n_len:0) / S_n_len
        }
    }
    par(op)
}
