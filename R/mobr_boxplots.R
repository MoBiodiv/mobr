#' Estimation of species richness
#' 
#' \code{calc_chao1}: estimation of species richness based on the methods
#' proposed in Chao (1984, 1987). 
#' 
#' This function is a trimmed version of \code{iNext::ChaoRichess} found at 
#' \url{https://github.com/JohnsonHsieh/iNEXT}. T. C. Hsieh, K. H. Ma and Anne
#' Chao are the orginal authors of the \code{iNEXT} package. 
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
#' Chao, A. (1987) Estimating the population size for capture-recapture data with
#' unequal catchability. Biometrics, 43, 783-791.
#' 
#' @export
calc_chao1 = function(x){
    if(!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
        stop("invalid data structure")
    if(is.matrix(x) | is.data.frame(x)){
        S_Chao1= apply(x, 1, calc_chao1)
    } else {
        n <- sum(x)
        D <- sum(x > 0)
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        if (f1 > 0 & f2 > 0)
            S_Chao1 <- D + (n - 1) / n * f1^2 / (2 * f2)
        else if (f1 > 1 & f2 == 0) # bias corrected form
            S_Chao1 <- D + (n - 1) / n * f1 * (f1 - 1) / (2 * (f2 + 1))
        else
            S_Chao1 <- D
    }
    return(S_Chao1)
}

#' Calculate probability of interspecific encounter (PIE)
#' 
#' PIE is also known as Simpson's evenness index and the Gini-Simpson index.
#' The formula of Hurlbert (1971) is used to calculate PIE:
#' 
#' \eqn{PIE = N /(N - 1) * (1 - p_i^2)}
#' 
#' where N is the total number of individuals and \eqn{p_i} is the relative abundance
#' of species i. This formulation uses sampling without replacement and it is
#' sometimes referred to as the bias corrected formulation of PIE. 
#' 
#' The code in this function borrows heavily from the function vegan::diversity()
#' but computes a different quantitiy. The function vegan::diversity() computes
#' PIE when sampling with replacement is assumed. The difference between the two 
#' formulations will decrease as N becomes large. Jari Oksanen and Bob O'Hara are
#' the original authors of the function vegan::diversity().
#' 
#' @inheritParams rarefaction
#' @param ENS boolean that determines if the effective number of species should
#' be returned or the raw PIE value. Defaults to FALSE
#'
#' @author Dan McGlinn
#' 
#' @references 
#' Hurlbert, S. H. 1971. The noncept of species diversity: a critique and
#'  alternative parameters. Ecology 52: 577–586.
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
        x = sweep(x, 1, total, "/")
    }
    else {
        total = sum(x)
        x = x / total
    }
    x = x * x
    if (length(dim(x)) > 1) 
        H = rowSums(x, na.rm = TRUE)
    else 
        H = sum(x, na.rm = TRUE)
    if (ENS) {
        # convert to effective number of species
        H = 1 / H
    } else {
        # correct for sampling with replacement
        H = total / (total - 1) * (1 - H)
        # check NA in data
        if (any(NAS = is.na(total)))
            H[NAS] = NA
    }        
    H[!is.finite(H) | total == 0] <- NA # set NA, when total == 0 or total == 1
    return(H)
}

# generate a single bootstrap sample of group-level biodiversity indices
boot_sample_groups <- function(abund_mat, index, n_rare)
{
   # sample rows and calculate abundance vector
   sample_dat <- by(abund_mat, INDICES = abund_mat$group_id, FUN = sample_frac, replace = T)
   class(sample_dat) <- "list"
   sample_dat <- bind_rows(sample_dat)
   
   # abundance distribution pooled in groups
   abund_group = aggregate(sample_dat[,-1], by = list(sample_dat[,1]), FUN = "sum")
   
   dat_groups <- calc_biodiv(abund_mat = abund_group[,-1],
                             groups = abund_group[,1],
                             index = index,
                             n_rare = n_rare)
   return(dat_groups)
}

# Get group-level p-values from permutation test
get_pval <- function(rand, obs, n_samples)
{
   n_extremes = sum(rand < -abs(obs) | rand > abs(obs))
   p_val = n_extremes/n_samples
   return(p_val)
}

# Get biodiversity indices
calc_biodiv <- function(abund_mat, groups, index, n_rare)
{
   out <- expand.grid(group = groups,
                      index = index[index != "S_rare"],
                      n_rare = NA, value = NA)
   N <- rowSums(abund_mat) 
   
   # Number of individuals -----------------------------------------------------
   if (any(index == "N")){
      out$value[out$index == "N"] <- N
   } 
   
   # Number of species ---------------------------------------------------------
   if (any(index == "S")){
      out$value[out$index == "S"] <-  rowSums(abund_mat > 0)
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if (any(index == "S_rare")){
      
         dat_S_rare <- expand.grid(group = groups,
                                   index = 'S_rare',
                                   n_rare = n_rare, value = NA)
         out <- rbind(out, dat_S_rare)
      
         out$value[out$index == "S_rare"] <- apply(abund_mat, 1, rarefaction,
                                                   method = 'indiv', effort= n_rare,
                                                   force = T)
   } # end rarefied richness

   # Asymptotic estimates species richness -------------------------------------
   if (any(index == "S_asymp")){
      
      S_asymp <- try(calc_chao1(abund_mat))
      if (class(S_asymp) == "try_error"){
         warning("The Chao richness estimator cannot be calculated for all groups.")
      } else {
         S_obs = rowSums(abund_mat > 0)
         S_asymp[!is.finite(S_asymp)] <- NA
      }
    
      out$value[out$index == "S_asymp"] <- S_asymp - S_obs
   }
   
   # Probability of Interspecific Encounter (PIE)-------------------------------
   if (any(index == "PIE")){
      
      plots_n01 = N <= 1 # Hurlbert's PIE can only be calculated for two or more individuals
      
      if (sum(plots_n01) == 1){
         warning(paste("There is", sum(plots_n01), "plot with less than two individuals.
                       This is removed for the calculation of PIE."))
         
      }
      
      if (sum(plots_n01) > 1){
         warning(paste("There are", sum(plots_n01), "plots with less than two individuals.
                       These are removed for the calculation of PIE."))
         
      }
      
      out$value[out$index == "PIE"] <- calc_PIE(abund_mat)
   }
   
   # Effective number of species based on PIE ----------------------------------
   if (any(index == "S_PIE")){
      plots_n01 = N <= 1
      S_PIE = calc_PIE(abund_mat, ENS=TRUE)
      S_PIE[plots_n01] <- NA
      out$value[out$index == "S_PIE"] <- S_PIE
   }
   
   return(out)
   
}

#Get F statistics from diversity indices and grouping vector
get_F_values <- function(div_dat, permute = F)
{
   if (permute)
      div_dat <- div_dat %>%
         group_by(index, n_rare) %>%
         mutate(group = sample(group))
   
   models <- div_dat %>%
      group_by(index, n_rare) %>%
      do(mod = lm(value ~ group, data = .))
   
   models <- models %>% mutate(F_val = anova(mod)$F[1],
                               mod = NULL) %>% ungroup()
   
   models$F_val[!is.finite(models$F_val)] <- NA
   
   return(models)
}

# Get group-level differences 
get_group_diff <- function(abund_mat, group_bin, index, n_rare,
                           permute = F)
{
   if (permute)
      group_bin <- sample(group_bin)
   
   abund_group <- aggregate(abund_mat, by = list(group_bin), FUN = "sum")
   
   dat_groups <- calc_biodiv(abund_mat = abund_group[,-1],
                             groups = abund_group[,1],
                             index = index,
                             n_rare = n_rare)
   delta_groups <- dat_groups %>%
      group_by(index, n_rare) %>%
      summarise(delta = diff(value)) %>%
      ungroup()
   
   return(delta_groups)
}


#' Calculate sample based and group based biodiversity statistics.
#' @inheritParams get_delta_stats
#' 
#' @param group_var A string that specifies which field in \code{mob_in$env} the 
#' data should be grouped by
#' 
#' @param ref_group The level in \code{group_var} that is used as control group for the
#' control vs. treatment test at the group level.
#' 
#' @param index The calculated biodiversity statistics. The options are
#' \itemize{
#'    \item \code{N} ... Number of individuals (total abundance)
#'    \item \code{S} ... Number of species
#'    \item \code{S_rare} ... Rarefied number of species
#'    \item \code{S_asymp} ... Estimated asymptotic species richness
#'    \item \code{PIE} ... Hurlbert's PIE (Probability of Interspecific Encounter)
#'    \item \code{S_PIE} ... Effective number of species based on PIE
#' }
#' See \emph{Details} for additional information on the biodiverstiy statistics.
#' 
#' @param n_rare_samples The standardized number of individuals used for the
#' calculation of rarefied species richness at the sample level. This can a be
#' single value or an integer vector. As default the minimum number of individuals
#' found across the samples is used, when this is not smaller than \code{n_rare_min}.
#' 
#' @param n_rare_min The minimum number of individuals considered for the calculation
#' of rarefied richness. Samples with less individuals then \code{n_rare_min} are
#' excluded from the analysis with a warning. Accordingly, when \code{n_rare_samples}
#' is set by the user it has to be higher than \code{n_rare_min}.
#' 
#' @param n_perm The number of permutations to use for testing for treatment effects.
#' 
#' @param boot_groups Use bootstrap resampling within groups to derive group-level
#' confidence intervals for all biodiversity indices. Default is \code{FALSE}. 
#' See \emph{Details} for information on the bootstrap approach.
#' 
#' @param conf_level Confidence level used for the calculation of group-level 
#' bootstrapped confidence intervals. Only used when \code{boot_groups = TRUE}.
#' 
#' @details 
#' 
#' \strong{BIODIVERSITY INDICES}
#' 
#' \strong{Rarefied species richness} is the expected number of species, given
#' a defined number of sampled individuals (Gotelli & Colwell 2001). Rarefied richness
#' at the sample level is calculated for the values provided in
#' \code{n_rare_samples} as long as these values are not smaller than the 
#' user-defined minimum value \code{n_rare_min}. In this case the minimum value
#' is used and samples with less individuals are discarded.
#' When no values for \code{n_rare_samples} are provided the observed 
#' minimum number of individuals of the samples is used, which is the standard
#' in rarefaction analysis (Gotelli & Colwell 2001). Because the number
#' of individuals is expected to scale linearly with sample area or effort,
#' at the group level the number of individuals for rarefaction is calculated as the minimum
#' number of samples within groups times \code{n_rare_samples}. For example,
#' when there are 10 samples within each group, \code{n_rare_groups} equals
#' \code{10 * n_rare_samples}. 
#' 
#' \strong{Asymptotic species richness} is calculated using the bias-corrected Chao1
#'  estimator (Chiu et al. 2014) provided in the function \code{\link[vegan]{estimateR}}.
#' 
#' \strong{PIE} represents the probability that two randomly drawn individuals
#' belong to the same species. Here we use the definition of Hurlbert (1971), which considers
#' sampling without replacement. PIE is closely related to the well-known Simpson 
#' diversity index, but the latter assumes sampling with replacement.
#' 
#' \strong{S_PIE} represents the Effective Number of Species derived from the PIE.
#' This corresponds to the the number of equally abundant species (i.e. a perfectly
#' even community), there would need to be to achieve the observed PIE (Jost 2006).
#' This means the higher the difference between S and S_PIE the more
#' uneven is the observed community. An intuitive interpretation of S_PIE is that
#' it corresponds to the number of dominant (highly abundant) species in the community.
#' 
#' For species richness \code{S}, rarefied richness \code{S_rare}, asymptotic
#' richness \code{S_asymp}, and the Effective Number of Species \code{S_PIE}
#' we also calculate beta-diversity using multiplicative partitioning
#' (Whittaker 1972, Jost 2007). That means for these indices we estimate beta-diversity 
#' as the ratio of group-level diversity (gamma) divided by sample-level diversity (alpha). 
#' 
#' \strong{PERMUTATION TESTS AND BOOTSTRAP}
#' 
#' We used permutation tests for assessing differences of the biodiversity statistics
#' among the groups (Legendre & Legendre 1998). At the sample level, one-way ANOVA 
#' (i.e. F-test) is implemented by shuffling treatment group labels across samples.
#' 
#'  At the group-level there is test whether the difference between treatment and control
#'  is significantly different from zero. For this purpose we permuted treatment group labels
#'  across samples, pooled the groups, and calculated the treatment-control difference for each
#'  permutation. This test only works for the comparison of two groups, i.e. treatment vs. control.
#'  The means the reference level specified by argument \code{ref_group} is compared
#'  with all other groups pooled.
#'  
#'  Especially, when there are more then two groups a bootstrap approachj can be used
#'  to test differences at the group level. When \code{boot_groups = T} instead
#'  of the group-level permutation test, there will be resampling of samples within
#'  groups to derive group-level confidence intervals for all biodiversity indices.
#'  The function output includes lower and upper confidence bounds and the median
#'  of the bootstrap samples. Please note that for the richness indices sampling
#'  with replacement corresponds to rarefaction to ca. 2/3 of the individuals, because
#'  the same samples occur several times in the resampled data sets.
#'  
#' 
#' @return A list of class \code{mob_stats} that contains sample-scale and 
#' group-scale biodiversity statistics, as well as the p-values for permutation tests
#' at both scales.
#' 
#' When \code{boot_groups = TRUE} there are no p-values at the group level. Instead
#' there is lower bound, median, and upper bound for each biodiversity index derived
#' from the bootstrap within groups.
#
#' @author Felix May and Dan McGlinn
#' 
#' @references 
#' 
#' Chiu, C.-H., Wang, Y.-T., Walther, B.A. & Chao, A. (2014). An improved nonparametric lower bound of species richness via a modified good-turing frequency formula. Biometrics, 70, 671-682.
#' 
#' Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology letters, 4, 379-391.
#' 
#' Hurlbert, S.H. (1971). The Nonconcept of Species Diversity: A Critique and Alternative Parameters. Ecology, 52, 577-586.
#' 
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.
#' 
#' Jost, L. (2007). Partitioning Diversity into Independent Alpha and Beta Components. Ecology, 88, 2427–2439.
#' 
#' Legendre, P. & Legendre, L.F.J. (1998). Numerical Ecology, Volume 24, 2nd Edition Elsevier, Amsterdam; Boston.
#' 
#' Whittaker, R.H. (1972). Evolution and Measurement of Species Diversity. Taxon, 21, 213-251.
#' 
#' @export
#' 
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, group_var = "group", ref_group = "uninvaded", n_perm = 20, n_rare_samples = c(5,10))
#' plot(inv_stats)
get_mob_stats = function(mob_in,
                         group_var,
                         ref_group = NULL,
                         index = c("N","S","S_rare","S_asymp","S_PIE"),
                         n_rare_samples = NULL,
                         n_rare_min = 5,
                         n_perm = 200,
                         boot_groups = F,
                         conf_level = 0.95
                         )
{
   if (n_perm < 1) 
       stop('Set nperm to a value greater than 1') 
   
   INDICES <- c("N", "S", "S_rare","S_asymp","PIE","S_PIE")
   index <- match.arg(index, INDICES, several.ok = TRUE)
   
   group_id  = factor(mob_in$env[, group_var])
   
   if (is.null(ref_group))
      ref_group = levels(group_id)[1]
   
   if (!ref_group %in% levels(group_id))
      stop("ref_group has to be one level in group_var!")
   
   group_id = relevel(group_id, ref_group)
   
   # Create factor with just two levels: treatment / control
   # for calculation of group-level differences
   group_bin = factor(rep("control", times = length(group_id)),
                      levels = c("control","treatment"))
   group_bin[group_id != ref_group] <- "treatment"
  
   # # Add labels to groups
   # group_type = c("(ctrl)",rep("(trt)", length(levels(group_id)) - 1))
   # group_labels <- paste(levels(group_id), group_type)
   # group_id <- factor(group_id, labels = group_labels) 
   
   print(index)
 
   # Get rarefaction level
   samples_N = rowSums(mob_in$comm) 
   samples_per_group <- table(group_id)
   
   if (any(is.null(n_rare_samples)) | !is.numeric(n_rare_samples)){ 
      N_min_sample = min(samples_N)
      n_rare_samples = N_min_sample
   } else {
      n_rare_samples <- floor(n_rare_samples)
   }
   
   if (any(n_rare_samples < n_rare_min)){
      warning(paste("The number of individuals for rarefaction analysis is too low and is set to the minimum of", n_rare_min,"individuals."))
      
      n_rare_samples[n_rare_samples < n_rare_min] <- n_rare_min
      n_rare_samples <- unique(n_rare_samples)
      
      print(paste("Number of individuals for rarefaction:",
                  paste(n_rare_samples, collapse = ", ")))
   }
  
   # Group-level indices
   n_rare_groups <- n_rare_samples*min(samples_per_group)
   
   # Abundance distribution pooled in groups
   abund_group = aggregate(mob_in$comm, by = list(group_id), FUN = "sum")
   
   dat_groups <- calc_biodiv(abund_mat = abund_group[,-1],
                             groups = abund_group[,1],
                             index = index,
                             n_rare = n_rare_groups)
   
   dat_samples <- calc_biodiv(abund_mat = mob_in$comm,
                              groups = group_id,
                              index = index,
                              n_rare = n_rare_samples)
   
   # beta-diversity
   
   # Number of species ---------------------------------------------------------
   if (any(index == "S")){
      gamma <- with(dat_groups, value[index == "S"])
      alpha <- with(dat_samples,  value[index == "S"])
      
      beta_S <- gamma[group_id]/alpha
      beta_S[!is.finite(beta_S)] <- NA
      
      dat_betaS <- data.frame(group = group_id,
                              index = "beta_S",
                              n_rare = NA,
                              value = beta_S)
      dat_samples <- rbind(dat_samples, dat_betaS)
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if ("S_rare" %in% index){  
      for (i in 1:length(n_rare_samples)){
         gamma <- with(dat_groups, value[index == "S_rare" & n_rare == n_rare_groups[i]])
         alpha <- with(dat_samples,  value[index == "S_rare" & n_rare == n_rare_samples[i]])
         
         beta_S_rare <- gamma[group_id]/alpha
         beta_S_rare[!is.finite(beta_S_rare)] <- NA
         
         dat_beta_S_rare <- data.frame(group = group_id,
                                       index = "beta_S_rare",
                                       n_rare = n_rare_samples[i],
                                       value = beta_S_rare)
         dat_samples <- rbind(dat_samples, dat_beta_S_rare)
      }
   } # end rarefied richness

   # Asymptotic estimates species richness -------------------------------------
   if ("S_asymp" %in% index){
      gamma <- with(dat_groups, value[index == "S_asymp"])
      alpha <- with(dat_samples,  value[index == "S_asymp"])
      
      beta_S_asymp <- gamma[group_id]/alpha
      beta_S_asymp[!is.finite(beta_S_asymp)] <- NA
      
      dat_beta_S_asymp <- data.frame(group = group_id,
                                     index = "beta_S_asymp",
                                     n_rare = NA,
                                     value = beta_S_asymp)
      dat_samples <- rbind(dat_samples, dat_beta_S_asymp)
   }
   
   # Effective number of species based on PIE ----------------------------------
   if ("S_PIE" %in% index){
      gamma <- with(dat_groups, value[index == "S_PIE"])
      alpha <- with(dat_samples,  value[index == "S_PIE"])
      
      beta_S_PIE <- gamma[group_id]/alpha
      beta_S_PIE[!is.finite(beta_S_PIE)] <- NA
      
      dat_beta_S_PIE <- data.frame(group = group_id,
                                     index = "beta_S_PIE",
                                     n_rare = NA,
                                     value = beta_S_PIE)
      dat_samples <- rbind(dat_samples, dat_beta_S_PIE)
   }
   
   # Significance tests
   
   # sample level
   F_obs <- get_F_values(dat_samples, permute = F)
   F_rand <- dplyr::bind_rows(replicate(n_perm, get_F_values(dat_samples, permute = T), simplify = F)) %>% ungroup()
   F_obs <- F_obs %>% mutate(F_val_obs = F_val,
                             F_val = NULL)
   F_rand <- left_join(F_rand, F_obs)
   
   p_val_samples <- F_rand %>% 
      group_by(index, n_rare) %>%
      summarise(p_val = sum(F_val_obs <= F_val)/n_perm) %>%
      ungroup()
   
   # group level
   if (!boot_groups){
      diff_obs <- get_group_diff(mob_in$comm, group_bin, index, n_rare = n_rare_groups,                                                                 permute = F)
      diff_rand <- bind_rows(replicate(n_perm, get_group_diff(mob_in$comm, group_bin,
                                                              index,
                                                              n_rare = n_rare_groups,
                                                              permute = T),
                                       simplify = F))
      diff_obs <- diff_obs %>% mutate(d_obs = delta,
                                      delta = NULL)
      diff_rand <- left_join(diff_rand, diff_obs)
      
      p_val_groups <- diff_rand %>% 
         group_by(index, n_rare) %>%
         summarise(p_val = get_pval(rand = delta, obs = first(d_obs),
                                    n_samples = n_perm)) %>% ungroup()
   } else {
      # bootstrap sampling within groups
      
      abund_dat <- cbind(group_id, mob_in$comm)
      
      boot_repl_groups <- replicate(n_perm,
                                    boot_sample_groups(abund_dat,
                                                       index = index,
                                                       n_rare = n_rare_groups),
                                    simplify = F)
      boot_repl_groups <- bind_rows(boot_repl_groups)
      
      alpha <- 1 - conf_level
      p <- c(alpha/2, 0.5, 1 - alpha/2)
     
      dat_groups <- boot_repl_groups %>% 
         group_by(group, index, n_rare) %>%
         do(setNames(data.frame(t(quantile(.$value, p, na.rm = T))),
                     c("lower","median","upper")))
   }

   # order output data frames by indices
   dat_samples$index <- factor(dat_samples$index,
                               levels = c("N",
                                          "S","beta_S",
                                          "S_rare","beta_S_rare",
                                          "S_asymp","beta_S_asymp",
                                          "PIE",
                                          "S_PIE","beta_S_PIE"))
   dat_samples <- dat_samples[order(dat_samples$index, dat_samples$n_rare, dat_samples$group),]
   
   dat_groups$index <- factor(dat_groups$index, levels = index)
   dat_groups <- dat_groups[order(dat_groups$index, dat_groups$n_rare, dat_groups$group),]
   
   #remove unused factor levels
   dat_samples$index <- factor(dat_samples$index)
   
   if (!boot_groups){
      
      p_val_groups$index <- factor(p_val_groups$index, levels = index)
      p_val_groups <- p_val_groups[order(p_val_groups$index),]
       
      out <- list(samples_stats = dat_samples,
                  groups_stats  = dat_groups,
                  samples_pval  = p_val_samples,
                  groups_pval  = p_val_groups,
                  p_min        = 1/n_perm)
   } else {
      out <- list(samples_stats = dat_samples,
                  groups_stats  = dat_groups,
                  samples_pval  = p_val_samples,
                  p_min        = 1/n_perm)
   }
  
   class(out) = 'mob_stats'
   return(out)
}

# Panel function for sample level results
samples_panel1 <- function(sample_dat, p_val, p_min, col, ylab = "",
                           main = "Sample scale", ...)
{
   if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
   else                          p_label <- bquote(p <= .(p_min))
   
   boxplot(value ~ group, data = sample_dat, main = main,
           ylab =  ylab, ylim = c(0, 1.1*max(sample_dat$value, na.rm = T)), 
           col = col, ...)
   mtext(p_label, side = 3, line = 0)  
}

# Panel function for group level results
groups_panel1 <- function(group_dat, p_val, p_min, col, ylab = "", main = "Group scale",
                          ...)
{
   if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
   else                          p_label <- bquote(p <= .(p_min))
   
   boxplot(value ~ group, data = group_dat, main = main,
           ylab = ylab, boxwex = 0, ylim = c(0, 1.1*max(group_dat$value, na.rm = T)),
           col = col, ...)
   points(value ~ group, data = group_dat, pch = 8, cex = 1.5, lwd = 2, col = col,
          ...)
   mtext(p_label, side = 3, line = 0)
}

# Panel function for group level results with confidence intervals
groups_panel2 <- function(group_dat, col, ylab = "", main = "Group scale", ...)
{
   boxplot(median ~ group, data = group_dat, main = main,
           ylab = ylab, boxwex = 0, ylim = c(0, 1.1*max(group_dat$upper)),
           col = col, ...)
   plotrix::plotCI(1:nrow(group_dat), group_dat$median, li = group_dat$lower,
                   ui = group_dat$upper, add = T, pch = 19, cex = 1.5, sfrac = 0.02,
                   col = col, ...)
}

#' Plot sample-level and group-level biodiversity statistics for a MoB analysis
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
#' is one figure for each index, with panels for sample-level and group-level results
#' as well as for beta-diversity when applicable. 
#' 
#' @param multi_panel A logical variable. If \code{multi_panel = TRUE} then a 
#' multipanel plot is produced, which shows observed, rarefied, and asymptotic 
#' species richness and S_PIE at the sample-level and the group-level.
#' This set of variables conveys a comprehensive picture of the underlying 
#' biodiversity changes. 
#' 
#' @param col a vector of colors for the two groups (control and treatment), set
#' to NA if no color is preferred
#' 
#' @param ... additional arguments to provide to boxplot, points, and confidence
#' interval functions
#' 
#' @author Felix May, Xiao Xiao, and Dan McGlinn 
#' 
#' @export
#' 
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' # without bootstrap CI for group scale
#' inv_stats = get_mob_stats(inv_mob_in, group_var = "group",
#'  ref_group = "uninvaded", n_perm = 20)
#' plot(inv_stats) 
#' 
#' windows(15,20)
#' plot(inv_stats, multi_panel = T)
#' # with bootstrap CI for group scale
#' inv_stats_boot = get_mob_stats(inv_mob_in, group_var = "group",
#'  ref_group = "uninvaded", n_perm = 20, boot_groups=T)
#' plot(inv_stats_boot)

plot.mob_stats = function(mob_stats, index = c("N","S","S_rare","S_asymp","S_PIE"),
                          multi_panel = FALSE, col=c("#2B83BA", "#FFC000"), ...)
{
   INDICES <- c("N", "S", "S_rare","S_asymp","PIE","S_PIE")
   
   if (multi_panel) index <- c("S","S_rare","S_asymp","S_PIE")
   index <- match.arg(index, INDICES, several.ok = TRUE)
   
   var_names <- levels(mob_stats$samples_stats$index)
   var_names2 <- var_names[var_names != "beta_S" & var_names != "beta_S_PIE"]
   
   index_match <- intersect(index, var_names)
   if (length(index_match) == 0)
      stop(paste("The indices", paste(index, collapse = ", "),"are missing in the input. Please choose other indices or re-run get_mob_stats with the indices of interest."))
   
   index_missing <- setdiff(index, var_names2)
   if (length(index_missing) > 0)
      warning(paste("The indices", paste(index, collapse = ", "),"are missing in the input and cannot be plotted."))
   
   if ("S_rare" %in% index_match){
      S_rare_samples <- filter(mob_stats$samples_stats, index == "S_rare")
      S_rare_groups <- filter(mob_stats$groups_stats, index == "S_rare")
      S_rare_len <- max(length(unique(S_rare_samples$n_rare)),
                        length(unique(S_rare_groups$n_rare)))
   } else {S_rare_len <- 0}
   
   if (multi_panel){
      n_rows <- 3 + S_rare_len
      op <- par(mfrow = c(n_rows,3), las = 1, cex.lab = 1.4, oma = c(0,1,0,0), xpd = NA)
   } 
   
   for (var in index_match){
      
      if (var %in% c("N","PIE")){
      
         if (!multi_panel)
            op <- par(mfrow = c(1,2), las = 1, cex.lab = 1.3,
                      oma = c(0,2,0,0), mar = c(4,3,5,1), xpd = NA)
         
         y_label <- switch(var,
                           "N" = expression(N),
                           "PIE" = expression(PIE))
         
         dat_samples <- filter(mob_stats$samples_stats, index == var)
         p_val <- with(mob_stats$samples_pval, p_val[index == var])
         samples_panel1(dat_samples, p_val = p_val, p_min = mob_stats$p_min,
                        ylab =  y_label, main = "Sample scale", col = col, ...)
         
         dat_groups <- filter(mob_stats$groups_stats, index == var)
         
         if (!is.null(mob_stats$groups_pval)){
            p_val <- with(mob_stats$groups_pval, p_val[index == var])
            groups_panel1(dat_groups, p_val = p_val, p_min = mob_stats$p_min, col = col, ...) 
         } else {
            groups_panel2(dat_groups, col = col, ...) 
         }
      }
      
      if (var %in% c("S", "S_asymp", "S_PIE")){
         
         if (!multi_panel)
            op <- par(mfrow = c(1,3), cex.lab = 1.6,
                      oma = c(0,2,0,0), mar = c(4,3,5,1), xpd = NA)
         
         if (multi_panel){
            if (var == "S_asymp") par(fig = c(0, 0.33, 1/n_rows, 2/n_rows), new = T)
            if (var == "S_PIE") par(fig = c(0, 0.33, 0       , 1/n_rows), new = T)
         }
         
         y_label <- switch(var,
                           "S" = expression(S),
                           "S_asymp" = expression(S[asymp]),
                           "S_PIE" = expression(S[PIE]))
         
         dat_samples <- filter(mob_stats$samples_stats, index == var)
         p_val <- with(mob_stats$samples_pval, p_val[index == var])
         samples_panel1(dat_samples, p_val = p_val, p_min = mob_stats$p_min,
                        ylab =  y_label, main = "Sample scale", col = col, ...)
         
         if (multi_panel){
            if (var == "S_asymp") par(fig = c(0.33, 0.67, 1/n_rows, 2/n_rows), new = T)
            if (var == "S_PIE") par(fig = c(0.33, 0.67, 0       , 1/n_rows), new = T)
         }
         
         beta_var <- paste("beta", var, sep = "_")
         dat_samples <- filter(mob_stats$samples_stats, index == beta_var)
         p_val <- with(mob_stats$samples_pval, p_val[index == beta_var])
         samples_panel1(dat_samples, p_val = p_val, p_min = mob_stats$p_min,
                        ylab =  "", main = "Beta-diversity across scales", col = col, ...)
         
         if (multi_panel){
            if (var == "S_asymp") par(fig = c(0.67, 1.0, 1/n_rows, 2/n_rows), new = T)
            if (var == "S_PIE") par(fig = c(0.67, 1.0, 0       , 1/n_rows), new = T)
         }
         
         dat_groups <- filter(mob_stats$groups_stats, index == var)
         
         if (!is.null(mob_stats$groups_pval)){
            p_val <- with(mob_stats$groups_pval, p_val[index == var])
            groups_panel1(dat_groups, p_val = p_val, p_min = mob_stats$p_min, col = col, ...) 
         } else {
            groups_panel2(dat_groups, col = col, ...) 
         }
      }
      
      if (var == "S_rare"){
         
         if (!multi_panel){
            op <- par(mfrow = c(S_rare_len, 3), las = 1, cex.lab = 1.6,
                      oma = c(0,2,0,0), mar = c(4,3,5,1), xpd = NA)
         }
         
         y_label <- expression(S[rare])
            
         n_rare_samples <- unique(S_rare_samples$n_rare)
         n_rare_groups <- unique(S_rare_groups$n_rare)
         
         for (i in 1:length(n_rare_samples)){
            
            if (multi_panel)
               par(fig = c(0, 0.33, (1+i)/n_rows, (2+i)/n_rows), new = T)
            
            dat_samples <- filter(S_rare_samples, n_rare == n_rare_samples[i])
            p_val <- with(mob_stats$samples_pval,
                          p_val[index == var & n_rare == n_rare_samples[i]])
            
            fig_title <- paste("Sample scale, n = ",n_rare_samples[i])
            
            samples_panel1(dat_samples, p_val = p_val, p_min = mob_stats$p_min,
                           ylab = y_label,
                           main = fig_title, col = col, ...)
            
            if (multi_panel)
               par(fig = c(0.33, 0.67, (1+i)/n_rows, (2+i)/n_rows), new = T)
            
            beta_var <- paste("beta", var, sep = "_")
            dat_samples <- filter(mob_stats$samples_stats, index == beta_var)
            p_val <- with(mob_stats$samples_pval,
                          p_val[index == beta_var & n_rare == n_rare_samples[i]])
            samples_panel1(dat_samples, p_val = p_val, p_min = mob_stats$p_min,
                           ylab = "", main = "Beta-diversity across scales", col = col, ...)
            
            if (multi_panel)
               par(fig = c(0.67, 1.0, (1+i)/n_rows, (2+i)/n_rows), new = T)
            
            dat_groups <- filter(S_rare_groups, n_rare == n_rare_groups[i])
            fig_title <- paste("Group scale, n = ", n_rare_groups[i])
            
            if (!is.null(mob_stats$groups_pval)){
               pval <- with(mob_stats$groups_pval,
                            p_val[index == var & n_rare == n_rare_groups[i]])
               groups_panel1(dat_groups, p_val = p_val, p_min = mob_stats$p_min,
                             ylab = "", main = fig_title, col = col, ...)
            } else {
               groups_panel2(dat_groups, main = fig_title, col = col, ...) 
            }
         }
         
         y_coords <- (S_rare_len:0)/S_rare_len
         
      }
   }
   
   par(op)
}