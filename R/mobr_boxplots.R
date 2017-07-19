#' Calculate probability of interspecific encounter (PIE)
#' 
#' PIE is also known as Simpson's evenness index and this function is 
#' a reduced form of the function vegan::diversity(). Jari Oksanen and Bob O'Hara
#' are the original authors of the function vegan::diversity().
#' 
#' In this function the formulate of Hurlbert (1971) is used to calculate PIE:
#' 
#' PIE = N/(N-1)*(1 - p_i^2),
#' 
#' where N is the total number of individuals and p_i is the relative abundance
#' of species i. 
#' 
#' 
#' @inheritParams rarefaction
#' @author Dan McGlinn
#' @keywords internal
#' 
#' @references 
#' Hurlbert, S. H. 1971. The Nonconcept of Species Diversity: A Critique and Alternative Parameters. - Ecology 52: 577–586.

#' 
#' @examples 
#' data(inv_comm)
#' calc_PIE(inv_comm)
calc_PIE = function(x) {
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
        H = apply(x, 1, sum, na.rm = TRUE)
    else 
        H = sum(x, na.rm = TRUE)
    H = total/(total - 1)* (1 - H)
    if (any(is.na(total))) 
        is.na(H) = is.na(total)
    H[!is.finite(H) | total == 0] <- NA # set NA, when total == 0 or total == 1
    return(H)
}

# Calculate biodiversity indices of single group
calc_biodiv_single_group <- function(abund_vec, index, n_rare){
   
   abund_vec <- as.numeric(abund_vec)
   abund_vec <- abund_vec[abund_vec > 0]
   
   out <- data.frame(var = factor(index[index != "S_rare"]), n_rare = NA, estimate = NA)
   if ("S_rare" %in% index){
      out <- rbind(out, expand.grid(var = "S_rare", n_rare = n_rare, estimate = NA)) 
   }   

      
   # Number of individuals -----------------------------------------------------
   if ("N" %in% index){
      out$estimate[out$var == "N"]  =  sum(abund_vec)   
   } 
   
   # Number of species ---------------------------------------------------------
   if ("S" %in% index){
      out$estimate[out$var == "S"]  = sum(abund_vec > 0)  
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if ("S_rare" %in% index){  
      inext1 <- iNEXT::iNEXT(abund_vec, q = 0, datatype = "abundance",
                                size = n_rare, se = F)
      out$estimate[out$var == "S_rare"] <- inext1$iNextEst$qD
   } # end rarefied richness
   
   
   # Asymptotic estimates species richness -------------------------------------
   if ("S_asymp" %in% index){
      
      S_asymp_group <- try(vegan::estimateR(abund_vec))
      if (class(S_asymp_group) == "try_error"){
         warning("The Chao richness estimator cannot be calculated for all groups.")
      } else {
         S_asymp_group = S_asymp_group["S.chao1"]
         S_asymp_group[!is.finite(S_asymp_group)] <- NA
      }
      
      dat_S_asymp <- data.frame(group = groups,
                                index = "S_asymp",
                                n_rare = NA,
                                value = S_asymp_group)
      dat_groups <- rbind(dat_groups, dat_S_asymp)
      
      out$estimate[out$var == "S_asymp"] <- S_asymp_group
   }
   
   # Probability of Interspecific Encounter (PIE)-------------------------------
   if ("PIE" %in% index){ 
      out$estimate[out$var == "PIE"]  = calc_PIE(abund_vec)
   }
   
   # Effective number of species based on PIE ----------------------------------
   if ("ENS_PIE" %in% index){
      ENS_PIE_groups = vegan::diversity(abund_vec, index = "invsimpson")
      ENS_PIE_groups[!is.finite(ENS_PIE_groups)] <- NA
      out$estimate[out$var == "ENS_PIE"] = ENS_PIE_groups
   }
   
   return(out)
}

# generate a single bootstrap sample of group-level biodiversity indices
boot_sample_group <- function(comm_dat, index, n_rare)
{
   n_sample <- nrow(comm_dat)
   sample_rows <- sample(1:n_sample, size = n_sample, replace = T)
   abund <- colSums(comm_dat[sample_rows, ])
   out <- calc_biodiv_single_group(abund, n_rare = n_rare, index = index)
   
   return(out)
}

# replicated bootstrap samples of group-level indices
repl_boot_samples <- function(comm_dat, index, n_rare, n_perm = 100)
{
   repl_biodiv <- replicate(n_perm, boot_sample_group(comm_dat, index = index,
                                                      n_rare = n_rare),
                            simplify = FALSE)
   return(repl_biodiv)
}

# get mean and confidence interval from bootstrap samples
get_mean_CI <- function(x, level = 0.95)
{
   alpha <- 1 - level
   p <- c(alpha/2, 0.5, 1 - alpha/2)
   mean_val <- rowMeans(x)
   CI_tab <- apply(x, MARGIN = 1, quantile, probs = p)
   CI_vec <- c(mean_val, CI_tab[1,], CI_tab[2,])
   names(CI_vec) <- paste(rep(colnames(CI_tab), 3),
                          rep(c("mean","low", "up"), each = ncol(CI_tab)), sep = "_")
   return(CI_vec)
}

# Get group-level p-values from permutation test
get_pval <- function(rand, obs, n_samples)
{
   n_extremes = sum(rand < -abs(obs) | rand > abs(obs))
   p_val = n_extremes/n_samples
   return(p_val)
}

# Get group-level biodiversity indices
calc_biodiv_groups <- function(abund_mat, groups, index, n_rare)
{
   dat_groups <- data.frame(group = factor(),
                            index = factor(),
                            n_rare = numeric(),
                            val = numeric())
   
   groups_N <- rowSums(abund_mat) 
   
   # Number of individuals -----------------------------------------------------
   if ("N" %in% index){
      dat_N <- data.frame(group = groups,
                          index = "N",
                          n_rare = NA,
                          value = groups_N)
      dat_groups <- rbind(dat_groups, dat_N)
   } 
   
   # Number of species ---------------------------------------------------------
   if ("S" %in% index){
      dat_S <- data.frame(group = groups,
                          index = "S",
                          n_rare = NA,
                          value = rowSums(abund_mat > 0) )
      dat_groups <- rbind(dat_groups, dat_S)
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if ("S_rare" %in% index){  
   
      for (i in 1:length(n_rare)){
         
         groups_low_n = groups_N < n_rare[i]
         
         if (sum(groups_low_n) == 1){
            warning(paste("There is",sum(groups_low_n),"group with less then",  n_rare[i],"individuals. This is removed for the calculation of rarefied richness."))
         }
         
         if (sum(groups_low_n) > 1){
            warning(paste("There are",sum(groups_low_n),"groups with less then",  n_rare[i],"individuals. These are removed for the calculation of rarefied richness."))
         }
         
         dat_S_rare <- data.frame(group = groups,
                                  index = "S_rare",
                                  n_rare = n_rare[i],
                                  value = rep(NA, nrow(abund_mat)))
         
         dat_S_rare$value[!groups_low_n] <- apply(abund_mat[!groups_low_n], MARGIN = 1,
                                                  rarefaction, method = "indiv",
                                                  effort = n_rare[i])
         dat_groups <- rbind(dat_groups, dat_S_rare)
      }
   } # end rarefied richness
   
   
   # Asymptotic estimates species richness -------------------------------------
   if ("S_asymp" %in% index){
      
      S_asymp_group <- try(vegan::estimateR(abund_mat))
      if (class(S_asymp_group) == "try_error"){
         warning("The Chao richness estimator cannot be calculated for all groups.")
      } else {
         S_asymp_group = S_asymp_group["S.chao1",]
         S_asymp_group[!is.finite(S_asymp_group)] <- NA
      }
    
      dat_S_asymp <- data.frame(group = groups,
                                index = "S_asymp",
                                n_rare = NA,
                                value = S_asymp_group)
      dat_groups <- rbind(dat_groups, dat_S_asymp)
   }
   
   # Probability of Interspecific Encounter (PIE)-------------------------------
   if ("PIE" %in% index){ 

      dat_PIE <- data.frame(group = groups,
                            index = "PIE",
                            n_rare = NA,
                            value = calc_PIE(abund_mat))
      dat_groups <- rbind(dat_groups, dat_PIE)
   }
   
   # Effective number of species based on PIE ----------------------------------
   if ("ENS_PIE" %in% index){
      
      ENS_PIE_groups = vegan::diversity(abund_mat, index = "invsimpson")
      ENS_PIE_groups[!is.finite(ENS_PIE_groups)] <- NA
      
      dat_ENS_PIE <- data.frame(group = groups,
                                index = "ENS_PIE",
                                n_rare = NA,
                                value = ENS_PIE_groups)
      dat_groups <- rbind(dat_groups, dat_ENS_PIE)
   }
   
   return(dat_groups)
   
}

#Get F statistics from diversity indices and grouping vector
get_F_values <- function(div_dat, permute = F)
{
   if (permute)
      div_dat <- mutate(div_dat, group = sample(group))
   
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
   
   dat_groups <- calc_biodiv_groups(abund_mat = abund_group[,-1],
                                    groups = abund_group[,1],
                                    index = index,
                                    n_rare = n_rare)
   dat_groups <- dat_groups %>%
      group_by(index, n_rare) %>%
      summarise(diff = value[group == "treatment"] - value[group == "control"]) %>%
      ungroup()
   
   return(dat_groups)
}


#' Calculate sample based and group based biodiversity statistics.
#' @inheritParams get_delta_stats
#' 
#' @param group_var A string that specifies which field in \code{mob_in$env} the data
#'   should be grouped by
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
#'    \item \code{ENS_PIE} ... Effective number of species based on PIE
#' }
#' See \emph{Details} for additional information on the biodiverstiy statistics.
#' 
#' @n_rare_samples The standardized number of individuals used for the calculation
#' of rarefied species richness at the sample level. This can a single value
#' or an integer vector.
#' 
#' @n_rare_groups The standardized number of individuals used for the calculation
#' of rarefied species richness at the group level. This can a single value
#' or an integer vector.
#' 
#' @n_perm The number of permutations to use for testing for treatment effects
#' 
#' @details 
#' \strong{Rarefied species richness} is the expected number of species, given
#' a defined number of sampled individuals (Gotelli & Colwell 2001). Rarefied richness
#' at the sample and group level is calculated for the values provided in
#' \code{n_rare_samples} and \code{n_rare_groups}, respectively.
#' When no values are provided the minimum number of individuals of the samples 
#' or groups is used. \code{S_rare} is calculated using the function \code{\link{rarefaction}}.
#' 
#' \strong{Asymptotic species richness} is calculated using the bias-corrected Chao1
#'  estimator (Chiu et al. 2014) provided in the function \code{\link[vegan]{estimateR}}.
#' 
#' \strong{PIE} represents the probility that two randomly drawn individuals
#' belong to the same species. Here we use the definition of Hurlbert (1971), which considers
#' sampling without replacement. PIE is closely related to the well-known Simpson 
#' diversity index, but the latter assumes sampling with replacement.
#' 
#' \strong{ENS_PIE} represents the Effective Number of Species derived from PIE.
#' This corresponds to the the number of equally abundant species, which together result
#' in the same PIE as the observed community, with - most likely - uneven abundances
#' (Jost 2006). That means the higher the difference between S and ENS_PIE the more
#' uneven is the observed community. An intuitive interpretation of ENS_PIE is that
#' it corresponds to the number of dominant species in the community.
#' 
#' For species richness \code{S} and the Effective Number of Species \code{ENS_PIE}
#' we also calculate beta-diversity using multiplicative partitioning
#' (Whittaker 1972, Jost 2007). That means for both indices we estimate beta-diversity as the ratio of group-level
#' diversity (gamma) divided by sample-level diversity (alpha). 
#' 
#' We used permutation tests for assessing differences of the biodiversity statistics
#' among the groups (Legendre & Legendre 1998). At the sample level, one-way ANOVA 
#' (i.e. F-test) is implemented by shuffling treatment group labels across samples.
#' 
#'  At the group-level there is test whether the difference between treatment and control is significantly different from zero. For this purpose we permuted treatment group labels across samples, pooled the groups, and calculated the treatment-control difference for each permutation.
#'  
#' 
#' @return A list of class \code{mob_stats} that contains sample-scale and 
#' group-scale biodiversity statistics, as well as the p-values for permutation tests
#' at both scales.
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
#' 
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, group_var = "group", ref_group = "uninvaded", n_perm = 100, n_rare_samples = c(5,10))
#' plot(inv_stats)
get_mob_stats = function(mob_in,
                         group_var,
                         ref_group = NULL,
                         index = c("N","S","S_rare","S_asymp","ENS_PIE"),
                         n_rare_samples = NULL,
                         n_rare_groups = NULL,
                         n_perm = 100
                         )
{
   if (n_perm < 1) 
       stop('Set nperm to a value greater than 1') 
   
   INDICES <- c("N", "S", "S_rare","S_asymp","PIE","ENS_PIE")
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
  
   # Add labels to groups
   group_type = c("(ctrl)",rep("(trt)", length(levels(group_id)) - 1))
   group_labels <- paste(levels(group_id), group_type)
   group_id <- factor(group_id, labels = group_labels) 
   
   print(index)
 
   # Abundance distribution pooled in groups
   abund_group = aggregate(mob_in$comm, by = list(group_id), FUN = "sum")
  
   # Group-level indices
   
   # Get rarefaction level
   groups_N <- rowSums(abund_group[,-1]) 
   if (any(is.null(n_rare_groups)) | !is.numeric(n_rare_groups)){   
      N_min_group = min(groups_N)
      n_rare_groups = N_min_group
   } else {
      n_rare_groups <- floor(n_rare_groups)
   }
   
   dat_groups <- calc_biodiv_groups(abund_mat = abund_group[,-1],
                                    groups = abund_group[,1],
                                    index = index,
                                    n_rare = n_rare_groups)
   
   # Sample-level indices
   dat_samples <- data.frame(group = factor(),
                             index = factor(),
                             n_rare = numeric(),
                             val = numeric()
                             )
   
   samples_N = rowSums(mob_in$comm) 
   
   # Number of individuals -----------------------------------------------------
   if ("N" %in% index){
      dat_N <- data.frame(group = group_id,
                          index = "N",
                          n_rare = NA,
                          value = samples_N)
      dat_samples <- rbind(dat_samples, dat_N)
   } 
   
   # Number of species ---------------------------------------------------------
   if ("S" %in% index){
      dat_S <- data.frame(group = group_id,
                          index = "S",
                          n_rare = NA,
                          value = rowSums(mob_in$comm > 0))
      
      beta_S <- dat_groups$value[dat_groups$index == "S"][group_id]/dat_S$value
      dat_betaS <- data.frame(group = group_id,
                              index = "beta_S",
                              n_rare = NA,
                              value = beta_S)
      
      dat_betaS$value[!is.finite(dat_betaS$value)] <- NA
      dat_samples <- rbind(dat_samples, dat_S, dat_betaS)
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if ("S_rare" %in% index){  
      
      if (any(is.null(n_rare_samples)) | !is.numeric(n_rare_samples)){ 
         N_min_sample = min(samples_N)
         n_rare_samples = N_min_sample
      }
      else {
         n_rare_samples <- floor(n_rare_samples)
      }
      
      if (any(n_rare_samples <= 1)){
         warning(paste("Comparisons of rarefied richness do not make sense for",
                       n_rare_samples,"individuals. Please choose higher value for \"n_rare_samples\"."))
      }
      
      for (i in 1:length(n_rare_samples)){
         
         plots_low_n = samples_N < n_rare_samples[i]
         
         if (sum(plots_low_n) == 1){
            warning(paste("There is",sum(plots_low_n),"plot with less then", n_rare_samples[i],"individuals. This plot are removed for the calculation of rarefied richness."))
         }
         
         if (sum(plots_low_n) > 1){
            warning(paste("There are",sum(plots_low_n),"plots with less then", n_rare_samples[i],"individuals. These plots are removed for the calculation of rarefied richness."))
         }
         
         dat_S_rare <- data.frame(group = group_id,
                                  index = "S_rare",
                                  n_rare = n_rare_samples[i],
                                  value = rep(NA, nrow(mob_in$comm)))
         
         dat_S_rare$value[!plots_low_n] = apply(mob_in$comm[!plots_low_n,], MARGIN = 1,
                                                rarefaction, method = "indiv",
                                                effort = n_rare_samples[i])
         dat_samples <- rbind(dat_samples, dat_S_rare)
      }
   } # end rarefied richness
   
   
  # Asymptotic estimates species richness -------------------------------------
  if ("S_asymp" %in% index){
     
      S_asymp_sample = try(vegan::estimateR(mob_in$comm))
      if (class(S_asymp_sample) == "try_error"){
         warning("The Chao richness estimator cannot be calculated for all samples.")
         S_asymp_sample <- rep(NA, nrow(mob_in$comm))
      } else {
         S_asymp_sample = S_asymp_sample["S.chao1",]
         S_asymp_sample[!is.finite(S_asymp_sample)] <- NA
      }
      
      dat_S_asymp <- data.frame(group = group_id,
                                index = "S_asymp",
                                n_rare = NA,
                                value = S_asymp_sample)
      dat_samples <- rbind(dat_samples, dat_S_asymp)
  }

   # Probability of Interspecific Encounter (PIE)-------------------------------
   if ("PIE" %in% index){ 
      
      plots_n01 = samples_N <= 1 # Hurlbert's PIE can only be calculated for two or more individuals
      
      if (sum(plots_n01) == 1){
         warning(paste("There is",sum(plots_n01), "plot with less than two individuals.
   This is removed for the calculation of PIE."))
      }
      
      if (sum(plots_n01) > 1){
         warning(paste("There are",sum(plots_n01), "plots with less than two individuals.
   These are removed for the calculation of PIE."))
      }
      
      dat_PIE <- data.frame(group = group_id,
                            index = "PIE",
                            n_rare = NA,
                            value = calc_PIE(mob_in$comm))
      dat_samples <- rbind(dat_samples, dat_PIE)
   }
   
   # Effective number of species based on PIE ----------------------------------
   if ("ENS_PIE" %in% index){
      
      plots_n01 = samples_N <= 1

      ENS_PIE_samples = vegan::diversity(mob_in$comm, index = "invsimpson")
      ENS_PIE_samples[plots_n01 | !is.finite(ENS_PIE_samples)] <- NA
      
      dat_ENS_PIE <- data.frame(group = group_id,
                                index = "ENS_PIE",
                                n_rare = NA,
                                value = ENS_PIE_samples)
      
      beta_ENS_PIE <-
         dat_groups$value[dat_groups$index == "ENS_PIE"][group_id]/dat_ENS_PIE$value
      beta_ENS_PIE[!is.finite(beta_ENS_PIE)] <- NA
      
      dat_beta_ENS_PIE <- data.frame(group = group_id,
                                     index = "beta_ENS_PIE",
                                     n_rare = NA,
                                     value = beta_ENS_PIE)
      dat_samples <- rbind(dat_samples, dat_ENS_PIE, dat_beta_ENS_PIE)
   }
   
   # Significance tests
   
   # sample level
   F_obs <- get_F_values(dat_samples, permute = F)
   F_rand <- dplyr::bind_rows(replicate(n_perm, get_F_values(dat_samples, permute = T),                                           simplify = F)) %>% ungroup()
   F_obs <- F_obs %>% mutate(F_val_obs = F_val,
                             F_val = NULL)
   F_rand <- left_join(F_rand, F_obs)
   
   p_val_samples <- F_rand %>% 
      group_by(index, n_rare) %>%
      summarise(p_val = sum(F_val_obs <= F_val)/n_perm) %>%
      ungroup()
   
   # group level
   diff_obs <- get_group_diff(mob_in$comm, group_bin, index, n_rare = n_rare_groups,                                     permute = F)
   diff_rand <- bind_rows(replicate(n_perm, get_group_diff(mob_in$comm, group_bin,
                                                           index,
                                                           n_rare = n_rare_groups,
                                                           permute = T),
                                    simplify = F))
   diff_obs <- diff_obs %>% mutate(d_obs = diff,
                                   diff = NULL)
   diff_rand <- left_join(diff_rand, diff_obs)
   
   p_val_groups <- diff_rand %>% 
      group_by(index, n_rare) %>%
      summarise(p_val = get_pval(rand = diff, obs = first(d_obs), n_samples = n_perm)) %>%
      ungroup()

   # ---------------------------------------------------------------------------
   # # old code - might be re-used later
   # 
   # # Group-based statistics
   # group_stats <- apply(abund_group[ ,-1], MARGIN = 1,
   #                      FUN = calc_biodiv_single_group, n_rare = rarefy_groups)
   # row.names(group_stats) <- paste(row.names(group_stats), "_obs", sep = "")
   # 
   # # bootstrap replicates
      # boot_samples <- by(mob_in$comm, INDICES = group_id, FUN = repl_boot_samples,
   #                    n_rare = n_rare_groups, n_perm = n_perm, index = index,
   #                    simplify = TRUE)

   # boot_samples <- by(mob_in$comm, INDICES = group_id, FUN = repl_boot_samples,
   #                    n_rare = rarefy_groups)
   # 
   # boot_CI <- sapply(boot_samples, get_mean_CI)
   # group_dat <- as.data.frame(t(rbind(group_stats, boot_CI)))
   # group_dat <- group_dat[,c("N_obs","N_mean","N_low","N_up",
   #                           "S_obs","S_mean","S_low","S_up",
   #                           "S_rare1_obs","S_rare1_mean","S_rare1_low","S_rare1_up",
   #                           "S_rare2_obs","S_rare2_mean","S_rare2_low","S_rare2_up",
   #                           "S_rare3_obs","S_rare3_mean","S_rare3_low","S_rare3_up",
   #                           "S_asymp_obs","S_asymp_mean","S_asymp_low","S_asymp_up",
   #                           "PIE_obs","PIE_mean","PIE_low","PIE_up",
   #                           "ENS_PIE_obs","ENS_PIE_mean","ENS_PIE_low","ENS_PIE_up")]
   # group <- factor(levels(group_id))
   # group_dat <- cbind(group, group_dat)
   # rownames(group_dat) <- NULL
   # 
   # # permutation tests
   # 
   
   out <- list(samples_stats = dat_samples,
               groups_stats  = dat_groups,
               samples_pval  = p_val_samples,
               groups_pval  = p_val_groups,
               n_perm        = n_perm)
  
   class(out) = 'mob_stats'
   return(out)
}

#' Plot sample-level and group-level biodiversity statistics for a MoB analysis
#' 
#' Plots a \code{mob_stats} object which is produced by the 
#' function \code{get_mob_stats}. The p value for each statistic
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
#' as well as beta-diversity for species richness \code{S} and the effective
#' number of species \code{ENS_PIE}. 
#' 
#' @param multi_panel A logical variable. If \code{multi_panel = TRUE) then a 
#' multipanel plot is produced, which shows observed, rarefied, and asymptotic 
#' species richness and ENS_PIE at the sample-level and the group-level.
#' This set of variables conveys a comprehensive picture of the underlying 
#' biodiversity changes. 
#' 
#' @author Felix May, Xiao Xiao, and Dan McGlinn 
#' 
#' @export
#' 
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, group_var = "group", ref_group = "uninvaded", n_perm = 100)
#' plot(inv_stats) 

plot.mob_stats = function(mob_stats, index = c("N","S","S_rare","S_asymp","ENS_PIE"),
                          multi_panel = FALSE)
{
   INDICES <- c("N", "S", "S_rare","S_asymp","PIE","ENS_PIE")
   
   if (multi_panel) index <- c("S","S_rare","S_asymp","ENS_PIE")
   index <- match.arg(index, INDICES, several.ok = TRUE)
   
   var_names <- levels(mob_stats$samples_stats$index)
   var_names2 <- var_names[var_names != "beta_S" & var_names != "beta_ENS_PIE"]
   
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
      op <- par(mfrow = c(n_rows,3), las = 1)
   } 
   
   for (var in index_match){
      
      if (var %in% c("N", "S_asymp","PIE")){
         
         if (!multi_panel)
            op <- par(mfrow = c(1,2), las = 1, cex.lab = 1.2)
         
         dat_samples <- filter(mob_stats$samples_stats, index == var)
         
         p_val <- with(mob_stats$samples_pval, p_val[index == var])
         if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
         else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
         
         if (multi_panel)
            par(fig = c(0, 0.33, 1/n_rows, 2/n_rows),new = T)
         boxplot(value ~ group, data = dat_samples, main = "Sample scale",
                 ylab =  var, ylim = c(0, 1.1*max(dat_samples$value, na.rm = T)))
         mtext(p_label, side = 3, line = 0)
         
         dat_groups <- filter(mob_stats$groups_stats, index == var)
         
         p_val <- with(mob_stats$groups_pval, p_val[index == var])
         if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
         else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
         
         if (multi_panel)
            par(fig = c(0.67, 1.0, 1/n_rows, 2/n_rows),new = T)
         boxplot(value ~ group, data = dat_groups, main = "Group scale",
                 ylab = "", boxwex = 0, ylim = c(0, 1.1*max(dat_groups$value, na.rm = T)))
         points(value ~ group, data = dat_groups, pch = 8, cex = 1.5, lwd = 2)
         mtext(p_label, side = 3, line = 0)
      }
      
      if (var %in% c("S", "ENS_PIE")){
         if (!multi_panel)
            op <- par(mfrow = c(1,3), las = 1, cex.lab = 1.2)
         
         dat_samples <- filter(mob_stats$samples_stats, index == var)
         
         p_val <- with(mob_stats$samples_pval, p_val[index == var])
         if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
         else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
        
         if (multi_panel & var == "ENS_PIE")
            par(fig = c(0, 0.33, 0, 1/n_rows),new = T)
         boxplot(value ~ group, data = dat_samples, main = "Sample scale",
                 ylab =  var, ylim = c(0, 1.1*max(dat_samples$value, na.rm = T)))
         mtext(p_label, side = 3, line = 0)
         
         beta_var <- paste("beta", var, sep = "_")
         dat_samples <- filter(mob_stats$samples_stats, index == beta_var)
         
         p_val <- with(mob_stats$samples_pval, p_val[index == var])
         if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
         else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
         
         if (multi_panel & var == "ENS_PIE")
            par(fig = c(0.33, 0.67, 0, 1/n_rows),new = T)
         boxplot(value ~ group, data = dat_samples, main = "Beta-diversity across scales",
                 ylim = c(0, 1.1*max(dat_samples$value, na.rm = T)))
         mtext(p_label, side = 3, line = 0)
         
         dat_groups <- filter(mob_stats$groups_stats, index == var)
         
         p_val <- with(mob_stats$groups_pval, p_val[index == var])
         if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
         else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
        
         if (multi_panel & var == "ENS_PIE")
            par(fig = c(0.67, 1.0, 0, 1/n_rows), new = T)
         boxplot(value ~ group, data = dat_groups, main = "Group scale",
                 ylab = "", boxwex = 0, ylim = c(0, 1.1*max(dat_groups$value, na.rm = T)))
         points(value ~ group, data = dat_groups, pch = 8, cex = 1.5, lwd = 2)
         mtext(p_label, side = 3, line = 0)
      }
      
      if (var == "S_rare"){
         
         if (!multi_panel){
            op <- par(mfcol = c(S_rare_len, 2), las = 1, cex.lab = 1.2)
         }
            
         n_rare_samples <- unique(S_rare_samples$n_rare)
         n_rare_groups <- unique(S_rare_groups$n_rare)
         
         for (j in 1:length(n_rare_samples)){
            
            dat_samples <- filter(S_rare_samples, n_rare == n_rare_samples[j])
            
            pval <- with(mob_stats$samples_pval,
                         p_val[index == var & n_rare == n_rare_samples[j]])
            if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
            else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
            
            fig_title <- paste("Sample scale, n = ",n_rare_samples[j])
            if (multi_panel)
               par(fig = c(0, 0.33,  (1+j)/n_rows, (2+j)/n_rows),new = T)
            boxplot(value ~ group, data = dat_samples, main = fig_title,
                    ylab =  "S_rare", ylim = c(0, 1.1*max(dat_samples$value, na.rm = T)))  
            mtext(p_label, side = 3, line = 0)
         }
         
         y_coords <- (S_rare_len:0)/S_rare_len
         
         for (j in 1:length(n_rare_groups)){
            
            dat_groups <- filter(S_rare_groups, n_rare == n_rare_groups[j])
            
            pval <- with(mob_stats$groups_pval,
                         p_val[index == var & n_rare == n_rare_groups[j]])
            if (p_val > 0 | is.na(p_val)) p_label <- bquote(p == .(p_val))
            else                          p_label <- bquote(p <= .(1/mob_stats$n_perm))
            
            fig_title <- paste("Group scale, n = ", n_rare_groups[j])
            
            if (!multi_panel) par(fig = c(0.5, 1.0, y_coords[j+1], y_coords[j]), new = T)
            else              par(fig = c(0.67, 1.0, (1+j)/n_rows, (2+j)/n_rows),new = T)
            boxplot(value ~ group, data = dat_groups, main = fig_title,
                    ylab =  "", boxwex = 0,
                    ylim = c(0, 1.1*max(dat_groups$value, na.rm = T)))
            points(value ~ group, data = dat_groups, pch = 8, cex = 1.5, lwd = 2)
            mtext(p_label, side = 3, line = 0)
         }
      }
   }
   
   par(op)
}


# old code ---------------------------------------------------------------------
# still kept for potential later use
# plot.mob_stats = function(mob_stats, 
#                           display=c('all','N', 'S', 'S_rare', 'S_asymp',
#                                     'PIE', 'ENS_PIE'),
#                           multipanel=FALSE, ...) {
#       if (multipanel) {
#           if ('all' %in% display) {
#               display = c('PIE', 'N', 'S', 'S asymp')
#               op = par(mfcol = c(3,5), las = 1, font.main = 1,
#                        mar = c(2, 2, 3, 2) + 0.1)
#           }
#           else
#               warning('The multipanel plots settings only make sense when display is set to "all"')
#       }# else if ('all' %in% display)
#        #   display = c('PIE', 'N', 'S', 'S asymp')
#       # PIE
#       if ('PIE' %in% display) {
#           if (multipanel)
#               par(fig = c(0.2, 0.4, 0.66, 1))
#           else
#              par(mfrow = c(1,2))
#           # panel_title = paste("PIE (p = ", mob_stats$pvalues$PIE,")\nplot scale",
#           #                     sep = "")
#           boxplot(PIE ~ group, data=mob_stats$samples, main = "PIE samples")
#           
#           if (multipanel)
#               par(fig = c(0.6, 0.8, 0.66, 1), new = T)
#           boxplot(PIE ~ group, data=mob_stats$groups, main = "PIE groups",boxwex = 0)
#           points(PIE ~ group, data=mob_stats$groups, pch = 19)
#       }
#       
#       # ENS PIE
#       if ('ENS_PIE' %in% display) {
#          if (multipanel)
#             par(fig = c(0.2, 0.4, 0.66, 1))
#          else
#             par(mfrow = c(1,3))
#          # panel_title = paste("ENS PIE (p = ", mob_stats$pvalues$ENS_PIE,")\nplot scale",
#          #                     sep = "")
#          boxplot(ENS_PIE ~ group, data=mob_stats$samples,
#                  main = "ENS PIE samples")
#          
#          # beta ENS PIE
#          # panel_title = paste("beta ENS PIE (p = ", mob_stats$pvalues$beta_ENS_PIE,")\n between scales",
#          #                     sep = "")
#          boxplot(beta_ENS_PIE ~ group, data=mob_stats$samples,
#                  main = "beta ENS PIE")
#          
#          if (multipanel)
#             par(fig = c(0.6, 0.8, 0.66, 1), new = T)
#          boxplot(ENS_PIE ~ group, data=mob_stats$groups, main = "ENS PIE groups",boxwex = 0)
#          points(ENS_PIE ~ group, data=mob_stats$groups, pch = 19)
#         
#       }
#       
#       if ('S' %in% display) {
#           if (multipanel)
#               par(fig = c(0, 0.2, 0.33, 0.66), new = T)
#           else
#               par(mfrow = c(1,3))
# 
#           # S observed samples
#           #panel_title = paste("S (p = ", mob_stats$pvalues$S,")\n plot scale",
#           #                    sep = "")
#           boxplot(S ~ group, data=mob_stats$samples,
#                   main = "S samples")
#           
#           # S beta
#           #panel_title = paste("beta S (p = ", mob_stats$pvalues$beta_S,")\n between scales",
#           #                    sep = "")
#           boxplot(beta_S ~ group, data = mob_stats$samples,
#                   main = "beta S")
#           
#           # S groups
#           if (multipanel)
#               par(fig = c(0.8, 1.0, 0.33, 0.66), new = T)
#           boxplot(S ~ group, data=mob_stats$groups, main = "S groups",boxwex = 0)
#           points(S ~ group, data=mob_stats$groups, pch = 19)
#           
# 
#       }
#       # if ('S rare' %in% display) {
#       #    
#       #    n_rows <- max(length(mob_stats$samples$S_rare), length(mob_stats$samples$S_rare)) 
#       #    
#       #    if (multipanel)
#       #         par(fig = c(0, 0.2, 0.33, 0.66), new = T)
#       #     else
#       #         par(mfrow = c(1,2))
#       # 
#       #     # S rarefied samples
#       #     panel_title = paste("S rarefied (p = ", mob_stats$pvalues$S_rare,
#       #                         ")\n plot scale", sep = "")
#       #     boxplot(S_rare ~ group, data=mob_stats$samples,
#       #             main = panel_title, ...)
#       #     
#       #     # S groups
#       #     if (multipanel)
#       #         par(fig = c(0.8, 1.0, 0.33, 0.66), new = T)
#       #     boxplot(S ~ group, data=mob_stats$groups, main = "S groups",boxwex = 0)
#       #     points(S ~ group, data=mob_stats$groups, pch = 19)
#       # 
#       # }
#       if ('N' %in% display) {
#           if (multipanel)
#               par(fig = c(0.2, 0.4, 0.33, 0.66), new = T)
#           else
#               par(mfrow = c(1,2))
#           # N samples
#           #panel_title = paste("N (p = ", mob_stats$pvalues$N,")\n plot scale",
#           #                    sep = "")
#           boxplot(N ~ group, data=mob_stats$samples, 
#                   main = "N samples")
#       
#           # N groups
#           if (multipanel)
#               par(fig = c(0.6, 0.8, 0.33, 0.66), new = T)
#           # y_limits <- with(mob_stats$groups, range(c(N_mean, N_low, N_obs, N_up)))
#           # boxplot(N_mean ~ levels(group), data=mob_stats$groups,
#           #         main = "N\ntreatment scale", boxwex = 0, ylim = y_limits)
#           # points(N_obs ~ group, data = mob_stats$groups, pch = 19)
#           # plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$N_mean,
#           #        li = mob_stats$groups$N_low, ui = mob_stats$groups$N_up,
#           #        add = T, pch = 1, ...)
#           boxplot(N ~ group, data=mob_stats$groups, main = "N groups",boxwex = 0)
#           points(N ~ group, data=mob_stats$groups, pch = 19)
#       }
#       if ('S asymp' %in% display) {
#           if (multipanel)
#               par(fig = c(0.2, 0.4, 0, 0.33), new = T)
#           else
#               par(mfrow = c(1,2))
#            # S asymptotic
#           if (all(is.na(mob_stats$samples$S_asymp))){
#               warning("Cannot plot asymptotic richness for the samples")
#           } else {
#               # panel_title = paste("S asympotic (p = ", mob_stats$pvalues$S_asymp,
#               #                     ")\n plot scale", sep = "")
#               boxplot(S_asymp ~ group, data=mob_stats$samples,
#                       main = "S asymptotic")
#           }
#           
#           if (multipanel)
#               par(fig = c(0.6, 0.8, 0, 0.33), new = T)
#           if (all(is.na(mob_stats$groups$S_asymp))){
#               warning("Cannot plot asymptotic richness for the groups")
#           } else {
#              boxplot(S_asymp ~ group, data=mob_stats$groups, main = "S asymptotic groups",boxwex = 0)
#              points(S_asymp ~ group, data=mob_stats$groups, pch = 19)
#           }
#        
#       }
#       if (multipanel)
#           par(op)
# }
# 

