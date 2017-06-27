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
    H[!is.finite(H)] <- NA
    return(H)
}

# Calculate biodiversity indices of single group

calc_biodiv_single_group <- function(abund_vec, n_rare){
   
   abund_vec <- abund_vec[abund_vec > 0]
   
   out_vec <- c("N"       = NA,
                "S"       = NA,
                "S_rare1"  = NA,
                "S_rare2"  = NA,
                "S_rare3"  = NA,
                "S_asymp" = NA,
                "PIE"     = NA,
                "ENS_PIE" = NA
                )

   if (length(abund_vec) > 0){
      out_vec["N"] <- sum(abund_vec)
      out_vec["S"] <- length(abund_vec)
      
      inext1 <- iNEXT::iNEXT(abund_vec, q = 0, datatype = "abundance", size = n_rare, se = F)
      out_vec["S_rare1"] <- inext1$iNextEst$qD[1] 
      out_vec["S_rare2"] <- inext1$iNextEst$qD[2] 
      out_vec["S_rare3"] <- inext1$iNextEst$qD[3] 
      
      if (sum(abund_vec) > 1){
         out_vec["PIE"] <- calc_PIE(abund_vec)
         out_vec["ENS_PIE"] <- 1/(1 - out_vec["PIE"])
      }
      S_asymp <- try(vegan::estimateR(abund_vec))
      if (class(S_asymp) == "try_error"){
         warning("The Chao richness estimator cannot be calculated for all groups.")
      } else {
         out_vec["S_asymp"] = S_asymp["S.chao1"]
      }
   }
   
   return(out_vec)
}

# generate a single bootstrap sample of group-level biodiversity indices
boot_sample_group <- function(comm_dat, n_rare)
{
   n_sample <- nrow(comm_dat)
   sample_rows <- sample(1:n_sample, size = n_sample, replace = T)
   abund <- colSums(comm_dat[sample_rows, ])
   out <- calc_biodiv_single_group(abund, n_rare = n_rare)
   
   return(out)
}

# replicated bootstrap samples of group-level indices
repl_boot_samples <- function(comm_dat, n_rare, n_perm = 100)
{
   repl_biodiv <- replicate(n_perm, boot_sample_group(comm_dat, n_rare = n_rare))
   return(repl_biodiv)
}

# get mean and confidence interval from bootstrap samples
get_mean_CI <- function(x, level = 0.95)
{
   alpha <- 1 - level
   p <- c(alpha/2, 1 - alpha/2)
   mean_val <- rowMeans(x)
   CI_tab <- apply(x, MARGIN = 1, quantile, probs = p)
   CI_vec <- c(mean_val, CI_tab[1,], CI_tab[2,])
   names(CI_vec) <- paste(rep(colnames(CI_tab), 3),
                          rep(c("mean","low", "up"), each = ncol(CI_tab)), sep = "_")
   return(CI_vec)
}


#' Calculate sample based and group based biodiversity statistics
#' @inheritParams get_delta_stats
#' @param group_var a string that specifies which field in mob_in$env the data
#'   should be grouped by
#' @param n_min the minimum number of individuals a plot must have to be included
#' in the estimate of rarefied richness
#' @param nperm the number of permutations to use for testing for treatment effects
#' @return a list of class \code{mob_stats} that contains p-values,
#'  sample-scale (i.e., plot) statistics, and treatment scale statistics
#' @author Felix May and Dan McGlinn
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, 'group', nperm=100)
#' plot(inv_stats, multipanel=TRUE)
get_mob_stats = function(mob_in, group_var, ref_group = NULL, n_min = 5, nperm = 100) {
   if (nperm < 1) 
       stop('Set nperm to a value greater than 1') 
   group_id  = factor(mob_in$env[, group_var])
   
   # Sample-based statistics
   N_sample = rowSums(mob_in$comm)      # individuals in each sample
   S_sample = rowSums(mob_in$comm > 0)  # species in each sample
   
   # rarefied richness
   N_min_sample = min(N_sample)
   rarefy_samples = c(n_min, floor(N_min_sample)/2, N_min_sample)
   rarefy_samples[rarefy_levels < n_min] <- n_min
   plots_low_n = N_sample < n_min
   
   if (sum(plots_low_n) > 0){
      warning(paste("There are",sum(plots_low_n),"plots with less then", n_min,"individuals. These plots are removed for the calculation of rarefied richness."))
   }

   S_rare = data.frame(S_rare1 = rep(NA, nrow(mob_in$comm)),
                       S_rare2 = rep(NA, nrow(mob_in$comm)),
                       S_rare3 = rep(NA, nrow(mob_in$comm)))
   
   S_rare[!plots_low_n,] <-  t(apply(mob_in$comm[!plots_low_n,], MARGIN = 1,
                                      FUN = rarefaction, method = "indiv",
                                      effort = rarefy_levels))

   
   # Probability of Interspecific Encounter
   plots_n01 = N_sample <= 1 # Hurlbert's PIE can only be calculated for two or more individuals
   
   if (sum(plots_n01) == 1){
      warning(paste("There is",sum(plots_n01), "plot with less than two individuals.
This is removed for the calculation of PIE."))
   }
   
   if (sum(plots_n01) > 1){
      warning(paste("There are",sum(plots_n01), "plots with less than two individuals.
These are removed for the calculation of PIE."))
   }
   
   PIE_sample = calc_PIE(mob_in$comm)
   ENS_PIE_sample = 1/(1 - PIE_sample)
   ENS_PIE_sample[!is.finite(ENS_PIE_sample)] <- NA
   
   PIE_sample[plots_n01] = NA
   ENS_PIE_sample[plots_n01] = NA
   
   # bias corrected Chao estimator
   S_asymp_sample = try(vegan::estimateR(mob_in$comm))
   if (class(S_asymp_sample) == "try_error"){
      warning("The Chao richness estimator cannot be calculated for all samples.")
      S_asymp_sample <- rep(NA, nrow(mob_in$comm))
   } else {
      if (is.matrix(S_asymp_sample)) S_asymp_sample = S_asymp_sample["S.chao1",]
      else                           S_asymp_sample = S_asymp_sample["S.chao1"]
      S_asymp_sample[!is.finite(S_asymp_sample)] <- NA
   }
   
   # ---------------------------------------------------------------------------
   # abundance distribution pooled in group
   abund_group = aggregate(mob_in$comm, by = list(group_id), FUN = "sum")
   
   # Group-based statistics
   N_group = rowSums(abund_group[ ,-1])      

   N_min_group = min(N_group)
   rarefy_groups = c(n_min, floor(N_min_group)/2, N_min_group)
   rarefy_groups[rarefy_groups < n_min] <- n_min
   
   groups_low_n = N_group < n_min
   
   if (sum(groups_low_n) == 1){
      warning(paste("There is",sum(plots_low_n),"group with less then", n_min,"individuals.
These are removed for the calculation of rarefied richness."))
   }
   
   if (sum(groups_low_n) > 1){
      warning(paste("There are",sum(plots_low_n),"groups with less then", n_min,"individuals.
These are removed for the calculation of rarefied richness."))
   }
   
   group_stats <- apply(abund_group[ ,-1], MARGIN = 1,
                        FUN = calc_biodiv_single_group, n_rare = rarefy_groups)
   row.names(group_stats) <- paste(row.names(group_stats), "_obs", sep = "")
  
   # bootstrap replicates
   boot_samples <- by(mob_in$comm, INDICES = group_id, FUN = repl_boot_samples,
                      n_rare = rarefy_groups)
   
   boot_CI <- sapply(boot_samples, get_mean_CI)
   group_dat <- as.data.frame(t(rbind(group_stats, boot_CI)))
   group_dat <- group_dat[,c("N_obs","N_mean","N_low","N_up",
                             "S_obs","S_mean","S_low","S_up",
                             "S_rare1_obs","S_rare1_mean","S_rare1_low","S_rare1_up",
                             "S_rare2_obs","S_rare2_mean","S_rare2_low","S_rare2_up",
                             "S_rare3_obs","S_rare3_mean","S_rare3_low","S_rare3_up",
                             "S_asymp_obs","S_asymp_mean","S_asymp_low","S_asymp_up",
                             "PIE_obs","PIE_mean","PIE_low","PIE_up",
                             "ENS_PIE_obs","ENS_PIE_mean","ENS_PIE_low","ENS_PIE_up")]
   group <- factor(levels(group_id))
   group_dat <- cbind(group, group_dat)
   rownames(group_dat) <- NULL
   
   #beta diversities
   beta_S       = group_dat$S_obs[group_id] / S_sample
   beta_ENS_PIE = group_dat$ENS_PIE_obs[group_id] / ENS_PIE_sample
   beta_S[!is.finite(beta_S)] <- NA
   beta_ENS_PIE[!is.finite(beta_ENS_PIE)] <- NA
   
   # permutation test for differences among samples
   F_obs <- data.frame(S             = anova(lm(S_sample ~ group_id))$F[1],
                       N             = anova(lm(N_sample ~ group_id))$F[1],
                       S_rare1       = anova(lm(S_rare$S_rare1 ~ group_id))$F[1],
                       S_rare2       = anova(lm(S_rare$S_rare2 ~ group_id))$F[1],
                       S_rare3       = anova(lm(S_rare$S_rare3 ~ group_id))$F[1],
                       S_asymp       = anova(lm(S_asymp_sample ~ group_id))$F[1],
                       PIE           = anova(lm(PIE_sample ~ group_id))$F[1],
                       ENS_PIE       = anova(lm(ENS_PIE_sample ~ group_id))$F[1],
                       beta_S        = anova(lm(beta_S ~ group_id))$F[1],
                       beta_ENS_PIE  = anova(lm(beta_S ~ group_id))$F[1]
                      )
   
   F_rand <- data.frame(S             = numeric(nperm),
                        N             = numeric(nperm),
                        S_rare1       = numeric(nperm),
                        S_rare2       = numeric(nperm),
                        S_rare3       = numeric(nperm),
                        S_asymp       = numeric(nperm),
                        PIE           = numeric(nperm),
                        ENS_PIE       = numeric(nperm),
                        beta_S        = numeric(nperm),
                        beta_ENS_PIE  = numeric(nperm)
                        )
   
   for (i in 1:nperm){
      group_id_rand     = sample(group_id)
      F_rand$S[i]       = anova(lm(S_sample ~ group_id_rand))$F[1]
      F_rand$N[i]       = anova(lm(N_sample ~ group_id_rand))$F[1]
      F_rand$S_rare1[i] = anova(lm(S_rare$S_rare1 ~ group_id_rand))$F[1]
      F_rand$S_rare2[i] = anova(lm(S_rare$S_rare2 ~ group_id_rand))$F[1]
      F_rand$S_rare3[i] = anova(lm(S_rare$S_rare3 ~ group_id_rand))$F[1]
      F_rand$S_asymp[i] = anova(lm(S_asymp_sample ~ group_id_rand))$F[1]
      F_rand$PIE[i]     = anova(lm(PIE_sample ~ group_id_rand))$F[1]
      F_rand$ENS_PIE[i] = anova(lm(ENS_PIE_sample ~ group_id_rand))$F[1]
      F_rand$beta_S[i]  = anova(lm(beta_S ~ group_id_rand))$F[1]
      F_rand$beta_ENS_PIE[i]  = anova(lm(beta_ENS_PIE ~ group_id_rand))$F[1]
   }
   
   p_S       = sum(F_obs$S       <= F_rand$S) / nperm
   p_N       = sum(F_obs$N       <= F_rand$N) / nperm
   p_S_rare1 = sum(F_obs$S_rare1 <= F_rand$S_rare1) / nperm
   p_S_rare2 = sum(F_obs$S_rare2 <= F_rand$S_rare2) / nperm
   p_S_rare3 = sum(F_obs$S_rare3 <= F_rand$S_rare3) / nperm
   p_S_asymp = sum(F_obs$S_asymp <= F_rand$S_asymp) / nperm
   p_PIE     = sum(F_obs$PIE     <= F_rand$PIE) / nperm
   p_ENS_PIE = sum(F_obs$ENS_PIE <= F_rand$ENS_PIE) / nperm
   p_beta_S  = sum(F_obs$beta_S <= F_rand$beta_S) / nperm
   p_beta_ENS_PIE  = sum(F_obs$beta_ENS_PIE <= F_rand$beta_ENS_PIE) / nperm
   
   stats_pvalues = data.frame(S       = p_S,
                              N       = p_N,
                              S_rare1 = p_S_rare1,
                              S_rare2 = p_S_rare2,
                              S_rare3 = p_S_rare3,
                              S_asymp = p_S_asymp,
                              PIE     = p_PIE,
                              ENS_PIE = p_ENS_PIE,
                              beta_S  = p_beta_S,
                              beta_ENS_PIE  = p_beta_ENS_PIE
                              )

   stats_samples = data.frame(group   = group_id,
                              N       = N_sample,
                              S       = S_sample,
                              S_rare1 = S_rare$S_rare1,
                              S_rare2 = S_rare$S_rare2,
                              S_rare3 = S_rare$S_rare3,
                              S_asymp = S_asymp_sample,
                              PIE     = PIE_sample,
                              ENS_PIE = ENS_PIE_sample,
                              beta_S  = beta_S,
                              beta_ENS_PIE  = beta_ENS_PIE
                              )
   
   stats_groups = group_dat
   
   out = list(n_rarefy_samples = rarefy_samples,
              n_rarefy_groups = rarefy_groups,
              pvalues = stats_pvalues,
              samples = stats_samples,
              groups  = stats_groups)
   class(out) = 'mob_stats'
   return(out)
   
}

#' Plot the scale agnostic statistics for a MoB analysis
#' 
#' Plots a \code{mob_stats} object which is produced by the 
#' function \code{get_mob_stats}. The p value for each statistic
#' is displayed in the plot title if applicable.
#' 
#' The user may specify which results to plot or simply to plot 
#' all the results. Note however that the rarefied species richness
#' results are not plotted by all
#' @param mob_stats a \code{mob_stats} object that has the samples and 
#' treatment level statistics
#' @param display eith a single character value of a vector of character
#' values. If set to \code{all} then plots of PIE, beta PIE, S, N, and asymptic
#' S will be produced, but not a plot of rarefied S. 
#' @param multipanel a boolean if TRUE then a pretty multipanel plot is produced 
#' setting this argument to TRUE only makes sense when displaying all the 
#' results (i.e., \code{display = 'all'})
#' @param ... additional graphical arguments to be supplied to the boxplot and
#' points functions
#' @author Felix May and Dan McGlinn 
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, 'group', nperm=100)
#' plot(inv_stats, multipanel=TRUE)
#' # with colors
#' plot(inv_stats, multipanel=TRUE, col=c('red', 'blue'))
#' # only display PIE results
#' par(mfrow=c(1,3))
#' plot(inv_stats, 'PIE', col=c('red', 'blue'))
plot.mob_stats = function(mob_stats, 
                          display=c('all','N', 'S', 'S rare', 'S asymp',
                                    'PIE', 'ENS PIE'),
                          multipanel=FALSE, ...) {
      if (multipanel) {
          if ('all' %in% display) {
              display = c('PIE', 'N', 'S', 'S asymp')
              op = par(mfcol = c(3,5), las = 1, font.main = 1,
                       mar = c(2, 2, 3, 2) + 0.1)
          }
          else
              warning('The multipanel plots settings only make sense when display is set to "all"')
      } else if ('all' %in% display)
          display = c('PIE', 'N', 'S', 'S asymp')
      # PIE
      if ('PIE' %in% display) {
          if (multipanel)
              par(fig = c(0.2, 0.4, 0.66, 1))
          else
             par(mfrow = c(1,2))
          panel_title = paste("PIE (p = ", mob_stats$pvalues$PIE,")\nplot scale",
                              sep = "")
          boxplot(PIE ~ group, data=mob_stats$samples,
                  main = panel_title, ...)
          
          if (multipanel)
              par(fig = c(0.6, 0.8, 0.66, 1), new = T)
          y_limits <- with(mob_stats$groups, range(c(PIE_mean, PIE_low, PIE_obs, PIE_up)))
          boxplot(PIE_obs ~ group, data=mob_stats$groups,
                  main = "PIE\ntreatment scale", boxwex = 0, ylim = y_limits)
          points(PIE_obs ~ group, data = mob_stats$groups, pch = 19)
          plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$PIE_mean,
                 li = mob_stats$groups$PIE_low, ui = mob_stats$groups$PIE_up,
                 add = T, pch = 1, ...)
      }
      
      # ENS PIE
      if ('ENS PIE' %in% display) {
         if (multipanel)
            par(fig = c(0.2, 0.4, 0.66, 1))
         else
            par(mfrow = c(1,3))
         panel_title = paste("ENS PIE (p = ", mob_stats$pvalues$ENS_PIE,")\nplot scale",
                             sep = "")
         boxplot(ENS_PIE ~ group, data=mob_stats$samples,
                 main = panel_title, ...)
         
         # beta ENS PIE
         panel_title = paste("beta ENS PIE (p = ", mob_stats$pvalues$beta_ENS_PIE,")\n between scales",
                             sep = "")
         boxplot(beta_ENS_PIE ~ group, data=mob_stats$samples,
                 main = panel_title, ...)
         
         if (multipanel)
            par(fig = c(0.6, 0.8, 0.66, 1), new = T)
         y_limits <- with(mob_stats$groups, range(c(ENS_PIE_mean, ENS_PIE_low, ENS_PIE_obs, ENS_PIE_up)))
         boxplot(ENS_PIE_obs ~ group, data=mob_stats$groups,
                 main = "ENS_PIE\ntreatment scale", boxwex = 0, ylim = y_limits)
         points(ENS_PIE_obs ~ group, data = mob_stats$groups, pch = 19)
         plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$ENS_PIE_mean,
                li = mob_stats$groups$ENS_PIE_low, ui = mob_stats$groups$ENS_PIE_up,
                add = T, pch = 1, ...)
      }
      
      if ('S' %in% display) {
          if (multipanel)
              par(fig = c(0, 0.2, 0.33, 0.66), new = T)
          else
              par(mfrow = c(1,3))

          # S observed samples
          panel_title = paste("S (p = ", mob_stats$pvalues$S,")\n plot scale",
                              sep = "")
          boxplot(S ~ group, data=mob_stats$samples,
                  main = panel_title, ...)
          
          # S beta
          panel_title = paste("beta S (p = ", mob_stats$pvalues$beta_S,")\n between scales",
                              sep = "")
          boxplot(beta_S ~ group, data=mob_stats$samples,
                  main = panel_title, ...)
          
          # S groups
          if (multipanel)
              par(fig = c(0.8, 1.0, 0.33, 0.66), new = T)
          y_limits <- with(mob_stats$groups, range(c(S_mean, S_low, S_obs, S_up)))
          boxplot(S_obs ~ levels(group), data=mob_stats$groups,
                  main = "S\ntreatment scale", boxwex = 0, ylim = y_limits)
          points(S_obs ~ group, data = mob_stats$groups, pch = 19)
          plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$S_mean,
                 li = mob_stats$groups$S_low, ui = mob_stats$groups$S_up,
                 add = T, pch = 1, ...)

      }
      if ('S rare' %in% display) {
          if (multipanel)
              par(fig = c(0, 0.2, 0.33, 0.66), new = T)
          else
              par(mfrow = c(1,2))

          # S rarefied samples
          panel_title = paste("S rarefied (p = ", mob_stats$pvalues$S_rare,
                              ")\n plot scale", sep = "")
          boxplot(S_rare ~ group, data=mob_stats$samples,
                  main = panel_title, ...)
          
          # S groups
          if (multipanel)
              par(fig = c(0.8, 1.0, 0.33, 0.66), new = T)
          y_limits <- with(mob_stats$groups, range(c(S_rare_mean, S_rare_low, S_rare_obs, S_rare_up)))
          boxplot(S_rare_mean ~ levels(group), data=mob_stats$groups,
                  main = "S rarefied\ntreatment scale", boxwex = 0, ylim = y_limits)
          points(S_rare_obs ~ group, data = mob_stats$groups, pch = 19)
          plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$S_rare_mean,
                 li = mob_stats$groups$S_rare_low, ui = mob_stats$groups$S_rare_up,
                 add = T, pch = 1, ...)

      }
      if ('N' %in% display) {
          if (multipanel)
              par(fig = c(0.2, 0.4, 0.33, 0.66), new = T)
          else
              par(mfrow = c(1,2))
          # N samples
          panel_title = paste("N (p = ", mob_stats$pvalues$N,")\n plot scale",
                              sep = "")
          boxplot(N ~ group, data=mob_stats$samples, 
                  main = panel_title, ...)
      
          # N groups
          if (multipanel)
              par(fig = c(0.6, 0.8, 0.33, 0.66), new = T)
          y_limits <- with(mob_stats$groups, range(c(N_mean, N_low, N_obs, N_up)))
          boxplot(N_mean ~ levels(group), data=mob_stats$groups,
                  main = "N\ntreatment scale", boxwex = 0, ylim = y_limits)
          points(N_obs ~ group, data = mob_stats$groups, pch = 19)
          plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$N_mean,
                 li = mob_stats$groups$N_low, ui = mob_stats$groups$N_up,
                 add = T, pch = 1, ...)
      }
      if ('S asymp' %in% display) {
          if (multipanel)
              par(fig = c(0.2, 0.4, 0, 0.33), new = T)
          else
              par(mfrow = c(1,2))
           # S asymptotic
          if (all(is.na(mob_stats$samples$S_asymp))){
              warning("Cannot plot asymptotic richness for the samples")
          } else {
              panel_title = paste("S asympotic (p = ", mob_stats$pvalues$S_asymp,
                                  ")\n plot scale", sep = "")
              boxplot(S_asymp ~ group, data=mob_stats$samples,
                      main = panel_title, ...)
          }
          
          if (multipanel)
              par(fig = c(0.6, 0.8, 0, 0.33), new = T)
          if (all(is.na(mob_stats$groups$S_asymp_obs))){
              warning("Cannot plot asymptotic richness for the groups")
          } else {
             y_limits <- with(mob_stats$groups, range(c(S_asymp_mean, S_asymp_low, S_asymp_obs, S_asymp_up))) 
             boxplot(S_asymp_obs ~ levels(group), data=mob_stats$groups,
                      main = "S asympotic\ntreatment scale", boxwex = 0, ylim = y_limits)
             points(S_asymp_obs ~ group, data = mob_stats$groups, pch = 19)
             plotCI(x = as.numeric(mob_stats$groups$group), y = mob_stats$groups$S_asymp_mean,
                    li = mob_stats$groups$S_asymp_low, ui = mob_stats$groups$S_asymp_up,
                    add = T, pch = 1, ...)
          }
          
      }
      if (multipanel)
          par(op)
}


#' Additional plotting function for betaPIE statistics 
#' @inheritParams plot.mob_stats
#' @author Felix May and Dan McGlinn
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' inv_stats = get_mob_stats(inv_mob_in, 'group', nperm=100)
#' plot_betaPIE(inv_stats)
plot_betaPIE = function(mob_stats)
{
   ngroups = nrow(mob_stats$groups)
   
   op = par(mfrow = c(1,2), las = 1, font.main = 1)
   pie_range = range(c(mob_stats$samples$PIE, mob_stats$groups$PIE),
                     na.rm=T)
   boxplot(PIE ~ group, data = mob_stats$sample, ylim = pie_range, ylab = "PIE",
           notch = T, main = "Sample and group PIE")
   points((1:ngroups)+0.1, mob_stats$groups$PIE, pch = 23, bg = "grey", cex = 1.5)
   
   boxplot(betaPIE ~ group, data = mob_stats$samples,
           ylab = " betaPIE", notch = T, main = "Group PIE - sample PIE")
   par(op)
}
