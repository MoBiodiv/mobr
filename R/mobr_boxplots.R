#' Calculate probability of interspecific encounter (PIE)
#' 
#' PIE is also known as Simpson's evenness index and this function is 
#' a reduced form of the function vegan::diversity(). Jari Oksanen and Bob O'Hara
#' are the original authors of the function vegan::diversity()
#' 
#' @inheritParams rarefaction
#' @author Dan McGlinn
#' @keywords internal
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
        H = 1 - H
    if (any(is.na(total))) 
        is.na(H) = is.na(total)
    return(H)
}


#' Calculate sample based and group based biodiversity statistics
#' @inheritParams get_delta_stats
#' @param group_var a string that specifies which field in mob_in$env the data
#'   should be grouped by
#' @importFrom SpadeR ChaoSpecies
#' @author Felix May and Dan McGlinn
#' @export
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' mob_stats(inv_mob_in)
mob_stats = function(mob_in, group_var, plot = T, n_min = 10, nperm = 100) {
   
   group_id  = factor(mob_in$env[, group_var])
   
   # Sample-based statistics
   N_sample = rowSums(mob_in$comm)      # individuals in each sample
   S_sample = rowSums(mob_in$comm > 0)  # species in each sample
   
   # rarefied richness
   N_min_sample = max(n_min, min(N_sample))
   plots_low_n = N_sample < n_min
   
   if (sum(plots_low_n) > 0){
      warning(paste("There are",sum(plots_low_n),"plots with less then", n_min,"individuals.
These are removed for the calculation of rarefied richness."))
   }

   S_rare_sample = rep(NA, nrow(mob_in$comm))
   S_rare_sample[!plots_low_n] = apply(mob_in$comm[!plots_low_n,], MARGIN = 1,
                                       FUN = rarefaction, method = "indiv",
                                       effort = N_min_sample)
   
   # Probability of Interspecific Encounter
   plots_n0 = N_sample == 0
   
   if (sum(plots_n0) > 0){
      warning(paste("There are",sum(plots_n0), "plots without any individuals.
These are removed for the calculation of PIE."))
   }
   
   PIE_sample = diversity(mob_in$comm, index = "simpson")
   ENS_PIE_sample = diversity(mob_in$comm, index = "invsimpson")
   PIE_sample[plots_n0] = NA
   ENS_PIE_sample[plots_n0] = NA
   
   # bias corrected Chao estimator
   S_asymp_sample = estimateR(mob_in$comm)["S.chao1",]
 
   # ---------------------------------------------------------------------------
   # abundance distribution pooled in group
   abund_group = aggregate(mob_in$comm, by = list(group_id), FUN = "sum")[ ,-1]
   
   # Group-based statistics
   N_group = rowSums(abund_group)      
   S_group = rowSums(abund_group > 0) 
   
   N_min_group = max(n_min, min(N_group))
   groups_low_n = N_group < n_min
   
   if (sum(groups_low_n) > 0){
      warning(paste("There are",sum(plots_low_n),"groups with less then", n_min,"individuals.
These are removed for the calculation of rarefied richness."))
   }
   
   S_rare_group = rep(NA, length(N_group))
   S_rare_group[!groups_low_n] = apply(abund_group[!groups_low_n], MARGIN = 1,
                                       FUN = rarefaction, method = "indiv",
                                       effort = N_min_group)
   
   # Probability of Interspecific Encounter
   groups_n0 = N_group == 0
   
   if (sum(groups_n0) > 0){
      warning(paste("There are",sum(groups_n0), "groups without any individuals.
                    These are removed for the calculation of PIE."))
   }
   
   PIE_group = diversity(abund_group, index = "simpson")
   ENS_PIE_group = diversity(abund_group, index = "invsimpson")
   PIE_group[groups_n0] = NA
   ENS_PIE_group[groups_n0] = NA
   
   # bias corrected Chao estimator
   S_asymp_group = estimateR(abund_group)["S.chao1",]
   
   #beta PIE
   betaPIE_sample = PIE_group[group_id] - PIE_sample
   
   # permutation test for differences among samples
   F_obs <- data.frame(S       = anova(lm(S_sample ~ group_id))$F[1],
                       N       = anova(lm(N_sample ~ group_id))$F[1],
                       S_rare  = anova(lm(S_rare_sample ~ group_id))$F[1],
                       PIE     = anova(lm(PIE_sample ~ group_id))$F[1],
                       ENS_PIE = anova(lm(ENS_PIE_sample ~ group_id))$F[1],
                       S_asymp = anova(lm(S_asymp_sample ~ group_id))$F[1],
                       betaPIE = anova(lm(betaPIE_sample ~ group_id))$F[1]
                      )
   
   F_rand <- data.frame(S       = numeric(nperm),
                        N       = numeric(nperm),
                        S_rare  = numeric(nperm),
                        PIE     = numeric(nperm),
                        ENS_PIE = numeric(nperm),
                        S_asymp = numeric(nperm),
                        betaPIE = numeric(nperm)
                        )
   
   for (i in 1:nperm){
      group_id_rand     = sample(group_id)
      F_rand$S[i]       = anova(lm(S_sample ~ group_id_rand))$F[1]
      F_rand$N[i]       = anova(lm(N_sample ~ group_id_rand))$F[1]
      F_rand$S_rare[i]  = anova(lm(S_rare_sample ~ group_id_rand))$F[1]
      F_rand$PIE[i]     = anova(lm(PIE_sample ~ group_id_rand))$F[1]
      F_rand$ENS_PIE[i] = anova(lm(ENS_PIE_sample ~ group_id_rand))$F[1]
      F_rand$S_asymp[i] = anova(lm(S_asymp_sample ~ group_id_rand))$F[1]
      F_rand$betaPIE[i] = anova(lm(betaPIE_sample ~ group_id_rand))$F[1]
   }
   
   p_S       = sum(F_obs$S       < F_rand$S)/nperm
   p_N       = sum(F_obs$N       < F_rand$N)/nperm
   p_S_rare  = sum(F_obs$S_rare  < F_rand$S_rare)/nperm
   p_PIE     = sum(F_obs$PIE     < F_rand$PIE)/nperm
   p_ENS_PIE = sum(F_obs$ENS_PIE < F_rand$ENS_PIE)/nperm
   p_S_asymp = sum(F_obs$S_asymp < F_rand$S_asymp)/nperm
   p_betaPIE = sum(F_obs$betaPIE < F_rand$betaPIE)/nperm
   
   if (plot == T){
      #windows(10,6)
      op = par(mfcol = c(3,5), las = 1, font.main = 1)
      
      # PIE
      par(fig = c(0.2, 0.4, 0.66, 1))
      panel_title = paste("PIE (p = ", p_PIE,")", sep = "")
      boxplot(PIE_sample ~ group_id, main = panel_title)
      
      par(fig = c(0.4, 0.6, 0.66, 1), new = T)
      panel_title = paste("beta PIE (p = ", p_betaPIE,")", sep = "")
      boxplot(betaPIE_sample ~ group_id, main = panel_title)
      
      par(fig = c(0.6, 0.8, 0.66, 1), new = T)
      boxplot(PIE_group ~ levels(group_id), main = "PIE", boxwex = 0)
      points(PIE_group, pch =19)
      
      # S observed samples
      par(fig = c(0, 0.2, 0.33, 0.66), new = T)
      panel_title = paste("S (p = ", p_S,")", sep = "")
      boxplot(S_sample ~ group_id, main = panel_title)
      
      # N samples
      par(fig = c(0.2, 0.4, 0.33, 0.66), new = T)
      panel_title = paste("N (p = ", p_N,")", sep = "")
      boxplot(N_sample ~ group_id, main = panel_title)
      
      # N groups
      par(fig = c(0.6, 0.8, 0.33, 0.66), new = T)
      boxplot(N_group ~ levels(group_id), main = "N", boxwex = 0)
      points(N_group, pch = 19)
      
      # S groups
      par(fig = c(0.8, 1.0, 0.33, 0.66), new = T)
      boxplot(S_group ~ levels(group_id), main = "S obs", boxwex = 0)
      points(S_group, pch = 19)
      
      # S asymptotic
      par(fig = c(0.2, 0.4, 0, 0.33), new = T)
      panel_title = paste("S asymp (p = ", p_S_asymp,")", sep = "")
      boxplot(S_asymp_sample ~ group_id, main = panel_title)
      
      par(fig = c(0.6, 0.8, 0, 0.33), new = T)
      boxplot(S_asymp_group ~ levels(group_id), main = "asymptotic S", boxwex = 0)
      points(S_asymp_group, pch =19)
   }
   
   stats_samples = data.frame(group   = group_id,
                              N       = N_sample,
                              S       = S_sample,
                              S_rare  = S_rare_sample,
                              S_asymp = S_asymp_sample,
                              PIE     = PIE_sample,
                              ENS_PIE = ENS_PIE_sample,
                              betaPIE = betaPIE_sample
                              )
   
   stats_groups = data.frame(group  = factor(levels(group_id), ordered = is.ordered(group_id)),
                             N      = N_group,
                             S      = S_group,
                             S_rare = S_rare_group,
                             S_asymp = S_asymp_group,
                             PIE    = PIE_group,
                             ENS_PIE = ENS_PIE_group
                             )
   
   return(list(samples = stats_samples,
               groups  = stats_groups))
   
}

#' Boxplots for sample based comparisons
#' @author Felix May and Dan McGlinn
#' @export
plot_samples = function(sample_stats, tukey = F, col=NULL)
{
   # create local diversity boxplots ------------------------------
   op = par(mfcol = c(2,3), las = 1, font.main = 1)
   
   #if (is.null(col))
   #  col = 
   
   # raw species richness
   par(fig = c(0.32, 0.68, 0.5, 1))
   test = kruskal.test(S ~ group, data = sample_stats)
   title = paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(S ~ group, data = sample_stats, ylab = "No. of species", notch = T, main = title)
   
   # number of individuals
   par(fig = c(0, 0.33, 0, 0.5), new = T)
   test = kruskal.test(N ~ group, data = sample_stats)
   title = paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(N ~ group, data = sample_stats, ylab = "No. of individuals", notch = T, main = title) 
   
   # number of species rarefied
   par(fig = c(0.33, 0.66, 0, 0.5), new = T)
   test = kruskal.test(S_rare ~ group, data = sample_stats)
   title = paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   ylabel = paste("No. of species (n = ",min(sample_stats$N),")",sep = "")
   boxplot(S_rare ~ group, data = sample_stats, ylab = ylabel , notch = T, main = title)
   
   # number of species extrapolated
   #par(fig = c(0.5, 0.75, 0, 0.5), new = T)
   #test = kruskal.test(S_ext ~ group, data = sample_stats)
   #title = paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   #ylabel = "No. of species extrapolated"
   #boxplot(S_ext ~ group, data = sample_stats, ylab = ylabel , notch = T, main = title)
   
   # PIE
   par(fig = c(0.66, 1, 0, 0.5), new = T)
   test = kruskal.test(PIE ~ group, data = sample_stats)
   title = paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(PIE ~ group, data = sample_stats, ylab = "PIE", notch = T, main = title)
   
   par(op)
   
   if (tukey == T){
      # create Tukey plot ------------------------------
      op = par(mfcol = c(2,4), las = 0, font.main = 1)
      
      par(fig = c(0.375, 0.625, 0.5, 1))
      plot(TukeyHSD(aov(S ~ group, data = sample_stats)))
      legend("topright",c("No. of species"), cex = 1, bty = "n")
      
      par(fig = c(0, 0.33, 0, 0.5), new = T)
      plot(TukeyHSD(aov(N ~ group, data = sample_stats)))
      legend("topright",c("No. of individuals"), cex = 1, bty = "n")
      
      par(fig = c(0.33, 0.66, 0, 0.5), new = T)
      plot(TukeyHSD(aov(S_rare ~ group, data = sample_stats)))
      legend("topright", paste("Species (n = ", min(sample_stats$N),")",sep = ""),
             cex = 1, bty = "n")
      
      #par(fig = c(0.5, 0.75, 0, 0.5), new = T)
      #plot(TukeyHSD(aov(S_ext ~ group, data = sample_stats)))
      #legend("topright",c("Species extrapolated"), cex = 1, bty = "n")
      
      par(fig = c(0.66, 1, 0, 0.5), new = T)
      plot(TukeyHSD(aov(PIE ~ group,  data = sample_stats)))
      legend("topright",c("PIE"), cex = 1, bty = "n")
      
      par(op)
   }
   
}


#' Plot for group-based comparisons
#' @importFrom plotrix plotCI
#' @export
plot_groups = function(group_stats)
{
   op = par(mfrow = c(2,2), las = 1, font.main = 1)
   
   minS = min(group_stats[,c("S", "S_rare")])
   maxS = max(group_stats[,c("S", "S_rare")])
   
   minN = min(group_stats[,c("N")])
   maxN = max(group_stats[,c("N")])
   
   minPIE = min(group_stats[,c("PIE")])
   maxPIE = max(group_stats[,c("PIE")])
   
   ngroups = nrow(group_stats)
   
   plot(S ~ group, data = group_stats, boxwex = 0, ylim = c(0.9*minS, 1.1*maxS),
        ylab = "", main = "Observed species richness")
   points(S ~ group, data = group_stats, pch = 19, cex = 1.5)
   #points((1:ngroups)-0.2, group_stats$S_rare, pch = 3, cex = 1.5 )
#   plotrix::plotCI((1:ngroups)+0.2, group_stats$S_ext_mean, pch = 1, cex = 1.5,
#          li = group_stats$S_ext_CIlow, ui = group_stats$S_ext_CIup, add = T)
   #legend("top",c("Observed S", "Rarefied S", "Extrapolated S"), pch = c(19,3,1), cex = 1)
   
   plot(S_rare ~ group, data = group_stats, boxwex = 0, ylim = c(0.9*minS, 1.1*maxS),
        ylab = "", main = "Rarefied species richness")
   points(S_rare ~ group, data = group_stats, pch = 19, cex = 1.5)
   #points((1:ngroups)-0.2, group_stats$S_rare, pch = 3, cex = 1.5 )
   #   plotrix::plotCI((1:ngroups)+0.2, group_stats$S_ext_mean, pch = 1, cex = 1.5,
   #          li = group_stats$S_ext_CIlow, ui = group_stats$S_ext_CIup, add = T)
   #legend("top",c("Observed S", "Rarefied S", "Extrapolated S"), pch = c(19,3,1), cex = 1)
   
   plot(N ~ group, data = group_stats, boxwex = 0, ylab = "", ylim = c(0.9*minN, 1.1*maxN),
        main = "Individuals")
   points(N ~ group, data =group_stats, pch = 19, cex = 1.5)
   
   plot(PIE ~ group, data = group_stats, boxwex = 0, ylab = "", ylim = c(0.9*minPIE, 1.1*maxPIE),
        main = "PIE")
   points(PIE ~ group, data =group_stats, pch = 19, cex = 1.5)
   
   par(op)
}

#' Plot for betaPIE
#' @export
plot_betaPIE = function(mob_stats)
{
   ngroups = nrow(mob_stats$groups)
   
   op = par(mfrow = c(1,2), las = 1, font.main = 1)
   pie_range = range(c(mob_stats$samples$PIE, mob_stats$groups$PIE))
   boxplot(PIE ~ group, data = mob_stats$sample, ylim = pie_range, ylab = "PIE",
           notch = T, main = "Sample and group PIE")
   points((1:ngroups)+0.1, mob_stats$groups$PIE, pch = 23, bg = "grey", cex = 1.5)
   
   boxplot(betaPIE ~ group, data = mob_stats$samples,
           ylab = " betaPIE", notch = T, main = "Group PIE - sample PIE")
   par(op)
}
   
# # Example 
# library(vegan)
# data(mite)
# data(mite.env)
# data(mite.xy)
# mite_comm = make_comm_obj(mite, cbind(mite.env, mite.xy))
# 
# mite_stats = mob_stats(mite_comm, "Topo")
# plot_samples(mite_stats$samples, tukey = T)
# plot_groups(mite_stats$groups)
# plot_betaPIE(mite_stats)
