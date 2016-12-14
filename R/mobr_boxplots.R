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
mob_stats = function(mob_in, group_var) {
   group_id  = factor(mob_in$env[, group_var])
   
   # Sample based statistics
   N_sample = rowSums(mob_in$comm)      # individuals in each sample
   S_sample = rowSums(mob_in$comm > 0) # species in each sample
   
   # rarefied species richness
   Nmin_sample = min(N_sample)              # lowest number of individuals in sample
   S_rare_sample = apply(mob_in$comm, MARGIN = 1,
                          FUN=function(x){rarefaction(as.numeric(x), method = "indiv", effort = Nmin_sample)})
   
   PIE_sample = calc_PIE(mob_in$comm)
   
   # extrapolated species richness - iChao1 estimator for each sample
   S_ext_sample = rep(NA, nrow(mob_in$comm))    
   #for (i in 1:nrow(mob_in$comm)) {
   #    if (sum(mob_in$comm[i, ]) > 4) {
   #        S = SpadeR::ChaoSpecies(mob_in$comm[i,], datatype = "abundance")
   #        S_ext_sample[i] = S$Species_table[5,1]
   #    } 
   #}
   # Group based statistics
   
   # abundance distribution pooled in group
   abund_group = by(mob_in$comm, group_id, function(x) colSums(x) )
   
   N_group     = sapply(abund_group, sum)
   S_group     = sapply(abund_group, function(x) sum(x > 0))
   
   Nmin_group = min(N_group) 
   S_rare_group = sapply(abund_group, function(x)
                         rarefaction(as.numeric(x), method = "indiv",
                                     effort = Nmin_group))
   
   PIE_group   = sapply(abund_group, calc_PIE)
   
   # calculate Chao estimator of species richness using JADE
   S_ext_group = NA
   #S_ext_group = sapply(abund_group, function(x) 
   #                      ChaoSpecies(x, datatype = "abundance")$Species_table[5,])
   betaPIE_sample = PIE_group[group_id] - PIE_sample
   
   stats_samples = data.frame(group   = group_id,
                               N       = N_sample,
                               S       = S_sample,
                               S_rare  = S_rare_sample,
#                               S_ext   = S_ext_sample,
                               PIE     = PIE_sample,
                               betaPIE = betaPIE_sample
                               )
   
   stats_groups = data.frame(group  = factor(levels(group_id), levels = levels(group_id), ordered = is.ordered(group_id)),
                              N      = N_group,
                              S      = S_group,
                              S_rare = S_rare_group,
#                              S_ext_mean  = S_ext_group[1,],
#                              S_ext_CIlow = S_ext_group[3,],
#                              S_ext_CIup  = S_ext_group[4,],
                              PIE    = PIE_group)
   
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
