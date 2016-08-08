# install Spade R package
# see http://chao.stat.nthu.edu.tw/wordpress/software_download/
#library(devtools)
#install_github('AnneChao/SpadeR')


# ------------------------------------------------------------------------------
# Calculate sample based and group based biodiversity statistics
mob_stats <- function(comm, group_var)
{
   require(vegan)
   require(SpadeR)
   
   group_id  = comm$env[, group_var]
   
   # Sample based statistics
   N_sample <- rowSums(comm$comm)      # individuals in each sample
   S_sample <- rowSums(comm$comm > 0) # species in each sample
   
   # rarefied species richness
   Nmin_sample <- min(N_sample)              # lowest number of individuals in sample
   S_rare_sample <- apply(comm$comm, MARGIN = 1,
                          FUN=function(x){rarefaction(as.numeric(x), method = "indiv", effort = Nmin_sample)})
   
   PIE_sample <- diversity(comm$comm, index = "simpson")
   
   # extrapolated species richness - iChao1 estimator for each sample
   S_ext_sample <- apply(comm$comm, MARGIN = 1,
                         FUN=function(x){ChaoSpecies(x, datatype = "abundance")$Species.Table[5,1]})
   
   # Group based statistics
   
   # abundance distribution pooled in group
   abund_group <- by(comm$comm, group_id, FUN = function(x){colSums(x)})
   
   N_group     <- sapply(abund_group, sum)
   S_group     <- sapply(abund_group, function(x){sum(x > 0)})
   
   Nmin_group <- min(N_group) 
   S_rare_group <- sapply(abund_group,
                          FUN=function(x){rarefaction(as.numeric(x), method = "indiv", effort = Nmin_group)})
   
   PIE_group   <- sapply(abund_group, diversity, index = "simpson")
   
   # calculate Chao estimator of species richness using JADE
   S_ext_group <- sapply(abund_group, 
                         FUN = function(x){ChaoSpecies(x, datatype = "abundance")$Species.Table[5, ]})
   
   betaPIE_sample <- PIE_group[group_id] - PIE_sample
   
   stats_samples <- data.frame(group   = group_id,
                               N       = N_sample,
                               S       = S_sample,
                               S_rare  = S_rare_sample,
                               S_ext   = S_ext_sample,
                               PIE     = PIE_sample,
                               betaPIE = betaPIE_sample
                               )
   
   stats_groups <- data.frame(group  = factor(levels(group_id), levels = levels(group_id), ordered = is.ordered(group_id)),
                              N      = N_group,
                              S      = S_group,
                              S_rare = S_rare_group,
                              S_ext_mean  = S_ext_group[1,],
                              S_ext_CIlow = S_ext_group[3,],
                              S_ext_CIup  = S_ext_group[4,],
                              PIE    = PIE_group)
   
   return(list(samples = stats_samples,
               groups  = stats_groups))
   
}

# ------------------------------------------------------------------------------
# Boxplots for sample based comparisons
plot_samples <- function(sample_stats, tukey = F)
{
   # create local diversity boxplots ------------------------------
   op <- par(mfcol = c(2,4), las = 1, font.main = 1)
   
   # raw species richness
   par(fig = c(0.32, 0.68, 0.5, 1))
   test <- kruskal.test(S ~ group, data = sample_stats)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(S ~ group, data = sample_stats, ylab = "No. of species", notch = T, main = title)
   
   # number of individuals
   par(fig = c(0, 0.25, 0, 0.5), new = T)
   test <- kruskal.test(N ~ group, data = sample_stats)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(N ~ group, data = sample_stats, ylab = "No. of individuals", notch = T, main = title) 
   
   # number of species rarefied
   par(fig = c(0.25, 0.5, 0, 0.5), new = T)
   test <- kruskal.test(S_rare ~ group, data = sample_stats)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   ylabel <- paste("No. of species (n = ",min(sample_stats$N),")",sep = "")
   boxplot(S_rare ~ group, data = sample_stats, ylab = ylabel , notch = T, main = title)
   
   # number of species extrapolated
   par(fig = c(0.5, 0.75, 0, 0.5), new = T)
   test <- kruskal.test(S_ext ~ group, data = sample_stats)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   ylabel <- "No. of species extrapolated"
   boxplot(S_ext ~ group, data = sample_stats, ylab = ylabel , notch = T, main = title)
   
   # PIE
   par(fig = c(0.75, 1, 0, 0.5), new = T)
   test <- kruskal.test(PIE ~ group, data = sample_stats)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(PIE ~ group, data = sample_stats, ylab = "PIE", notch = T, main = title)
   
   par(op)
   
   if (tukey == T){
      # create Tukey plot ------------------------------
      op <- par(mfcol = c(2,4), las = 0, font.main = 1)
      
      par(fig = c(0.375, 0.625, 0.5, 1))
      plot(TukeyHSD(aov(S ~ group, data = sample_stats)))
      legend("topright",c("No. of species"), cex = 1, bty = "n")
      
      par(fig = c(0, 0.25, 0, 0.5), new = T)
      plot(TukeyHSD(aov(N ~ group, data = sample_stats)))
      legend("topright",c("No. of individuals"), cex = 1, bty = "n")
      
      par(fig = c(0.25, 0.5, 0, 0.5), new = T)
      plot(TukeyHSD(aov(S_rare ~ group, data = sample_stats)))
      legend("topright", paste("Species (n = ", min(sample_stats$N),")",sep = ""),
             cex = 1, bty = "n")
      
      par(fig = c(0.5, 0.75, 0, 0.5), new = T)
      plot(TukeyHSD(aov(S_ext ~ group, data = sample_stats)))
      legend("topright",c("Species extrapolated"), cex = 1, bty = "n")
      
      par(fig = c(0.75, 1, 0, 0.5), new = T)
      plot(TukeyHSD(aov(PIE ~ group,  data = sample_stats)))
      legend("topright",c("PIE"), cex = 1, bty = "n")
      
      par(op)
   }
   
}


# ------------------------------------------------------------------------------
# Plot for group-based comparisons
plot_groups <- function(group_stats)
{
   require(plotrix)
   op <- par(mfrow = c(1,3), las = 1, font.main = 1)
   
   minS <- min(group_stats[, 3:7])
   maxS <- max(group_stats[, 3:7])
   
   ngroups <- nrow(group_stats)
   
   plot(S ~ group, data = group_stats, boxwex = 0, ylim = c(0.95*minS, 1.05*maxS),
        ylab = "", main = "Species")
   points(S ~ group, data = group_stats, pch = 19, cex = 1.5)
   points((1:ngroups)-0.2, group_stats$S_rare, pch = 3, cex = 1.5 )
   plotCI((1:ngroups)+0.2, group_stats$S_ext_mean, pch = 1, cex = 1.5,
          li = group_stats$S_ext_CIlow, ui = group_stats$S_ext_CIup, add = T)
   legend("top",c("Observed S", "Rarefied S", "Extrapolated S"), pch = c(19,3,1), cex = 1)
   
   plot(N ~ group, data = group_stats, boxwex = 0, ylab = "", main = "Individuals")
   points(N ~ group, data =group_stats, pch = 19)
   
   plot(PIE ~ group, data = group_stats, boxwex = 0, ylab = "", main = "PIE")
   points(PIE ~ group, data =group_stats, pch = 19)
   
   par(op)
}

# ------------------------------------------------------------------------------
# Plot for betaPIE
plot_betaPIE <- function(mob_stats)
{
   ngroups <- nrow(mob_stats$groups)
   
   op <- par(mfrow = c(1,2), las = 1, font.main = 1)
   pie_range <- range(c(mob_stats$samples$PIE, mob_stats$groups$PIE))
   boxplot(PIE ~ group, data = mob_stats$sample, ylim = pie_range, ylab = "PIE",
           notch = T, main = "Sample and group PIE")
   points((1:ngroups)+0.1, mob_stats$groups$PIE, pch = 23, bg = "grey", cex = 1.5)
   
   boxplot(betaPIE ~ group, data = mob_stats$samples,
           ylab = " betaPIE", notch = T, main = "Group PIE - sample PIE")
   par(op)
}
   
# ------------------------------------------------------------------------------
# # Example 
# library(vegan)
# data(mite)
# data(mite.env)
# data(mite.xy)
# mite_comm = make_comm_obj(mite, cbind(mite.env, mite.xy))
# 
# mite_stats <- mob_stats(mite_comm, "Topo")
# plot_samples(mite_stats$samples, tukey = T)
# plot_groups(mite_stats$groups)
# plot_betaPIE(mite_stats)
