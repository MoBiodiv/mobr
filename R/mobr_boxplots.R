#' Boxplots for initial comparison of diversity in different groups of samples
#'
#' @param comm Community object
#' @param env_var String with the name of the grouping variable
#'
#' @examples
#' 
#' library(vegan)
#' data(mite)
#' data(mite.env)
#' data(mite.xy)
#' mite_comm = make_comm_obj(mite, cbind(mite.env, mite.xy))
#' boxplot(mite_comm, "Topo")
#' 
boxplot.comm <- function(comm, env_var){
   
   require(vegan)
   
   env_data = comm$env[, env_var]
   groups = unique(env_data)
   
   nInd_sample <- rowSums(comm$comm)      # inidividuals in each sample
   nSpec_sample <- rowSums(comm$comm > 0) # species in each sample
   
   n.min <- min(nInd_sample)  # lowest number of individuals observed
   nSpec_rare <- apply(comm$comm, MARGIN = 1, FUN=function(x){rarefaction(as.numeric(x), method = "indiv", effort = n.min)})
   
   PIE_sample <- diversity(comm$comm, index = "simpson")
   PIE_group <- by(comm$comm, env_data, FUN = function(x){diversity(colSums(x), index = "simpson")})
   
   # calculate Chao estimator of species richness using JADE
   nSpec_jade <- apply(comm$comm, MARGIN = 1, FUN = function(x){nrow(SpecDist(x))})
   
   # create diversity boxplots ------------------------------
   windows(width = 10, height = 8, title = "Boxplots")
   op <- par(mfcol = c(2,4), las = 1, cex.main = 1.5, cex.lab = 1.5,
             cex.axis = 1.1, font.main = 1)
   
   # raw species richness
   par(fig = c(0.32, 0.68, 0.5, 1), new = T)
   test <- kruskal.test(nSpec_sample ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(nSpec_sample ~ env_data, ylab = "No. of species", notch = T, main = title)
   
   # number of individuals
   par(fig = c(0, 0.25, 0, 0.5), new = T)
   test <- kruskal.test(nInd_sample ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(nInd_sample ~ env_data, ylab = "No. of individuals", notch = T, main = title) 
   
   # number of species rarefied
   par(fig = c(0.25, 0.5, 0, 0.5), new = T)
   test <- kruskal.test(nSpec_rare ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   ylabel <- paste("No. of species (n = ",n.min,")",sep = "")
   boxplot(nSpec_rare ~ env_data, ylab = ylabel , notch = T, main = title)
   
   # number of species extrapolated
   par(fig = c(0.5, 0.75, 0, 0.5), new = T)
   test <- kruskal.test(nSpec_jade ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   ylabel <- "No. of species extrapolated"
   boxplot(nSpec_jade ~ env_data, ylab = ylabel , notch = T, main = title)
   
   # PIE
   par(fig = c(0.75, 1, 0, 0.5), new = T)
   test <- kruskal.test(PIE_sample ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 3), sep = "")
   boxplot(PIE_sample ~ env_data, ylab = "PIE", notch = T, main = title)
   
   
   # create Tukey plot ------------------------------
   windows(width = 15, height = 8, title = "Tukey HSD")
   op <- par(mfcol = c(2,4), las = 0, cex.main = 1.5, cex.lab = 1.5,
             cex.axis = 1.1, font.main = 1)
   
   par(fig = c(0.375, 0.625, 0.5, 1), new = T)
   plot(TukeyHSD(aov(nSpec_sample ~ env_data)))
   legend("topright",c("No. of species"), cex = 1.5, bty = "n")
   
   par(fig = c(0, 0.25, 0, 0.5), new = T)
   plot(TukeyHSD(aov(nInd_sample ~ env_data)))
   legend("topright",c("No. of individuals"), cex = 1.5, bty = "n")
   
   par(fig = c(0.25, 0.5, 0, 0.5), new = T)
   plot(TukeyHSD(aov(nSpec_rare ~ env_data)))
   legend("topright", paste("No. of species (n = ",n.min,")",sep = ""),
          cex = 1.5, bty = "n")
   
   par(fig = c(0.5, 0.75, 0, 0.5), new = T)
   plot(TukeyHSD(aov(nSpec_jade ~ env_data)))
   legend("topright",c("No. of species extrapolated"), cex = 1.5, bty = "n")
   
   par(fig = c(0.75, 1, 0, 0.5), new = T)
   plot(TukeyHSD(aov(PIE_sample ~ env_data)))
   legend("topright",c("PIE"), cex = 1.5, bty = "n")

   
   # 
   # # PIE plot and beta PIE
   # 
   # boxplot(PIE_sample ~ env_data, ylab = "PIE", notch = T)
   # 
   # beta_PIE <- PIE_plot - PIE_sample
   # test <- kruskal.test(beta_PIE ~ env_data)
   # 
   # title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 2), sep = "")
   # boxplot(beta_PIE ~ env_data, ylab = "beta PIE", notch = T, main = title)
   # plot(TukeyHSD(aov(beta_PIE ~ env_data)))
   
   outdat <- data.frame(Groups       = env_data,
                        nInd_sample  = nInd_sample,
                        nSpec_sample = nSpec_sample,
                        nSpec_rare   = nSpec_rare,
                        nSpec_jade   = nSpec_jade,
                        PIE_sample   = PIE_sample)
   
   return(outdat)
}

