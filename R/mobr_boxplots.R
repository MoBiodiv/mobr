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
   PIE_plot <- diversity(colSums(comm$comm), index = "simpson")
   
   par(mfcol = c(2,5))
   
   test <- kruskal.test(nInd_sample ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 2), sep = "")
   boxplot(nInd_sample ~ env_data, ylab = "No. of individuals", notch = T, main = title) 
   plot(TukeyHSD(aov(nInd_sample ~ env_data)))
   
   test <- kruskal.test(nSpec_sample ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 2), sep = "")
   boxplot(nSpec_sample ~ env_data, ylab = "No. of species", notch = T, main = title)
   plot(TukeyHSD(aov(nSpec_sample ~ env_data)))
   
   test <- kruskal.test(nSpec_rare ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 2), sep = "")
   boxplot(nSpec_rare ~ env_data, ylab = "No. of species - rarefied", notch = T, main = title)
   plot(TukeyHSD(aov(nSpec_rare ~ env_data)))
   
   test <- kruskal.test(PIE_sample ~ env_data)
   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 2), sep = "")
   boxplot(PIE_sample ~ env_data, ylab = "PIE of samples", notch = T, main = title)
   plot(TukeyHSD(aov(PIE_sample ~ env_data)))
   
   beta_PIE <- PIE_plot - PIE_sample
   test <- kruskal.test(beta_PIE ~ env_data)

   title <- paste("Kruskal-Test: p = ", round(test$p.value, digits = 2), sep = "")
   boxplot(beta_PIE ~ env_data, ylab = "beta PIE", notch = T, main = title)
   plot(TukeyHSD(aov(beta_PIE ~ env_data)))
   
   outdat <- data.frame(Groups       = env_data,
                        nInd_sample  = nInd_sample,
                        nSpec_sample = nSpec_sample,
                        nSpec_rare   = nSpec_rare,
                        PIE_sample   = PIE_sample,
                        beta_PIE      = beta_PIE)
   
   return(outdat)
}

