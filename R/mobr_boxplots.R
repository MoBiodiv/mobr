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
   
   par(mfrow = c(2,3))
   boxplot(nInd_sample ~ comm$env[, env_var], ylab = "No. of individuals") 
   boxplot(nSpec_sample ~ comm$env[, env_var], ylab = "No. of species") 
   boxplot(nSpec_rare ~ comm$env[, env_var], ylab = "No. of species - rarefied")
   
   boxplot(PIE_sample ~ comm$env[, env_var], ylab = "PIE of samples")
   
   betaPIE <- PIE_plot - PIE_sample
   boxplot(betaPIE ~ comm$env[, env_var], ylab = "beta PIE")
}

