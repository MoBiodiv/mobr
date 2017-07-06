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

get_pval <- function(rand, obs, n_samples)
{
   n_extremes = sum(rand < -abs(obs) | rand > abs(obs))
   p_val = n_extremes/n_samples
   return(p_val)
}

# Get group-level biodiversity indices
calc_biodiv_groups <- function(abund_mat, groups, index, n_rare)
{
   out <- list("group" = groups)
   
   groups_N <- rowSums(abund_mat) 
   
   # Number of individuals -----------------------------------------------------
   if ("N" %in% index){
      out$N  =  groups_N   
   } 
   
   # Number of species ---------------------------------------------------------
   if ("S" %in% index){
      out$S  = rowSums(abund_mat > 0)  
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if ("S_rare" %in% index){  
      
      out$S_rare <- list()
      
      n_rare <- floor(n_rare)
      if (is.null(n_rare) | !is.numeric(n_rare)){   
         N_min_group = min(groups_N)
         n_rare = N_min_group
      }
      
      for (i in 1:length(n_rare)){
         
         groups_low_n = groups_N < n_rare[i]
         
         if (sum(groups_low_n) == 1){
            warning(paste("There is",sum(groups_low_n),"group with less then",  n_rare[i],"individuals. This is removed for the calculation of rarefied richness."))
         }
         
         if (sum(groups_low_n) > 1){
            warning(paste("There are",sum(groups_low_n),"groups with less then",  n_rare[i],"individuals. These are removed for the calculation of rarefied richness."))
         }
         
         S_rare = rep(NA, nrow(abund_mat))
         S_rare[!groups_low_n] = apply(abund_mat[!groups_low_n], MARGIN = 1,
                                       rarefaction, method = "indiv",
                                       effort = n_rare[i])
         out$S_rare[[paste("n =", n_rare[i])]] = S_rare
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
      
      out$S_asymp <- S_asymp_group
   }
   
   # Probability of Interspecific Encounter (PIE)-------------------------------
   if ("PIE" %in% index){ 

      PIE_groups  = calc_PIE(abund_mat)
      out$PIE  = PIE_groups
      }
   
   # Effective number of species based on PIE ----------------------------------
   if ("ENS_PIE" %in% index){
      
      ENS_PIE_groups = vegan::diversity(abund_mat, index = "invsimpson")
      ENS_PIE_groups[!is.finite(ENS_PIE_groups)] <- NA
      
      out$ENS_PIE = ENS_PIE_groups
   }
   
   return(out)
   
}

#Get F statistics from diversity indices and grouping vector
get_F_values <- function(div_list, permute = F)
{
   group_id <- div_list$group
   if (permute)
      group_id <- sample(group_id)
   
   F_list <- list() # F values for sample level 
   
   for (i in 2:length(div_list)){
      if (names(div_list[i]) != "S_rare"){
         lm1 <- lm(div_list[[i]] ~ group_id)
         F_list[[names(div_list[i])]] <- anova(lm1)$F[1]
      } else {
         for (j in 1:length(div_list[[i]])){
            lm1 <- lm(div_list[[i]][[j]] ~ group_id)  
            F1 <- anova(lm1)$F[1]
            names(F1) <- names(div_list[["S_rare"]][j])
            F_list[["S_rare"]][j] <- list(F1)
         }
            
      }
   }
   
   return(F_list)
}

# Get group-level differences 
get_group_diff <- function(abund_mat, group_bin, index, n_rare,
                           permute = F)
{
   if (permute)
      group_bin <- sample(group_bin)
   
   abund_group = aggregate(abund_mat, by = list(group_bin), FUN = "sum")
   
   group_list <- calc_biodiv_groups(abund_mat = abund_group[,-1],
                                    groups = abund_group[,1],
                                    index = index,
                                    n_rare = n_rare)
   
   diff_list <- list() # differences for group-level 
   
   for (i in 2:length(group_list)){
      var <- names(group_list[i])
      if (var != "S_rare"){
         diff_list[[var]] <- group_list[[var]][2] - group_list[[var]][1]
      } else {
         for (j in 1:length(group_list[[i]])){
            d1 <-  group_list$S_rare[[j]][2] - group_list$S_rare[[j]][1]
            names(d1) <- names(group_list[["S_rare"]][j])
            diff_list[["S_rare"]][j] <- list(d1)
         }
      }
   }
   
   return(diff_list)
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
#' inv_stats = get_mob_stats(inv_mob_in, 'group', ref_group = "uninvaded", nperm=100, n_rare_samples = c(5,10))
#' plot(inv_stats)
get_mob_stats = function(mob_in,
                         group_var,
                         ref_group = NULL,
                         index = c("N","S","S_rare","S_asymp","PIE","ENS_PIE"),
                         n_rare_samples = NULL,
                         n_rare_groups = NULL,
                         nperm = 100)
{
   if (nperm < 1) 
       stop('Set nperm to a value greater than 1') 
   
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
   
   index <- match.arg(index, several.ok = TRUE)
   print(index)
 
   # Abundance distribution pooled in groups
   abund_group = aggregate(mob_in$comm, by = list(group_id), FUN = "sum")
   
   #Define output list
   out <- list("samples" = list("group" = group_id))
   
   # Group-level indices
   out$groups <- calc_biodiv_groups(abund_mat = abund_group[,-1],
                                    groups = abund_group[,1],
                                    index = index,
                                    n_rare = n_rare_groups)
   
   samples_N = rowSums(mob_in$comm) 
   
   # Number of individuals -----------------------------------------------------
   if ("N" %in% index){
      out$samples$N <- samples_N
   } 
   
   # Number of species ---------------------------------------------------------
   if ("S" %in% index){
      out$samples$S = rowSums(mob_in$comm > 0) 
       
      beta_S = out$group$S[group_id] / out$samples$S 
      beta_S[!is.finite(beta_S)] <- NA
      
      out$samples$beta_S <- beta_S
   }  
   
   # Rarefied richness ---------------------------------------------------------
   if ("S_rare" %in% index){  
      
      # sample level
      out$samples$S_rare <- list()
      
      n_rare_samples <- floor(n_rare_samples)
      if (is.null(n_rare_sample) | !is.numeric(n_rare_samples)){   
         N_min_sample = min(samples_N)
         n_rare_samples = N_min_sample
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
         
         S_rare = rep(NA, nrow(mob_in$comm))
         S_rare[!plots_low_n] = apply(mob_in$comm[!plots_low_n,], MARGIN = 1,
                                      rarefaction, method = "indiv",
                                      effort = n_rare_samples[i])
         out$samples$S_rare[[paste("n =", n_rare_samples[i])]] = S_rare
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
      
      out$samples$S_asymp <- S_asymp_sample
  }

   # Probability of Interspecific Encounter (PIE)-------------------------------
   if ("PIE" %in% index){ 
      
      plots_n01 = out$samples$N <= 1 # Hurlbert's PIE can only be calculated for two or more individuals
      
      if (sum(plots_n01) == 1){
         warning(paste("There is",sum(plots_n01), "plot with less than two individuals.
   This is removed for the calculation of PIE."))
      }
      
      if (sum(plots_n01) > 1){
         warning(paste("There are",sum(plots_n01), "plots with less than two individuals.
   These are removed for the calculation of PIE."))
      }
      
      PIE_samples = calc_PIE(mob_in$comm)
      out$samples$PIE = PIE_samples
   }
   
   # Effective number of species based on PIE ----------------------------------
   if ("ENS_PIE" %in% index){
      
      plots_n01 = out$samples$N <= 1

      ENS_PIE_samples = vegan::diversity(mob_in$comm, index = "invsimpson")
      ENS_PIE_samples[plots_n01] <- NA
      
      out$samples$ENS_PIE = ENS_PIE_samples
      
      beta_ENS_PIE = out$groups$ENS_PIE[group_id] / ENS_PIE_samples
      beta_ENS_PIE[!is.finite(beta_ENS_PIE)] <- NA
      out$samples$beta_ENS_PIE = beta_ENS_PIE
   }
   
   # Significance tests
   F_obs <- get_F_values(out$samples, permute = F)
   F_rand_list <- replicate(nperm, get_F_values(out$samples, permute = T), simplify = F)
   
   diff_obs <- get_group_diff(mob_in$comm, group_bin, index,n_rare = n_rare_groups, permute = F)
   diff_rand_list <- replicate(nperm, get_group_diff(mob_in$comm, group_bin,
                                                     index, n_rare = n_rare_groups,
                                                     permute = T), simplify = F)
   
   out$p_values$samples <- list()
   var_names <- names(out$samples)[-1]
   for (var in var_names){
      if (var != "S_rare"){
         F_rand <- sapply(F_rand_list,"[[",var)
         p_val <- sum(F_obs[[var]] <= F_rand) / nperm
         out$p_values$samples[[var]] <- p_val 
      } else {
         for (j in 1:length(F_obs$S_rare)){
            F_rand <- sapply(F_rand_list,function(list){list$S_rare[[j]]})
            p_val <- sum(F_obs$S_rare[[j]] <= F_rand) / nperm
            names(p_val) <- names(F_obs$S_rare[[j]])
            out$p_values$samples$S_rare[j] <- list(p_val)
         }
      }
   }
   
   out$p_values$groups <- list()
   var_names <- names(out$groups)[-1]
   for (var in var_names){
      if (var != "S_rare"){
         diff_rand <- sapply(diff_rand_list,"[[",var)
         p_val <- get_pval(diff_rand, diff_obs[[var]], nperm)
         out$p_values$groups[[var]] <- p_val 
      } else {
         for (j in 1:length(diff_obs$S_rare)){
            diff_rand <- sapply(diff_rand_list,function(list){list$S_rare[[j]]})
            p_val <- get_pval(diff_rand, diff_obs$S_rare[[j]], nperm)
            names(p_val) <- names(diff_obs$S_rare[[j]])
            out$p_values$groups$S_rare[j] <- list(p_val)
         }
      }
   }
   
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
#' plot(inv_stats) 

plot.mob_stats = function(mob_stats)
{
   var_names <- names(mob_stats$samples)[-1]
   var_names2 <- var_names[var_names != "beta_S" & var_names != "beta_ENS_PIE"]
   
   for (var in var_names2){
      
      if (var %in% c("N", "S_asymp","PIE")){
         
         par(mfrow = c(1,2), las = 1, cex.lab = 1.2)
         
         y_sample <- mob_stats$samples[[var]]
         p_val <- mob_stats$p_values$samples[[var]]
         fig_title <- paste("Sample scale\np =",p_val)
         
         boxplot(y_sample ~ group, data=mob_stats$samples, main = fig_title,
                 ylab =  var)
         
         y_group <- mob_stats$groups[[var]]
         p_val <- mob_stats$p_values$groups[[var]]
         fig_title <- paste("Group scale\np =",p_val)
         boxplot(y_group ~ group, data=mob_stats$groups, main = fig_title,
                 ylab = "", boxwex = 0)
         points(y_group ~ group, data=mob_stats$groups, pch = 19)
      }
      
      if (var %in% c("S","ENS_PIE")){
         par(mfrow = c(1,3), las = 1, cex.lab = 1.2)
         
         y_sample <- mob_stats$samples[[var]]
         p_val <- mob_stats$p_values$samples[[var]]
         fig_title <- paste("Sample scale\np =",p_val)
         boxplot(y_sample ~ group, data=mob_stats$samples, main = fig_title,
                 ylab =  var)
         
         beta_var <- paste("beta",var,sep = "_")
         y_beta <- mob_stats$samples[[beta_var]]
         p_val <- mob_stats$p_values$samples[[beta_var]]
         fig_title <- paste("Beta-diversity\np =",p_val)
         boxplot(y_beta ~ group, data=mob_stats$samples, main = fig_title)
         
         y_group <- mob_stats$groups[[var]]
         p_val <- mob_stats$p_values$groups[[var]]
         fig_title <- paste("Group scale\np =",p_val)
         boxplot(y_group ~ group, data=mob_stats$groups, main = fig_title,
                 ylab = "", boxwex = 0)
         points(y_group ~ group, data=mob_stats$groups, pch = 19)
      }
      
      if (var == "S_rare"){
         n_rows <- max(length(mob_stats$samples$S_rare),
                       length(mob_stats$groups$S_rare))   
         
         par(mfcol = c(n_rows,2), las = 1, cex.lab = 1.2)
         
         for (j in 1:length(mob_stats$samples$S_rare)){
            y_sample <- mob_stats$samples$S_rare[[j]]
            p_val <- mob_stats$p_values$samples$S_rare[[j]]
            fig_title <- paste("Sample scale",names(p_val),"\np =",p_val)
            boxplot(y_sample ~ group, data=mob_stats$samples, main = fig_title,
                    ylab =  "S_rare")   
         }
         
         y_coords <- (n_rows:0)/n_rows
         for (j in 1:length(mob_stats$groups$S_rare)){
            y_group <- mob_stats$groups$S_rare[[j]]
            
            p_val <- mob_stats$p_values$groups$S_rare[[j]]
            fig_title <- paste("Group scale",names(mob_stats$groups$S_rare[j]),
                               "\np =",p_val)
            
            par(fig = c(0.5, 1.0, y_coords[j+1], y_coords[j]), new = T)
            boxplot(y_group ~ group, data = mob_stats$group, main = fig_title,
                    ylab =  "S_rare", boxwex = 0)
            points(y_group ~ group, data=mob_stats$groups, pch = 19)
         }
      }
   }
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

