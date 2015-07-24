## A rough sketch of what we need in the module
## Self reminder: need two pieces of inputs: 
## 1. a data frame of plot characteristics, with columns: plot ID, treatment, x, y, area
## 2. a data frame where the 1st column is plot ID, and subsequent columns are species abundances

## Additional Functions:
## 3. Plot of diff curves
## 4. pair-wise t-test (plot, table)

## Functions
require(rareNMtests)
require(pracma)
## Functions to obtain the three kinds of curves: sample-based explicit, sample-based implicit, 
## and individual-based rescaled to sample
get_avg_dens = function(dat_sp, dat_plot, ScaleBy){
  # Auxillary function to obtain average density within plot for rescaling
  if(!is.na(ScaleBy)){ # Ask user to specify which treatment is used as "standard" for rescaling. 
    dat_scale = dat_sp[match(dat_plot[which(dat_plot[, 2] == ScaleBy), 1], dat_sp[, 1]), 2:ncol(dat_sp)]
    avg_dens = mean(apply(dat_scale, 1, sum))
  }
  else{ # If unspecified, uses min density among treatments.
    trmts = unique(dat_plot[, 2])
    avg_dens = min(sapply(trmts, function(x) 
      mean(apply(dat_sp[match(dat_plot[which(dat_plot[, 2] == x), 1], dat_sp[, 1]), 2:ncol(dat_sp)], 1, sum))))
  }
  return(avg_dens)
}

rarefy_sample_explicit = function(dat_sp, dat_plot, treatment){
  plot_trmt = dat_plot[dat_plot[, 2] == treatment, ]
  sp_trmt = dat_sp[match(plot_trmt[, 1], dat_sp[, 1]), ]
  explic_loop = matrix(0, nrow(sp_trmt), nrow(sp_trmt))
  pair_dist = as.matrix(dist(plot_trmt[, 3:4]))
  for (i in 1:nrow(plot_trmt)){
    focal_site = plot_trmt[i, 1]
    dist_to_site = pair_dist[i, ]
    new_order = sample(1:nrow(plot_trmt)) # Shuffle plots, so that tied grouping is not biased by original order.
    #new_order = 1:nrow(plot_trmt) # Alternative: no shuffles
    plots_new = plot_trmt[new_order, 1]
    dist_new = dist_to_site[new_order]
    plots_new_ordered = plots_new[order(dist_new)]
    plots_new_ordered = c(focal_site, plots_new_ordered[plots_new_ordered != focal_site]) # Move focal site to the front
    sp_ordered = sp_trmt[match(sp_trmt[, 1], plots_new_ordered), 2:ncol(sp_trmt)]
    sp_bool = as.data.frame(ifelse(sp_ordered[ ,1:ncol(sp_ordered)] == 0, 1, 0)) # 1 for absence, 0 for presence
    rich = cumprod(sp_bool)
    explic_loop[ , i] = as.numeric(ncol(dat_sp) - 1 - rowSums(rich))
  }
  explic_S = apply(explic_loop, 1, mean)
  return(explic_S)
}

rarefy_sample_implicit = function(dat_sp, dat_plot, treatment){
  plot_trmt = dat_plot[dat_plot[, 2] == treatment, ]
  sp_trmt = dat_sp[match(plot_trmt[, 1], dat_sp[, 1]), ]
  sample_S = rarefaction.sample(sp_trmt[, 2:ncol(sp_trmt)])[, 2]
  return(sample_S)
}

rarefy_individual_rescaled = function(dat_sp, dat_plot, treatment, ScaleBy = NA){
  avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
  sp_trmt = dat_sp[match(dat_plot[dat_plot[, 2] == treatment, 1], dat_sp[, 1]), ]
  trmt_sad = as.numeric(apply(sp_trmt[, 2:ncol(sp_trmt)], 2, sum))
  trmt_sad = trmt_sad[trmt_sad != 0]
  ind_S = rarefaction.individual(trmt_sad)[, 2]
  adjusted_n = 1:length(ind_S) / avg_dens
  ind_S_plot = pchip(adjusted_n, ind_S, 1:floor(max(adjusted_n))) # interpolate to integer plot counts
  return(ind_S_plot)
}

get_deltaSsad = function(dat_sp, dat_plot, treatment1, treatment2, ScaleBy = NA){
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  rescaled_ind = sapply(c(treatment1, treatment2), function(x) 
    rarefy_individual_rescaled(dat_sp, dat_plot, x, ScaleBy)[1:min(nplots)])
  deltaSsad = as.numeric(na.omit(rescaled_ind[, 2] - rescaled_ind[, 1]))
  return(deltaSsad)
}

get_deltaSN = function(dat_sp, dat_plot, treatment1, treatment2){
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  implicit_sample = sapply(c(treatment1, treatment2), function(x) rarefy_sample_implicit(dat_sp, dat_plot, x)[1:min(nplots)])
  deltaSsad = get_deltaSsad(dat_sp, dat_plot, treatment1, treatment2)[1:min(nplots)]
  deltaSN = as.numeric(na.omit(implicit_sample[, 2] - implicit_sample[, 1] - deltaSsad))
  return(deltaSN)
}

get_deltaSagg = function(dat_sp, dat_plot, treatment1, treatment2){
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  implicit_sample = sapply(c(treatment1, treatment2), function(x) rarefy_sample_implicit(dat_sp, dat_plot, x)[1:min(nplots)])
  explicit_sample = sapply(c(treatment1, treatment2), function(x) rarefy_sample_explicit(dat_sp, dat_plot, x)[1:min(nplots)])
  deltaSagg = as.numeric(na.omit(explicit_sample[, 2] - explicit_sample[, 1] - (implicit_sample[, 2] - implicit_sample[, 1])))
  return(deltaSagg)
}

null_sad = function(dat_sp, dat_plot, treatment1, treatment2, nperm = 1000, CI = 0.95, ScaleBy = NA){
  # This function generates null confidence intervals for the deltaSsad curve assuming 
  # that samples of the two treatments come from the same SAD.
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
  dat_plot_trmts = dat_plot[dat_plot[, 2] %in% c(treatment1, treatment2), ]
  dat_sp = dat_sp[match(dat_plot_trmts[, 1], dat_sp[, 1]), ] # Only keep the plots for the two given treatments
  sp_extend = unlist(sapply(1:nrow(dat_sp), function(x) rep(2:ncol(dat_sp), dat_sp[x, 2:ncol(dat_sp)])))
  trmt_extend = unlist(sapply(1:nrow(dat_sp), 
                              function(x) rep(dat_plot[dat_plot[, 1] == dat_sp[x, 1], 2], sum(dat_sp[x, 2:ncol(dat_sp)]))))
  deltaSsad_perm = matrix(NA, nperm, min(nplots))
  for (i in 1:nperm){
    trmt_S = matrix(NA, 2, min(nplots))
    trmt_shuffle = sample(trmt_extend) # Shuffle treatment label of each individual  
    sad_shuffle = sapply(c(treatment1, treatment2), function(x) 
      as.numeric(table(factor(sp_extend[trmt_shuffle == x], levels = 2:ncol(dat_sp)))))
    ind_S = sapply(1:2, function(x) rarefaction.individual(sad_shuffle[, x])[, 2])
    trmt_S = sapply(1:2, function(x) 
      pchip(1:length(unlist(ind_S[x])) / avg_dens, unlist(ind_S[x]), 1:floor(length(unlist(ind_S[x])) / avg_dens)))
    deltaSsad_perm[i, ] = (as.numeric(na.omit(unlist(trmt_S[2])[1:min(nplots)] - unlist(trmt_S[1])[1:min(nplots)])))[1:ncol(deltaSsad_perm)]
  }
  quant_lower = as.numeric(na.omit(apply(deltaSsad_perm, 2, function(x) quantile(x, (1-CI)/2, na.rm = T))))
  quant_higher =  as.numeric(na.omit(apply(deltaSsad_perm, 2, function(x) quantile(x, (1+CI)/2, na.rm = T))))
  return(list(lowerCI = quant_lower, upperCI = quant_higher))
}

null_N = function(dat_sp, dat_plot, treatment1, treatment2, nperm = 1000, CI = 0.95){
  # This function generates null confidence intervals for the deltaSN curve assuming 
  # that the two treatments do not differe systematically in N.
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  dat_plot_trmts = dat_plot[dat_plot[, 2] %in% c(treatment1, treatment2), ]
  dat_sp = dat_sp[match(dat_plot_trmts[, 1], dat_sp[, 1]), ] # Only keep the plots for the two given treatments
  n_plot = apply(dat_sp[, 2:ncol(dat_sp)], 1, sum) # Abundance within each plot
  sad_row = sapply(1:nrow(dat_sp), function(x) rep(2:ncol(dat_sp), dat_sp[x, 2:ncol(dat_sp)]))
  deltaSN_perm = matrix(NA, nperm, min(nplots))
  for (i in 1:nperm){
    n_plot_shuffle = sample(n_plot)
    dat_sp_perm = as.data.frame(matrix(0, nrow(dat_sp), ncol(dat_sp)))
    dat_sp_perm[, 1] = dat_sp[, 1]
    new_counts = sapply(1:nrow(dat_sp), function(x) as.numeric(table(factor(sample(unlist(sad_row[x]), 
      n_plot_shuffle[x], replace = T), levels = 2:ncol(dat_sp)))))
    dat_sp_perm[, 2:ncol(dat_sp_perm)] = as.data.frame(t(new_counts))
    deltaSN_perm[i, ] = (get_deltaSN(dat_sp_perm, dat_plot, treatment1, treatment2))[1:min(nplots)]  
  }
  quant_lower = as.numeric(na.omit(apply(deltaSN_perm, 2, function(x) quantile(x, (1-CI)/2, na.rm = T))))
  quant_higher =  as.numeric(na.omit(apply(deltaSN_perm, 2, function(x) quantile(x, (1+CI)/2, na.rm = T))))
  return(list(lowerCI = quant_lower, upperCI = quant_higher))
}

null_agg = function(dat_sp, dat_plot, treatment1, treatment2, nperm = 1000, CI = 0.95){
  # This function generates null confidence intervals for the deltaSagg curve assuming 
  # that plots within each treatment have no spatial structure.
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  trmt_index = list(which(dat_plot[, 2] == treatment1), which(dat_plot[, 2] == treatment2))
  deltaSagg_perm = matrix(NA, nperm, min(nplots))
  for (i in 1:nperm){
    index_shuffle = list(sample(unlist(trmt_index[1])), sample(unlist(trmt_index[2])))
    new_order = unlist(index_shuffle)[order(unlist(trmt_index))]
    plot_shuffle_loc = dat_plot
    plot_shuffle_loc[unlist(trmt_index), 3:4] = plot_shuffle_loc[new_order, 3:4]
    deltaSagg_perm[i, ] = get_deltaSagg(dat_sp, plot_shuffle_loc, treatment1, treatment2)
  }
  quant_lower = as.numeric(na.omit(apply(deltaSagg_perm, 2, function(x) quantile(x, (1-CI)/2, na.rm = T))))
  quant_higher =  as.numeric(na.omit(apply(deltaSagg_perm, 2, function(x) quantile(x, (1+CI)/2, na.rm = T))))
  return(list(lowerCI = quant_lower, upperCI = quant_higher))
}

table_effect_on_S = function(dat_sp, dat_plot, treatment1, treatment2, ScaleBy = NA){
  # Returns a data frame with the effects of SAD, N, and aggregation on diversity across scales
  nplots = c(nrow(dat_plot[dat_plot[, 2] == treatment1, ]), nrow(dat_plot[dat_plot[, 2] == treatment2, ]))
  explicit_sample = sapply(c(treatment1, treatment2), function(x) rarefy_sample_explicit(dat_sp, dat_plot, x)[1:min(nplots)])
  overall = c(0, as.numeric(na.omit(explicit_sample[, 2] - explicit_sample[, 1])))
  deltaSsad = c(0, get_deltaSsad(dat_sp, dat_plot, treatment1, treatment2, ScaleBy))
  deltaSN = c(0, get_deltaSN(dat_sp, dat_plot, treatment1, treatment2))
  deltaSagg = c(0, get_deltaSagg(dat_sp, dat_plot, treatment1, treatment2))
  # Rarefy to desired abundances
  avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
  max_level = floor(log10(avg_dens * min(nplots)))
  out = sapply(list(overall, deltaSsad, deltaSN, deltaSagg), function(x)
    pchip(0:min(nplots) * avg_dens, x, 10 ^ (1:max_level)))
  out = as.data.frame(t(out))
  row.names(out) = c('overall', 'SAD', 'N', 'aggregation')
  names(out) = as.character(10 ^ (1:max_level))
  out$maxN = c(overall[length(overall)], deltaSsad[length(deltaSsad)], deltaSN[length(deltaSN)], deltaSagg[length(deltaSagg)])
  return(out)
}

## Functions for plotting
plotSADs = function(dat_sp, dat_plot, col = NA){
  # TO DO: add check to ensure that col is the same length as treatments
  require(scales)
  trmts = unique(dat_plot[, 2])
  if (is.na(col)) col = rainbow(length(trmts))
  plot(1, type="n", xlab='% abundance (log scale)', ylab='% species', xlim=c(0.01, 1), ylim=c(0, 1))
  for(i in 1:length(trmts)){
    col_trmt = col[i]
    plots_trmt = dat_plot[dat_plot[, 2] == trmts[i], 1]
    dat_trmt = dat_sp[match(plots_trmt, dat_sp[, 1]), 2:ncol(dat_sp)]
    for (j in 1:nrow(dat_trmt)){
      sad_row = as.numeric(sort(dat_trmt[j, dat_trmt[j, ] != 0]))
      s_cul = 1:length(sad_row) / length(sad_row)
      n_cul = sapply(1:length(sad_row), function(x) sum(sad_row[1:x]) / sum(sad_row))
      lines(n_cul, s_cul, col = alpha(col_trmt, 0.5), lwd = 1, type = 'l')
    }
  }
  legend('bottomright', trmts, col = col, lwd = 2)
}

plotSNpie = function(dat_sp, dat_plot, col = NA){
  # TO DO: add check to ensure that col is the same length as treatments
  require(rgl)
  trmts = unique(dat_plot[, 2])
  if (is.na(col)) col = rainbow(length(trmts))
  S_list = sapply(1:nrow(dat_sp), function(x) length(which(dat_sp[x, 2:ncol(dat_sp)] != 0)))
  N_list = apply(dat_sp[, 2:ncol(dat_sp)], 1, sum)
  PIE_list = sapply(1:nrow(dat_sp), function(x) 
    N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, 2:ncol(dat_sp)] / N_list[x])^2)))
  trmt_list = as.character(dat_plot[match(dat_plot[, 1], dat_sp[, 1]), 2])
  col_list = sapply(trmt_list, function(x) col[which(trmts == x)])
  plot3d(S_list, N_list, PIE_list, 'S', 'N', 'PIE', col = col_list, size = 8)
}