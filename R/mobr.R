## A rough sketch of what we need in the module Self reminder: need two pieces of
## inputs: 1. a data frame of plot characteristics, with columns: plot,
## group, x, y 2. a data frame where rownames are plot ids, and
## subsequent columns are species abundances

## Functions
require(pracma)


## define mobr object
make_comm_obj = function(comm, plot_attr, binary=FALSE) {
    # possibly make group_var and ref_group mandatory arguments
    out = list()
    out$tests = list(indiv=T, sampl=T, spat= T)
    # carry out some basic checks
    if (nrow(comm) < 5) {
        stop("Number of plots in community is less than five therefore only individual rarefaction will be computed")
        out$tests$samp = FALSE
        out$tests$spat = FALSE
    }
    if (nrow(comm) != nrow(plot_attr))
        stop("Number of plots in community does not equal number of plots in plot attribute table")
    if (any(row.names(comm) != row.names(plot_attr)))
        warning("Row names of community and plot attributes tables do not match")
    if (binary)  {
        warning("Only spatially-explict sampled based forms of rarefaction can be computed on binary data")
        out$tests$indiv = FALSE
        out$tests$samp = FALSE
    } else {
        if (max(comm) == 1)
            warning("Maximum abundance is 1 which suggests data is binary, change the binary argument to TRUE")
    }
    if (any(colSums(comm) == 0)) {
        warning("Some species have zero occurrences and will be dropped from the community table")
        comm = comm[colSums(comm) != 0, ]
    }
    out$comm = comm
    spat_cols = which(names(plot_attr) %in% c('x', 'y'))
    if (length(spat_cols) > 0) {
        out$env = plot_attr[ , -spat_cols]
        out$spat = plot_attr[ , spat_cols]
    }
    else {
        out$sampling$spat = FALSE
        out$env = plot_attr
        out$spat = NULL
    }
    class(out) = 'comm'
    return(out)
}

print.comm = function(x) {
    cat('Printing the head of each attribute in the object\n')
    print(lapply(x, head))
}

print.mobr = function(...) {
   # print rarefaction and delta rarefaction summaries
}

plot.mobr = function(...) {
  # plot rarefation and delta rarefaction curves
}

summary.mobr = function(...) {
   #  print summary anova style table
}

rarefaction = function(x, method, effort=NULL) {
    # analytical formulations from Cayuela et al. 2015. Ecological and biogeographic null hypotheses for
    # comparing rarefaction curves. Ecological Monographs 85:437–454.
    # Appendix A: http://esapubs.org/archive/mono/M085/017/appendix-A.php
    # possible inputs: sad or sp x site
    # community matrix examine input properties to determine if input is
    # compatible with method
    if (!any(method %in% c('indiv', 'samp')))
        stop('method must be "indiv" or "samp" for individual or sample based rarefaction, respectively')
    if (method == 'samp') {
        if (is.null(dim(x)))
            stop('For sample based rarefaction "x" must be a site x species matrix as the input')
        else {
            x = (x > 0) * 1             
            n = nrow(x) # empty sites are counted as samples 
            x = colSums(x)
        }
    }
    if (method == 'indiv') {
        if (!is.null(dim(x)))
            x = colSums(x)
        n = sum(x)
    }
    x = x[x > 0] # drop species with no observations
    S = length(x)
    if (is.null(effort))
        effort = 1:n
    if (any(effort > n))
        warning('"effort" larger than total number of samples')
    ldiv = lchoose(n, effort)
    p = matrix(0, length(effort), S)
    for (i in seq_along(effort)) {
        p[i, ] = ifelse(n - x < effort[i], 0, 
                        exp(lchoose(n - x, effort[i]) - ldiv[i]))
    }
    out = rowSums(1 - p)
    names(out) = effort
    class(out) = 'rarefaction'
    return(out)
}

table_effect_on_S = function(dat_sp, dat_plot, groups, ScaleBy = NA) {
    # Returns a data frame with the effects of SAD, N, and aggregation on diversity
    # across scales
    # not b/c of spatial ties the values will change every time 
    # this is calculated therefore best pratice may be
    # tst = replicate(20, table_effect_on_S(dat_sp, dat_plot, groups, ScaleBy), simplify=FALSE)
    # plyr::aaply(plyr::laply(tst, as.matrix), c(2, 3), mean)
    nplots = table(dat_plot$group)
    explicit_sample = sapply(groups, function(x) 
                             rarefy_sample_explicit(dat_sp, dat_plot, x, 1:min(nplots)))
    overall = as.numeric(na.omit(explicit_sample[ , 2] - explicit_sample[ , 1]))
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, groups)
    deltaSN = get_deltaSN(dat_sp, dat_plot, groups, ScaleBy) # why is this call diff
    deltaSagg = get_deltaSagg(dat_sp, dat_plot, groups)
    # Rarefy to desired abundances
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    max_level = floor(log10(avg_dens * min(nplots)))
    out = as.data.frame(matrix(NA, 4, max_level))
    row.names(out) = c("overall", "SAD", "N", "aggregation")
    names(out) = as.character(10^(1:max_level))
    for (row in c(2, 3)) {
        deltaS = unlist(list(overall, deltaSsad, deltaSN, deltaSagg)[row])
        out[row, ] = sapply(10^(1:max_level), function(x) 
                            ifelse(length(deltaS) >= x, deltaS[x], NA))
    }
    for (row in c(1, 4)) {
        deltaS = unlist(list(overall, deltaSsad, deltaSN, deltaSagg)[row])
        out_row = pchip(xi=(0:length(deltaS)) * avg_dens, yi=c(0, deltaS), 
                        x=10^(1:min(max_level, floor(log10(length(deltaS) * avg_dens)))))
        out[row, 1:length(out_row)] = out_row
    }
    out = cbind(out, c(overall[length(overall)], deltaSsad[length(deltaSsad)], 
                deltaSN[length(deltaSN)], deltaSagg[length(deltaSagg)]))
    names(out)[max_level + 1] = length(deltaSsad)
    return(out)
}

pairwise_t = function(dat_sp, dat_plot, groups, lower_N = NA) {
    dat_plot_grps = dat_plot[dat_plot$group %in% groups, ]
    dat_sp = dat_sp[match(dat_plot_grps$plot, row.names(dat_sp)), ]
    S_list = rowSums(dat_sp > 0)
    N_list = rowSums(dat_sp)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, ]/N_list[x])^2)))
    if (is.na(lower_N)) {
        rarefied_S_list = apply(dat_sp, 1, function(x) 
                                rarefaction(x, 'indiv', effort = 1:min(N_list)))
    } else {
        # Remove plots with abundance below lower_N in the analysis of rarefied S
        rarefied_S_list = apply(dat_sp, 1, function(x) 
                                if (sum(x) < lower_N)
                                    rep(NA, lower_N)
                                else 
                                    rarefaction(x, 'indiv', effort = 1:lower_N))
        if (any(is.na(rarefied_S_list))) 
            print("Warning: some plots are removed in rarefaction.")
    }
    out = as.data.frame(matrix(NA, 5, 4))
    stats_list = list(rarefied_S_list, N_list, PIE_list, S_list)
    for (i in 1:length(stats_list)) {
        stat = unlist(stats_list[i])
        stat_1 = stat[dat_plot$group == groups[1]]
        stat_2 = stat[dat_plot$group == groups[2]]
        stat_1 = stat_1[!is.na(stat_1)]
        stat_2 = stat_2[!is.na(stat_2)]
        out[ , i] = c(mean(stat_1), sd(stat_1), 
                      mean(stat_2), sd(stat_2), 
                      t.test(stat_1, stat_2)$p.val)
    }
    names(out) = c("S_rarefied", "N", "PIE", "S_raw")
    row.names(out) = c(paste(groups[1], "(mean)", sep = ""), 
                       paste(groups[1], "(sd)", sep = ""), 
                       paste(groups[2], "(mean)", sep = ""), 
                       paste(groups[2], "(sd)", sep = ""), "p_value")
    # Boxplots
    par(mfrow = c(2, 2))  # This is not ideal but I cannot get layout to work in Rstudio
    plot_names = c(paste("Rarified S at N=", 
                         ifelse(is.na(lower_N), min(N_list), lower_N), sep = ""),
                   "N", "PIE", "Raw S")
    plot_names = sapply(1:4, function(x) 
                        paste(plot_names[x], " (p=", round(out[5, x], 6), ")", sep = ""))
    for (i in 1:length(stats_list)) {
        stat = unlist(stats_list[i])
        stat_1 = stat[dat_plot$group == groups[1]]
        stat_2 = stat[dat_plot$group == groups[2]]
        stat_1 = stat_1[!is.na(stat_1)]
        stat_2 = stat_2[!is.na(stat_2)]
        boxplot(stat_1, stat_2, names = c(groups[1], groups[2]), main = plot_names[i])
    }
    return(out)
}

## Functions for plotting
plotEffectS = function(dat_sp, dat_plot, groups, Nperm = 1000,
                       CI = 0.95, ScaleBy = NA) {
    # get individual densities
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    #avg_A = mean(dat_plot$area)
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, groups)
    deltaSN = get_deltaSN(dat_sp, dat_plot, groups, ScaleBy)
    deltaSagg = get_deltaSagg(dat_sp, dat_plot, groups)
    deltaS_list = list(deltaSsad, deltaSN, deltaSagg)
    
    sad_CI = null_sad(dat_sp, dat_plot, groups, Nperm, CI)
    N_CI = null_N(dat_sp, dat_plot, groups, Nperm, CI, ScaleBy)
    agg_CI = null_agg(dat_sp, dat_plot, groups, Nperm, CI)
    null_list = list(sad_CI, N_CI, agg_CI)
    
    par(mfrow = c(1, 3), mar = c(5, 4, 6.5, 2) + 0.1)
    main_list = c("SAD", "N", "Aggregation")
    for (i in 1:3) {
        deltaS = unlist(deltaS_list[i])
        null_CI = unlist(null_list[i], recursive = F)
        if (i %in% c(1, 2)) {
            plot(1:length(null_CI$lowerCI), null_CI$lowerCI, type = "l",
                 col = "grey84", ylab = expression(Delta ~ "S"), 
                 xlab = "Number of Individuals", 
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(1, length(deltaS)), cex.lab = 1.5, cex.axis = 1.5)
            polygon(c(1:length(null_CI$lowerCI), length(null_CI$lowerCI):1), 
                    c(null_CI$lowerCI, rev(null_CI$upperCI)), col = "grey84", border = NA)
            lines(1:length(null_CI$mean), null_CI$mean, type = "l", lwd = 2)
            par(new = T)
            plot((1:length(deltaS))/avg_dens * 1, deltaS, col = "red", lwd = 2, 
                 type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(1, length(deltaS))/avg_dens * 1)
            axis(3, cex.axis = 1.5)
            mtext("Number of Plots", side = 3, line = 2.5)
            title(main = main_list[i], line = 4.5, cex.main = 1.8)
        } else {
            plot(1 * 1:length(null_CI$lowerCI), null_CI$lowerCI, type = "l", 
                 col = "grey84", ylab = expression(Delta ~ "S"), xlab = "Number of Plots",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(0, 1 * length(deltaS)), cex.lab = 1.5, cex.axis = 1.5)
            polygon(c(1 * 1:length(null_CI$lowerCI), 1 * length(null_CI$lowerCI):1), 
                    c(null_CI$lowerCI, rev(null_CI$upperCI)), col = "grey84", border = NA)
            lines(1 * 1:length(null_CI$mean), null_CI$mean, type = "l", lwd = 2)
            par(new = T)
            plot((1:length(deltaS)) * avg_dens, deltaS, col = "red", lwd = 2,
                 type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)), 
                 xlim = c(0, avg_dens * length(deltaS)))
            axis(3, cex.axis = 1.5)
            mtext("Number of Individuals", side = 3, line = 2.5)
            title(main = main_list[i], line = 4.5, cex.main = 1.8)
        }
    }
}

plotSADs = function(dat_sp, dat_plot, col = NA) {
    # TO DO: add check to ensure that col is the same length as treatments
    require(scales)
    par(mfrow = c(1, 1))
    grps = unique(dat_plot$group)
    if (is.na(col)) 
        col = rainbow(length(grps))
    plot(1, type = "n", xlab = "% abundance (log scale)", ylab = "% species", 
         xlim = c(0.01, 1), ylim = c(0, 1), log = "x")
    for (i in 1:length(grps)) {
        col_grp = col[i]
        plots_grp = dat_plot[dat_plot$group == grps[i], 1]
        dat_grp = dat_sp[match(plots_grp, row.names(dat_sp)), ]
        for (j in 1:nrow(dat_grp)) {
            sad_row = as.numeric(sort(dat_grp[j, dat_grp[j, ] != 0]))
            s_cul = 1:length(sad_row)/length(sad_row)
            n_cul = sapply(1:length(sad_row), function(x) sum(sad_row[1:x]) / sum(sad_row))
            lines(n_cul, s_cul, col = alpha(col_grp, 0.5), lwd = 1, type = "l")
        }
    }
    legend("bottomright", grps, col = col, lwd = 2)
}

plotSNpie = function(dat_sp, dat_plot, col = NA) {
    # TO DO: add check to ensure that col is the same length as treatments
    require(rgl)
    grps = unique(dat_plot$group)
    if (is.na(col)) 
        col = rainbow(length(grps))
    S_list = rowSums(dat_sp > 0)
    N_list = rowSums(dat_sp)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, ]/N_list[x])^2)))
    grp_list = as.character(dat_plot$group[match(dat_plot$plot, row.names(dat_sp))])
    col_list = sapply(grp_list, function(x) col[which(grps == x)])
    plot3d(S_list, N_list, PIE_list, "S", "N", "PIE", col = col_list, size = 8)
} 

# Auxillary function: difference between the ind-based rarefaction and the sample-based rarefaction for one group
#   with the evaluation sample size (number of individuals) defined by ref_dens
# Output: a two-column data frame, with sample size
#   and deltaS (effect of N)
effect_of_N = function(comm_group, ref_dens){
  group_sad = colSums(comm_group)
  S_samp_rare = as.numeric(rarefaction(comm_group, 'samp'))
  plot_to_ind = round(ref_dens * (1:nrow(comm_group)))
  plot_to_ind = plot_to_ind[which(plot_to_ind <= sum(comm_group))]
  interp_S_samp_rare = pchip(c(1, ref_dens * (1:nrow(comm_group))), c(1, S_samp_rare), plot_to_ind)
  S_indiv_rare = rarefaction(group_sad, 'indiv', effort = plot_to_ind)
  deltaS = as.numeric(S_indiv_rare) - interp_S_samp_rare
  out = data.frame(cbind(plot_to_ind, deltaS))
  names(out) = c('effort', 'deltaS')
  return(out)
}

# enforce_min_group_size = function(comm, group_data, min_group_size) {
#     group_cts = table(group_data)
#     if (any(group_cts < min_group_size)) {
#         small_groups = names(group_cts)[group_cts < min_group_size]
#         large_groups = names(group_cts)[group_cts >= min_group_size]
#         if (length(large_groups) == 0)
#             stop(paste('No groups have at least', min_group_size, 'replicates'))
#         if (length(large_groups) == 1) 
#             stop(paste('Only the group', large_groups, 'has at least',
#                        min_group_size, 'replicates'))
#         warning(paste('The groups', paste(small_groups, collapse=', '),
#                       'have less than', min_group_size, 
#                       'plots and therefore will be dropped'))
#         row_indices = which(env_data == small_groups)
#         comm$comm = comm$comm[-row_indices, ]
#         comm$env = comm$env[-row_indices, ]
#         comm$spat = comm$spat[-row_indices, ]
#     }
#     return(comm)
# }
# 
# Auxillary function: spatially-explicit sample-based rarefaction 
# 
rarefy_sample_explicit = function(comm_one_group, xy_one_group) {
  #plot_grp = dat_plot[dat_plot$group == group, ] 
  #sp_grp = dat_sp[dat_plot$group == group, ]
  row.names(comm_one_group) = as.character(seq(nrow(comm_one_group)))
  explicit_loop = matrix(0, nrow(comm_one_group), nrow(comm_one_group))
  pair_dist = as.matrix(dist(xy_one_group))
  for (i in 1:nrow(comm_one_group)) {
    focal_site = row.names(comm_one_group)[i]
    dist_to_site = pair_dist[i, ]
    # Shuffle plots, so that tied grouping is not biased by original order.
    new_order = sample(1:nrow(comm_one_group))  
    plots_new = row.names(comm_one_group)[new_order]
    dist_new = dist_to_site[new_order]
    plots_new_ordered = plots_new[order(dist_new)]
    # Move focal site to the front
    plots_new_ordered = c(focal_site, 
                          plots_new_ordered[plots_new_ordered != focal_site])  
    comm_ordered = comm_one_group[match(row.names(comm_one_group), plots_new_ordered), ]
    # 1 for absence, 0 for presence
    comm_bool = as.data.frame((comm_ordered == 0) * 1) 
    rich = cumprod(comm_bool)
    explicit_loop[ , i] = as.numeric(ncol(comm_one_group) - rowSums(rich))
  }
  explicit_S = apply(explicit_loop, 1, mean)
  return(explicit_S)
}

get_delta_stats = function(comm, env_var, ref_group=NULL, 
                           tests=c('indiv', 'sampl', 'spat'),
                           type='discrete', log_scale=FALSE, inds=NULL, min_plot = NULL, 
                           density_stat ='mean', corr='spearman', 
                           nperm=1000) {
  # Inputs:
  # comm - a 'comm' type object, with attributes ...
  # type - if the envronmental variable is 'discrete' vs 'continous'
  # env_var - string, name of the envronmental variable
  # ref_group - the group that will be used as a reference in pair-wise comparisons (ie., all other groups will be compared with it).
  #   This argument is only needed then type == 'discrete'.
  # test - tests to be included. A single value in c('indiv', 'sampl', 'spat'), or a combination of them as a list. 
  #   The default is test = c('indiv', 'sampl', 'spat') (all three tests).
  # type - pair-wise comparisons ('discrete') or regression ('continuous'). For the discrete case, ref_group has to be specified as the baseline.
  # inds - argument for individual-based rarefaction. If given a single number, will be taken as the number of points for
  #   for rarefaction. log_scale will then be employed (see below). If given a list of numbers, will be taken as the levels of N 
  #   for rarefaction. Default (NULL) will result in rarefaction at all possible N.
  # min_plots - minimal number of plots for each group for test 'spat', where plots are randomized within groups for null test.
  #   If a value is given, all groups with less than min_plots plots will be removed. If NULL (default), all groups are kept.
  #   Warnings will be issued if 1. there is only one group left in discrete case, or 2. there are less than three groups left in continuous case, 
  #   or 3. ref_group is removed in discrete case. In these three cases, the test will not be carried out.
  #   Another warning will be issued if any of the remaining groups have less than five plots, but the test will be carried out in this case.
  # log_scale - If number of points for rarefaction is given, log_scale = TRUE leads to values evenly spaced on log scale,
  #   log_scale = FAlSE (default) leads to values evenly spaced on arithmetic scale.
  # density_stat - reference density used in converting number of plots to number of individuals. Can take one of three values:
  #   'mean', 'max', 'min'. If 'mean', the average plot-level abundance across all plots are used. If 'min' or 'max', 
  #   the minimum/maximum plot-level abundance is used.
  # corr - which correlation is used for S/delta S vs envronmental variable in the continuous case. Can be 'spearman' (default, rank correlation)
  #   or 'pearson' (less recommended because of potential nonlinearity)
  # nperm - number of iterations for null models
  # Output:
  # out - a 'mobr' type object, with attributes...
  # check tests
    env_data = comm$env[, env_var]
    groups = unique(env_data)
    if (!(type %in% c('continuous', 'discrete')))
      stop('"type" has to be "discrete" or "continuous".')
    if (type == 'continuous' & !(corr %in% c('spearman', 'pearson')))
      stop('"corr" has to be "spearman" or "pearson".')
    if ('factor' %in% class(env_data)) {
        if (type == 'continuous') {
            group_vals = data.frame(groups=groups, 
                                    values=as.integer(env_data)[match(groups, env_data)])
            warning(paste(env_var, 'is a factor but will be treated as a continous variable for the analysis which the following values'))
            print(group_vals)
        }
        else if (type == 'discrete') {
            if (is.null(ref_group))
                stop('For a discrete analysis you must specify a ref_group to compare groups to')
            else if (!(ref_group %in% env_data))
                stop(paste('Reference group is not present in', env_var))
        }
    } else if (type == 'discrete') {
        warning(paste(env_var, 'is not a factor and each unique value will be treated as a grouping variable'))
        if (is.null(ref_group))
            stop('For a discrete analysis you must specify a ref_group to compare groups to')
        else if (!(ref_group %in% env_data))
            stop(paste('Reference group is not present in', env_var))
    }
    test_status = sapply(tests, function(x) 
                         eval(parse(text=paste('comm$tests$', x, sep=''))))
    approved_tests = tests[which(test_status == TRUE)]
    if (any(test_status == FALSE)) {
        tests_string = paste(tests[which(tests %in% approved_tests)], collapse=' and ')
        warning(paste('Based upon the attributes of the community object only the following tests will be performed:',
                      tests_string))
    }
    # TODO: groups as row names in all outputs?
    out = list()  # This is the object with the outputs
    out$type = type

    group_sad = aggregate(comm$comm, by=list(env_data), sum)
    group_levels = group_sad[ , 1]
    group_sad = group_sad[ , -1]
    group_minN = min(rowSums(group_sad))
    group_plots = data.frame(table(env_data)) # Number of plots within each group
    plot_abd = apply(comm$comm, 1, sum)
    #rows_keep_group = which(env_data %in% group_keep)
    #comm_group = comm$comm[rows_keep_group, ]
    #env_var_keep = comm$env[rows_keep_group, env_var]
    #keep_group_sad = aggregate(comm_group, by = list(env_var_keep), FUN = sum)
    #row.names(keep_group_sad) = keep_group_sad[, 1]
    #keep_group_sad = keep_group_sad[, -1]
    if (density_stat == 'mean')
        ref_dens = sum(comm$comm) / nrow(comm$comm)
    else if (density_stat == 'max')
        ref_dens = max(rowSums(comm$comm))
    else if (density_stat == 'min')
       ref_dens = min(rowSums(comm$comm))
    else 
       stop('The argument ref must be set to mean, min or max')
    
    # 1. Individual-based rarefaction (effect of SAD) vs env_var vs N
    if ('indiv' %in% approved_tests) {
      if (is.null(inds))
        ind_sample_size = seq(group_minN)
      else if (length(inds) > 1)
        ind_sample_size = inds
      else {
        if (log_scale == T)
          ind_sample_size = floor(exp(seq(inds) * log(group_minN) / inds))
        else
          ind_sample_size = floor(seq(inds) * group_minN / inds)
      }
      ind_rare = data.frame(apply(group_sad, 1, function(x) 
        rarefaction(x, 'indiv', ind_sample_size)))
      row.names(ind_rare) = NULL
      out$indiv_rare = cbind(ind_sample_size, ind_rare)
      names(out$indiv_rare) = c('sample', as.character(group_levels))
      
      if (type == 'continuous'){
        ind_cor = apply(ind_rare, 1, function(x) 
          cor(x, as.numeric(group_levels), method = corr))
        
        # Null test 
        sp_extent = unlist(apply(group_sad, 1, function(x) rep(1:ncol(group_sad), x)))
        env_extent = rep(group_levels, time = apply(group_sad, 1, sum))
        null_ind_r_mat = matrix(NA, nperm, length(ind_sample_size))
        for (i in 1:nperm){
          env_shuffle = sample(env_extent)
          sad_shuffle = sapply(group_levels, function(x) 
            as.integer(table(factor(sp_extent[env_shuffle == x], levels=1:ncol(group_sad)))))
          perm_ind_rare = apply(sad_shuffle, MARGIN = 2, function(x)
            rarefaction(x, 'indiv', ind_sample_size))
          null_ind_r_mat[i, ] = apply(perm_ind_rare, 1, function(x){cor(x, as.numeric(group_levels), method = corr)})
        }
        ind_r_null_CI = apply(null_ind_r_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975))) # 95% CI
        out$continuous$indiv = data.frame(cbind(ind_sample_size, ind_cor, t(ind_r_null_CI)))
        names(out$continuous$indiv) = c('effort_ind', 'r_emp', 'r_null_low', 'r_null_median', 'r_null_high')
      }
      else { # discrete case
        ref_sad = group_sad[which(as.character(group_levels) == as.character(ref_group)), ]
        out$discrete$indiv = data.frame(sample = numeric(), group = character(), deltaS_emp = numeric(), deltaS_null_low = numeric(), 
                                        deltaS_null_median = numeric(), deltaS_null_high = numeric(), stringsAsFactors = F)
        for (group in group_levels){
          if (as.character(group) != as.character(ref_group)){
            deltaS = out$indiv_rare[, as.character(group)] - out$indiv_rare[, as.character(ref_group)]
            group_levels_pairwise = c(as.character(ref_group), as.character(group))
            level_sad = group_sad[which(as.character(group_levels) == as.character(group)), ]
            comp_sad = rbind(ref_sad, level_sad)
            sp_extent = unlist(apply(comp_sad, 1, function(x) rep(1:ncol(comp_sad), x)))
            env_extent = rep(group_levels_pairwise, time = apply(comp_sad, 1, sum))
            null_ind_deltaS_mat = matrix(NA, nperm, length(ind_sample_size))
            for (i in 1:nperm){
              env_shuffle = sample(env_extent)
              sad_shuffle = sapply(group_levels_pairwise, function(x) 
                as.integer(table(factor(sp_extent[env_shuffle == x], levels=1:ncol(comp_sad)))))
              perm_ind_rare = apply(sad_shuffle, MARGIN = 2, function(x)
                rarefaction(x, 'indiv', ind_sample_size))
              null_ind_deltaS_mat[i, ] = perm_ind_rare[, as.character(group)] - perm_ind_rare[, as.character(ref_group)]
            }
            ind_deltaS_null_CI = apply(null_ind_deltaS_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            ind_group = data.frame(cbind(rep(as.character(group), length(ind_sample_size)),ind_sample_size,  
                                         deltaS, t(ind_deltaS_null_CI)))
            out$discrete$indiv = rbind(out$discrete$indiv, ind_group, stringsAsFactors = F)
          }
        }
        names(out$discrete$indiv) = c('group', 'effort_ind', 'deltaS_emp', 'deltaS_null_low', 'deltaS_null_median', 'deltaS_null_high')
      }
    }
    
    # Sample-based spatially-implicit and -explicit rarefaction 
    if ('sampl' %in% approved_tests | 'spat' %in% approved_tests){
      if ('spat' %in% approved_tests)
        out$sample_rare = data.frame(group = character(), sample_plot = numeric(), impl_S = numeric(), expl_S = numeric())
      else
        out$sample_rare = data.frame(group = character(), sample_plot = numeric(), impl_S = numeric())
      for (group in group_levels){
        comm_group = comm$comm[as.character(env_data) == as.character(group), ]
        impl_S = as.numeric(rarefaction(comm_group, 'samp'))
        sample_rare_group = data.frame(cbind(rep(as.character(group), length(impl_S)), seq(length(impl_S)), impl_S))
        if ('spat' %in% approved_tests){
          xy_group = comm$spat[as.character(env_data) == as.character(group), ]
          expl_S = rarefy_sample_explicit(comm_group, xy_group)
          sample_rare_group = cbind(sample_rare_group, expl_S)
        }
        out$sample_rare = rbind(out$sample_rare, sample_rare_group)
      }
      if ('spat' %in% approved_tests)
        names(out$sample_rare) = c('group', 'sample_plot', 'impl_S', 'expl_S')
      else
        names(out$sample_rare) = c('group', 'sample_plot', 'impl_S')
    }
    
    # 2. Sample-based rarefaction (effect of density) vs env_var vs N
    if ('sampl' %in% approved_tests){
      # TODO: Checks?
      if (type == 'continuous'){
        if (density_stat == 'mean')
          ref_dens = sum(comm$comm) / nrow(comm$comm)
        else if (density_stat == 'max')
          ref_dens = max(rowSums(comm$comm))
        else if (density_stat == 'min')
          ref_dens = min(rowSums(comm$comm))
        else 
          stop('The argument ref must be set to mean, min or max')
        
        effect_N_by_group = data.frame(matrix(NA, ncol = 4, nrow = max(group_plots$Freq)))
        for (i in 1:length(group_levels)){
          group = group_levels[i]
          comm_group = comm$comm[which(env_data == group), ]
          group_effect_N = effect_of_N(comm_group, ref_dens)
          if (i == 1)
            effect_N_by_group[, i] = group_effect_N$effort[1:nrow(effect_N_by_group)]
          effect_N_by_group[, i + 1] = group_effect_N$deltaS[1:nrow(effect_N_by_group)]
        }
        effect_N_by_group = effect_N_by_group[complete.cases(effect_N_by_group), ]
        r_emp = apply(effect_N_by_group[, 2:4], 1, function(x)
          cor(x, as.numeric(group_levels), method = corr))
        
        # Null model
        null_N_r_mat = matrix(NA, nperm, length(r_emp))
        for (i in 1:nperm){
          plot_abd_perm = as.numeric(sample(plot_abd))
          sp_draws = sapply(1:nrow(comm$comm), function(x)
            sample(rep(1:ncol(comm$comm), as.numeric(group_sad[which(group_levels == env_data[x]), ])),
                   size = plot_abd_perm[x], replace = T))
          comm_perm = t(sapply(1:nrow(comm$comm), function(x)
            table(c(1:ncol(comm$comm), sp_draws[[x]])) - 1 ))
          effect_N_perm = data.frame(matrix(NA, ncol = 4, nrow = max(group_plots$Freq)))
          for (j in 1:length(group_levels)){
            group = group_levels[j]
            comm_group = comm_perm[which(env_data == group), ]
            group_N_perm = effect_of_N(comm_group, ref_dens)
            if (j == 1)
              effect_N_perm[, j] = group_N_perm$effort[1:nrow(effect_N_perm)]
            effect_N_perm[, j + 1] = group_N_perm$deltaS[1:nrow(effect_N_perm)]
          }
          effect_N_perm = effect_N_perm[complete.cases(effect_N_perm), ]
          null_N_r_mat[i, ] = apply(effect_N_perm[, 2:4], 1, function(x)
            cor(x, as.numeric(group_levels), method = corr))[1:ncol(null_N_r_mat)]
        }
        N_r_null_CI = apply(null_N_r_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
        out$continuous$N = data.frame(cbind(effect_N_by_group[, 1], r_emp, t(N_r_null_CI)))
        names(out$continuous$N) = c('effort_ind', 'r_emp', 'r_null_low', 'r_null_median', 'r_null_high')
      }
      else {
        for (group in group_levels){
          if (as.character(group) != as.character(ref_group)){
            comm_2groups = comm$comm[which(as.character(env_data) %in% 
                                       c(as.character(group), as.character(ref_group))), ]
            env_2groups = env_data[which(as.character(env_data) %in% 
                                           c(as.character(group), as.character(ref_group)))]
            if (density_stat == 'mean')
              ref_dens = sum(comm_2groups) / nrow(comm_2groups)
            else if (density_stat == 'max')
              ref_dens = max(rowSums(comm_2groups))
            else if (density_stat == 'min')
              ref_dens = min(rowSums(comm_2groups))
            else 
              stop('The argument ref must be set to mean, min or max')
            
            effect_N_ref = effect_of_N(comm$comm[which(as.character(env_data) == as.character(ref_group)), ], ref_dens)
            effect_N_group = effect_of_N(comm$comm[which(env_data == group), ], ref_dens)
            ddeltaS_group = effect_N_group$deltaS[1:min(nrow(effect_N_ref), nrow(effect_N_group))] - 
              effect_N_ref$deltaS[1:min(nrow(effect_N_ref), nrow(effect_N_group))]
            
            null_N_deltaS_mat = matrix(NA, nperm, length(ddeltaS_group))
            for (i in 1:nperm){
              plot_abd_perm = sample(plot_abd[which(as.character(env_data) %in% c(as.character(group), as.character(ref_group)))])
              sp_draws = sapply(1:nrow(comm_2groups), function(x)
                sample(rep(1:ncol(comm_2groups), as.numeric(group_sad[which(group_levels == env_2groups[x]), ])),
                       size = plot_abd_perm[x], replace = T))
              comm_perm = t(sapply(1:nrow(comm_2groups), function(x)
                table(c(1:ncol(comm_2groups), sp_draws[[x]])) - 1 ))
              
              ID_ref = as.character(which(as.character(env_data) == as.character(ref_group)))  
              N_ref_perm = effect_of_N(comm_perm[which(row.names(comm_2groups) %in% ID_ref), ], ref_dens)
              ID_group = as.character(which(as.character(env_data) == as.character(group)))
              N_group_perm = effect_of_N(comm_perm[which(row.names(comm_2groups) %in% ID_group), ], ref_dens)
              ddeltaS_perm = N_ref_perm$deltaS[1:min(nrow(N_ref_perm), nrow(N_group_perm))] - 
                N_group_perm$deltaS[1:min(nrow(N_ref_perm), nrow(N_group_perm))]
              null_N_deltaS_mat[i, ] = ddeltaS_perm[1:ncol(null_N_deltaS_mat)]  
            }
            N_deltaS_null_CI = apply(null_N_deltaS_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
            N_group = data.frame(cbind(rep(as.character(group), length(ddeltaS_group)),effect_N_group$effort[1:length(ddeltaS_group)],  
                                         ddeltaS_group, t(N_deltaS_null_CI)))
            out$discrete$N = rbind(out$discrete$N, N_group, stringsAsFactors = F)
          }
        }
        names(out$discrete$N) = c('group', 'effort_sample', 'ddeltaS_emp', 'ddeltaS_null_low', 
                                  'ddeltaS_null_median', 'ddeltaS_null_high')
      }
    }
    # 3. Sample-based spatially-explicit rarecation (effect of aggregation) vs env_var vs N
    if ('spat' %in% approved_tests){
      if (!is.null(min_plot))
        group_keep = group_plots[which(group_plots$Freq>=min_plot), 1]
      else
        group_keep = group_levels
      
      if (length(group_keep) == 1 & type == 'discrete')
        warning('Error: pair-wise comparison cannot be conducted on one group.')
      else if (length(group_keep) < 3 & type == 'continuous')
        warning('Error: correlation analysis cannot be conducted with less than three groups.')
      else if (!(as.character(ref_group) %in% as.character(group_keep)) & type == 'discrete')
        warning('Error: reference group does not have enough plots and have been dropped.')
      else {
        sample_rare_keep = out$sample_rare[which(out$sample_rare$group %in% as.character(group_keep)), ]
        sample_rare_keep$deltaS = as.numeric(as.character(sample_rare_keep$impl_S)) - 
          as.numeric(as.character(sample_rare_keep$expl_S))
        if (min(group_plots$Freq[group_plots[, 1] %in% group_keep]) < 5)
          warning('Warning: some groups have less than 5 plots. The results of the null model are not very informative.')
        
        if (type == 'continuous'){
          min_plot_group = min(group_plots$Freq[which(as.character(group_plots[, 1]) %in% as.character(group_keep))])
          r_emp = sapply(seq(min_plot_group), function(x)
            cor(sample_rare_keep$deltaS[which(sample_rare_keep$sample_plot == x)], 
                as.numeric(sample_rare_keep$group[which(sample_rare_keep$sample_plot == x)]), method = corr))
          # Null test
          null_agg_r_mat = matrix(NA, nperm, min_plot_group)
          for (i in 1:nperm){
            xy_perm = comm$spat[sample(nrow(comm$spat)), ]
            deltaS_perm = c()
            for (group in unique(sample_rare_keep$group)){
              comm_group = comm$comm[as.character(env_data) == as.character(group), ]
              xy_perm_group = xy_perm[as.character(env_data) == as.character(group), ]
              expl_S_perm = rarefy_sample_explicit(comm_group, xy_perm_group)
              deltaS_perm = c(deltaS_perm, as.numeric(as.character(sample_rare_keep$impl_S[sample_rare_keep$group == group])) - expl_S_perm)
            }
            null_agg_r_mat[i, ] = sapply(seq(min_plot_group), function(x)
              cor(deltaS_perm[which(sample_rare_keep$sample_plot == x)], 
                  as.numeric(sample_rare_keep$group[which(sample_rare_keep$sample_plot == x)]), method = corr))
          }
          agg_r_null_CI = apply(null_agg_r_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
          out$continuous$agg = data.frame(cbind(seq(min_plot_group), r_emp, t(agg_r_null_CI)))
          names(out$continuous$agg) = c('effort_sample', 'r_emp', 'r_null_low', 'r_null_median', 'r_null_high')
        }
        
        else {
          ref_sample = sample_rare_keep[which(as.character(sample_rare_keep$group) == as.character(ref_group)), ]
          out$discrete$agg = data.frame(group = character(), sample_plot = numeric(), ddeltaS_emp = numeric(), 
                                        ddeltaS_null_low = numeric(), ddeltaS_median = numeric(), 
                                        ddeltaS_high = numeric())
          ref_comm = comm$comm[as.character(env_data) == as.character(ref_group), ]
          impl_S_ref = sample_rare_keep$impl_S[which(as.character(sample_rare_keep$group) == as.character(ref_group))]
          impl_S_ref = as.numeric(as.character(impl_S_ref))
          for (group in unique(group_keep)){
            if ((as.character(group) != as.character(ref_group))){
              min_plot_group = min(group_plots$Freq[which(as.character(group_plots[, 1]) %in% 
                                                            c(as.character(group), as.character(ref_group)))])
              ddeltaS_group = sample_rare_keep$deltaS[sample_rare_keep$group == group][1:min_plot_group] - 
                sample_rare_keep$deltaS[as.character(sample_rare_keep$group) == as.character(ref_group)][1:min_plot_group]
              
              comm_group = comm$comm[as.character(env_data) == as.character(group), ]
              impl_S_group = sample_rare_keep$impl_S[which(as.character(sample_rare_keep$group) == as.character(group))]
              impl_S_group = as.numeric(as.character(impl_S_group))
              
              null_agg_deltaS_mat = matrix(NA, nperm, min_plot_group)
              for (i in 1:nperm){
                xy_perm = comm$spat[sample(nrow(comm$spat)), ]
                xy_perm_group = xy_perm[as.character(env_data) == as.character(group), ]
                xy_perm_ref = xy_perm[as.character(env_data) == as.character(ref_group), ]
                expl_S_perm_group = rarefy_sample_explicit(comm_group, xy_perm_group)
                expl_S_perm_ref = rarefy_sample_explicit(ref_comm, xy_perm_ref)
                null_agg_deltaS_mat[i, ] = impl_S_group[1:min_plot_group] - expl_S_perm_group[1:min_plot_group] - 
                  (impl_S_ref[1:min_plot_group] - expl_S_perm_ref[1:min_plot_group])
              }
              agg_deltaS_null_CI = apply(null_agg_deltaS_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
              agg_group = data.frame(cbind(rep(as.character(group), min_plot_group),1:min_plot_group,  
                                           ddeltaS_group, t(agg_deltaS_null_CI)))
              out$discrete$agg = rbind(out$discrete$agg, agg_group, stringsAsFactors = F)
            }
          }
          names(out$discrete$agg) = c('group', 'effort_sample', 'ddeltaS_emp', 'ddeltaS_null_low', 
                                      'ddeltaS_null_median', 'ddeltaS_null_high')
        }
      }
    }
    class(out) = 'mobr'
    return(out)
}
