library(devtools)
install_github('JohnsonHsieh/Jade')
library(Jade)

require(pracma)


#' Create a community ('comm') object.
#' 
#' The 'comm' object will be passed on for analyses of biodiversity across scales.
#' 
#'  @param comm plot (rows) by species (columns) matrix. Values can be species abundances
#'  or presence/absence (1/0).
#'  @param plot_attr matrix which includes the environmental attributes and spatial 
#'  coordinates of the plots. Environmental attributes are mandatory, while spatial
#'  coordinates are not. If spatial coordinates are provided, the column(s) has to have
#'  names "x" and/or "y". 
#'  @param binary whether the plot by species matrix "comm" is in abundances or presence/absence.
#'  @return a "comm" object with four attributes. "comm" is the plot by species matrix. 
#'  "env" is the environmental attribute matrix, without the spatial coordinates. "spat" 
#'  contains the spatial coordinates (1-D or 2-D). "tests" specifies whether each of the 
#'  three tests in the biodiversity analyses is allowed by data.
#'  @export
#'  @examples
#'  {
#'  library(vegan)
#'  data(mite)
#'  data(mite.env)
#'  data(mite.xy)
#'  mite_comm = make_comm_obj(mite, cbind(mite.env, mite.xy))
#'  }
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
        comm = comm[, colSums(comm) != 0]
    }
    out$comm = data.frame(comm)
    spat_cols = which(names(plot_attr) %in% c('x', 'y'))
    if (length(spat_cols) > 0) {
        out$env = data.frame(plot_attr[ , -spat_cols])
        out$spat = data.frame(plot_attr[ , spat_cols])
    }
    else {
        out$tests$spat = FALSE
        out$env = data.frame(plot_attr)
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

plot.mobr = function(mobr, group = NULL) {
  # plot rarefation and delta rarefaction curves
  # Input: 
  # mobr object
  # type: 'discrete' or 'continuous'
  # group: which group to plot. Only required for type = 'discrete' and there are more than one 
  #   pair-wise comparison
  type = mobr$type
  tests = c('indiv', 'N', 'agg')
  names = c('Effect of SAD', 'Effect of N', 'Effect of Aggregation')
  par(mfrow = c(1, 3))
  xlabs = c('number of individuals', 'number of individuals', 'number of plots')
  if (type == 'discrete'){
    ylabs = rep('delta-S', 3)
    if (is.null(group) & length(unique(mobr[[type]][[tests[1]]][, 1])) > 1)
      stop("Error: 'group' has to be specified.")
    for (i in 1:3){
      if (is.null(group))
        mobr_group_test = mobr[[type]][[tests[i]]]
      else {
        mobr_group_test = mobr[[type]][[tests[i]]]
        mobr_group_test = mobr_group_test[which(as.character(mobr_group_test$group) == as.character(group)), ]
      }
      mobr_group_test = mobr_group_test[complete.cases(mobr_group_test), ]
      for (icol in 2:ncol(mobr_group_test))
        mobr_group_test[, icol] = as.numeric(as.character(mobr_group_test[, icol]))
      plot(mobr_group_test[, 2], mobr_group_test[, 3], lwd = 2, type = 'l', col = 'red', 
            xlab = xlabs[i], ylab = ylabs[i], xlim = c(0, max(mobr_group_test[, 2])), main = names[i],
            ylim = c(min(mobr_group_test[, 3:ncol(mobr_group_test)], na.rm = T), max(mobr_group_test[, 3:ncol(mobr_group_test)], na.rm = T)))
      polygon(c(mobr_group_test[, 2], rev(mobr_group_test[, 2])), 
              c(mobr_group_test[, 4], rev(mobr_group_test[, 6])), col = '#C1CDCD', border = NA)
      lines(mobr_group_test[, 2], mobr_group_test[, 3], lwd = 2, type = 'l', col = 'red')
      lines(mobr_group_test[, 2], mobr_group_test[, 5], lwd = 2, type = 'l')
    }
  }
  else{
    ylabs = rep('r', 3)
    for (i in 1:3){
      mobr_group_test = mobr[[type]][[tests[i]]]
      mobr_group_test = mobr_group_test[complete.cases(mobr_group_test), ]
      for (icol in 1:ncol(mobr_group_test))
        mobr_group_test[, icol] = as.numeric(as.character(mobr_group_test[, icol]))
      plot(mobr_group_test[, 1], mobr_group_test[, 2], lwd = 2, type = 'l', col = 'red', 
            xlab = xlabs[i], ylab = ylabs[i], xlim = c(0, max(mobr_group_test[, 1])), main = names[i],
            ylim = c(-1, 1))
      polygon(c(mobr_group_test[, 1], rev(mobr_group_test[, 1])), 
              c(mobr_group_test[, 3], rev(mobr_group_test[, 5])), col = '#C1CDCD', border = NA)
      lines(mobr_group_test[, 1], mobr_group_test[, 2], lwd = 2, type = 'l', col = 'red')
      lines(mobr_group_test[, 1], mobr_group_test[, 4], lwd = 2, type = 'l')
    }
  }
}

summary.mobr = function(...) {
   #  print summary anova style table
}

plot_rarefy = function(mobr){
  # Plot the three curves for an mobr project, 
  # separated by groups
  # Output is a 1*3 figure with the three curves (of each group) separated into subplots
  
  cols = rainbow(ncol(mobr$indiv_rare) - 1)
  par(mfrow = c(1, 3), oma=c(0,0,2,0))
  for (icol in 2:ncol(mobr$indiv_rare)){
    if (icol == 2)
      plot(mobr$indiv_rare$sample, mobr$indiv_rare[, icol], lwd = 2, type = 'l', 
           col = cols[icol - 1], xlab = 'N individuals', ylab = 'Rarefied S',
           main = 'Individual-based Rarefaction', xlim = c(0, max(mobr$indiv_rare$sample)),
           ylim = c(min(mobr$indiv_rare[, -1]), max(mobr$indiv_rare[, -1])))
    else
      lines(mobr$indiv_rare$sample, mobr$indiv_rare[, icol], lwd = 2, col = cols[icol - 1])
  }
  
  groups = unique(mobr$sample_rare$group)
  for (i in 1:length(groups)){
    group = groups[i]
    dat_group = mobr$sample_rare[mobr$sample_rare$group == group, ]
    if (i == 1)
      plot(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$impl_S)), lwd = 2, type = 'l',
           xlab = 'N samples', ylab = 'Rarefied S', col = cols[i], ylim = c(0, mobr$sample_rare$impl_S),
           main = 'Sample-based Rarefaction')
    else
      lines(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$impl_S)), lwd = 2, col = cols[i])
  }
  
  for (i in 1:length(groups)){
    group = groups[i]
    dat_group = mobr$sample_rare[mobr$sample_rare$group == group, ]
    if (i == 1)
      plot(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$expl_S)), lwd = 2, type = 'l',
           xlab = 'N samples', ylab = 'Rarefied S', col = cols[i],ylim = c(0, mobr$sample_rare$expl_S),
           main = 'Accumulation Curve')
    else
      lines(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$expl_S)), lwd = 2, col = cols[i])
  }
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

# Auxillary function: difference between the ind-based rarefaction and the sample-based rarefaction for one group
#   with the evaluation sample size (number of individuals) defined by ref_dens,
#   evaluated at specified points (given by "effort")
# Output: a two-column data frame, with sample size (effort)
#   and deltaS (effect of N)
effect_of_N = function(comm_group, ref_dens, effort){
  group_sad = colSums(comm_group)
  S_samp_rare = as.numeric(rarefaction(comm_group, 'samp'))
  effort = effort[which(effort <= min(sum(comm_group), ref_dens * nrow(comm_group)))]
  interp_S_samp_rare = pchip(c(1, ref_dens * (1:nrow(comm_group))), c(1, S_samp_rare), effort)
  S_indiv_rare = rarefaction(group_sad, 'indiv', effort = effort)
  deltaS = interp_S_samp_rare - as.numeric(S_indiv_rare)
  out = data.frame(cbind(effort, deltaS))
  names(out) = c('effort', 'deltaS')
  return(out)
}

# Auxillary function: spatially-explicit sample-based rarefaction 
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

#' Conduct the MOBR tests on drivers of biodiversity across scales.
#' 
#' There are three tests, on effects of 1. the shape of the SAD, 2. treatment/group-level density,
#' 3. degree of aggregation. The user can specificy to conduct one or more of these tests.
#' 
#'  @param comm "comm" object created by make_comm_obj()
#'  @param env_var a character string specifying the environmental variable in comm$env used
#'  to separate plots into distinct groups. This is the explanatory variable.
#'  @param ref_group one value of env_var, used to define the reference group to which all other groups are compared with
#'  when "type" is discrete. It is not needed when "test" is continuous.
#'  @param tests specifies which one or more of the three tests ('indiv', 'sampl', 'spat') are to be performed. 
#'  Default is to include all three tests.
#'  @param type "discrete" or "continuous". If "discrete", pair-wise comparisons are conducted between all other groups and
#'  the reference group. If "continuous", a correlation analysis is conducted between the response variables and env_var.
#'  @param inds effort size at which the individual-based rarefaction curves are to be evaluated, and to which the sample-based
#'  rarefaction curves are to be interpolated. It can take three types of values, a single integer, a vector of 
#'  intergers, and NULL. If inds = NULL (default), the curves are evaluated at every possible effort size, from 1 to 
#'  the total number of individuals within the group (slow). If inds is a single integer, it is taken as the number 
#'  of points at which the curves are evaluated; the positions of the points are determined by the "log_scale" argument.
#'  If inds is a vector of integers, it is taken as the exact points at which the curves are evaluated.
#'  @param log_scale if "inds" is given a single integer, "log_scale" determines the position of the points. If log_scale is TRUE,
#'  the points are equally spaced on logarithmic scale. If it is FALSE (default), the points are equally spaced on arithmetic scale.
#'  @param min_plot minimal number of plots for test 'spat', where plots are randomized within groups as null test. If it is given
#'  a value, all groups with fewer plots than min_plot are removed for this test. If it is NULL (default), all groups are kept. Warnings
#'  are issued if 1. there is only one group left and "type" is discrete, or 2. there are less than three groups left and "type" is continuous,
#'  or 3. reference group ("ref_group") is removed and "type" is discrete. In these three scenarios, the function will terminate.
#'  A different warning is issued if any of the remaining groups have less than five plots (which have less than 120 permutations), but the 
#'  test will be carried out.
#'  @param density_stat reference density used in converting number of plots to numbers of individuals, a step in test "sampl". It can take
#'  one of the three values: "mean", "max", or "min". If it is "mean", the average plot-level abundance across plots (all plots when "type"
#'  is "continuous, all plots within the two groups for each pair-wise comparison when "type" is "discrete") are used. If it is "min" or "max",
#'  the minimum/maximul plot-level density is used.
#'  @param corr which kind of correlation to use when "type" is "continuous". It can take two values,
#'  "spearman" or "pearson". "spearman" (default) is generally recommended because the relationship between 
#'  the response and "env_var" may not be linear.
#'  @param nperm number of iterations to run for null tests.
#'
#'  @return a "mobr" object with attributes...
#'  @export
#'  @examples
#'  {
#'  library(vegan)
#'  data(mite)
#'  data(mite.env)
#'  data(mite.xy)
#'  mite_comm = make_comm_obj(mite, cbind(mite.env, mite.xy))
#'  mite_comm_discrete = get_delta_stats(mite_comm, 'Shrub', ref_group = 'None', inds = 20)
#'  }

get_delta_stats = function(comm, env_var, group_var=NULL, ref_group=NULL, 
                           tests=c('indiv', 'sampl', 'spat'),
                           type='discrete', inds=NULL, log_scale=FALSE, min_plot = NULL, 
                           density_stat ='mean', corr='spearman', 
                           nperm=1000) {
  # Inputs:
  # comm - a 'comm' type object, with attributes ...
  # type - if the envronmental variable is 'discrete' vs 'continous'
  # env_var - string, name of the envronmental variable
  # group_var - string, optional name of field in comm$env that grouping should
  #   carried out on.
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
    env_data = comm$env[ , env_var]
    if (is.null(group_var)) 
        groups = env_data
    else {
        groups = unique(comm$env[ , group_var])
    }
    if (!(type %in% c('continuous', 'discrete')))
        stop('"type" has to be "discrete" or "continuous".')
    if (type == 'continuous' & !(corr %in% c('spearman', 'pearson')))
      stop('"corr" has to be "spearman" or "pearson".')
    if (type == 'discrete') {
      if (is.null(ref_group))
        stop('For a discrete analysis you must specify a ref_group to compare groups to')
      else if (!(ref_group %in% env_data))
        stop(paste('Reference group is not present in', env_var))
    }
    if ('factor' %in% class(env_data) & type == 'continuous') {
        group_vals = data.frame(groups=groups, 
                                values=as.integer(env_data)[match(groups, env_data)])
        warning(paste(env_var, 'is a factor but will be treated as a continous variable for the analysis which the following values'))
        print(group_vals)
    } else if (!('factor' %in% class(env_data)) & type == 'discrete') 
        warning(paste(env_var, 'is not a factor and each unique value will be treated as a grouping variable'))

    test_status = sapply(tests, function(x) 
                         eval(parse(text=paste('comm$tests$', x, sep=''))))
    approved_tests = tests[which(test_status == TRUE)]
    if (any(test_status == FALSE)) {
        tests_string = paste(tests[which(tests %in% approved_tests)], collapse=' and ')
        cat(paste('Based upon the attributes of the community object only the following tests will be performed:',
                  tests_string))
    }
    out = list()  # This is the object with the outputs
    out$type = type

    group_sad = aggregate(comm$comm, by=list(groups), sum)
    if (is.null(group_var))
        group_levels = group_sad[ , 1]
    else
        group_levels = tapply(env_data, list(groups), mean)
    group_sad = group_sad[ , -1]
    group_minN = min(rowSums(group_sad))
    group_plots = data.frame(table(groups)) # Number of plots within each group
    plot_abd = apply(comm$comm, 1, sum)
    if (density_stat == 'mean')
        ref_dens = sum(comm$comm) / nrow(comm$comm)
    else if (density_stat == 'max')
        ref_dens = max(rowSums(comm$comm))
    else if (density_stat == 'min')
       ref_dens = min(rowSums(comm$comm))
    else 
       stop('The argument ref must be set to mean, min or max')
    
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
    ind_sample_size = unique(c(1, ind_sample_size)) # Force (1, 1) to be included
    
    # 1. Individual-based rarefaction (effect of SAD) vs env_var vs N
    if ('indiv' %in% approved_tests) {
      ind_rare = data.frame(apply(group_sad, 1, function(x) 
                            rarefaction(x, 'indiv', ind_sample_size)))
      row.names(ind_rare) = NULL
      out$indiv_rare = cbind(ind_sample_size, ind_rare)
      names(out$indiv_rare) = c('sample', as.character(group_levels))
      
      if (type == 'continuous'){
        ind_cor = apply(ind_rare, 1, function(x) 
                        cor(x, as.numeric(group_levels), method=corr))
        # Null test 
        sp_extent = unlist(apply(group_sad, 1, function(x) rep(1:ncol(group_sad), x)))
        env_extent = rep(group_levels, times=rowSums(group_sad))
        null_ind_r_mat = matrix(NA, nperm, length(ind_sample_size))
        for (i in 1:nperm){
          overall_sad_lumped = as.numeric(colSums(group_sad))
          meta_freq = SpecDist(overall_sad_lumped)$probability
          sad_perm = sapply(as.numeric(rowSums(group_sad)), function(x)
            data.frame(table(sample(1:length(meta_freq), x, replace = T, prob = meta_freq)))[, 2])
          
          perm_ind_rare = apply(sad_perm, MARGIN = 2, function(x)
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
            level_sad = group_sad[which(as.character(group_levels) == as.character(group)), ]
            comp_sad = rbind(ref_sad, level_sad)

            null_ind_deltaS_mat = matrix(NA, nperm, length(ind_sample_size))
            for (i in 1:nperm){
              comp_sad_lumped = as.numeric(colSums(comp_sad))
              meta_freq = SpecDist(comp_sad_lumped)$probability
              sad_perm = sapply(c(sum(level_sad), sum(ref_sad)), function(x)
                data.frame(table(sample(1:length(meta_freq), x, replace = T, prob = meta_freq)))[, 2])
              perm_ind_rare = sapply(sad_perm, function(x)
                rarefaction(x, 'indiv', ind_sample_size))
              null_ind_deltaS_mat[i, ] = perm_ind_rare[, 1] - perm_ind_rare[, 2]
            }
            ind_deltaS_null_CI = apply(null_ind_deltaS_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
            ind_group = data.frame(cbind(rep(as.character(group), length(ind_sample_size)),ind_sample_size,  
                                         deltaS, t(ind_deltaS_null_CI)))
            out$discrete$indiv = rbind(out$discrete$indiv, ind_group)
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
        effect_N_by_group = data.frame(matrix(NA, ncol=length(group_levels)+1,
                                              nrow=max(group_plots$Freq)))
        for (i in 1:length(group_levels)){
          if (is.null(group_var)) {
              group = group_levels[i]
              comm_group = comm$comm[which(env_data == group), ]
          } else {
              group = names(group_levels)[i]
              comm_group = comm$comm[which(groups == group), ]
          }
          group_effect_N = effect_of_N(comm_group, ref_dens, ind_sample_size)
          if (i == 1)
            effect_N_by_group[, i] = group_effect_N$effort[1:nrow(effect_N_by_group)]
          effect_N_by_group[, i + 1] = group_effect_N$deltaS[1:nrow(effect_N_by_group)]
        }
        effect_N_by_group = effect_N_by_group[complete.cases(effect_N_by_group), ]

        r_emp = apply(effect_N_by_group[ , -1], 1, function(x)
                      cor(x, as.numeric(group_levels), method = corr))
        
        # Null model
        null_N_r_mat = matrix(NA, nperm, length(r_emp))
        for (i in 1:nperm){
          plot_abd_perm = as.numeric(sample(plot_abd))
          ## TODO: clean up this shuffling code so it is easier to read!
          if (is.null(group_var))
              sp_draws = sapply(1:nrow(comm$comm), function(x)
                                sample(rep(1:ncol(comm$comm), 
                                       as.numeric(group_sad[which(group_levels == env_data[x]), ])),
                                       size = plot_abd_perm[x], replace = T))
          else
              sp_draws = sapply(1:nrow(comm$comm), function(x)
                                sample(rep(1:ncol(comm$comm), 
                                       as.numeric(group_sad[which(names(group_levels) == groups[x]), ])),
                                       size = plot_abd_perm[x], replace = T))
          comm_perm = t(sapply(1:nrow(comm$comm), function(x)
                               table(c(1:ncol(comm$comm), sp_draws[[x]])) - 1 ))
          effect_N_perm = data.frame(matrix(NA, ncol=length(group_levels)+1,
                                            nrow=max(group_plots$Freq)))
          for (j in 1:length(group_levels)){
              if (is.null(group_var)) {
                  group = group_levels[j]
                  comm_group = comm_perm[which(env_data == group), ]
              } else {
                  group = names(group_levels)[j]
                  comm_group = comm_perm[which(groups == group), ]
              }
              group_N_perm = effect_of_N(comm_group, ref_dens, ind_sample_size)
              if (j == 1)
                  effect_N_perm[, j] = group_N_perm$effort[1:nrow(effect_N_perm)]
              effect_N_perm[, j + 1] = group_N_perm$deltaS[1:nrow(effect_N_perm)]
          }
          effect_N_perm = effect_N_perm[complete.cases(effect_N_perm), ]
          # If the output is not long enough, fill it with NA's
          null_N_r_mat[i, ] = apply(effect_N_perm[, -1], 1, function(x)
                                    cor(x, as.numeric(group_levels), 
                                        method = corr))[1:ncol(null_N_r_mat)]

        }
        N_r_null_CI = apply(null_N_r_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
        out$continuous$N = data.frame(cbind(effect_N_by_group[, 1], r_emp, t(N_r_null_CI)))
        names(out$continuous$N) = c('effort_ind', 'r_emp', 'r_null_low', 'r_null_median', 'r_null_high')
      }
      else {
        out$discrete$N = data.frame(matrix(0, 0, 6))
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
            
            effect_N_ref = effect_of_N(comm$comm[which(as.character(env_data) == as.character(ref_group)), ], ref_dens, ind_sample_size)
            effect_N_group = effect_of_N(comm$comm[which(env_data == group), ], ref_dens, ind_sample_size)
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
              
              N_ref_perm = effect_of_N(comm_perm[which(as.character(env_2groups) == as.character(ref_group)), ], 
                                       ref_dens, ind_sample_size)
              N_group_perm = effect_of_N(comm_perm[which(as.character(env_2groups) == as.character(group)), ], 
                                         ref_dens, ind_sample_size)
              ddeltaS_perm = N_group_perm$deltaS[1:min(nrow(N_ref_perm), nrow(N_group_perm))] - 
                              N_ref_perm$deltaS[1:min(nrow(N_ref_perm), nrow(N_group_perm))]
              null_N_deltaS_mat[i, ] = ddeltaS_perm[1:ncol(null_N_deltaS_mat)]  
            }
            N_deltaS_null_CI = apply(null_N_deltaS_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
            N_group = data.frame(cbind(rep(as.character(group), length(ddeltaS_group)),effect_N_group$effort[1:length(ddeltaS_group)],  
                                         ddeltaS_group, t(N_deltaS_null_CI)))
            out$discrete$N = rbind(out$discrete$N, N_group)
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
        stop('Error: pair-wise comparison cannot be conducted on one group.')
      else if (length(group_keep) < 3 & type == 'continuous')
        stop('Error: correlation analysis cannot be conducted with less than three groups.')
      else if (!(as.character(ref_group) %in% as.character(group_keep)) & type == 'discrete')
        stop('Error: reference group does not have enough plots and have been dropped.')
      else {
        sample_rare_keep = out$sample_rare[which(out$sample_rare$group %in% as.character(group_keep)), ]
        sample_rare_keep$deltaS = as.numeric(as.character(sample_rare_keep$expl_S)) - 
                                  as.numeric(as.character(sample_rare_keep$impl_S))
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
              deltaS_perm = c(deltaS_perm, as.numeric(expl_S_perm - as.character(sample_rare_keep$impl_S[sample_rare_keep$group == group])))
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
                null_agg_deltaS_mat[i, ] = expl_S_perm_group[1:min_plot_group] - impl_S_group[1:min_plot_group] - 
                  (expl_S_perm_ref[1:min_plot_group] - impl_S_ref[1:min_plot_group])
              }
              agg_deltaS_null_CI = apply(null_agg_deltaS_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
              agg_group = data.frame(cbind(rep(as.character(group), min_plot_group),1:min_plot_group,  
                                           ddeltaS_group, t(agg_deltaS_null_CI)))
              out$discrete$agg = rbind(out$discrete$agg, agg_group)
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

plot_9_panels(ombr)