library(Jade)
library(pracma)
library(RandomFields)

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

print.mobr = function(x) {
    cat('Printing the head of each attribute in the object\n')
    print(lapply(x, head))
}

plot.mobr = function(mobr, group = NULL, par_args=NULL, 
                     same_scale=FALSE) {
  # plot rarefation and delta rarefaction curves
  # Input: 
  # mobr object
  # type: 'discrete' or 'continuous'
  # group: which group to plot. Only required for type = 'discrete' and there are more than one 
  #   pair-wise comparison
  type = mobr$type
  tests = c('indiv', 'N', 'agg')
  names = c('Effect of SAD', 'Effect of N', 'Effect of Aggregation')
  if(!is.null(par_args))
    eval(parse(text=paste('par(', par_args, ')')))
  else
    par(mfrow = c(1, 3))
  xlabs = c('number of individuals', 'number of individuals', 'number of plots')
  if (type == 'discrete'){
    ylabs = c('delta-S', rep('delta-delta-S', 2))
    if (is.null(group) & length(unique(mobr[[type]][[tests[1]]][, 1])) > 1)
      stop("Error: 'group' has to be specified.")
    if(mobr$log_scale) {
      plot_log='x'
      xmin = 1
    }
    else {
      plot_log=''
      xmin = 0
    }
    if (same_scale)
      ylim = range(lapply(mobr$discrete, function(x)
                          lapply(x[ , -(1:2)], function(y)
                                 as.numeric(as.character(y)))))
    for (i in 1:3){
      if (i == 3)
        plot_log=''
      if (is.null(group))
        mobr_group_test = mobr[[type]][[tests[i]]]
      else {
        mobr_group_test = mobr[[type]][[tests[i]]]
        mobr_group_test = mobr_group_test[which(as.character(mobr_group_test$group) == as.character(group)), ]
      }
      mobr_group_test = mobr_group_test[complete.cases(mobr_group_test), ]
      for (icol in 2:ncol(mobr_group_test))
        mobr_group_test[, icol] = as.numeric(as.character(mobr_group_test[, icol]))
      if (!same_scale)
        ylim = c(min(mobr_group_test[, 3:ncol(mobr_group_test)], na.rm = T),
                 max(mobr_group_test[, 3:ncol(mobr_group_test)], na.rm = T))
      plot(mobr_group_test[, 2], mobr_group_test[, 3], lwd = 2, type = 'l', col = 'red', 
            xlab = xlabs[i], ylab = ylabs[i], xlim = c(xmin, max(mobr_group_test[, 2])), main = names[i],
            ylim = ylim, log=plot_log)
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

plot_rarefy = function(mobr, col=NULL){
  # Plot the three curves for an mobr project, 
  # separated by groups
  # Output is a 1*3 figure with the three curves (of each group) separated into subplots
  if (is.null(col[1]))
    col = rainbow(ncol(mobr$indiv_rare) - 1)
  par(mfrow = c(1, 3), oma=c(0,0,2,0))
  groups = unique(mobr$sample_rare$group)
  
  for (i in 1:length(groups)){
    group = groups[i]
    dat_group = mobr$sample_rare[mobr$sample_rare$group == group, ]
    if (i == 1)
      plot(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$expl_S)), lwd = 2, type = 'l',
           xlab = 'N samples', ylab = 'Rarefied S', col = col[i],ylim = c(0, max(as.numeric(as.character(mobr$sample_rare$expl_S)))),
           main = 'Accumulation Curve')
    else
      lines(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$expl_S)), lwd = 2, col = col[i])
  }
  
  for (i in 1:length(groups)){
    group = groups[i]
    dat_group = mobr$sample_rare[mobr$sample_rare$group == group, ]
    if (i == 1)
      plot(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$impl_S)), lwd = 2, type = 'l',
           xlab = 'N samples', ylab = 'Rarefied S', col = col[i], ylim = c(0, max(as.numeric(as.character(mobr$sample_rare$impl_S)))),
           main = 'Sample-based Rarefaction')
    else
      lines(as.numeric(as.character(dat_group$sample_plot)), as.numeric(as.character(dat_group$impl_S)), lwd = 2, col = col[i])
  }
  
  for (icol in 2:ncol(mobr$indiv_rare)){
    if (icol == 2)
      plot(mobr$indiv_rare$sample, mobr$indiv_rare[, icol], lwd = 2, type = 'l', 
           col = col[icol - 1], xlab = 'N individuals', ylab = 'Rarefied S',
           main = 'Individual-based Rarefaction', xlim = c(0, max(mobr$indiv_rare$sample)),
           ylim = c(min(mobr$indiv_rare[, -1]), max(mobr$indiv_rare[, -1])))
    else
      lines(mobr$indiv_rare$sample, mobr$indiv_rare[, icol], lwd = 2, col = col[icol - 1])
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
    #class(out) = 'rarefaction'
    return(out)
}

# Auxillary function: difference between the ind-based rarefaction and the
# sample-based rarefaction for one group with the evaluation sample size (number
# of individuals) defined by ref_dens, evaluated at specified points (given by
# inds) Output: a two-column data frame, with sample size
# (effort) and deltaS (effect of N)
effect_of_N = function(comm_group, ref_dens, inds){
  nplots = nrow(comm_group)
  group_dens = sum(comm_group) / nplots
  # calcualte the number of individuals sampled
  # as each of the nplots are collected
  samp_effort = round((1:nplots) * group_dens)
  #effort = effort[which(effort <= min(sum(comm_group), ref_dens * nrow(comm_group)))]
  # use individual based rarefaction evaluated at the 
  # number of individuals sampled for each plot
  S_samp = rarefaction(comm_group, 'indiv', samp_effort)
  # rescale and interpolate this plot based S to individual based 
  # using the ref_density (i.e., not the observed density)
  rescaled_effort = round(1:nplots * ref_dens)
  if (max(rescaled_effort) > max(inds)) 
    warning('Extrapolating the rarefaction curve because the number of rescaled individuals is larger than the inds argument')
  interp_S_samp = pchip(c(1, rescaled_effort),
                        c(1, S_samp), inds)
  S_indiv = rarefaction(comm_group, 'indiv', inds)
  deltaS = interp_S_samp - S_indiv
  out = data.frame(inds=inds, deltaS=deltaS)
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

#' Permute community matrix within groups
#' 
#' Two types of permutation can be carried out: 
#' 1) each individual of each species is reassigned a plot randomly which removes 
#' any patterns due to within and between plot spatial aggregation, but maintains
#' species group abundance and therefore observed group N,
#' 2) the total number of individuals in a plot is shuffled and then that many
#' individuals are drawn randomly from the group specific species-abundance
#' distribution for each plot which provides a mean of removing group differenes
#' in the total number of individuals. 
#' 
#' @param comm a site-by-species matrix
#' @param swap either 'noagg' or 'swapN' 
#' @param groups optional argument that is a vector of group ids which specify
#'   which group each site is associated with. If is NULL then all rows of the
#'   community matrix are assumed to be members of the same group
#'   
#' @return a permuted site-by-species matrix
#' 
#' @examples 
#' S = 3
#' N = 20
#' nplots = 4
#' comm = matrix(rpois(S*nplots, 1),
#'               ncol=S, nrow=nplots)
#' comm
#' groups = rep(1:2, each=2)
#' groups
#' permute_comm(comm, 'noagg')
#' permute_comm(comm, 'noagg', groups)
#' permute_comm(comm, 'swapN')
#' permute_comm(comm, 'swapN', groups)
permute_comm = function(comm, swap, groups=NULL) {
    if (!(is.matrix(comm) | is.data.frame(comm)))
        stop('comm must be a matrix or data.frame')
    if (is.null(groups))
        groups = rep(1, nrow(comm)) 
    group_levels = unique(groups)
    S = ncol(comm)
    comm_group_perm = matrix(NA, ncol=S, nrow=nrow(comm))
    if (swap == 'swapN')
        Nperm = sample(rowSums(comm))
    for(i in seq_along(group_levels)) {
        row_indices = groups == group_levels[i]
        group_comm = comm[row_indices, ]
        sp_abu = colSums(group_comm)
        plot_ids = 1:nrow(group_comm)
        if (swap == 'noagg') {
            tmp_comm = sapply(sp_abu, function(x) 
                              table(c(sample(plot_ids, x,
                                             replace=T),
                              plot_ids)) - 1)
        } 
        else if (swap == 'swapN') {
            sp_draws = sapply(plot_ids, function(x)
                              sample(rep(1:S, sp_abu),
                                     size=Nperm[x],
                                     replace = T))
            tmp_comm = t(sapply(plot_ids, function(x)
                         table(c(sp_draws[[x]], 1:S)) - 1 ))
        }
        else 
            stop('The argument swap must be either "noagg" or "swapN"')
        comm_group_perm[row_indices, ] = tmp_comm
    }  
    return(comm_group_perm)
}
  

avg_perm_rare = function(comm, swap, groups=NULL, nperm=1000, effort=NULL){
    S = replicate(nperm, 
                  rarefaction(permute_comm(comm, swap, groups),
                              'samp', effort))
    Savg = apply(S, 1, mean)
    return(Savg)
}

swap_binary_species = function(comm, groups){
  ###ToDO incoperate into permute_comm() function if needed
  # This function converts the plot by sp matrix into binary,
  #   then swap the presences among the plots.
  #   In this way the (overall, across-all-plots) intraspecific 
  #   aggregation pattern is maintained, and equalized among the 
  #   treatments
  # comm is the plot by species matrix
  # groups is the grouping factor, the same length as nrow(comm)
  comm_binary = (comm > 0) * 1
  pa_group = aggregate(comm_binary, by = list(groups), sum)
  pa_group = (pa_group[, -1] > 0)
  comm_out = matrix(nrow = nrow(comm_binary), ncol = ncol(comm_binary))
  
  for (sp in 1:ncol(comm_binary)){
    sp_swap = sample(comm_binary[, sp])
    sp_pa_group = aggregate(sp_swap, by = list(groups), sum)
    while(any((sp_pa_group[, 2] > 0) != pa_group[, sp])){
      sp_swap = sample(comm_binary[, sp])
      sp_pa_group = aggregate(sp_swap, by = list(groups), sum)
    }
    comm_out[, sp] = sp_swap
  }
  return(comm_out)
}

# Attempt to maintain some spatial autcorrelation to improve
# type 1 error of spatial null model
# Note: function not complete 
samp_ssad = function(comm, groups){
  ords = apply(aggregate(comm, list(groups), sum)[ , -1], 1,
               order, decreasing=TRUE)
  group_levels = unique(groups)
  comm_rank = comm
  #sapply(seq_along(group_levels), function(x)
  #        comm[groups == group_levels[x], ords[ , x]],
  #       simplify='array')
  for(i in seq_along(group_levels)) {
    row_bool = groups == group_levels[i]
    comm_rank[row_bool, ] = comm[row_bool, ords[ , i]]
  }
  group_sad = aggregate(comm_rank, list(groups), sum)[ , -1]
  comm_perm = comm_rank
  for (sp in 1:ncol(comm_rank)){
    coin = ifelse(runif(1) < 0.5, 1, 2)
    row_bool = groups == group_levels[coin]
    if (all(group_sad[ , sp] > 0)) {
#        comm_perm[row_bool, sp] = sample(
    }
 
  }
  return(comm_out)
}



#' Conduct the MOBR tests on drivers of biodiversity across scales.
#' 
#' There are three tests, on effects of 1. the shape of the SAD, 2.
#' treatment/group-level density, 3. degree of aggregation. The user can
#' specificy to conduct one or more of these tests.
#' 
#' @param comm "comm" object created by make_comm_obj()
#' @param env_var a character string specifying the environmental variable to in
#'   comm$env used for a comparison
#' @param group_var an optional character string to specify group membership. If
#'   it is NULL then each unique value of env_var is used as the group variable
#' @param ref_group a character string used to define the reference group to
#'   which all other groups are compared with when "type" is discrete. It is not
#'   needed when "type" is continuous.
#' @param tests specifies which one or more of the three tests ('indiv',
#'   'sampl', 'spat') are to be performed. Default is to include all three
#'   tests.
#' @param type "discrete" or "continuous". If "discrete", pair-wise comparisons
#'   are conducted between all other groups and the reference group. If
#'   "continuous", a correlation analysis is conducted between the response
#'   variables and env_var.
#' @param inds effort size at which the individual-based rarefaction curves are
#'   to be evaluated, and to which the sample-based rarefaction curves are to be
#'   interpolated. It can take three types of values, a single integer, a vector
#'   of intergers, and NULL. If inds = NULL (default), the curves are evaluated
#'   at every possible effort size, from 1 to the total number of individuals
#'   within the group (slow). If inds is a single integer, it is taken as the
#'   number of points at which the curves are evaluated; the positions of the
#'   points are determined by the "log_scale" argument. If inds is a vector of
#'   integers, it is taken as the exact points at which the curves are
#'   evaluated.
#' @param log_scale if "inds" is given a single integer, "log_scale" determines
#'   the position of the points. If log_scale is TRUE, the points are equally
#'   spaced on logarithmic scale. If it is FALSE (default), the points are
#'   equally spaced on arithmetic scale.
#' @param min_plot minimal number of plots for test 'spat', where plots are
#'   randomized within groups as null test. If it is given a value, all groups
#'   with fewer plots than min_plot are removed for this test. If it is NULL
#'   (default), all groups are kept. Warnings are issued if 1. there is only one
#'   group left and "type" is discrete, or 2. there are less than three groups
#'   left and "type" is continuous, or 3. reference group ("ref_group") is
#'   removed and "type" is discrete. In these three scenarios, the function will
#'   terminate. A different warning is issued if any of the remaining groups
#'   have less than five plots (which have less than 120 permutations), but the 
#'   test will be carried out.
#' @param density_stat reference density used in converting number of plots to
#'   numbers of individuals, a step in test "sampl". It can take one of the
#'   three values: "mean", "max", or "min". If it is "mean", the average
#'   plot-level abundance across plots (all plots when "type" is "continuous,
#'   all plots within the two groups for each pair-wise comparison when "type"
#'   is "discrete") are used. If it is "min" or "max", the minimum/maximul
#'   plot-level density is used.
#' @param corr which kind of correlation to use when "type" is "continuous". It
#'   can take two values, "spearman" or "pearson". "spearman" (default) is
#'   generally recommended because the relationship between the response and
#'   "env_var" may not be linear.
#' @param nperm number of iterations to run for null tests.
#'   
#' @return a "mobr" object with attributes...
#' @export
#'  @examples
#'  {
#'  library(vegan)
#'  data(mite)
#'  data(mite.env)
#'  data(mite.xy)
#'  mite_comm = make_comm_obj(mite, data.frame(mite.env, mite.xy))
#'  mite_comm_discrete = get_delta_stats(mite_comm, 'Shrub',
#'                                       ref_group = 'None', inds = 20)
#'  }

get_delta_stats = function(comm, env_var, group_var=NULL, ref_group=NULL, 
                           tests=c('indiv', 'sampl', 'spat'),
                           type='discrete', inds=NULL, log_scale=FALSE,
                           min_plot = NULL, density_stat ='mean',
                           corr='spearman', nperm=1000) {
    env_data = comm$env[ , env_var]
    if (is.null(group_var)) 
        groups = as.character(env_data)
    else {
        groups = as.character(unique(comm$env[ , group_var]))
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
        tests_string = paste(tests[which(tests %in% approved_tests)],
                             collapse=' and ')
        cat(paste('Based upon the attributes of the community object only the following tests will be performed:',
                  tests_string))
    }
    out = list()  
    out$type = type
    out$log_scale = log_scale
    S = ncol(comm$comm)
    group_sad = aggregate(comm$comm, by=list(groups), sum)
    if (is.null(group_var))
        group_levels = group_sad[ , 1]
    else
        group_levels = tapply(env_data, list(groups), mean)
    group_sad = group_sad[ , -1]
    group_minN = min(rowSums(group_sad))
    group_plots = data.frame(table(groups)) # Number of plots within each group
    plot_abd = rowSums(comm$comm)
    if (density_stat == 'mean')
        indiv_dens = sum(comm$comm) / nrow(comm$comm)
    else if (density_stat == 'max')
        indiv_dens = max(rowSums(comm$comm))
    else if (density_stat == 'min')
        indiv_dens = min(rowSums(comm$comm))
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
                            data.frame(table(sample(1:length(meta_freq),
                                                    x, replace = T,
                                                    prob = meta_freq)))[, 2])
          perm_ind_rare = apply(sad_perm, MARGIN = 2, function(x)
                                rarefaction(x, 'indiv', ind_sample_size))
          null_ind_r_mat[i, ] = apply(perm_ind_rare, 1, function(x){
                                      cor(x, as.numeric(group_levels),
                                          method = corr)})
        }
        ind_r_null_CI = apply(null_ind_r_mat, 2, function(x)
                              quantile(x, c(0.025, 0.5, 0.975))) # 95% CI
        out$continuous$indiv = data.frame(cbind(ind_sample_size, ind_cor,
                                                t(ind_r_null_CI)))
        names(out$continuous$indiv) = c('effort_ind', 'r_emp', 'r_null_low',
                                        'r_null_median', 'r_null_high')
      }
      else { # discrete case
        ref_sad = group_sad[which(as.character(group_levels) ==
                                  as.character(ref_group)), ]
        out$discrete$indiv = data.frame(sample = numeric(), group = character(),
                                        deltaS_emp = numeric(),
                                        deltaS_null_low = numeric(), 
                                        deltaS_null_median = numeric(),
                                        deltaS_null_high = numeric(),
                                        stringsAsFactors = F)
        for (group in group_levels){
          if (as.character(group) != as.character(ref_group)){
            deltaS = out$indiv_rare[, as.character(group)] - out$indiv_rare[, as.character(ref_group)]
            level_sad = group_sad[which(as.character(group_levels) == as.character(group)), ]
            comp_sad = rbind(ref_sad, level_sad)

            null_ind_deltaS_mat = matrix(NA, nperm, length(ind_sample_size))
            cat('\nComputing null model for SAD effect\n')
            pb <- txtProgressBar(min = 0, max = nperm, style = 3)
            for (i in 1:nperm){
              setTxtProgressBar(pb, i)
              comp_sad_lumped = as.numeric(colSums(comp_sad))
              meta_freq = SpecDist(comp_sad_lumped)$probability
              sad_perm = lapply(c(sum(level_sad), sum(ref_sad)), function(x)
                data.frame(table(sample(1:length(meta_freq), x, replace = T,
                                        prob = meta_freq)))[, 2])
              perm_ind_rare = sapply(sad_perm, function(x)
                                    rarefaction(x, 'indiv', ind_sample_size))
              null_ind_deltaS_mat[i, ] = perm_ind_rare[, 1] - perm_ind_rare[, 2]
              setTxtProgressBar(pb, i)
            }
            close(pb)
            ind_deltaS_null_CI = apply(null_ind_deltaS_mat, 2, function(x)
                                       quantile(x, c(0.025, 0.5, 0.975)))
            ind_group = data.frame(cbind(rep(as.character(group),
                                             length(ind_sample_size)),
                                         ind_sample_size,  
                                         deltaS, t(ind_deltaS_null_CI)))
            out$discrete$indiv = rbind(out$discrete$indiv, ind_group)
          }
        }
        names(out$discrete$indiv) = c('group', 'effort_ind', 'deltaS_emp',
                                      'deltaS_null_low', 'deltaS_null_median',
                                      'deltaS_null_high')
      }
    }
    
    # Sample-based spatially-implicit and -explicit rarefaction -----------
    if ('sampl' %in% approved_tests | 'spat' %in% approved_tests){
      if ('spat' %in% approved_tests)
        out$sample_rare = data.frame(group = character(), sample_plot = numeric(), impl_S = numeric(), expl_S = numeric())
      else
        out$sample_rare = data.frame(group = character(), sample_plot = numeric(), impl_S = numeric())
      for (group in group_levels){
        comm_group = comm$comm[as.character(env_data) == as.character(group), ]
        #impl_S = as.numeric(rarefaction(comm_group, 'samp'))
        # replace sample-based rarefaction curve with the new curve
        #nplots = nrow(comm_group)
        #group_dens = sum(comm_group) / nplots
        #samp_effort = round((1:nplots) * group_dens)
        #S_samp = rarefaction(comm_group, 'indiv', samp_effort)
        #rescaled_effort = round(1:nplots * ref_dens)
        #interp_S_samp = pchip(c(1, rescaled_effort),
        #                      c(1, S_samp), ind_sample_size)
        impl_S = avg_perm_rare(comm_group, 'noagg', nperm = 100)
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
    
    # 2. Sample-based rarefaction (effect of density) vs env_var vs N----------
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
          group_effect_N = effect_of_N(comm_group, indiv_dens, ind_sample_size, 
                                       permute=T)
          if (i == 1)
            effect_N_by_group[, i] = group_effect_N$effort[1:nrow(effect_N_by_group)]
          effect_N_by_group[, i + 1] = group_effect_N$deltaS[1:nrow(effect_N_by_group)]
        }
        effect_N_by_group = effect_N_by_group[complete.cases(effect_N_by_group), ]

        r_emp = apply(effect_N_by_group[ , -1], 1, function(x)
                      cor(x, as.numeric(group_levels), method = corr))
        
        # Null model -------
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
              group_N_perm = effect_of_N(comm_group, indiv_dens, ind_sample_size)
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
      else if (type == 'discrete') {
          for (group in group_levels) {
              if (group != ref_group){
                  row_bool = group == groups | ref_group == groups
                  comm_groups = comm$comm[row_bool, ]
                  env_groups = env_data[row_bool]
              if (density_stat == 'mean')
                  ref_dens = sum(comm_groups) / nrow(comm_groups)
              else if (density_stat == 'max')
                  ref_dens = max(rowSums(comm_groups))
              else if (density_stat == 'min')
                  ref_dens = min(rowSums(comm_groups))
              else 
                  stop('The argument ref must be set to mean, min or max')
              
              N_eff = sapply(c(ref_group, group), function(x) 
                             effect_of_N(comm_groups[env_groups == x, ], ref_dens,
                                         ind_sample_size)$deltaS)
              ddeltaS_group = N_eff[ , 2] - N_eff[ , 1]
              null_N_deltaS_mat = matrix(NA, nperm, length(ddeltaS_group))
              cat('\nComputing null model for N effect\n')
              pb <- txtProgressBar(min = 0, max = nperm, style = 3)
              for (i in 1:nperm){
                  setTxtProgressBar(pb, i)
                  # swap plot abu between group 1 and each other group
                  comm_perm = permute_comm(comm_groups, 'swapN', env_groups)  
                  N_eff_perm = sapply(c(ref_group, group), function(x) 
                                      effect_of_N(comm_perm[env_groups == x, ],
                                                  ref_dens, ind_sample_size)$deltaS)
                  null_N_deltaS_mat[i, ] = N_eff_perm[ , 2] - N_eff_perm[ , 1]
              }
              close(pb)
              N_deltaS_null_CI = apply(null_N_deltaS_mat, 2, function(x)
                                       quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
              N_group = data.frame(group, ind_sample_size, ddeltaS_group,
                                   t(N_deltaS_null_CI))
              out$discrete$N = rbind(out$discrete$N, N_group)
          }
        }
        names(out$discrete$N) = c('group', 'effort_sample', 'ddeltaS_emp', 'ddeltaS_null_low', 
                                  'ddeltaS_null_median', 'ddeltaS_null_high')
      }
    }
    # 3. Sample-based spatially-explicit rarefaction (effect of aggregation) vs env_var vs N
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
            deltaS_perm = c()
            comm_perm = swap_binary_species(comm$comm, groups)
            for (group in unique(sample_rare_keep$group)){
              comm_group = comm_perm[as.character(env_data) == as.character(group), ]
              xy_group = comm$spat[as.character(env_data) == as.character(group), ]
              expl_S_perm = rarefy_sample_explicit(comm_group, xy_group)
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
          xy_ref = comm$spat[as.character(env_data) == as.character(ref_group), ]
          impl_S_ref = sample_rare_keep$impl_S[which(as.character(sample_rare_keep$group) == as.character(ref_group))]
          impl_S_ref = as.numeric(as.character(impl_S_ref))
          for (group in unique(group_keep)){
            if ((as.character(group) != as.character(ref_group))){
              min_plot_group = min(group_plots$Freq[which(as.character(group_plots[, 1]) %in% 
                                                            c(as.character(group), as.character(ref_group)))])
              ddeltaS_group = sample_rare_keep$deltaS[sample_rare_keep$group == group][1:min_plot_group] - 
                sample_rare_keep$deltaS[as.character(sample_rare_keep$group) == as.character(ref_group)][1:min_plot_group]
              
              group_for_2 = env_data[which(as.character(env_data) %in% c(as.character(ref_group), as.character(group)))]
              comm_group = comm$comm[as.character(env_data) == as.character(group), ]
              impl_S_group = sample_rare_keep$impl_S[which(as.character(sample_rare_keep$group) == as.character(group))]
              impl_S_group = as.numeric(as.character(impl_S_group))
              
              null_agg_deltaS_mat = matrix(NA, nperm, min_plot_group)
              cat('\nComputing null model for aggregation effect\n')
              pb <- txtProgressBar(min = 0, max = nperm, style = 3)
              for (i in 1:nperm){
                setTxtProgressBar(pb, i)
                xy_perm = comm$spat[sample(nrow(comm$spat)), ]
                xy_perm_group = xy_perm[as.character(env_data) == as.character(group), ]
                xy_perm_ref = xy_perm[as.character(env_data) == as.character(ref_group), ]
                expl_S_perm_group = rarefy_sample_explicit(comm_group, xy_perm_group)
                expl_S_perm_ref = rarefy_sample_explicit(ref_comm, xy_perm_ref)
                null_agg_deltaS_mat[i, ] = expl_S_perm_group[1:min_plot_group] - impl_S_group[1:min_plot_group] - 
                  (expl_S_perm_ref[1:min_plot_group] - impl_S_ref[1:min_plot_group])
              }
              close(pb)
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

plot_grp_rads = function(comm_obj, env_var, col=NA,
                         log='') {
  #plot group pooled rank species abundance distribution
  par(mfrow = c(1, 1))
  env_data = comm_obj$env[ , env_var]
  grps = unique(env_data)
  if (is.na(col[1])) 
    col = rainbow(length(grps))
  sads = aggregate(comm_obj$comm, by=list(comm_obj$env[ , env_var]), 
                   sum)
  grps = as.character(sads[,1])
  sads = as.matrix(sads[,-1])
  sads = ifelse(sads == 0, NA, sads)
  plot(1:10, 1:10, type='n', 
       xlab='rank', ylab='abundance',
       log=log, xlim=c(1, ncol(comm_obj$comm)), 
       ylim=range(sads, na.rm=T), 
       cex.lab = 1.5, cex.axis = 1.5)
  for(i in 1:nrow(sads)) 
    lines(1:sum(!is.na(sads[i, ])), sort(sads[i, ], dec=T),
          col=col[i], lwd=2)
  legend('topright', grps, col=col, bty='n', lty=1,
         lwd=3, cex=2)
}


plotSADs = function(comm_obj, env_var, col = NA) {
  # TO DO: add check to ensure that col is the same length as treatments
  require(scales)
  par(mfrow = c(1, 1))
  env_data = comm_obj$env[ , env_var]
  grps = unique(env_data)
  if (is.na(col[1])) 
    col = rainbow(length(grps))
  plot(1, type = "n", xlab = "% abundance (log scale)", ylab = "% species", 
       xlim = c(0.01, 1), ylim = c(0, 1), log = "x")
  for (i in 1:length(grps)) {
    col_grp = col[i]
    comm_grp = comm_obj$comm[env_data == grps[i], ]
    for (j in 1:nrow(comm_grp)) {
      sad_row = as.numeric(sort(comm_grp[j, comm_grp[j, ] != 0]))
      s_cul = 1:length(sad_row)/length(sad_row)
      n_cul = sapply(1:length(sad_row), function(x) sum(sad_row[1:x]) / sum(sad_row))
      lines(n_cul, s_cul, col = alpha(col_grp, 0.5), lwd = 1, type = "l")
    }
  }
  legend("topleft", legend=grps, col = col, lwd = 2, bty='n')
}

plotSNpie = function(comm_obj, env_var, col = NA) {
  # TO DO: add check to ensure that col is the same length as treatments
  require(rgl)
  env_data = comm_obj$env[ , env_var]
  grps = unique(env_data)
  if (is.na(col[1])) 
    col = rainbow(length(grps))
  S_list = rowSums(comm_obj$comm > 0)
  N_list = rowSums(comm_obj$comm)
  PIE_list = sapply(1:nrow(comm$comm), function(x) 
    N_list[x]/(N_list[x] - 1) * (1 - sum((comm_obj$comm[x, ]/N_list[x])^2)))
  col_list = sapply(env_data, function(x) col[which(grps == x)])
  plot3d(S_list, N_list, PIE_list, "S", "N", "PIE", col = col_list, size = 8)
} 

plot_9_panels = function(mobr, trt_group, ref_group,
                         par_args=NULL, same_scale=FALSE){
  type = mobr$type
  if (type == 'continuous')
    stop("Currently this plot only works for mobr object with type discrete.")
  else{
    cols = c('red', 'blue')
    deltaS_col = 'turquoise'
    ddeltaS_col = 'magenta'
    if (!is.null(par_args))
      eval(parse(text=paste('par(', par_args, ')')))
    else 
      par(mfrow = c(3, 3))
    if (same_scale)
      ylim = range(lapply(mobr$discrete, function(x)
                          lapply(x[ , -(1:2)], function(y)
                                 as.numeric(as.character(y)))))
    if(mobr$log_scale) {
      plot_log = 'x'
      xmin = 1
    }
    else {
      plot_log = ''
      xmin = 0
    }
    # Create the three sets of curves
    plot(mobr$indiv_rare$sample, mobr$indiv_rare[[trt_group]], xlab = 'Number of individuals', 
         ylab = 'Richness(S)', xlim = c(xmin, max(mobr$indiv_rare$sample)), 
         ylim = c(0, max(mobr$indiv_rare[, -1])), type = 'l', lwd = 2, 
         col = cols[1], cex.lab = 1.5, cex.axis = 1.5, main = 'Individual', cex.main = 2,
         log=plot_log)
    lines(mobr$indiv_rare$sample, mobr$indiv_rare[[ref_group]], 
          type = 'l', lwd = 2, col = cols[2])

    mobr$sample_rare[, -1] = lapply(mobr$sample_rare[, -1], function(x)
                                    as.numeric(as.character(x)))
    sample_rare_group = mobr$sample_rare[mobr$sample_rare == trt_group, ]
    sample_rare_ref = mobr$sample_rare[mobr$sample_rare == ref_group, ]
    
    plot(1:nrow(sample_rare_group), sample_rare_group$impl_S, 
         xlab = 'Number of plots', ylab = 'Richness (S)', 
         xlim = c(xmin, max(nrow(sample_rare_ref), nrow(sample_rare_group))),
         ylim = c(min(sample_rare_ref$impl_S, sample_rare_group$impl_S), 
                  max(sample_rare_ref$impl_S, sample_rare_group$impl_S)), 
         type = 'l', lwd = 2, col = cols[1], cex.lab = 1.5, cex.axis = 1.5, 
         main = 'Sample', cex.main = 2)
    lines(1:nrow(sample_rare_ref), sample_rare_ref$impl_S, type = 'l', 
          lwd = 2, col = cols[2])
    
    plot(1:nrow(sample_rare_group), sample_rare_group$expl_S, 
         xlab = 'Number of plots', ylab='Richness (S)',
         xlim = c(xmin, max(nrow(sample_rare_ref), nrow(sample_rare_group))),
         ylim = c(min(sample_rare_ref$expl_S, sample_rare_group$expl_S), 
                  max(sample_rare_ref$expl_S, sample_rare_group$expl_S)), 
         type = 'l', lwd = 2, col = cols[1], cex.lab = 1.5, cex.axis = 1.5, 
         main = 'Spatial', cex.main = 2)
    lines(1:nrow(sample_rare_ref), sample_rare_ref$expl_S, type = 'l', 
          lwd = 2, col = cols[2])
    
    # Create the plots for the three delta-S between groups
    deltaS_Sind = mobr$indiv_rare[[trt_group]] - mobr$indiv_rare[[ref_group]]
    plot(mobr$indiv_rare$sample, deltaS_Sind, ylim = c(min(deltaS_Sind, 0), max(deltaS_Sind, 0)),
         cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = 2, col = deltaS_col,
         xlab = 'Number of individuals', ylab = 'delta S',
         log=plot_log)
    abline(h = 0, lwd = 2, lty = 2)
    
    minN = min(nrow(sample_rare_group), nrow(sample_rare_ref))
    delta_Ssample = sample_rare_group$impl_S[1:minN] - sample_rare_ref$impl_S[1:minN]
    plot(seq(minN), delta_Ssample, ylim = c(min(delta_Ssample, 0), max(delta_Ssample, 0)),
         cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = 2, col = deltaS_col,
         xlab = 'Number of plots', ylab = 'delta S' )
    abline(h = 0, lwd = 2, lty = 2)

    delta_Sspat = sample_rare_group$expl_S[1:minN] - sample_rare_ref$expl_S[1:minN]
    plot(seq(minN), delta_Sspat, ylim = c(min(delta_Sspat, 0), max(delta_Sspat, 0)),
         cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = 2, col = deltaS_col,
         xlab = 'Number of plots', ylab = 'delta S' )
    abline(h = 0, lwd = 2, lty = 2)
    
    # Create the plots for the three d-delta S
    mobr$discrete$ind[, -1] = lapply(mobr$discrete$ind[, -1], function(x)
      as.numeric(as.character(x))) 
    delta_Sind = mobr$discrete$ind[which(as.character(mobr$discrete$ind$group) == as.character(trt_group)), ]
    if (!same_scale)
      ylim = range(delta_Sind[ , -(1:2)])
    plot(delta_Sind$effort_ind, delta_Sind$deltaS_emp, 
         ylim = ylim, log=plot_log,
         cex.axis = 1.5, cex.lab = 1.5, type = 'n',
         xlab = 'Number of individuals', ylab = 'delta S')
    polygon(c(delta_Sind$effort_ind, rev(delta_Sind$effort_ind)), 
            c(delta_Sind$deltaS_null_low, rev(delta_Sind$deltaS_null_high)),
            col = '#C1CDCD', border = NA)
    abline(h = 0, lwd = 2, lty = 2)
    lines(delta_Sind$effort_ind, delta_Sind$deltaS_emp,
          lwd = 2, col = ddeltaS_col)
        
    mobr$discrete$N[, -1] = lapply(mobr$discrete$N[, -1], function(x)
      as.numeric(as.character(x))) 
    ddelta_Ssample = mobr$discrete$N[which(as.character(mobr$discrete$N$group) == as.character(trt_group)), ]
    if (!same_scale)
      ylim = range(ddelta_Ssample[ , -(1:2)])
    plot(ddelta_Ssample$effort_sample, ddelta_Ssample$ddeltaS_emp,
         ylim = ylim, log=plot_log,
         cex.axis = 1.5, cex.lab = 1.5, type = 'n', 
         xlab = 'Number of individuals', ylab = 'delta-delta S')
    polygon(c(ddelta_Ssample$effort_sample, rev(ddelta_Ssample$effort_sample)), 
            c(ddelta_Ssample$ddeltaS_null_low, rev(ddelta_Ssample$ddeltaS_null_high)),
            col = '#C1CDCD', border = NA)
    abline(h = 0, lwd = 2, lty = 2)
    lines(ddelta_Ssample$effort_sample, ddelta_Ssample$ddeltaS_emp,
          lwd = 2, col = ddeltaS_col)
    
    mobr$discrete$agg[, -1] = lapply(mobr$discrete$agg[, -1], function(x)
      as.numeric(as.character(x))) 
    ddelta_Sspat = mobr$discrete$agg[which(as.character(mobr$discrete$agg$group) == as.character(trt_group)), ]
    if (!same_scale)
      ylim = range(ddelta_Sspat[ , -(1:2)])
    plot(ddelta_Sspat$effort_sample, ddelta_Sspat$ddeltaS_emp,
         ylim = ylim, log='',
         cex.axis = 1.5, cex.lab = 1.5, type = 'n', 
         xlab = 'Number of plots', ylab = 'delta-delta S')
    polygon(c(ddelta_Sspat$effort_sample, rev(ddelta_Sspat$effort_sample)), 
            c(ddelta_Sspat$ddeltaS_null_low, rev(ddelta_Sspat$ddeltaS_null_high)),
            col = '#C1CDCD', border = NA)
    abline(h = 0, lwd = 2, lty = 2)
    lines(ddelta_Sspat$effort_sample, ddelta_Sspat$ddeltaS_emp, 
          lwd = 2, col = ddeltaS_col)
  }
}
