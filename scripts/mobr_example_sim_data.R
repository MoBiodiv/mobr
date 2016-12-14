# Sensitivity test for Mobr functions
# using Felix May's MoBspatial package
.libPaths(c(.libPaths(), "C:/Users/Xiao/Documents/R/win-library/3.2"))

library(MoBspatial)
source('C:\\Users\\Xiao\\Documents\\GitHub\\mobr\\R\\mobr.R')
#library(mobr) # This doesn't work at the moment; will have to wait until the mobr package is built


# Simulate a community using Sim.Thomas.Community() and wrangle it into a list
#   which can later be combined to create the comm object for mobr.
# Inputs:
# S, N, cv, sigma - parameters to simulate the community. xmax and ymax are fixed at 1.
# Note that now S is S_pool (total richness for the metacommunity), while 
# N is N_local (number of individuals for the local group)
# Because of the sampling process, the local group may not have all S species.
# sqrt_numplots - sqrt(number of plots) within the domain. This is how the domain will be divided.
# Ouput:
# A list with three components:
# $comm - plot by species matrix
# $coords - coordinates of the plots
# $params - parameters used to create the simulation, including S, N, cv (shape of the SAD),
#   and sigma (spatial aggregation)
sim_comm_single_pars = function(S, N, cv, sigma, sqrt_numplots){
  sim_comm = Sim.Thomas.Community(S, N, cv, sigma)
  out = list()
  out$comm = matrix(0, sqrt_numplots^2, S)
  names(out$comm) = sapply(1:S, function(x) paste('species', as.character(x), sep = ''))
  out$coords = matrix(NA, sqrt_numplots^2, 2)
  x_span = 1 / sqrt_numplots 
  y_span = 1 / sqrt_numplots 
  
  nrow = 1
  for (i in 1:sqrt_numplots){
    xlim = c(i-1, i) * x_span
    for (j in 1:sqrt_numplots){
      ylim = c(j-1, j) * y_span
      cond = which(sim_comm$X >= xlim[1] & sim_comm$X <= xlim[2] &
                     sim_comm$Y >= ylim[1] & sim_comm$Y <= ylim[2])
      comm_plot = sim_comm[cond, ]
      plot_sad = data.frame(table(comm_plot$SpecID))
      for (k in 1:nrow(plot_sad)){
        out$comm[nrow, which(names(out$comm) == as.character(plot_sad[k, 1]))] = 
          plot_sad[k, 2]
      }
      out$coords[nrow, ] = c(mean(xlim), mean(ylim))
      nrow = nrow + 1
    }
  }
  
  out$coords = data.frame(out$coords)
  names(out$coords) = c('x', 'y')
  out$params = data.frame(t(c(S, N, cv, sigma)))
  names(out$params) = c('S', 'N', 'cv', 'sigma')
  return(out)
}

# Combine multiple simulations with different parameters into a comm object
#   to be passed on for MOBR analysis
# The inputs should either be a single value (which will be taken as constant across treatments)
#   or vectors of the same length
sim_comm_multi_pars = function(S, N, cv, sigma, sqrt_numplots){
  # check that the input parameters have the same dimension
  # If a parameter only has one level, it would be taken as constant
  #   across all groups
  lengths = c(length(S), length(N), length(cv), length(sigma), length(sqrt_numplots))
  max_lengths = max(lengths)
  for (i in 1:length(lengths)){
    if (lengths[i] > 1 & lengths[i] < max_lengths)
      stop("Error: the lengths of the input parameters need to match, or equal to 1.")
  }
  comm = data.frame(matrix(NA, nrow = 0, ncol = max(S)))
  coords = data.frame(matrix(NA, nrow = 0, ncol = 2))
  names(coords) = c('x', 'y')
  env = data.frame(matrix(NA, nrow = 0, ncol = 4))
  
  for (j in 1:max_lengths){
    out = sim_comm_single_pars(ifelse(is.na(S[j]), S[1], S[j]),
                               ifelse(is.na(N[j]), N[1], N[j]),
                               ifelse(is.na(cv[j]), cv[1], cv[j]),
                               ifelse(is.na(sigma[j]), sigma[1], sigma[j]),
                               ifelse(is.na(sqrt_numplots[j]), sqrt_numplots[1], sqrt_numplots[j]))
    if (ncol(out$comm) < max(S)){
      comm_extent = data.frame(matrix(0, nrow(out$comm), max(S)))
      sampl_S = sample(1:max(S), ncol(out$comm))
      comm_extent[, sampl_S] = out$comm
    }
    else
      comm_extent = out$comm
    comm = rbind(comm, comm_extent)
    numplots = (ifelse(is.na(sqrt_numplots[j]), sqrt_numplots[1], sqrt_numplots[j])) ** 2
    coords = rbind(coords, out$coords)
    env = rbind(env, data.frame(matrix(out$params, nrow = numplots, ncol = length(out$params), byrow = T)))
  }
  names(env) = c('S', 'N', 'cv', 'sigma')
  for (icol in 1:ncol(env))
    env[, icol] = as.numeric(env[, icol])
  comm_obj = make_comm_obj(comm, cbind(coords, env))
  return(comm_obj)
}

const_SAD = 0; const_SAD_error = 0
const_N = 0; const_N_error = 0
const_agg = 0; const_agg_error = 0

N_SAD = 0; N_SAD_error = 0
N_N = 0; N_N_error = 0
N_agg = 0; N_agg_error = 0

SAD_SAD = 0; SAD_SAD_error = 0
SAD_N = 0; SAD_N_error = 0
SAD_agg = 0; SAD_agg_error = 0

agg_SAD = 0; agg_SAD_error = 0
agg_N = 0; agg_N_error = 0
agg_agg = 0; agg_agg_error = 0

for (i in 1:50){
  # (Almost) no change
  comm = sim_comm_multi_pars(100, c(401, 400), 1, 0.2, 4)
  mobr_const = get_delta_stats(comm, 'N', ref_group = 400, inds = 30, nperm = 100)
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\const\\', 
            '3curves_', as.character(i), '.png', sep = ''))
  plot_rarefy(mobr_const)
  dev.off()
  
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\N\\', 
            'Noeffect_', as.character(i), '.png', sep = ''))
  plot(mobr_const)
  dev.off()
  
  focal_dat = mobr_const$discrete$indiv[complete.cases(mobr_const$discrete$indiv),]
  for (irow in 2:nrow(focal_dat)){
    const_SAD = const_SAD + 1
    if (as.numeric(as.character(focal_dat$deltaS_emp[irow])) < as.numeric(as.character(focal_dat$deltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$deltaS_emp[irow])) > as.numeric(as.character(focal_dat$deltaS_null_high[irow])))
      const_SAD_error = const_SAD_error + 1
  }

  focal_dat = mobr_const$discrete$N[complete.cases(mobr_const$discrete$N),]
  for (irow in 2:nrow(focal_dat)){
    const_N = const_N + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      const_N_error = const_N_error + 1
  }
  
  focal_dat = mobr_const$discrete$agg[complete.cases(mobr_const$discrete$agg),]
  for (irow in 2:nrow(focal_dat)){
    const_agg = const_agg + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      const_agg_error = const_agg_error + 1
  }
  
  # Changing N 
  comm = sim_comm_multi_pars(100, c(1000, 400), 1, 0.2, 4)
  mobr_N = get_delta_stats(comm, 'N', ref_group = 400, inds = 30, nperm = 100)
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\N\\', 
            '3curvesN_', as.character(i), '.png', sep = ''))
  plot_rarefy(mobr_N)
  dev.off()
  
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\N\\', 
            'effectN_', as.character(i), '.png', sep = ''))
  plot(mobr_N)
  dev.off()
  
  focal_dat = mobr_N$discrete$indiv[complete.cases(mobr_N$discrete$indiv),]
  for (irow in 2:nrow(focal_dat)){
    N_SAD = N_SAD + 1
    if (as.numeric(as.character(focal_dat$deltaS_emp[irow])) < as.numeric(as.character(focal_dat$deltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$deltaS_emp[irow])) > as.numeric(as.character(focal_dat$deltaS_null_high[irow])))
      N_SAD_error = N_SAD_error + 1
  }
  
  focal_dat = mobr_N$discrete$N[complete.cases(mobr_N$discrete$N),]
  for (irow in 2:nrow(focal_dat)){
    N_N = N_N + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      N_N_error = N_N_error + 1
  }
  
  focal_dat = mobr_N$discrete$agg[complete.cases(mobr_N$discrete$agg),]
  for (irow in 2:nrow(focal_dat)){
    N_agg = N_agg + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      N_agg_error = N_agg_error + 1
  }
  
  # Changing SAD
  comm = sim_comm_multi_pars(100, 1000, c(1, 2), 0.2, 4)
  mobr_sad = get_delta_stats(comm, 'cv', ref_group = 1, inds = 30, nperm = 100)
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\SAD\\', 
            '3curvesSAD_', as.character(i), '.png', sep = ''))
  plot_rarefy(mobr_sad)
  dev.off()
  
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\SAD\\', 
            'effectSAD_', as.character(i), '.png', sep = ''))
  plot(mobr_sad)
  dev.off()

  focal_dat = mobr_sad$discrete$indiv[complete.cases(mobr_sad$discrete$indiv),]
  for (irow in 2:nrow(focal_dat)){
    SAD_SAD = SAD_SAD + 1
    if (as.numeric(as.character(focal_dat$deltaS_emp[irow])) < as.numeric(as.character(focal_dat$deltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$deltaS_emp[irow])) > as.numeric(as.character(focal_dat$deltaS_null_high[irow])))
      SAD_SAD_error = SAD_SAD_error + 1
  }
  
  focal_dat = mobr_sad$discrete$N[complete.cases(mobr_sad$discrete$N),]
  for (irow in 2:nrow(focal_dat)){
    SAD_N = SAD_N + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      SAD_N_error = SAD_N_error + 1
  }
  
  focal_dat = mobr_sad$discrete$agg[complete.cases(mobr_sad$discrete$agg),]
  for (irow in 2:nrow(focal_dat)){
    SAD_agg = SAD_agg + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      SAD_agg_error = SAD_agg_error + 1
  }
  
  # Changing aggregation
  comm = sim_comm_multi_pars(100, 1000, 1, c(0.1, 1), 4)
  mobr_aggr = get_delta_stats(comm, 'sigma', ref_group = 0.1, inds = 30, nper = 100)
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\aggr\\', 
            '3curvesaggr_', as.character(i), '.png', sep = ''))
  plot_rarefy(mobr_aggr)
  dev.off()
  
  png(paste('C:\\Users\\Xiao\\Dropbox\\projects\\mobr\\sensitivity\\aggr\\', 
            'effectaggr_', as.character(i), '.png', sep = ''))
  plot(mobr_aggr)
  dev.off()
  
  focal_dat = mobr_aggr$discrete$indiv[complete.cases(mobr_aggr$discrete$indiv),]
  for (irow in 2:nrow(mobr_aggr$discrete$indiv)){
    agg_SAD = agg_SAD + 1
    if (as.numeric(as.character(focal_dat$deltaS_emp[irow])) < as.numeric(as.character(focal_dat$deltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$deltaS_emp[irow])) > as.numeric(as.character(focal_dat$deltaS_null_high[irow])))
      agg_SAD_error = agg_SAD_error + 1
  }
  
  focal_dat = mobr_aggr$discrete$N[complete.cases(mobr_aggr$discrete$N),]
  for (irow in 2:nrow(focal_dat)){
    agg_N = agg_N + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      agg_N_error = agg_N_error + 1
  }
  
  focal_dat = mobr_aggr$discrete$agg[complete.cases(mobr_aggr$discrete$agg),]
  for (irow in 2:nrow(focal_dat)){
    agg_agg = agg_agg + 1
    if (as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) < as.numeric(as.character(focal_dat$ddeltaS_null_low[irow])) |
        as.numeric(as.character(focal_dat$ddeltaS_emp[irow])) > as.numeric(as.character(focal_dat$ddeltaS_null_high[irow])))
      agg_agg_error = agg_agg_error + 1
  }
}

print(paste("Constant parameters, detected differences in SAD, N, aggregation:",
            as.character(const_SAD_error / const_SAD), 
            as.character(const_N_error / const_N),
            as.character(const_agg_error / const_agg)), sep = '  ')

print(paste("Changing SAD, detected differences in SAD, N, aggregation:",
            as.character(SAD_SAD_error / SAD_SAD), 
            as.character(SAD_N_error / SAD_N),
            as.character(SAD_agg_error / SAD_agg)), sep = '  ')

print(paste("Changing N, detected differences in SAD, N, aggregation:",
            as.character(N_SAD_error / N_SAD), 
            as.character(N_N_error / N_N),
            as.character(N_agg_error / N_agg)), sep = '  ')

print(paste("Changing aggregation, detected differences in SAD, N, aggregation:",
            as.character(agg_SAD_error / agg_SAD), 
            as.character(agg_N_error / agg_N),
            as.character(agg_agg_error / agg_agg)), sep = '  ')

