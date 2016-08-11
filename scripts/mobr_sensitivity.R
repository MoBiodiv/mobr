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

Niter = 5 # Repeat simulation Niter times
sqrt_numplots = 4 # Each group divided into 4*4 = 16 plots
results = data.frame(matrix(0, nrow = 4, ncol = 6))
row.names(results) = c('const', 'N', 'SAD', 'agg')
names(results) = c('SADTot', 'SADErr', 'NTot', 'NErr', 'AggTot', 'AggErr')

par_table = data.frame(matrix(c(100, 1000, 1, 0.2, 
                                100, 1001, 1, 0.2,  
                                100, 400, 1, 0.2,  
                                100, 1000, 2, 0.2, 
                                100, 1000, 1, 1), byrow = T, ncol = 4))
names(par_table) = c('S', 'N', 'cv', 'sigma')
row.names(par_table) = c('ref', row.names(results))

ref_par = par_table['ref', ]
for (type in row.names(results)){
  type_par = par_table[type, ]
  par_list = c()
  for (j in 1:length(ref_par)){
    if (type_par[j] == ref_par[j]) par_list = c(par_list, as.numeric(type_par[j]))
    else {
      env_var = names(par_table)[j]
      par_list = c(par_list, list(c(as.numeric(ref_par[j]),
                                    as.numeric(type_par[j]))))
    }
  }
  # "Equivalent" groups are defined by tiny difference in N
  if (type == 'const') type_mobr = 'N' 
  else type_mobr = type
  
  for (i in 1:Niter){
    comm = sim_comm_multi_pars(par_list[[1]], par_list[[2]], par_list[[3]], 
                               par_list[[4]], sqrt_numplots)
    mobr = get_delta_stats(comm, env_var, inds = 30, nperm = 100, 
                           ref_group = as.numeric(ref_par[env_var]))
    # Optional plotting:
    # plot_rarefy(mobr)
    # plot(mobr)
    
    effects = c('indiv', 'N', 'agg')
    for (k in 1:length(effects)){
      effect = effects[k]
      focal_dat = mobr$discrete[[effect]]
      focal_dat = focal_dat[complete.cases(focal_dat), ]
      if (k == 3) start = 1 # Ignore first row (n = 1) for SAD and N
      else start = 2
      for (irow in 2:nrow(focal_dat)){
        results[type, 2 * k - 1] = results[type, 2 * k - 1] + 1
        if (as.numeric(as.character(focal_dat[irow, 3])) < as.numeric(as.character(focal_dat[irow, 4])) |
            as.numeric(as.character(focal_dat[irow, 3])) > as.numeric(as.character(focal_dat[irow, 6])))
          results[type, 2 * k] = results[type, 2 * k] + 1
      }
    }
  }
}

for (type in row.names(results)){
  print(paste(type, ", detected differences in SAD, N, aggregation:",
              as.character(results[type, 2] / results[type, 1]), 
              as.character(results[type, 4] / results[type, 3]),
              as.character(results[type, 6] / results[type, 5])), sep = '  ')
  
}

