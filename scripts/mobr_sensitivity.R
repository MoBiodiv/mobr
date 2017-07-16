# Sensitivity test for Mobr functions
# using Felix May's MoBspatial package
.libPaths("C://Users//Xiao//Documents//R//win-library//3.2")
library(mobsim)
library(mobr)
# Generate an SAD from a poisson lognormal distribution with fixed S, N, and parameters
# This function is copied and modified from the one in MoBspatial for our purpose
# The set.seed steps ensure that the underlying distribution remains completely
#     unchanged when the shape parameters are the same
SAD.lognorm <- function(S.pool, N.local, mean.abund = 100, cv.abund = 1,
                        rand.seed = 10){
    sd.abund <- mean.abund*cv.abund

    sigma1 <- sqrt(log(sd.abund^2/mean.abund^2 +1))
    mu1 <- log(mean.abund) - sigma1^2/2

    set.seed(rand.seed)
    abund.pool <- rlnorm(S.pool, meanlog=mu1, sdlog=sigma1)
    relabund.pool <- sort(abund.pool/sum(abund.pool), decreasing = T)
    set.seed(NULL)
    abund.local <- table(sample(1:S.pool, N.local, replace = T, prob = relabund.pool))
    names(abund.local) <- paste("species", names(abund.local), sep = "")

    return(list(abund = abund.local,
                mean.theor = N.local/S.pool,
                sd.theor   = N.local/S.pool * cv.abund,
                mean.sim   = mean(abund.local),
                sd.sim     = sd(abund.local))
    )
}

# Simulate a community using Sim.Thomas.Community() and wrangle it into a list
#   which can later be combined to create the comm object for mobr.
# Inputs:
# S, N, cv, sigma - parameters to simulate the community. xmax and ymax are fixed at 1.
# Note that now S is S_pool (total richness for the metacommunity), while
# N is N_local (number of individuals for the local group)
# Because of the sampling process, the local group may not have all S species.
# sqrt_numplots - sqrt(number of plots) within the domain. This is how the domain will be divided.
# rand.seed - random seed for SAD, which controls whether the underlying metacommunity
#               SAD remains unchanged or not.
# Ouput:
# A list with three components:
# $comm - plot by species matrix
# $coords - coordinates of the plots
# $params - parameters used to create the simulation, including S, N, cv (shape of the SAD),
#   and sigma (spatial aggregation)
sim_comm_single_pars = function(S, N, cv, sigma, sqrt_numplots, rand.seed = 10){
    sim_sad = SAD.lognorm(S, N, cv.abund = cv, rand.seed = rand.seed)
    sim_comm = sim_thomas_coords(sim_sad$abund, sigma)$census
    out = list()
    out$comm = matrix(0, sqrt_numplots^2, S)
    colnames(out$comm) = sapply(1:S, function(x) paste('species', as.character(x),
                                                    sep = ''))
    out$coords = matrix(NA, sqrt_numplots^2, 2)
    x_span = 1 / sqrt_numplots
    y_span = 1 / sqrt_numplots

    nrow = 1
    for (i in 1:sqrt_numplots){
        xlim = c(i-1, i) * x_span
        for (j in 1:sqrt_numplots){
            ylim = c(j-1, j) * y_span
            cond = which(sim_comm$x >= xlim[1] & sim_comm$x <= xlim[2] &
                             sim_comm$y >= ylim[1] & sim_comm$y <= ylim[2])
            comm_plot = sim_comm[cond, ]
            plot_sad = data.frame(table(comm_plot$species))
            for (k in 1:nrow(plot_sad)){
                out$comm[nrow, which(colnames(out$comm) == as.character(plot_sad[k, 1]))] =
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
# The inputs should either be a single value (which will be taken as constant
#   across treatments) or vectors of the same length
sim_comm_multi_pars = function(S, N, cv, sigma, sqrt_numplots){
    # check that the input parameters have the same dimension
    # If a parameter only has one level, it would be taken as constant
    #   across all groups
    lengths = c(length(S), length(N), length(cv), length(sigma),
                length(sqrt_numplots))
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
        } else {
            comm_extent = out$comm
        }
        comm = rbind(comm, comm_extent)
        numplots = (ifelse(is.na(sqrt_numplots[j]), sqrt_numplots[1],
                           sqrt_numplots[j])) ** 2
        coords = rbind(coords, out$coords)
        env = rbind(env, data.frame(matrix(out$params, nrow = numplots,
                                           ncol = length(out$params), byrow = T)))
    }
    names(env) = c('S', 'N', 'cv', 'sigma')
    for (icol in 1:ncol(env))
        env[, icol] = as.numeric(env[, icol])
    comm_obj = make_mob_in(comm, cbind(coords, env))
    return(comm_obj)
}

# Compare a group with one or more parameters changed to the reference group
# Inputs:
# ref_pars, comp_var - a vecotr of length 5, for S, N, cv, sigma, and sqrt_numplots
#     Currently we assume only N, cv, and sigma are changing
# output:
# A vector of length 9. The first three elements show if N, cv, and sigma have
#     changed (NA if no change, new value if changed). The last six elements
#     are the total number of points evaluated for SAD, N & aggregation, and the
#     number of times where the points are outside of the 95% CI.

comp_two_groups = function(ref_pars, comp_pars, sqrt_numplots, Niter){
    if (ref_pars[1] != comp_pars[1])
        stop('For now, assume S to be the same between groups.')
    new_par_val = rep(NA, 3)
    par_list = list()
    par_names = c('N', 'cv', 'sigma')
    for (j in 2:4){
        if (ref_pars[j] == comp_pars[j]){
            par_list[[j - 1]] = ref_pars[j]
        } else{
            par_list[[j - 1]] = c(ref_pars[j], comp_pars[j])
            new_par_val[j - 1] = comp_pars[j]
            env_var = par_names[j - 1]
            ref_group = ref_pars[j]
        }
    }
    eval_counts = rep(0, 6)
    for (i in 1:Niter){
        comm = sim_comm_multi_pars(ref_pars[1], par_list[[1]], par_list[[2]],
                                   par_list[[3]], sqrt_numplots)
        mobr = get_delta_stats(comm, env_var, inds = 30, nperm = 200,
                               ref_group = ref_group)

        effects = c('SAD', 'N', 'agg')
        for (k in 1:length(effects)){
            effect = effects[k]
            focal_dat = mobr[[effect]]
            focal_dat = focal_dat[complete.cases(focal_dat), ]
            if (k == 3) start = 1 # Ignore first row (n = 1) for SAD and N
            else start = 2
            for (irow in start:nrow(focal_dat)){
                eval_counts[2 * k - 1] = eval_counts[2 * k - 1] + 1
                if (as.numeric(as.character(focal_dat[irow, 3])) < as.numeric(as.character(focal_dat[irow, 4])) |
                    as.numeric(as.character(focal_dat[irow, 3])) > as.numeric(as.character(focal_dat[irow, 6])))
                   eval_counts[2 * k] = eval_counts[2 * k] + 1
            }
        }
    }
    return(c(new_par_val, eval_counts))
}

Niter = 500 # Repeat simulation Niter times
sqrt_numplots = 4 # Each group divided into 4*4 = 16 plots
results = data.frame(matrix(NA, nrow = 0, ncol = 9))
ref_pars = c(100, 1000, 2, 10) # This is the reference to which all other cases will be compared
ref_pars_almost_equal = c(100, 1001, 2, 10)
out = comp_two_groups(ref_pars, ref_pars_almost_equal, sqrt_numplots, Niter)
results = rbind(results, out)
names(results) = c('N', 'cv', 'sigma', 'SADTot', 'SADErr', 'NTot',
                   'NErr', 'AggTot', 'AggErr')
write.csv(results, 'C:\\Users\\Xiao\\Documents\\GitHub\\mobr\\scripts\\mobr_sensitivity.csv',
          row.names = F, quote = F)

# Scenarios where only one aspect has changed
par_list = list(N = c(700, 800, 900), SAD = c(0.5, 1, 1.5),
                agg = c(0.02, 0.05, 0.1))
for (i in 1:length(par_list)){
   pars = par_list[[i]]
   for (j in 1:length(pars)){
      comp_pars = ref_pars
      comp_pars[i + 1] = pars[j]
      out = comp_two_groups(ref_pars, comp_pars, sqrt_numplots, Niter)
      results = rbind(results, out)
      write.csv(results, 'C:\\Users\\Xiao\\Documents\\GitHub\\mobr\\scripts\\mobr_sensitivity.csv',
                row.names = F, quote = F)
   }
}

# Scenarios where two aspects have changed simultaneously
# Using only the intermediate levels from the previous test
comp_pars_full = c(ref_pars[1], par_list[[1]][2], par_list[[2]][2],
                  par_list[[3]][2])
for (i in 2:length(comp_pars_full)){
   comp_pars = comp_pars_full
   comp_pars[i] = ref_pars[i]
   out = comp_two_groups(ref_pars, comp_pars, sqrt_numplots, Niter)
   results = rbind(results, out)
   write.csv(results, 'C:\\Users\\Xiao\\Documents\\GitHub\\mobr\\scripts\\mobr_sensitivity.csv',
             row.names = F, quote = F)
}

# Scenario where all three aspects have changed
out = comp_two_groups(ref_pars, comp_pars_full, sqrt_numplots, Niter)
results = rbind(results, out)

write.csv(results, 'C:\\Users\\Xiao\\Documents\\GitHub\\mobr\\scripts\\mobr_sensitivity.csv',
          row.names = F, quote = F)
