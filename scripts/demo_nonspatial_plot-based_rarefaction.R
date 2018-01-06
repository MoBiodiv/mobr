## objectives of this script
# to demonstrate that the analytical approach of rescaling
# the individual-based rarefaction curve to compute the non-spatial, plot-based
# rarefaction curve
# 
library(mobr)
library(mobsim)
library(vegan)
set.seed(7)
# simulate a community 
S_pool = 100
N = c(1000, 500)
n_quadrats = 10
maps = lapply(N, function(n) mobsim::sim_thomas_community(S_pool, n, 'lnorm'))
names(maps) = c('hi', 'lo')

pdf('./figs/figS1_rescaling_rarefaction.pdf')
#maps of the two treatments
plot(maps$hi, axes=F, xlab='', ylab='', main='High density treatment')
plot(maps$lo, axes=F, xlab='', ylab='', main='Low density treatment')

# sample n_quadrats from the mapped communities
comms = lapply(maps, function(x) sample_quadrats(x, n_quadrats, plot = F))

density = lapply(comms, function(x) sum(x$spec_dat) / n_quadrats)
names(density) = c('hi', 'lo')
N = sapply(comms, function(x) sum(x$spec_dat))
N_max = max(N)
N_min = min(N)

# sampling effort for individual-based curve
effort_indiv = lapply(N, function(n) 1:n)

# compute effort for plot-based curve by computing the 
# number of individuals on average per number of samples using 
# the treatment specific density. Round this quantity to the nearest
# integer so it can be used in a formula for individual-based rarefaction:
effort_plot = lapply(density, function(x) x * 1:n_quadrats)

# also compute the number of individuals encountered on average with a plot
# from either treatment by computing mean density
density_avg = mean(unlist(density))
effort_plot_avg = 1:n_quadrats * density_avg

# compute non-spatial plot-based curve by computing the individual-based
# rarefaction curves at the number of individuals observed for a given
# number of plots
rare_plot = lapply(1:2, function(i) 
                   rarefaction(comms[[i]]$spec_dat, 'indiv', effort_plot[[i]],
                               dens_ratio = 1 + 1e-14)) 
names(rare_plot) = c('hi', 'lo')

# compute individual-based rarefaction
rare_indiv = lapply(1:2, function(i)
                   rarefaction(comms[[i]]$spec_dat, 'indiv', effort_indiv[[i]]))
names(rare_indiv) = c('hi', 'lo')

# compute the rescaled curves based on density at every integer so that 
# comparisons with the individual based curves can be made at each integer
rare_indiv_rescale = lapply(1:2, function(i)
                           rarefaction(comms[[i]]$spec_dat, 'indiv', 1:N_min, 
                                       dens_ratio = density_avg / density[[i]]))
names(rare_indiv_rescale) = c('hi', 'lo')
# plot the rarefaction curves
ylim = range(rare_indiv, rare_plot, rare_indiv_rescale)
lwd = 2

#par(mfrow=c(2,2))
# individual based curves
plot(rare_indiv$hi, col='red', type='p', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, lwd=lwd)
points(rare_indiv$lo, col='blue', lwd=lwd)
abline(v=N_min)
legend('bottomright', c('low N', 'high N'), col=c('blue', 'red'), pch=1, 
       bty='n')

# the non-spatial 
plot(effort_plot_avg, rare_plot$hi, col='red', type='p',
     xlab='Number of Individuals', ylab='Species richness',
     ylim=ylim, xlim=range(1, effort_plot_avg), lwd=lwd, axes=F, frame.plot=T)
lines(effort_plot_avg, rare_plot$lo, col='blue', type='p', lwd=lwd)
axis(side=1)
axis(side=2)
axis(side=3, at = effort_plot_avg, labels=1:10)
mtext(side=3, 'Number of plots', padj=-4)
abline(v=N_min)
# add the gamma analyatical form
lines(rare_indiv_rescale$hi, col='red', lty=1, lwd=lwd)
lines(rare_indiv_rescale$lo, col='blue', lty=1, lwd=lwd)

dev.off()

## older superfluous code. 

# From 1 to the minimum number of individuals between the two treatments (N_min)
# we need to compute the differences between the individual-based and the 
# plot-based rarefaction curves and then take the difference of those differences. 
# This is difficult because the plot-based curve is in units of plots rather
# than units of individuals and therefore it needs to be rescaled so that 
# the differences line up along the x-axis. 

# spline interpolate the plot-based curves from 1 to N_min
samp_rare_smooth = lapply(1:2, function(i) 
                          pracma::pchip(effort_samp_avg, rare_samp[[i]], 1:N_min))
names(samp_rare_smooth) = c('hi', 'lo')

lines(1:N_min, samp_rare_smooth$hi, col='red')
lines(1:N_min, samp_rare_smooth$lo, col='blue')

# compute difference of interpolated lines to get delta S
deltaS_smooth = samp_rare_smooth$hi - samp_rare_smooth$lo
deltaS_analyt = rarefaction(comms$hi$spec_dat, 'indiv', 1:N_min, 
                            dens_ratio = density_avg / density$hi) - 
                rarefaction(comms$lo$spec_dat, 'indiv', 1:N_min, 
                            dens_ratio = density_avg / density$lo)

plot(deltaS_smooth, ylim=range(deltaS_smooth, deltaS_analyt),
     xlab='Number of individuals', ylab='delta S (effect of N and SAD)', type='p', lwd=lwd)
lines(deltaS_analyt, col='purple', lwd=lwd)

# compute additional difference after considering diff in indiv_rare
deltaS_indiv = rarefaction(comms$hi$spec_dat, 'indiv', 1:N_min) - 
               rarefaction(comms$lo$spec_dat, 'indiv', 1:N_min)
ddeltaS_smooth = deltaS_smooth - deltaS_indiv
ddeltaS_analyt = deltaS_analyt - deltaS_indiv

plot(ddeltaS_smooth, ylim=range(ddeltaS_smooth, ddeltaS_analyt),
     xlab='Number of individuals', ylab='delta S (effect of N)',
     type='p', lwd=lwd)
lines(ddeltaS_analyt, col='purple', lwd=lwd)

# show how the curve is shifted
# how each curve is rescaled
plot(rare_indiv_rescale$hi[1:N_min], col='red', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, xlim=c(1, N_max), lty=2, lwd=lwd)
lines(rare_indiv$hi, col='red', type='l', lwd=lwd)
lines(rare_indiv_rescale$lo[1:N_min], col='blue', lty=2, lwd=lwd)
lines(rare_indiv$lo, col='blue', type='l', lwd=lwd)
abline(v=N_min)
## on seperate graphs
plot(rare_indiv_rescale$hi[1:N_min], col='red', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, xlim=c(1, N_max), lty=2, lwd=lwd)
lines(rare_indiv$hi, col='red', type='l', lwd=lwd)
plot(rare_indiv_rescale$lo[1:N_min], col='blue', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, xlim=c(1, N_max), lty=2, lwd=lwd)
lines(rare_indiv$lo, col='blue', type='l', lwd=lwd)
#rescaled comparison 
plot(rare_indiv_rescale$hi[1:N_min], col='red', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, lty=2, lwd=2)
lines(rare_indiv_rescale$lo, col='blue', lty=2, lwd=2)

dev.off()




plot(effort, S_sam_un, type='p', col='blue',
     ylim=c(0,85), xlim=c(0,Nmin), xlab='# of indiv', main='interpolation via smoother')
lines(S_sam_un2, col='blue', lwd=2)
points(effort, S_sam_in, col='red')
lines(S_sam_in2, col='red', lwd=2)




ylim = range(sam_ns_rare, sam_rare)
plot(sam_ns_rare[[1]], col='red', type='l', xlab='Number of plots',
     ylim=ylim, lwd=2, lty=3)
lines(sam_rare[[1]], col='red', lwd=2)
lines(sam_ns_rare[[2]], col='blue', lwd=2, lty=3)
lines(sam_rare[[2]], col='blue', lwd=2)




# need to rescale this curve 


# how each curve is rescaled
plot(ind_rescale_rare[[1]][1:nmin], col='red', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, xlim=c(1, nmax), lty=2, lwd=lwd)
lines(ind_rare[[1]], col='red', type='l', lwd=lwd)
plot(ind_rescale_rare[[2]][1:nmin], col='blue', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, xlim=c(1, nmax), lty=2, lwd=lwd)
lines(ind_rare[[2]], col='blue', type='l', lwd=lwd)
#rescaled comparison 
plot(ind_rescale_rare[[1]][1:nmin], col='red', type='l', xlab='Number of individuals',
     ylab='Species richness', ylim=ylim, lty=2, lwd=2)
lines(ind_rescale_rare[[2]], col='blue', lty=2, lwd=2)
# delta curve
plot(ind_rescale_rare[[2]] - ind_rescale_rare[[1]][1:nmin], xlab='Number of individuals', 
     ylab='delta S', type='l', lwd=2)
dev.off()

# need to demonstrate curves match
#1) match between an average abu indiv to sample curve (analytical)
plot(sam_ns_perm_rare[[1]], col='red', type='p', xlab='Number of plots',
     ylab='Species richness', ylim=ylim, lwd=lwd)
lines(sam_ns_rare[[1]], lwd=lwd, col='blue')
## the randomization is actually under the analytical because it has
# a soft constraint b/c sampling from SAD is with replacement therefore,
# on average Snull < Sobs. 

#2) match between randomization and analytical rescaling raw curve
replicate(10, rarefaction(


#2) the delta curves match
permute_comm(comms[[1]]$spec_dat



sam_ns_rare = lapply(comms, function(x) 
                rarefaction(x$spec_dat, 'indiv', 
                  effort= round(sum(x$spec_dat)/n_quadrats * 1:n_quadrats)))
                                                    
sam_ns_perm_rare = lapply(comms, function(x) 
                     apply(
                       replicate(100, 
                         rarefaction(
                           permute_comm(x$spec_dat, 'noagg'), 'samp')), 1, mean))
#sam_ns_rare = lapply(comms, function(x) 
#                apply(
#                  replicate(100, 
#                    rarefaction(
#                      simulate(
#                        nullmodel(x$spec_dat, 'r0_ind'), nsim=1)[ , , 1], 'samp')), 1, mean))

sam_rare = lapply(comms, function(x) rarefaction(x$spec_dat, 'samp'))
             
             
             

library(mobr)
library(vegan)

data(inv_plot_attr)
data(inv_comm)
sad = aggregate(inv_comm, list(inv_plot_attr$group), sum)[,-1]
Nmin = min(rowSums(sad))
nplots = sum(inv_plot_attr$group == 'invaded')

ref_density = mean(rowSums(sad)/nplots)
group_density_in = rowSums(sad)[1] / nplots
dens_ratio_in = ref_density / group_density_in
group_density_un = rowSums(sad)[2] / nplots
dens_ratio_un = ref_density / group_density_un

# individual rarefaction curves across all N indiv
S_ind_un = rarefaction(inv_comm[inv_plot_attr$group == 'uninvaded', ],
                         'indiv')
S_ind_in = rarefaction(inv_comm[inv_plot_attr$group == 'invaded', ],
                 'indiv')
effort = round(1:nplots * ref_density)
un_effort = round(1:nplots * group_density_un)
in_effort = round(1:nplots * group_density_in)

S_sam_un = rarefaction(inv_comm[inv_plot_attr$group == 'uninvaded', ],
                         'indiv', un_effort)
S_sam_in = rarefaction(inv_comm[inv_plot_attr$group == 'invaded', ],
                 'indiv', in_effort)


par(mfrow=c(2,2))
# individual based rarefaction to Nmax
plot(S_ind_un, type='p', col='blue', main='individual-based rarefaction',
      ylim=c(0,85), xlab='# of indiv')
lines(S_ind_in, type='p', col='red')

# interpolate the individual curves at 1:Nmin individuals
rescale_effort_in = 1:Nmin * (ref_density / group_density_in)
rescale_effort_un = 1:Nmin * (ref_density / group_density_un)

S_ind_un2 = pracma::pchip(rescale_effort_un, S_ind_un[1:Nmin], 1:Nmin)
S_ind_in2 = pracma::pchip(rescale_effort_in, S_ind_in[1:Nmin], 1:Nmin)

plot(rescale_effort_un, S_ind_un[1:Nmin], type='p', col='blue',
     ylim=c(0,85), xlim=c(0,Nmin), xlab='# of indiv', main='interpolation via smooter')
lines(S_ind_un2, col=1, lwd=2)
points(rescale_effort_in, S_ind_in[1:Nmin], col='red')
lines(S_ind_in2, type='l', col=1, lwd=2)


plot(1:Nmin, rarefaction(inv_comm[inv_plot_attr$group == 'uninvaded', ], 'indiv', 1:Nmin, 
                         dens_ratio = dens_ratio_un), col='blue')
lines(1:Nmin, rarefaction(inv_comm[inv_plot_attr$group == 'invaded', ], 'indiv', 1:Nmin, 
                         dens_ratio = dens_ratio_in), col='red')

# create sample based curves from individual curves

plot(S_sam_un, type='p', col='blue',
     ylim=c(0,85), xlab='# of plots')
lines(S_sam_in, type='p', col='red')

#interpolate these sample curves

S_sam_un2 = pracma::pchip(effort, S_sam_un, 1:Nmin)
S_sam_in2 = pracma::pchip(effort, S_sam_in, 1:Nmin)

plot(effort, S_sam_un, type='p', col='blue',
     ylim=c(0,85), xlim=c(0,Nmin), xlab='# of indiv', main='interpolation via smoother')
lines(S_sam_un2, col='blue', lwd=2)
points(effort, S_sam_in, col='red')
lines(S_sam_in2, col='red', lwd=2)

# compare deltaS to analytical solution
deltaS_N_smoother = S_sam_un2 - S_ind_un[1:Nmin] 
deltaS_N_analy = deltaS_N(inv_comm[inv_plot_attr$group == 'uninvaded', ], ref_density, 1:Nmin)
plot(deltaS_N_smoother, ylim=range(deltaS_N_analy$deltaS, deltaS_N_smoother))
lines(deltaS_N_analy, col='red')
               
                



##-------------------------------------------------



S_perm_in = apply(replicate(100, rarefaction(permute_comm(inv_comm, 'noagg',
                                                 inv_plot_attr$group)[inv_plot_attr$group == 'invaded', ], 
                                    'samp', 1:nplots)), 1, mean)
S_perm_un = apply(replicate(100, rarefaction(permute_comm(inv_comm, 'noagg',
                                                 inv_plot_attr$group)[inv_plot_attr$group == 'uninvaded', ], 
                                    'samp', 1:nplots)), 1, mean)


S_ind_in = rarefaction(inv_comm[inv_plot_attr$group == 'invaded', ],
                      'indiv', effort = round(group_density_in * 1:nplots))
S_rescaled_in = rarefaction(inv_comm[inv_plot_attr$group == 'invaded', ],
                      'indiv', round(group_density_in * 1:nplots), dens_ratio = dens_ratio_in)
S_ind_un = rarefaction(inv_comm[inv_plot_attr$group == 'uninvaded', ],
                      'indiv', round(group_density_in * 1:nplots))
S_rescaled_un = rarefaction(inv_comm[inv_plot_attr$group == 'uninvaded', ],
                      'indiv', round(group_density_in * 1:nplots), dens_ratio = dens_ratio_un)

plot(1:nplots * group_density_in, S_perm_in, col='red')
lines(1:Nmin, S_ind_in)
lines(1:Nmin, S_rescaled_in, col='red')

dS_perm = (S_perm_in - S_ind_in) - (S_perm_un - S_ind_un)
dS_anly = (S_rescaled_in - S_ind_in) - (S_rescaled_un - S_ind_un)

plot(rarefaction(inv_comm[inv_plot_attr$group == 'invaded', ],
                      'indiv', effort = 1:Nmin) - 
      rarefaction(inv_comm[inv_plot_attr$group == 'uninvaded', ],
                      'indiv', effort = 1:Nmin), log='x') 
plot(dS_perm, log='x', ylim=c(-50, 0))
lines(dS_anly)

##-------------
data(mite)
x = mite

x = dat[[5]]$comm

S = rarefaction(x, 'indiv')
Sr = rescale_rarefaction(x, 1, 1:sum(x))
time1 = system.time(replicate(5, rarefaction(x, 'indiv'))) / 5 
time2 = system.time(replicate(5, rescale_rarefaction(x, 1, 1:sum(x)))) / 5
time1
time2
time2[3] - time1[3]
time2[3] / time1[3]

#Sv = rarefy(colSums(x), 1:sum(x))
cor(S, Sr)


par(mfrow=c(1,1))  
plot(S, Sr)
abline(a=0, b=1, col='red')
#plot(S, Sv)
#abline(a=0, b=1, col='red')

# compare interpolation approach to analytical approach ---------------------

rescale_rarefaction = function(comm, ratio, effort){
    if (!is.null(dim(comm)))
        x = colSums(comm)
    n = sum(x)
    x = x[x > 0] 
    S = length(x)
    effort = effort[effort/ratio <= n]
    ldiv = lgamma(n - effort/ratio + 1) - lgamma(n + 1)
    p = matrix(0, length(effort), S)
    for (i in seq_along(effort)) {
        for (j in 1:S){
            p[i, j] = ifelse(n - x[j] < effort[i]/ratio, 0, 
                            exp(lgamma(n - x[j] + 1) - lgamma(n - x[j] - effort[i]/ratio + 1)
                                + ldiv[i]))
        }
    }
    out = rowSums(1 - p)
    return(out)
}

# uses the analytical approach
deltaS_N = function(comm, ref_dens, inds){
    nplots = nrow(comm)
    group_dens = sum(comm) / nplots
    ratio = ref_dens / group_dens
    interp_S_samp = rescale_rarefaction(comm, ratio, inds)
    # Ensure that interp_S_samp has the right length (filled with NA's if needed)
    interp_S_samp = interp_S_samp[1:length(inds)]
    S_indiv = rarefaction(comm, 'indiv', inds)
    deltaS = interp_S_samp - S_indiv
    out = data.frame(inds = inds, deltaS = deltaS)
    return(out)
}


# here is the intperolation approach
deltaS_N2 = function(comm, ref_dens, inds){
    nplots = nrow(comm)
    group_dens = sum(comm) / nplots
    # assign sampling effort such that if the density has higher
    # than the ref_density that the sampling effort goes high 
    # enough so that when it is eventually rescaled (i.e., shrunk)
    # it has spanned a large enough domain to be compared across inds
    # Note: the effort increment is one individual but for datasets
    # with lots of individuals (hundreds of thousands) this will be 
    # slow. 
    if (group_dens > ref_dens)
        effort = 1:(ref_dens * nplots)
    else
        effort = 1:max(inds)
    # use individual based rarefaction to generate richness values
    S_samp = rarefaction(comm, 'indiv', effort)
    # rescale the effort to reflect what it would be if no treatment
    # differences in density
    rescale_effort = effort * (ref_dens / group_dens)
    # No extrapolation of the rescaled rarefaction curve, only interpolation
    interp_S_samp = pracma::pchip(rescale_effort, S_samp, inds)    
    S_indiv = rarefaction(comm, 'indiv', inds)
    deltaS = interp_S_samp - S_indiv
    out = data.frame(inds = inds, deltaS = deltaS)
    return(out)
}

## --
comm = dat$mor$comm
env = dat$mor$env
comm_un = comm[env$groups == 'before', ]
comm_in = comm[env$groups == 'after', ]
inds = 1:min(sum(comm_un), sum(comm_in))
nplots = nrow(comm_un)
group_dens_un = sum(comm_un) / nplots
group_dens_in = sum(comm_in) / nplots
ref_dens = mean(c(group_dens_un, group_dens_in))
#if (group_dens > ref_dens)
    effort = 1:(ref_dens * nplots)
#else
    effort = 1:max(inds)
comm = comm_un
group_dens = group_dens_un
S_samp = rarefaction(comm, 'indiv', effort)
rescale_effort = effort * (ref_dens / group_dens)
interp_S_samp = pracma::pchip(rescale_effort, S_samp, inds)    
# ratio is the rescaling factor, e.g., mean_density/treatment_density
interp_S_samp2 = rescale_rarefaction(comm, ref_dens / group_dens, inds)
S_indiv = rarefaction(comm, 'indiv', inds)

plot(interp_S_samp, interp_S_samp2)
abline(a=0, b=1, col='red')

plot(S_indiv, ylim=range(S_indiv, interp_S_samp2))
points(interp_S_samp, col='green3')
lines(interp_S_samp2, col='red', lwd=2)
legend('bottomright', c('hypergeo', 'rescaled gamma', 'rescaled interpolated'), 
       col=c('black', 'red', 'green'), lty=c(NA, 1, NA), 
       lwd=c(NA, 2, NA), pch=c(19, NA, 19), bty='n')



##------------------------
comm = dat$inv$comm
nplots = nrow(comm)
group_dens = sum(comm) / nplots
ratio = 0.95
ref_dens = ratio * group_dens
inds = 1:sum(comm)
rs_inds = 1:(2*ref_dens - group_dens)
#if (group_dens > ref_dens)
    effort = 1:(ref_dens * nplots)
#else
    effort = 1:max(inds)
S_samp = rarefaction(comm, 'indiv', effort)
rescale_effort = effort * (ref_dens / group_dens)
interp_S_samp = pracma::pchip(rescale_effort, S_samp, inds)    
# ratio is the rescaling factor, e.g., mean_density/treatment_density
interp_S_samp2 = rescale_rarefaction(comm, ref_dens / group_dens, rs_inds)
S_indiv = rarefaction(comm, 'indiv', inds)

plot(S_indiv, ylim=range(S_indiv, interp_S_samp2))
points(interp_S_samp, col='green3')
lines(interp_S_samp2, col='red', lwd=2)
legend('bottomright', c('hypergeo', 'rescaled gamma', 'rescaled interpolated'), 
       col=c('black', 'red', 'green'), lty=c(NA, 1, NA), 
       lwd=c(NA, 2, NA), pch=c(19, NA, 19), bty='n')


#




