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

