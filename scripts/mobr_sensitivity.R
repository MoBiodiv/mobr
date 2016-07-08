# Sensitivity test for Mobr functions
# using Felix May's MoBspatial package

library(MoBspatial)
#library(mobr) # This doesn't work at the moment; will have to wait until the mobr package is built


# Simulate a community using Sim.Thomas.Community() and wrangle it into a list
#   which can later be combined to create the comm object for mobr.
# Inputs:
# S, N, cv, sigma - parameters to simulate the community. xmax and ymax are fixed at 1.
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
  out$comm = matrix(NA, sqrt_numplots^2, nlevels(sim_comm$SpecID))
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
      plot_sad = data.frame(table(comm_plot$SpecID))[, 2]
      out$comm[nrow, ] = plot_sad
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
