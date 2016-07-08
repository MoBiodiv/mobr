# Sensitivity test for Mobr functions
# using Felix May's MoBspatial package

library(MoBspatial)
#library(mobr) # This doesn't work at the moment; will have to wait until the mobr package is built

# This function draws samples (plots) from a simulated community by dividing the domain up into plots.
# Inputs:
# sim_comm - a community simulated by Sim.Thomas.Community() or Sim.Poisson.Community()
#   It has three columns: the X, Y coordinates of an individual and its SpecID
# sqrt_N - sqrt(number of samples). The x- and y- coordinates of the domain will be divided
#   by this number to create the plots.
# xmax, ymax - size of the simulated community in X & Y coordinates
# Output:
# A list with two components:
# $comm: A plot (sample) * species matrix. Species are allowed to have all zero values.
# $coords: (x, y) coordinates at the center of each plot

sample_comm = function(sim_comm, sqrt_N, xmax = 1, ymax = 1){
  out = list()
  out$comm = matrix(NA, sqrt_N^2, nlevels(sim_comm$SpecID))
  out$coords = matrix(NA, sqrt_N^2, 2)
  x_span = xmax / sqrt_N
  y_span = ymax / sqrt_N
  nrow = 1
  for (i in 1:sqrt_N){
    xlim = c(i-1, i) * x_span
    for (j in 1:sqrt_N){
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
  return(out)
}

