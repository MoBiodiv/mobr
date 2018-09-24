library(mobr)

data(inv_comm)
data(inv_plot_attr)
inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)

# run analyses
inv_mob_stats = get_mob_stats(inv_mob_in, 'group', n_perm=999)
inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_group='uninvaded',
                              type='discrete', log_scale=TRUE, n_perm=999,
                              overall_p=TRUE)

cols = c(rgb(0, 37, 112, maxColorValue = 255),    # dark navy
         rgb(0, 177, 240, maxColorValue = 255),   # royal blue
         rgb(0, 160, 73, maxColorValue = 255))    # green


# create graphics
pdf('./figs/inv_mob_stats.pdf', height = 7*0.5)
plot(inv_mob_stats, col = cols[2:3])
dev.off()

pdf('./figs/inv_delta_stats.pdf')
plot(inv_mob_out, 'invaded', 'uninvaded', leg_loc = 'bottomright')
dev.off()

pdf('./figs/inv_mob_stacked.pdf')
overlap_effects(inv_mob_out, 'invaded', leg_loc = NA,
                col = cols)
overlap_effects(inv_mob_out, 'invaded', 'stacked', prop=T, leg_loc = NA,
                col = cols)
dev.off()


## generate case study + conceptual fig ---------

plot_delta_con = function(mob_out, trt_group, ref_group, same_scale=FALSE, 
                          log='', display=c('rarefaction', 'delta S', 'ddelta S'),
                          lwd=3, par_args=NULL, ...) {
  type = mob_out$type
  tests = mob_out$tests
  if (type == 'continuous')
    stop("Currently this plot only works for mob_out object with type discrete.")
  cols = list()
  cols$trt = "#00B1F0"     # blue
  cols$ref = "#00A049"     # green
  cols$deltaS = rgb(179, 216, 160, maxColorValue = 255) # light green
  cols$ddeltaS = rgb(213, 239, 252, maxColorValue = 255) # light blue
  cols$ddeltaSline = rgb(197, 233, 251, maxColorValue = 255) #light blue but darker
  if (is.null(par_args)) {
    par_args = paste('mfrow = c(', length(display) + 1, ',',
                     length(tests), '), mgp = c(2.5, 1, 0)',  sep='')
    
  } 
  eval(parse(text=paste('par(', par_args, ')')))  
  if (same_scale) {
    # not currently implemented for the delta S plots
    if ('rarefaction' %in% display) {
      if ('agg' %in% tests) 
        S_cols = c('impl_S', 'expl_S')
      else
        S_cols = 'impl_S'
      ylim_rare = range(list(mob_out$indiv_rare[ , -1],
                             mob_out$sample_rare[ , S_cols]))
    }
  }
  mob_out$sample_rare[, -1] = lapply(mob_out$sample_rare[, -1], function(x)
    as.numeric(as.character(x)))
  sample_rare_trt = mob_out$sample_rare[mob_out$sample_rare == trt_group, ]
  sample_rare_ref = mob_out$sample_rare[mob_out$sample_rare == ref_group, ]
  if ('rarefaction' %in% display) {
    groups = c(trt_group, ref_group)
    if ('agg' %in% mob_out$tests) {
      if (!same_scale)
        ylim_rare = c(0, max(mob_out$sample_rare$expl_S))
      dat_group = mob_out$sample_rare
      N1 = sample_rare_trt$sample_plot
      N2 = sample_rare_ref$sample_plot
      S1 = sample_rare_trt$expl_S
      S2 = sample_rare_ref$expl_S
      plot(N1, S1, lwd = lwd,
           type = 'n', xlab = 'Number of plots',
           ylab = 'Richness (S)', col = cols$trt,
           ylim = ylim_rare,
           main = 'sSBR', cex.axis = 1.5, cex.lab = 1.5,
           log=log, frame.plot=F)
      polygon(c(N1, rev(N2)), c(S1, rev(S2)), border=NA,
              col = cols$deltaS)            
      lines(N1, S1, lwd = lwd, col = cols$trt)
      lines(N2, S2, lwd = lwd, col = cols$ref)
      
      #legend('topleft', as.character(groups), col=as.character(unlist(cols)), 
      #       lty=1, lwd=lwd, bty='n')          
    }
    if ('N' %in% mob_out$tests) {
      if (!same_scale)
        ylim_rare = c(0, max(mob_out$sample_rare$impl_S))
      S1 = sample_rare_trt$impl_S
      S2 = sample_rare_ref$impl_S
      plot(N1, S1,
           lwd = lwd, type = 'n', xlab = 'Number of plots',
           ylab = 'Richness (S)', col = cols$trt, 
           ylim = ylim_rare,
           main = 'nsSBR', cex.axis = 1.5, cex.lab = 1.5,
           log = log, frame.plot=F)
      polygon(c(N1, rev(N2)), c(S1, rev(S2)), border=NA,
              col = cols$deltaS)            
      lines(N1, S1, lwd = lwd, col = cols$trt)
      lines(N2, S2, lwd = lwd, col = cols$ref)
    }
    if ('SAD' %in% mob_out$tests) {
      if (!same_scale)
        ylim_rare = range(mob_out$indiv_rare[, -1])
      N1 = mob_out$indiv_rare$sample
      N2 = N1
      S1 = mob_out$indiv_rare$invaded
      S2 = mob_out$indiv_rare$uninvaded          
      plot(N1, S1,
           lwd = lwd, type = 'n', col = cols$trt, xlab = 'Number of individuals', 
           ylab = 'Richness (S)', main = 'IBR', 
           cex.axis = 1.5, cex.lab = 1.5, log=log, frame.plot=F)
      polygon(c(N1, rev(N2)), c(S1, rev(S2)), border=NA,
              col = cols$deltaS)            
      lines(N1, S1, lwd = lwd, col = cols$trt)
      lines(N2, S2, lwd = lwd, col = cols$ref)
    }        
  }    
  if ('delta S' %in% display) {
    minN = min(nrow(sample_rare_trt), nrow(sample_rare_ref))
    if ('agg' %in% mob_out$tests) {
      delta_Sspat = sample_rare_trt$expl_S[1:minN] - 
        sample_rare_ref$expl_S[1:minN]
      plot(seq(minN), delta_Sspat, 
           ylim = c(min(delta_Sspat, 0), max(delta_Sspat, 0)),
           cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = lwd,
           col = cols$deltaS, xlab = 'Number of plots',
           ylab = expression(Delta * 'S'),
           main = 'Aggr., N, & SAD effects', frame.plot=F, ...)
      abline(h = 0, lwd = 1, lty = 2)
    }
    if ('N' %in% mob_out$tests) {
      delta_Ssample = sample_rare_trt$impl_S[1:minN] - 
        sample_rare_ref$impl_S[1:minN]
      plot(seq(minN), delta_Ssample, 
           ylim = c(min(delta_Ssample, 0), max(delta_Ssample, 0)),
           cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = lwd, lty = 2,
           col = cols$deltaS, xlab = 'Number of plots', 
           ylab = expression(Delta * 'S'),
           main = 'N & SAD effects', frame.plot=F, ...)
      abline(h = 0, lwd = 1, lty = 2)
    }
    if ('SAD' %in% mob_out$tests) {
      # Create the plots for the three delta-S between groups
      deltaS_Sind = mob_out$indiv_rare[[trt_group]] - 
        mob_out$indiv_rare[[ref_group]]
      plot(mob_out$indiv_rare$sample, deltaS_Sind,
           ylim = c(min(deltaS_Sind, 0), max(deltaS_Sind, 0)),
           cex.axis = 1.5, cex.lab = 1.5, type = 'l', lwd = lwd, lty = 4,
           col = cols$deltaS, xlab = 'Number of individuals', 
           ylab = expression(Delta * 'S'),
           main = 'SAD effect', log=log,
           frame.plot=F, ...)
      abline(h = 0, lwd = 1, lty = 2)
      
    }       
    
  }
  
  # panels showing both deta S on one graph ---------
  minN = min(nrow(sample_rare_trt), nrow(sample_rare_ref))
  delta_Sspat = sample_rare_trt$expl_S[1:minN] - 
    sample_rare_ref$expl_S[1:minN]
  delta_Ssample = sample_rare_trt$impl_S[1:minN] - 
    sample_rare_ref$impl_S[1:minN]
  delta_Sind = mob_out$indiv_rare[[trt_group]] - 
    mob_out$indiv_rare[[ref_group]]
  # we need to interpolate a version of 
  # delta_Ssample so that differences with delta_Sind can be made
  plot_dens = mob_out$density_stat$plot_dens
  virt_effort = mob_out$indiv_rare$sample
  delta_Ssample_interp = pracma::pchip(1:minN * plot_dens, delta_Ssample,
                                       virt_effort)    
  
  ##
  plot(seq(minN), delta_Sspat, 
       ylim = range(delta_Sspat, 0),
       cex.axis = 1.5, cex.lab = 1.5, type = 'n', lwd = lwd,
       col = cols$deltaS, xlab = 'Number of plots',
       ylab = expression(Delta * 'S'), frame.plot=F)
  polygon(c(seq(minN), rev(seq(minN))), c(delta_Sspat, rev(delta_Ssample)),
          col = cols$ddeltaS, border=NA)
  abline(h = 0, lwd = 1, lty = 2)
  lines(seq(minN), delta_Sspat, lwd = lwd, col = cols$deltaS)
  lines(seq(minN), delta_Ssample, lwd = lwd, col = cols$deltaS, lty = 2)
  ##
  N_indiv = mob_out$indiv_rare$sample
  plot(N_indiv , delta_Sind, 
       ylim = range(0, delta_Ssample, delta_Sind),
       cex.axis = 1.5, cex.lab = 1.5, type = 'n', lwd = lwd,
       col = cols$deltaS, xlab = 'Number of plots',
       ylab = expression(Delta * 'S'), frame.plot=F)
  polygon(c(N_indiv , rev(N_indiv)), 
          c(delta_Ssample_interp, rev(delta_Sind)),
          col = cols$ddeltaS, border=NA)
  abline(h = 0, lwd = 1, lty = 2)
  lines(N_indiv, delta_Ssample_interp, lwd = lwd, col = cols$deltaS, lty = 2)
  lines(N_indiv, delta_Sind, lwd = lwd, col = cols$deltaS, lty = 4)
  ## null plot
  plot(1:10, 1:10, type='n', axes=F, frame.plot=F, xlab='', ylab='')
  
  if ('ddelta S' %in% display) {
    # Create the plots for the three ddelta S
    ylim_ddelta = range(lapply(mob_out[tests], function(x)
      lapply(x[ , -(1:2)], function(y)
        as.numeric(as.character(y)))), na.rm=T)
    if ('agg' %in% mob_out$tests) {
      mob_out$agg[, -1] = lapply(mob_out$agg[, -1], function(x)
        as.numeric(as.character(x))) 
      ddelta_Sspat = mob_out$agg[which(as.character(mob_out$agg$group) == as.character(trt_group)), ]
      if (!same_scale)
        ylim = range(ddelta_Sspat[ , -(1:2)])
      plot(ddelta_Sspat$effort_sample, ddelta_Sspat$ddeltaS_emp,
           ylim = ylim_ddelta, log='',
           cex.axis = 1.5, cex.lab = 1.5, type = 'n', 
           xlab = 'Number of plots', 
           ylab = expression(Delta * 'S'),
           main = 'Aggr. effect', frame.plot=F, ...)
      polygon(c(ddelta_Sspat$effort_sample,
                rev(ddelta_Sspat$effort_sample)), 
              c(ddelta_Sspat$ddeltaS_null_low,
                rev(ddelta_Sspat$ddeltaS_null_high)),
              col = '#C1CDCD', border = NA)
      abline(h = 0, lwd = 1, lty = 2)
      lines(ddelta_Sspat$effort_sample, ddelta_Sspat$ddeltaS_emp, 
            lwd = lwd, col = cols$ddeltaSline)
    }
    if ('N' %in% mob_out$tests) {
      mob_out$N[, -1] = lapply(mob_out$N[, -1], function(x)
        as.numeric(as.character(x))) 
      ddelta_Ssample = mob_out$N[which(as.character(mob_out$N$group) == as.character(trt_group)), ]
      if (!same_scale)
        ylim = range(ddelta_Ssample[ , -(1:2)])
      plot(ddelta_Ssample$effort_sample, ddelta_Ssample$ddeltaS_emp,
           ylim = ylim_ddelta, log=log,
           cex.axis = 1.5, cex.lab = 1.5, type = 'n', 
           xlab = 'Number of individuals', 
           ylab = expression(Delta * 'S'),
           main = 'N effect', frame.plot=F, ...)
      polygon(c(ddelta_Ssample$effort_sample, 
                rev(ddelta_Ssample$effort_sample)), 
              c(ddelta_Ssample$ddeltaS_null_low, 
                rev(ddelta_Ssample$ddeltaS_null_high)),
              col = '#C1CDCD', border = NA)
      abline(h = 0, lwd = 1, lty = 2)
      lines(ddelta_Ssample$effort_sample, ddelta_Ssample$ddeltaS_emp,
            lwd = lwd, col = cols$ddeltaSline)
    }
    if ('SAD' %in% mob_out$tests) {
      mob_out$ind[, -1] = lapply(mob_out$ind[, -1], function(x)
        as.numeric(as.character(x))) 
      delta_Sind = mob_out$SAD[which(as.character(mob_out$SAD$group) == as.character(trt_group)), ]
      if (!same_scale)
        ylim = range(delta_Sind[ , -(1:2)])
      plot(delta_Sind$effort_ind, delta_Sind$deltaS_emp, 
           ylim = ylim_ddelta, log=log,
           cex.axis = 1.5, cex.lab = 1.5, type = 'n',
           xlab = 'Number of individuals', ylab = expression(Delta * 'S'),
           main = 'SAD effect', frame.plot=F, ...)
      polygon(c(delta_Sind$effort_ind, rev(delta_Sind$effort_ind)), 
              c(delta_Sind$deltaS_null_low, rev(delta_Sind$deltaS_null_high)),
              col = '#C1CDCD', border = NA)
      abline(h = 0, lwd = 1, lty = 2)
      lines(delta_Sind$effort_ind, delta_Sind$deltaS_emp,
            lwd = lwd, col = cols$ddeltaSline)
    }         
    
  }
}
pdf('./figs/inv_delta_stats_conceptual.pdf')

plot_delta_con(inv_mob_out, 'invaded', 'uninvaded', same_scale = T)

dev.off()
