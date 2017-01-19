
avg_perm_rare = function(comm, swap, groups=NULL, nperm=1000, effort=NULL){
    S = replicate(nperm, 
                  rarefaction(permute_comm(comm, swap, groups),
                              'samp', effort))
    Savg = apply(S, 1, mean)
    return(Savg)
}



S = 10
nplots = 10

comm_trt = matrix(round(rlnorm(S*nplots, 0, 1)), nplots, S)
comm_cnt = matrix(round(rlnorm(S*nplots, 1, 1)), nplots, S)

comm = rbind(comm_trt, comm_cnt)

sad_trt = colSums(comm[1:nplots, ])
sad_cnt = colSums(comm[(nplots + 1):(2 * nplots), ])

inds = 1:min(sum(sad_trt), sum(sad_cnt))

group_dens_trt = sum(sad_trt) / nplots
group_dens_cnt = sum(sad_cnt) / nplots
samp_effort_trt = round((1:nplots) * group_dens_trt) #essentially plot sampling
samp_effort_cnt = round((1:nplots) * group_dens_cnt)

S_ns_trt = rarefaction(comm_trt, 'indiv', samp_effort_trt)
S_ns_cnt = rarefaction(comm_cnt, 'indiv', samp_effort_cnt)

S_ind_trt = rarefaction(comm_trt, 'indiv', inds)
S_ind_cnt = rarefaction(comm_cnt, 'indiv', inds)

S_samp_trt = rarefaction(comm_trt, 'samp')

S_nsp_trt = avg_perm_rare(comm_trt, 'noagg')

rescaled_effort = round(1:nplots * sum(sad_trt, sad_cnt) / (2 * nplots))

interp_S_trt = pracma::pchip(c(1, rescaled_effort), c(1, S_ns_trt), inds)
interp_S_cnt = pracma::pchip(c(1, rescaled_effort), c(1, S_ns_cnt), inds)

plot(samp_effort_trt, S_ns_trt, type='o', col='red', lwd=2, log='xy',
     xlim=range(1, samp_effort_trt, samp_effort_cnt),
     ylim=range(1, S_ns_trt, S_ns_cnt))
lines(samp_effort_cnt, S_ns_cnt, type='o', col='blue')
lines(S_ind_trt, col='darkred')
lines(S_ind_cnt, col='darkblue')
lines(inds, interp_S_trt, col='red', lty=2, lwd=2)
lines(inds, interp_S_cnt, col='blue', lty=2, lwd=2)

#
plot(samp_effort_trt, S_ns_trt, type='o', col='red', lwd=2, log='',
     xlim=range(1, samp_effort_trt),
     ylim=range(1, S_ns_trt))
lines(S_ind_trt, col='darkred')
lines(inds, interp_S_trt, col='red', lty=2, lwd=2)
lines(samp_effort_trt, S_samp_trt, col='blue')
lines(samp_effort_trt, S_nsp_trt, col='pink')
#
plot(1:nplots, S_samp_trt, col='blue', type='o')
lines(1:nplots, S_nsp_trt, col='red')

####-------------------------
S = 30
nplots = 40

# two diff comms so that there is a within plot agg effect 
# such that sample based rarefaction and nonspatial rarefaction
# will differ from one another 
comm1 = matrix(round(rlnorm(S * nplots, -4, 1)), nplots / 2, S)
comm2 = matrix(round(rlnorm(S * nplots, -1, 1)), nplots / 2, S)
comm = rbind(comm1, comm2)

sad = colSums(comm)

inds = 1:sum(sad)

group_dens = sum(sad) / nplots
samp_effort = round((1:nplots) * group_dens) #essentially plot sampling

S_samp = rarefaction(comm, 'samp')               # sample-based raref
S_nsp = avg_perm_rare(comm, 'noagg', nperm=1000) # sample-based perm raref
S_ns = rarefaction(comm, 'indiv', samp_effort)   # indiv-based raref plot scale
S_ind = rarefaction(comm, 'indiv', inds)         # indiv-based raref ind scale

# consider situation in which trt of focus is high density
ref_dens_hi = group_dens - 1
# consider situtaion in which trt of focus is low density
ref_dens_lo = group_dens + 1
rescaled_effort_hi = round(1:nplots * ref_dens_hi)
rescaled_effort_lo = round(1:nplots * ref_dens_lo)
interp_S_hi_den = pracma::pchip(c(1, rescaled_effort_hi), c(1, S_ns), inds)
interp_S_lo_den = pracma::pchip(c(1, rescaled_effort_lo), c(1, S_ns), inds)

# 

interp_S_hi_den2 = pracma::pchip(inds * (ref_dens_hi / group_dens), S_ind, inds)
interp_S_lo_den2 = pracma::pchip(inds * (ref_dens_lo / group_dens), S_ind, inds)


plot(inds, S_ind, type='l')
#points(inds, interp_S_hi_den, col='red')
#points(inds, interp_S_lo_den, col='blue')
#points(inds * (ref_dens_hi / group_dens) , S_ind, col='pink')
#points(inds * (ref_dens_lo / group_dens) , S_ind, col='purple')
lines(inds, interp_S_hi_den2, col='red')
lines(inds, interp_S_lo_den2, col='blue')

# demonstrate that the non-spatial rarefaction curve
# is different from sample based rarefaction and
# that it is identical to individual based rarefaction
# that is is on the sampling scale of plots 
plot(1:nplots, S_samp, col='blue', type='o')
points(1:nplots, S_nsp, col='red')
lines(1:nplots, S_ns, col='green3') 
# add the full indiv rarefaction curve to demonstrate
# that there shouldn't be a need for interpoloation because
# we have the non-spatial curve at any value we would want to
# evaulate it at. 
par(new=T)
plot(inds, S_ind, col='orange', axes=F, xlab='', ylab='')
        
# We need to compute the difference between the individual 
# based curve and the nonspatial curve 






plot(samp_effort, S_ns, type='o', col='red', lwd=2, log='',
     xlim=range(1, samp_effort), ylim=range(1, S_ns))
lines(S_ind, col='darkred')
lines(inds, interp_S, col='red', lty=2, lwd=2)
lines(samp_effort, S_samp, col='blue')
lines(samp_effort, S_nsp, col='pink')

