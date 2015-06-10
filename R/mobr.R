#' @title Calcualte site richness and abundance
#' 
#' @description 
#' This function calculates total abundance (N), richness (S), 
#' and rarified richness (S_rare to the smallest number of individuals 
#' across all sites
#' @param x a data.frame, a matrix, or a community matrix that has 
#' the class 'indiv_occ', 'site_occ', 'site_by_sp' associated with it
#' @return 
#' A data.frame with each sites summary values
#' @details 
#' rarefied richness is computed under an assumption of random sampling
#' without replacement (i.e., the hypergeometric distribution)
#' @examples 
#' nindiv= 100
#' nsites = 10
#' stems = data.frame(site = rep(1:nsites, length.out=nindiv),
#'                    sp = sample(nrow, replace = T))
#' class(stem_occ) = c('data.frame', 'indiv_occ')
#' site_summary(stem_occ)
#' @export
site_summary = function(x) {
    require(vegan)
    if (!any(class(x) %in% c('indiv_occ', 'site_occ', 'site_by_sp'))) {
        stop('x must be either of class indiv_occ, stem_occ, or site_by_sp')
    }
    if (any(class(x) == 'indiv_occ')) {
        sites = sort(unique(x$site))
        N = tapply(quad_abu$abu, list(quad_abu$site),sum)
        S = tapply(x$sp, list(x$site), function(y) length(unique(y)))
        sads = tapply(stem_occ$sp, list(stem_occ$site), get_sad)
        S_rare = sapply(sads, rarefy, min(N))
    }
    else if (any(class(x) == 'site_occ')) {
        sites = sort(unique(x$site))
        N = tapply(x$abu, list(x$site), sum)
        S = tapply(x$sp, list(x$site), function(y) length(y))
        S_rare = tapply(x$abu, list(x$site), rarefy, min(N))
    }
    else if (any(class(x) == 'site_by_sp')) {
        sites = row.names(x)
        N = apply(x, 1, sum)
        S = apply(x > 0, 1, sum)
        S_rare = rarefy(x, min(N))
    }
    out = data.frame(site=sites, N, S, S_rare)
    return(out)
}

#' @title Expected number of species from negative binomial
#' 
#' @description
#' assuming that there is a common aggregation parameter k
#' that does not change with scale.
#' 
#' @param rsad relative species abundance distribution (p_i's)
#' @param n_indiv how many individuals are sampled
#' @param k aggregation parameter (e.g., see Green & Plotkin 2007)
#' k cannot be between [-max(sad) * n_indiv, 0].
#' 
#' @return the average expected number of species under the
#' negative binomial distribution for a sample size of n_indiv
#' 
#' @examples
#' sad = rpois(100, 5)
#' rsad = sad / sum(sad)
#' expS_negbin(rsad, 1:10, 0.5)
#' @export
expS_negbin = function(rsad, n_indiv, k){
    S_0 = length(rsad)
    S_n = sapply(n_indiv, function(n) S_0 - sum((k / (rsad * n + k)) ^ k))
    return (S_n)
}


#' @title rank abundance distribution
#' 
#' @param indiv_ids a vector of species names were each element of the vector
#'        represents a different individual
#' @examples
#' indivs = c('A', 'A', 'B', 'C', 'C', 'C')
#' get_sad(indivs)
#' @export
get_sad = function(indiv_ids){
    # is there a good reason to return a ranked vector of abunances. my pref
    # is for non-ranked sad b/c its less likely to confuse a user 
    # I think we should rename this one b/c sad is short for species-abundance distribution
    # which I don't think many people assume is a rank abundance-distribution
    abu_list = as.numeric(table(indiv_ids))
    return(sort(abu_list, decreasing = T))
}

force_S = function(sad, newS){
  ## Force the richness to a new value newS, without changing the shape of the SAD.
  ## Here it is assumed that the SAD is a Poisson lognormal and the parameters are MLEs.
  ## Arguments:
  ## sad: a list of species abundances
  ## newS: desirable new level of richness
  ## Returns:
  ## a list of relative abundances of length newS coming from the same Poisson lognormal distribution.
  require(poilog)
  pars = as.numeric(poilogMLE(sad, startVals = c(mu = mean(log(sad)), sig = sd(log(sad))))$par)
  newsad = rpoilog(newS, pars[1], pars[2])
  while(length(newsad) < newS){
    newsp = rpoilog(1, pars[1], pars[2], keep0 = T)
    if (newsp != 0){
      newsad = c(newsad, newsp)
    }
  }
  return(newsad)
}

near_neigh_ind = function(data, nperm=20){
    # The input data has three columns: x, y, and species ID for each individual.
    N = nrow(data)
    S = matrix(0, nrow=nperm, ncol=N)
    for (i in 1:nperm) {
        data = data[sample(1:dim(data)[1]), ]
        focal_row = sample(dim(data)[1], 1)
        # Compute Euclidean distances
        x_diff = data[, 1] - as.numeric(data[focal_row, 1])
        y_diff = data[, 2] - as.numeric(data[focal_row, 2])
        dist_row = sqrt(x_diff^2 + y_diff^2)
        data_order = data[order(dist_row), ]
        #vec_list = lapply(1:dim(data_order)[1], seq)
        #lapply(vec_list, length(unique(data_order[vec_list, 3])))
        for (n in 1:N) {
            sp_id_list = data_order[1:n, 3]
            n_rich = length(unique(sp_id_list))
            S[i, n] = n_rich
        }
    }
    S = apply(S, 2, mean)
    return(S)
}

near_neigh_quadrat = function(data){
  # The input data is a data frame, with the first three columns being 
  # quadrat ID, x, and y. Column 4 and beyong are the abundances for 
  # species in each quadrat. 
  data = data[sample(1:dim(data)[1]), ]
  data_spec = data[, -(1:3)]
  pair_dist = as.matrix(dist(data[, 2:3]))
  focal_row = sample(dim(data)[1], 1)
  dist_row = pair_dist[focal_row, ]
  data_order = data_spec[order(dist_row), ]
  data_bool = as.data.frame(ifelse(data_order[, 1:dim(data_order)[2]] == 0, 1, 0))
  data_rich = cumprod(data_bool)
  S = as.numeric(dim(data_spec)[2] - rowSums(data_rich))
  N = as.numeric(cumsum(rowSums(data_order)))
  return(list(S = S, N = N))
}

diff_rarefy = function(sad_extend, label_extend){
  ## Sub-function for rare_ind() to create the rarefaction curves for the 
  ## two treatments, and return the difference.
  trtmt_list = sort(unique(label_extend))
  n1 = length(label_extend[label_extend == trtmt_list[1]])
  min_n = min(n1, length(label_extend) - n1)
  
  list_of_rarefied_sads = c()
  for (trtmt in trtmt_list){
    sad_trtmt = sad_extend[label_extend == trtmt]
    sad_numeric = as.numeric(table(sad_trtmt))
    sad_rarefy = as.numeric(rarefy(sad_numeric, 1:min_n))
    list_of_rarefied_sads = c(list_of_rarefied_sads, list(sad_rarefy))
  }
  delta_s = unlist(list_of_rarefied_sads[1]) - unlist(list_of_rarefied_sads[2])
  return(delta_s)
}

avg_rarefy = function(sad_extend, sample_extend, min_n){
  ## Sub-function for rare_ind_avg() to create average rarefaction curves among
  ## samples with the same treatment.
  sample_list = sort(unique(sample_extend))
  list_of_rarefied_sads = as.data.frame(matrix(nrow = length(sample_list), 
                                               ncol = min_n))
  for (i in 1:length(sample_list)){
    sample = sample_list[i]
    sad_sample = sad_extend[sample_extend == sample]
    sad_numeric = as.numeric(table(sad_sample))
    sad_rarefy = as.numeric(rarefy(sad_numeric, 1:min_n))
    list_of_rarefied_sads[i, ] = sad_rarefy
  }
  avg_s = apply(list_of_rarefied_sads, 2, mean)
  return(avg_s)
}

rare_ind = function(table_of_sads, nperm = 100){
  ## Rarefy on the individual basis for SADs from two treatment groups.
  ## Arguments:
  ## table_of_sads: a data frame, with each row being the SAD for one sample from a 
  ## treatment. The first column is treatment label, and the subsequent columns are 
  ## the abundances of specific species in the sample. Zeros are allowed.
  ## nperm: number of permutations
  library(vegan)
  list_of_labels = as.vector(table_of_sads[, 1])
  sad_extend = c()
  label_extend = c()
  for (label in unique(list_of_labels)){
    dat_label = table_of_sads[which(list_of_labels == label), 2:dim(table_of_sads)[2]]
    sad_label = as.numeric(apply(dat_label, 2, sum))
    sad_extend = c(sad_extend, rep(1:length(sad_label), sad_label))
    label_extend = c(label_extend, rep(label, sum(sad_label)))
  }
  # Do the rarefaction
  delta_s_orig = diff_rarefy(sad_extend, label_extend)
  delta_s_perm = as.data.frame(matrix(nrow = nperm, ncol = length(delta_s_orig)))
  for (i in 1:nperm){
    print(i)
    label_new = sample(label_extend, length(label_extend))
    delta_s_perm[i, ] = diff_rarefy(sad_extend, label_new)
  }
  
  quant95 = apply(delta_s_perm, 2, quantile, c(.025, .975))
  delta_s_comb = as.data.frame(matrix(nrow = length(delta_s_orig), ncol = 4))
  delta_s_comb[, 1] = 1:length(delta_s_orig)
  delta_s_comb[, 2] = delta_s_orig
  delta_s_comb[, 3] = as.numeric(quant95[1, ])
  delta_s_comb[, 4] = as.numeric(quant95[2, ])
  names(delta_s_comb) = c('N', 'delta_S', 'null_lo', 'null_hi')
  return(delta_s_comb)
}

rare_ind_avg = function(table_of_sads, nperm = 100){
  ## Rarefy on the individual basis for SADs from two treatment groups. 
  ## Unlike rare_ind(), here the rarefaction curves are averaged among samples 
  ## within treatment groups, thus delta-S only goes to Nmin, where Nmin is the 
  ## lowest totally abundance among samples across treatments.
  library(vegan)
  list_of_labels = as.vector(table_of_sads[, 1])
  trtmt_rows = list(which(list_of_labels==unique(list_of_labels)[1]), 
                    which(list_of_labels==unique(list_of_labels)[2]))
  sample_n = apply(table_of_sads[, 2:dim(table_of_sads)[2]], 1, sum)
  min_n = min(sample_n)
  
  avg_S = as.data.frame(matrix(0, nrow = 2, ncol = min_n))
  sample_extend_full = c()
  sad_extend_full = c()
  for (i in 1:2){
    sad_rows = table_of_sads[unlist(trtmt_rows[i]), 2:dim(table_of_sads)[2]]
    sad_extend = c()
    sample_extend = c()
    for (j in 1:dim(sad_rows)[1]){
      sad_row = sad_rows[j, ]
      sad_extend = c(sad_extend, rep(1:length(sad_row), sad_row))
      sample_extend = c(sample_extend, rep(unlist(trtmt_rows[i])[j], sum(sad_row)))
    }
    avg_S[i, ] = avg_rarefy(sad_extend, sample_extend, min_n)
    sample_extend_full = c(sample_extend_full, sample_extend)
    sad_extend_full = c(sad_extend_full, sad_extend)
  }
  delta_s_orig = avg_S[1, ] - avg_S[2, ]
  
  delta_s_perm = as.data.frame(matrix(nrow = nperm, ncol = length(delta_s_orig)))
  for (i in 1:nperm){
    print(i)
    avg_S = as.data.frame(matrix(0, nrow = 2, ncol = min_n))
    sample_full_new = sample(sample_extend_full, length(sample_extend_full))
    for (j in 1:2){
      sample_new = sample_full_new[which(sample_full_new %in% unlist(trtmt_rows[j]))]
      sad_new = sad_extend_full[which(sample_full_new %in% unlist(trtmt_rows[j]))]
      avg_S[j, ] = avg_rarefy(sad_new, sample_new, min_n)
    }
    delta_s_perm[i, ] = avg_S[1, ] - avg_S[2, ]
  }
  
  quant95 = apply(delta_s_perm, 2, quantile, c(.025, .975))
  delta_s_comb = as.data.frame(matrix(nrow = length(delta_s_orig), ncol = 4))
  delta_s_comb[, 1] = 1:length(delta_s_orig)
  delta_s_comb[, 2] = as.numeric(delta_s_orig)
  delta_s_comb[, 3] = as.numeric(quant95[1, ])
  delta_s_comb[, 4] = as.numeric(quant95[2, ])
  names(delta_s_comb) = c('N', 'delta_S', 'null_lo', 'null_hi')
  return(delta_s_comb)
}

get_acc_avg = function(pooled_sads) {
    # pooled_sads = matrix of sads each row is a site, each column a sp
    acc = apply(pooled_sads, 1, function(x) rarefy(x, 1:sum(x)))
    # average rarefaction curves within a treatment
    min_N = min(sapply(acc, length))
    acc_std = t(sapply(acc, function(x) x[1:min_N]))
    acc_avg = aggregate(acc_std, list(rownames(acc_std)), mean)
    return(acc_avg)
}

get_acc_delta = function(pooled_sads) {
    acc = get_acc_avg(pooled_sads)
    delta = acc[1, -1] - acc[2, -1] 
    return(as.numeric(delta))
}


perm_labels = function(pooled_sads, nperm=2, quantiles=c(.025, .975)) {
    # pooled_sads = list of sads each element of the list is a
    # different treatment
    obs_delta = get_acc_delta(pooled_sads)
    swap = function(x){ rownames(x) = sample(rownames(x)) ; return(x)}
    perms = replicate(nperm, get_acc_delta(swap(pooled_sads)))
    null_qt = apply(perms, 1, quantile, quantiles)
    out = data.frame(N = 1:length(obs_delta), delta_S = obs_delta,
                     null_lo = null_qt[1, ], null_hi = null_qt[2, ])
    return(out)
}

addCI = function(x, y_lo, y_hi, col, data=NULL) {
    ## if data is not null then all of the arguments
    ## including data itself must be text strings
    if (!is.null(data)) {
        x = eval(parse(text=paste(data, '$', x, sep='')))
        y_lo = eval(parse(text=paste(data, '$', y_lo, sep='')))
        y_hi = eval(parse(text=paste(data, '$', y_hi, sep='')))
    }  
    xvals = c(x, rev(x))
    yvals = c(y_lo, rev(y_hi))
    polygon(xvals, yvals, border=NA, col=col)
}

mat2psp = function(sp_mat, xy_coord, N=NULL, M=NULL)
{
    ##place site by species matrix (x) into an S x N x M array where N >= M
    ##a multidimensional array that eases computing 
    ##replaces the old function 'grid.pres'
    ##Note: if area rectangular the first dimension of xy_coord does NOT
    ## have to be larger than the second dimension
    if (is.null(N) )
        N = sqrt(nrow(sp_mat))
    if (is.null(M) )
        M = N
    if (nrow(sp_mat) != nrow(xy_coord))
        stop('Number of samples in species matrix must match number of samples in xy-coordinates')
    if (N < M)
        stop('N should be >= M')
    if (N * M != nrow(sp_mat))
        stop('Number of specified samples (N X M) must be equal to the number of samples in sp_mat')
    S = ncol(sp_mat)
    ## order rows of sp_mat based on xy_coords so that they are placed
    ## in the multidimensional array in the correct arrangement
    proper_order = order(xy_coord[ , 2], xy_coord[ , 1])
    sp_mat = sp_mat[proper_order, ]
    psp = array(sp_mat, dim=c(N, M, S))
    psp = aperm(psp, c(3, 1, 2))
    spSums = apply(psp, 1, sum)
    ## drop species that never occur
    if(any(spSums %in% 0))
        psp = psp[spSums > 0, , ]
    return(psp)
}

psp2mat = function(psp) {
    out = apply(psp, 1, as.vector)
    return(out)
}

coverage = function(abu) {
    # from Chao and Jost 2012 equ 4a
    f1 = sum(abu == 1)
    f2 = sum(abu == 2)
    n = sum(abu)
    1 - f1 / n * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
}

getDAR = function(psp, grains, mv_window=FALSE, avg=TRUE) {
    ## Purpose: to construct spatially explict diversity area based upon a
    ## mapped grid of occurances
    ## Arguments:
    ## psp: community array (i.e., S x N x M abundance array where N >= M)
    ## grains: the areas in pixels for which to compute the SAR
    ##         only grains that have integer log base 2 are considered
    ## mv_window: FALSE indicates that a non-moving window SAR will be calculated
    require(vegan)
    if (class(psp) != 'array')
        stop('psp must be a community array (S X N X M)')
    grains = grains[log2(grains) == round(log2(grains))]
    ## define the size of sampling units on each side
    lenN = n_pixels_long(log2(grains))
    lenM = n_pixels_wide(log2(grains))
    # the columns of div are S, D, and J'
    div = matrix(0, ncol=5, nrow=sum(max(grains)/grains))
    S = dim(psp)[1]
    N = dim(psp)[2]
    M = dim(psp)[3]
    if (M > N) {
        stop('The first spatial dimension of psp must be larger than or equal to the second 
             (i.e. psp[S,N,M] where N >= M)')
    }
    for (l in seq_along(grains)) {
        if (grains[l] == 1) {  # if area=1
            row_max = max(grains)/grains[1]
            # richness
            div[1:row_max, 1] = apply(psp > 0, 2:3, sum) # if S = 1 may not work
            # abundance
            div[1:row_max, 2] = apply(psp, 2:3, sum)
            # shannons
            div[1:row_max, 3] = diversity(psp2mat(psp), index = 'shannon')
            # simpsons
            div[1:row_max, 4] = 1 - diversity(psp2mat(psp), index = 'simpson')
            # coverage
            div[1:row_max, 5] = sapply(apply(psp2mat(psp), 1, sum), coverage)
            icount = row_max + 1
        }
        else{
            if (mv_window) {
                brksN = 1:(N - lenN[l] + 1)
                brksM = 1:(M - lenM[l] + 1)
            }
            else{
                brksN = seq(1, N, lenN[l])
                brksM = seq(1, M, lenM[l])
            }  
            for (n in brksN) {
                for (m in brksM) {
                    psp_tmp = psp[ , n:(n + (lenN[l] - 1)),
                                     m:(m + (lenM[l] - 1))]
                    #if (S == 1) {
                    #    div[row_indices, 1] = any(psp_tmp > 0) * 1
                    #    div[row_indices, 2] = sum(psp_tmp)
                    #    div[row_indices, 3] = 0
                    #    div[row_indices, 4] = 1
                    #    div[row_indices, 5] = NaN
                    #}
                    #else {
                        # richness
                        div[icount, 1] = sum(apply(psp_tmp > 0, 1, sum) > 0)
                        # abuncande
                        div[icount, 2] = sum(psp_tmp)
                        abu = apply(psp_tmp, 1, sum)
                        # shannonn
                        div[icount, 3] = diversity(abu, index = 'shannon')
                        # simpsons
                        div[icount, 4] = 1 - diversity(abu, index = 'simpson')
                        # coverage
                        div[icount, 5] = coverage(abu)
                        
                        #}
                    icount = icount + 1
                }
            }
        }
    }
    grains = rep(grains, times=max(grains)/grains)
    out = cbind(grains, div)  
    colnames(out) = c('grains', 'richness', 'indiv', 'shannon', 'simpson', 'coverage')
    if (avg) {
        out = aggregate(dat[ , -(1:2)], by=list(dat$grains, dat$site), mean, na.rm=T)
        names(out) = c('grain', 'site', names(out)[-(1:2)])
    }
    return(out)
}  

getEAR = function(psp, grains, mv_window=FALSE)
{
    ## Purpose: to construct spatially explict Endemics Area Relationship
    ## based upon a mapped grid of occurances
    ## Arguments:
    ## psp: community array (i.e., S x N x M abundance array where N >= M)
    ## grains: the areas in pixels for which to compute the SAR
    ##         only grains that have integer log base 2 are considered
    ## mv_window: FALSE indicates that a non-moving window SAR will be calculated
    stop('Error function has not been tested!')
    if (class(psp) != 'array')
        stop('psp must be a community array (S X N X M)')
    grains = grains[log2(grains) == round(log2(grains))]
    ## define the size of sampling units on each side
    lenN = n_pixels_long(log2(grains))
    lenM = n_pixels_wide(log2(grains))
    sr = rep(0, length(grains))
    ind = rep(0, length(grains))
    std = rep(0, length(grains))
    cs = rep(0, length(grains))
    S = dim(psp)[1]
    N = dim(psp)[2]
    M = dim(psp)[3]
    if (M > N) {
        stop('The first spatial dimension of psp must be larger than or equal to the second 
             (i.e. psp[S,N,M] where N >= M)')
    }       
    for (l in seq_along(grains)) {
        if (mv_window) {
            brksN = 1:(N - lenN[l] + 1)
            brksM = 1:(M - lenM[l] + 1)
        }
        else{
            brksN = seq(1, N, lenN[l])
            brksM = seq(1, M, lenM[l])
        }  
        sr_vec = NULL
        for (n in brksN) {
            for (m in brksM) {
                n_indices = n:(n + (lenN[l] - 1))
                m_indices = m:(m + (lenM[l] - 1))
                psp_in = psp[ , n_indices, m_indices]
                psp_out = psp[ , -n_indices, -m_indices]
                
                if (S == 1) {
                    sr_vec = c(sr_vec, any(psp_tmp > 0) * 1)
                    ind[l] = ind[l] + sum(psp_tmp)
                }
                else {
                    
                    sr_vec = c(sr_vec, sum(apply(psp_tmp > 0, 1, sum) > 0))
                    ind[l] = ind[l] + sum(apply(psp_tmp, 1, sum))
                }
                cs[l] = cs[l] + 1
            }
        }
        sr[l] = sum(sr_vec)
        std[l] = sd(sr_vec)
    }
    out = cbind(grains, sr / cs, ind / cs, cs, std)  
    colnames(out) = c('grains', 'richness', 'indiv', 'count', 'sr_std')
    return(out)
    }  


##3.9##
aggr_comm_matrix = function(mat, coords, bisect, grains, binary=FALSE){
    ## Purpose: This function generates aggregated community matrices for each
    ## spatial grain that is specified. This function is only approrpriate for 
    ## data from a regular square spatial grid.  Coordinates, bisect, and 
    ## grains will be generated if not supplied
    ## Inputs:
    ## mat: site x species matrix assumed to be sampled from a square grid
    ## coords: two column spatial x and y coordinates
    ## bisect: specifies the levels of bisection that the community should be 
    ##   aggregated at
    ## grains: 
    ## binary: boolean, if TRUE binary pres/abse matrix returned
    S = ncol(mat)
    N = nrow(mat)
    bisect_start = log2(N)
    if (missing(bisect)) {
        bisect = log2(bisect_start) : 2
    } 
    if (any(bisect > bisect_start)) {
        stop(paste('The number of bisections must be less than ',
                   bisect_start, sep=''))  
    }
    if (missing(coords)) {
        coords = expand.grid(1 : n_pixels_long(bisect_start),
                             1 : n_pixels_wide(bisect_start))
    }
    D = dist(coords)
    minD = min(D)
    domain = c(min(coords[ , 1]), max(coords[ , 1]) + minD,
               min(coords[ , 2]), max(coords[ , 2]) + minD)
    xdiff = abs(domain[1] - domain[2]) 
    ydiff = abs(domain[3] - domain[4]) 
    
    if (xdiff > ydiff) {
        xlengths = xdiff / n_pixels_long(bisect)
        ylengths = ydiff / n_pixels_wide(bisect)
    }
    else if (xdiff < ydiff) {
        xlengths = xdiff / n_pixels_wide(bisect)
        ylengths = ydiff / n_pixels_long(bisect)
    }
    else if (xdiff == ydiff) { ## extent is a sqr.
        xlengths = xdiff / sqrt(2^bisect)
        ylengths = ydiff / sqrt(2^bisect)
    }
    else
        stop('Function cannot figure out how to split up the area')
    n_quadrats = sum(2^bisect)
    if (missing(grains)) {
        grains = round(xlengths * ylengths, 2)
    }
    comms = matrix(NA, nrow=sum(n_quadrats), ncol=S + 3)
    if (is.null(colnames(mat))) {
        colnames(mat) = paste('sp', 1:ncol(mat), sep='')
    }  
    colnames(comms) = c('grain', 'x', 'y', colnames(mat))
    irow = 1
    for (i in seq_along(bisect)) {
        if (bisect[i] == bisect_start) {
            indices = irow : (irow + N - 1)
            comms[indices, 1:3] = cbind(grains[i], coords[ , 1], coords[ , 2])    
            comms[indices, -(1:3)] = as.matrix(mat)
            irow = max(indices) + 1
        }
        else {
            xbreaks = seq(domain[1], domain[2], xlengths[i])
            ybreaks = seq(domain[3], domain[4], ylengths[i]) 
            for (x in 1:(length(xbreaks) - 1)) {
                for (y in 1:(length(ybreaks) - 1)) {
                    inQuad =  xbreaks[x] <= coords[ , 1] & coords[ , 1] < xbreaks[x + 1] & 
                        ybreaks[y] <= coords[ , 2] & coords[ , 2] < ybreaks[y + 1]
                    comms[irow, 1:3] = c(grains[i], x, y)
                    comms[irow, -(1:3)] = apply(mat[inQuad, ], 2, sum)
                    irow = irow + 1 
                }
            }
        }  
    }
    if (binary)
        comms[ , -(1:3)] = as.numeric(comms[ , -(1:3)] > 0)
    return(comms)
}

n_pixels_long = function(i_bisect){
    ## returns the number of pixels on the side of a grid with more or equal pixels
    ## after i bisection events
    ## old function name: len
    2^floor((i_bisect + 1) / 2) 
}

n_pixels_wide = function(i_bisect){
    ## returns the number of pixels on the side of a grid with less or equal pixels
    ## after i bisecttion events
    ## old function name: wid
    2^floor(i_bisect / 2) 
} 

##3.19 ##
make_comm_matrix = function(spnum, S, coords, n_quadrats, domain, abu = NULL,
                            grainSuffix=NULL, rm_absent_sp=TRUE) {
    ## Output: 
    ## A community matrix where each row is a differnet pixel on a grid.  
    ## Arguments:
    ## spnum : an integer specifying species identities
    ## S : the size of the species pool may be larger than the number of unique 
    ##     spnum
    ## coords : two column matrix (x,y) specifying the spatial coordinates of each stem
    ## n_quadrats : the number of quadrats at each spatial grain
    ## domain : specifies the spatial domain of the area:  (xmin, xmax, ymin, ymax)
    ## abu: abundance associated with each record, if NULL then it is set to 1
    ##      individual per record
    ## grainSuffix : if supplied the grain column will have this appended to it
    ##               so that it is clear what community this corresponds with
    ## rm_absent_sp: boolean that defaults to TRUE to remove species
    ##  who no longer occur in the site x species matrix after subsetting base
    ## on the defined spatial domain (i.e., argument 'domain' specifies a smaller
    ## area than spatial coordiantes are provided for)
    xdiff = abs(domain[2] - domain[1])
    ydiff = abs(domain[4] - domain[3])
    if (xdiff > ydiff) {
        xlengths = xdiff / n_pixels_long(log2(n_quadrats))
        ylengths = ydiff / n_pixels_wide(log2(n_quadrats))
    }  
    else if (xdiff < ydiff) {
        xlengths = xdiff / n_pixels_wide(log2(n_quadrats))
        ylengths = ydiff / n_pixels_long(log2(n_quadrats))
    }
    else if (xdiff == ydiff) {
        xlengths = ylengths = rep(NA, length(n_quadrats))
        for (i in seq_along(n_quadrats)) {
            if (log2(n_quadrats[i]) %% 2 == 1) {
                ## if # of bisections is odd then arbitrarily 
                ## make x dimension have more cells
                xlengths[i] = xdiff / n_pixels_long(log2(n_quadrats[i]))
                ylengths[i] = ydiff / n_pixels_wide(log2(n_quadrats[i]))
            }  
            else {
                xlengths[i] = xdiff / sqrt(n_quadrats[i])
                ylengths[i] = ydiff / sqrt(n_quadrats[i])
            } 
        }  
    }
    else
        stop('Function cannot figure out how to split up the area')
    comms = matrix(NA, nrow=sum(n_quadrats), ncol=S + 3)
    colnames(comms) = c('grain', 'x', 'y', paste('sp', 1:S, sep=''))
    irow = 1
    for (i in seq_along(n_quadrats)) {
        xbreaks = seq(domain[1], domain[2], xlengths[i])
        ybreaks = seq(domain[3], domain[4], ylengths[i]) 
        for (x in 1:(length(xbreaks) - 1)) {
            for (y in 1:(length(ybreaks) - 1)) {
                inQuad =  xbreaks[x] <= coords[ , 1] & coords[ , 1] < xbreaks[x + 1] & 
                    ybreaks[y] <= coords[ , 2] & coords[ , 2] < ybreaks[y + 1]
                if (is.null(grainSuffix)) {
                    comms[irow, c(1:3)] = c(paste(round(xlengths[i] * ylengths[i], 2), sep=''),
                                            x, y)
                }
                else {
                    comms[irow, c(1:3)] = c(paste(round(xlengths[i] * ylengths[i], 2),
                                                  grainSuffix, sep=''), x, y)
                }
                if (is.null(abu) ){
                    comms[irow, -c(1:3)] = as.integer(table(c(spnum[inQuad],1:S)) - 1)
                }
                else {
                    comms[irow, -c(1:3)] =  as.integer(table(c(unlist(mapply(
                        rep, spnum[inQuad], abu[inQuad])), 1:S)) - 1)
                }  
                irow = irow + 1 
            }
        }
    }
    if (rm_absent_sp) {
        cols_to_rm = which(apply(comms[ , -(1:3)], 2, function(x) all(x == '0'))) + 3
        if (length(cols_to_rm) > 0)
            comms = comms[ , -cols_to_rm]
    }
    return(comms)
}


expS_binom = function(sad, n_indiv) {
    ## Expected number of species from Coleman (1981), Eq. 3.11
    ## Arguments:
    ## sad: species abundance distribution
    ## n_indiv: how many individuals 
    ## Returns:
    ## the average expected number of species under the binomial distr for
    ## a sample of area A out of A0 and species abundances n
    S = sapply(n_indiv, function(x) sum(1 - (1 - (sad/sum(sad))) ^ x))
    return(S)
}

dexpS_binom = function(sad, n_indiv) {
    pi = sad / sum(sad)
    S = length(sad)
    dS = sapply(n_indiv, function(x) 
        (x * sum((1 - pi)^x * log(1 - pi))) / (S - sum(1 - pi)^x))
    return(dS)
}

get_sample_stats = function(dat_in){
  ## This function calculates the characteristics of each sample.
  ## Input data have to take a specific format, where each row is the record of one individual,
  ##     and it must have the following columns (additional columns are allowed): 
  ##     sample - identifies in which sample the individual occurs
  ##     treatment - treatment associated with the sample
  ##     spcode - species code for the individual
  ## It returns a data frame where each row is a sample, with the following columns:
  ## sample: sample ID
  ## treatment: treatment associated with the sample
  ## N: total number of individuals
  ## obsS: total number of observed species
  ## PIE: 1 - sum(p_i^2) 
  ## rareS: rarefied S at min(N). This value depends on other samples.
  ## MIHS: difference between obsS and rareS. Thus the min value among the rows will always be zero. 
  library(vegan)
  samples = unique(dat_in$sample)
  dat_out = as.data.frame(matrix(nrow = length(samples), ncol = 7))
  names(dat_out) = c('sample', 'treatment', 'N', 'obsS', 'PIE', 'rareS', 'MIHS')
  dat_out$sample = samples
  dat_out$N = as.numeric(table(dat_in$sample))
  minN = min(dat_out$N)
  for (i in 1:dim(dat_out)[1]){
    dat_in_sample = dat_in[dat_in$sample == samples[i], ]
    dat_out$treatment[i] = dat_in_sample$treatment[1]
    dat_out$N[i] = dim(dat_in_sample)[1]
    dat_out$obsS[i] = length(unique(dat_in_sample$spcode))
    sp_counts = as.numeric(table(dat_in_sample$spcode))
    dat_out$PIE[i] = dat_out$N[i] / (dat_out$N[i] - 1) * (1 - sum((sp_counts / sum(sp_counts))^2))
    dat_out$rareS[i] = rarefy(sp_counts, minN)
  }
  dat_out$MIHS = dat_out$obsS - dat_out$rareS
  return(dat_out)
}

get_sample_stats_pair = function(dat_in){
  ## This function computes the characteristics of an input data frame with paired plots.
  ## Input data: a data frame with the following columns (in this order):
  ## treatment, plot, pair label (plots in the same pair have the same label), sample (within plots), 
  ##     species ID, abundance
  ## Output: a data frame with the following columns:
  ## treatment, plot, pair label, N (abundance within plot), obsS (observed richness within plot), 
  ##     PIE (probability of interspecific encounter), rareS (rarefied S at the lowest level of N), 
  ##     MIHS (obsS - rareS)
  library(vegan)
  # Remove factors in dat_in
  i = sapply(dat_in, is.factor)
  dat_in[i] = lapply(dat_in[i], as.character)
  
  unique_plots = unique(dat_in[, 1:3])
  dat_out = as.data.frame(matrix(nrow = dim(unique_plots)[1], ncol = 8))
  names(dat_out) = c('treatment', 'plot', 'pair_label', 'N', 'obsS', 'PIE', 'rareS', 'MIHS')
  Ns = aggregate(dat_in[, 6] ~ dat_in[, 1] + dat_in[, 2] + dat_in[, 3], FUN = sum)[, 4]
  Nmin = min(Ns)
  for (i in 1:dim(unique_plots)[1]){
    plot = unique_plots[i, ]
    dat_plot = merge(dat_in, plot)
    dat_out[i, 1:3] = plot
    dat_out$obsS[i] = length(unique(dat_plot[, 5]))
    dat_out$N[i] = sum(dat_plot[, 6])
    plot_sp_counts = aggregate(x = dat_plot[, 6], FUN = sum, by = list(sp = dat_plot[, 5]))[, 2]
    dat_out$PIE[i] = dat_out$N[i] / (dat_out$N[i] - 1) * (1 - sum((plot_sp_counts / sum(plot_sp_counts))^2))
    dat_out$rareS[i] = rarefy(plot_sp_counts, Nmin)
  }
  dat_out$MIHS = dat_out$obsS - dat_out$rareS
  return (dat_out)
}

initial_tests_S_N = function(dat_sample, parametric, plot){
  ## This function takes the output from get_sample_stats() and perform tests on whether samples differ 
  ## by treatments. The test will be ANoVA if treatment is categorical, linear regression if it is 
  ## continuous.
  ## Inputs: 
  ## dat_sample: output from get_sample_stats(), with 6 columns and one sample in each row.
  ## parametric: If TRUE, p-value is read directly from the output of lm(). If FALSE, p-value is 
  ##     computed by shuffling the treatment labels 1,000 times. Note that when the sample size is
  ##     low there may not be 1,000 unique permutations of the treatment labels.
  ## plot: If TRUE, a 2*3 plot is created to illustrate the comparisons of the five variables.
  ## Output:
  ## A list of the five p-values: p_N, p_obsS, p_PIE, p_rareS, p_MIHS
  p_vec = c()
  f_vec = c()
  for (i in 3:7){
    model = summary(lm(dat_sample[, i]~dat_sample[, 2]))
    f = model$fstat
    f_vec = c(f_vec, f[1])
    p_vec = c(p_vec, pf(f[1], f[2], f[3], lower = F))
  }
  if (parametric == F) { # Non-parametric test
    f_frame = as.data.frame(matrix(nrow = 1000, ncol = 5))
    for (j in 1:1000){
      dat_rand = dat_sample
      dat_rand$treatment = sample(dat_rand$treatment, length(dat_rand$treatment))
      for (i in 3:7){
        model = summary(lm(dat_rand[, i]~dat_rand[, 2]))
        f = model$fstat
        f_frame[j, i-2] = f[1]
      }
    }
    p_vec = c()
    for (i in 1:5){
      p_vec = c(p_vec, length(f_frame[f_frame[, i]>f_vec[i], i]) / 1000)
    }
  }
  if (plot == T){
    # Currently it is assumed that treatment has distinct levels.
    par(mfrow = c(2, 3))
    col_names = c('N', 'Observed S', 'PIE', 'Rarefied S', 'MIH delta-S')
    for (i in 1:5){
      boxplot(dat_sample[, i + 2] ~ dat_sample[, 2], main = col_names[i],
              las = 2)
      mtext(paste('p=', p_vec[i], sep = ''), cex = 0.8)
    }
  }
  return(list(p_N = p_vec[1], p_obs = p_vec[2], p_PIE = p_vec[3], p_rareS = p_vec[4], 
              p_MIHS = p_vec[5]))
}

initial_test_S_N_pair = function(dat_in, plot){
  ## This function takes the output from get_sample_stats_pair() and perform pairwise t tests on 
  ## whether plots differ by treatments. (I think) the treatment can only be categorical here.
  ## This function does not allow non-parametric test by shuffling the treatment labels, because the 
  ## number of possible combinations is likely to be extremely limited.
  ## Inputs: 
  ## dat_in: output from get_sample_stats_pair(), with 8 columns and one plot in each row.
  ## plot: If TRUE, a 2*3 plot is created to illustrate the comparisons of the five variables.
  ## Output:
  ## A list of the five p-values: p_N, p_obsS, p_PIE, p_rareS, p_MIHS
  dat_in = dat_in[with(dat_in, order(pair_label, treatment)), ]
  rows_trtmt1 = c()
  rows_trtmt2 = c()
  for (pair in unique(dat_in$pair_label)){
   rows_trtmt1 = c(rows_trtmt1, which(dat_in$pair_label == pair)[1])
   rows_trtmt2 = c(rows_trtmt2, which(dat_in$pair_label == pair)[2])
  }
  
  p_vec = c()
  t_vec = c()
  trtmt = unique(dat_in$treatment)
  for (i in 4:8){
    vals_1 = dat_in[rows_trtmt1, i]
    vals_2 = dat_in[rows_trtmt2, i]
    model = t.test(vals_1, vals_2, paired = T)
    t_vec = c(t_vec, as.numeric(model$stat))
    p_vec = c(p_vec, as.numeric(model$p.val))
  }
  if (plot == T){
    par(mfrow = c(2, 3))
    col_names = c('N', 'Observed S', 'PIE', 'Rarefied S', 'MIH delta-S')
    for (i in 1:5){
      boxplot(dat_in[, i + 3] ~ dat_in[, 1], main = col_names[i],
              las = 2)
      mtext(paste('p=',round(p_vec[i], digits = 6), sep = ''), cex = 0.8)
    }
  }
  return(list(p_N = p_vec[1], p_obs = p_vec[2], p_PIE = p_vec[3], p_rareS = p_vec[4], 
              p_MIHS = p_vec[5]))
}
