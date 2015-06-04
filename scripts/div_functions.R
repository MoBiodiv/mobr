expS_negbin = function(rsad, n_indiv, k){
  ## Expected number of species from negative binomial, 
  ## assuming that there is a common aggregation parameter k
  ## that does not change with scale.
  ## Arguments: 
  ## rsad: relative species abundance distribution (p_i's)
  ## n_indiv: how many individuals are sampled
  ## k: aggregation parameter (e.g., see Green & Plotkin 2007)
  ## k cannot be between [-max(sad) * n_indiv, 0].
  ## Returns: 
  ## the average expected number of species under the negative binomial
  ## distribution for a sample size of n_indiv
  S_0 = length(rsad)
  S_n = sapply(n_indiv, function(n) S_0 - sum((k / (rsad * n + k)) ^ k))
  return (S_n)
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

get_sad = function(list_sp){
  ## Return the list of abundances, ranked from the most abundant to the least abundant,
  ## given a list of individuals with species names.
  abd_list = as.numeric(table(list_sp))
  return(sort(abd_list, decreasing = T))
}

force_S = function(sad, newS){
  ## Force the richness to a new value newS, without changing the shape of the SAD.
  ## Here it is assumed that the SAD is a Poisson lognormal and the parameters are MLEs.
  ## Arguments:
  ## sad: a list of species abundances
  ## newS: desirable new level of richness
  ## Returns:
  ## a list of relative abundances of length newS coming from the same Poisson lognormal distribution.
  library(poilog)
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
