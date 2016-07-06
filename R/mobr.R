## A rough sketch of what we need in the module Self reminder: need two pieces of
## inputs: 1. a data frame of plot characteristics, with columns: plot,
## group, x, y 2. a data frame where rownames are plot ids, and
## subsequent columns are species abundances

## Functions
require(rareNMtests)
require(pracma)
require(plyr)
require(vegan)

# Global constants
min_ind = 5
min_plot = 5

## define mobr object
make_comm_obj = function(comm, plot_attr, binary=FALSE, ref_var=NULL, ref_group=NULL) {
    # possibly make ref_var and ref_group mandatory arguments
    out = list()
    out$analyses = list(indiv=T, sampl=T, spat= T)
    # carry out some basic checks
    if (nrow(comm) < 5) {
        stop("Number of plots in community is less than five therefore only individual rarefaction will be computed")
        out$analyses$samp = FALSE
        out$analyses$spat = FALSE
    }
    if (nrow(comm) != nrow(plot_attr))
        stop("Number of plots in community does not equal number of plots in plot attribute table")
    if (any(row.names(comm) != row.names(plot_attr)))
        warning("Row names of community and plot attributes tables do not match")
    if (binary)  {
        warning("Only spatially-explict sampled based forms of rarefaction can be computed on binary data")
        out$analyses$indiv = FALSE
        out$analyses$samp = FALSE
    } 
    else {
        if (max(comm) == 1)
            warning("Maximum abundance is 1 which suggests data is binary, change the binary argument to TRUE")
    }
    if (any(colSums(comm) == 0)) {
        warning("Some species have zero occurrences and will be dropped from the community table")
        comm = comm[colSums(comm) != 0, ]
    }
    if (!is.null(ref_var)) {
        if (is.null(ref_group)) {
            stop('The reference group (argument "ref_group") must also be provided')
        } else {
            ref_data = eval(parse(text=paste('plot_attr$', ref_var, sep='')))
            if (!(ref_group %in% ref_data))
                stop('Reference group is not present in the reference variable')
            group_cts = table(ref_data)
            if (any(group_cts < 5)) {
                warning('Some groups in the reference variable have fewer than 5 replicates and therefore will be dropped from future comparisions')
                row_indices = which(ref_data == names(group_cts)[group_cts < 5])
                comm = comm[-row_indices, ]
                plot_attr = plot_attr[-row_indices, ]
            }    
        }
    }
    out$comm = comm
    spat_cols = which(names(plot_attr) %in% c('x', 'y'))
    if (length(spat_cols) > 0) {
        out$envi = plot_attr[ , -spat_cols]
        out$spat = plot_attr[ , spat_cols]
    }
    else {
        out$analyses$spat = FALSE
        out$envi = plot_attr
        out$spat = NULL
    }
    out$ref_var = ref_var
    out$ref_group = ref_group
    class(out) = 'comm'
    return(out)
}

print.comm = function(x) {
    cat('Printing the head of each attribute in the object\n')
    print(lapply(x, head))
}

print.mobr = function(...) {
   # print rarefaction and delta rarefaction summaries
}

plot.mobr = function(...) {
  # plot rarefation and delta rarefaction curves
}

summary.mobr = function(...) {
   #  print summary anova style table
}

## Functions to obtain the three kinds of curves: sample-based explicit,
## sample-based implicit, and individual-based rescaled to sample
get_avg_dens = function(comm_obj, ref_var=NULL, ref_group=NULL) {
    # Auxillary function to obtain average density within plot for rescaling Ask user
    # to specify which group is used as 'standard' for rescaling.
    if (!is.null(ref_var)) {
        if (is.null(ref_group)) 
            stop('The reference group (argument "ref_group") must also be provided')
        else
            ref_data = eval(parse(text=paste('plot_attr$', ref_var, sep='')))
    }            
    if (!is.na(ScaleBy)) {
        avg_dens = mean(apply(dat_sp[dat_plot$group == ScaleBy, ], 1, sum))
    } else {
        # If unspecified, uses min density among groups.
        groups = unique(dat_plot$group)
        avg_dens = min(sapply(groups, function(x) 
                              mean(apply(dat_sp[dat_plot$group == x, ], 1, sum))))
    }
    return(avg_dens)
}

rarefaction = function(x, method, effort=NULL) {
    # analytical formulations from Cayuela et al. 2015. Ecological and biogeographic null hypotheses for
    # comparing rarefaction curves. Ecological Monographs 85:437–454.
    # Appendix A: http://esapubs.org/archive/mono/M085/017/appendix-A.php
    # possible inputs: sad or sp x site
    # community matrix examine input properties to determine if input is
    # compatible with method
    if (!any(method %in% c('indiv', 'samp')))
        stop('method must be "indiv" or "samp" for individual or sample based rarefaction, respectively')
    if (method == 'samp') {
        if (is.null(dim(x)))
            stop('For sample based rarefaction "x" must be a site x species matrix as the input')
        else {
            x = (x > 0) * 1             
            n = nrow(x) # empty sites are counted as samples 
            x = colSums(x)
        }
    }
    if (method == 'indiv') {
        if (!is.null(dim(x)))
            x = colSums(x)
        n = sum(x)
    }
    x = x[x > 0] # drop species with no observations
    S = length(x)
    if (is.null(effort))
        effort = 1:n
    if (any(effort > n))
        warning('"effort" larger than total number of samples')
    ldiv = lchoose(n, effort)
    p = matrix(0, length(effort), S)
    for (i in seq_along(effort)) {
        p[i, ] = ifelse(n - x < effort[i], 0, 
                        exp(lchoose(n - x, effort[i]) - ldiv[i]))
    }
    out = rowSums(1 - p)
    names(out) = effort
    class(out) = 'rarefaction'
    return(out)
}

rarefy_sample_explicit = function(dat_sp, dat_plot, group, effort) {
    plot_grp = dat_plot[dat_plot$group == group, ] 
    sp_grp = dat_sp[dat_plot$group == group, ]
    explicit_loop = matrix(0, nrow(sp_grp), nrow(sp_grp))
    pair_dist = as.matrix(dist(plot_grp[ , c('x', 'y')]))
    for (i in 1:nrow(plot_grp)) {
        focal_site = plot_grp$plot[i]
        dist_to_site = pair_dist[i, ]
        # Shuffle plots, so that tied grouping is not biased by original order.
        new_order = sample(1:nrow(plot_grp))  
        # new_order = 1:nrow(plot_grp) # Alternative: no shuffles
        plots_new = plot_grp$plot[new_order]
        dist_new = dist_to_site[new_order]
        plots_new_ordered = plots_new[order(dist_new)]
        # Move focal site to the front
        plots_new_ordered = c(focal_site, 
                              plots_new_ordered[plots_new_ordered != focal_site])  
        sp_ordered = sp_grp[match(row.names(sp_grp), plots_new_ordered), ]
        # 1 for absence, 0 for presence
        sp_bool = as.data.frame((sp_ordered == 0) * 1) 
        rich = cumprod(sp_bool)
        explicit_loop[ , i] = as.numeric(ncol(dat_sp) - 1 - rowSums(rich))
    }
    explicit_S = apply(explicit_loop, 1, mean)
    return(explicit_S[effort])
}

rarefy_sample_implicit = function(dat_sp, dat_plot, group, effort=NULL) {
    sp_grp = dat_sp[dat_plot$group == group, ]
    sample_S = rarefaction(sp_grp, 'samp', effort)
    return(sample_S)
}

rarefy_individual = function(dat_sp, dat_plot, group, effort=NULL) {
    sad = colSums(dat_sp[dat_plot$group == group, ])
    ind_S = rarefaction(sad, 'indiv', effort)
    return(ind_S)
}

get_deltaSsad = function(dat_sp, dat_plot, groups, Nind=NULL) {
    if (is.null(Nind))
        Nind = min(tapply(rowSums(dat_sp),  dat_plot$group, sum))
    # the first element of the argument groups must be the reference group
    rescaled_ind = sapply(groups, function(x) 
                          rarefy_individual(dat_sp, dat_plot, x, 1:Nind))
    deltaSsad = rescaled_ind[ , 2] - rescaled_ind[ , 1]
    return(deltaSsad)
}

get_deltaSN = function(dat_sp, dat_plot, groups, ScaleBy = NA, Nind=NULL) {
    if (is.null(Nind))
        Nind = min(tapply(rowSums(dat_sp),  dat_plot$group, sum))
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy) 
    # shouldn't the above line be:  'mean(rowSums(dat_sp))' ?
    effort = 1:min(table(dat_plot$group[dat_plot$group %in% groups]))
    implicit_sample = sapply(groups, function(x)
                             rarefy_sample_implicit(dat_sp, dat_plot, x, effort))
    # DJM: It may be more conservative to aggregate the individual rarefaction result
    # rather than interpolate the sample based differences to an artifical level of 
    # precision 
    # rescale sample based differences via intpolation to reflect numbers of individuals sampled
    # Use S_diff(1) = 0 for interpolation
    # here the differences are limited to comparisons of two groups
    deltaSsamp = pchip(xi=c(1, effort * avg_dens), 
                       yi=c(0, implicit_sample[ , 2] - implicit_sample[ , 1]),
                       x=1:Nind)
    # if the calculation of the sad rarefaction is slow then possibly it should
    # be allowed as an argument as well.
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, groups)
    min_nind = min(length(deltaSsamp), length(deltaSsad))
    deltaSN = deltaSsamp[1:min_nind] - deltaSsad[1:min_nind]
    return(deltaSN)
}

get_deltaSagg = function(dat_sp, dat_plot, groups) {
    min_nplots = min(table(dat_plot$group[dat_plot$group %in% groups]))
    implicit_sample = sapply(groups, function(x) 
                             rarefy_sample_implicit(dat_sp, dat_plot, x, 1:min_nplots))
    explicit_sample = sapply(groups, function(x) 
                             rarefy_sample_explicit(dat_sp, dat_plot, x, 1:min_nplots))
    deltaSagg = as.numeric(na.omit(explicit_sample[ , 2] - explicit_sample[ , 1] - 
                                   (implicit_sample[ , 2] - implicit_sample[ , 1])))
    return(deltaSagg)
}

null_sad = function(dat_sp, dat_plot, groups, nperm = 1000, CI = 0.95) {
    # This function generates null confidence intervals for the deltaSsad curve
    # assuming that samples of the two treatments come from the same SAD.
    # Only keep the plots for the two given treatments
    dat_sp = dat_sp[dat_plot$group %in% groups, ] 
    sp_extend = unlist(apply(dat_sp, 1, function(x) 
                             rep(1:ncol(dat_sp), x)))
    n_plot = apply(dat_sp, 1, sum)
    grp_extend = rep(dat_plot$group, times=n_plot)
    Nind = min(tapply(n_plot, dat_plot$group, sum)) 
    dS_perm = matrix(NA, nperm, Nind)
    for (i in 1:nperm) {
        # Shuffle treatment label of each individual  
        grp_shuffle = sample(grp_extend) 
        # see vegan::rrarefy for a possibly better way to shuffle indiv
        sad_shuffle = sapply(groups, function(x) 
                             as.integer(table(factor(
                                 sp_extend[grp_shuffle == x], levels=1:ncol(dat_sp)))))
        ind_S = apply(sad_shuffle, 2, rarefaction, 'indiv', 1:Nind)
        # this step is explicitly pairwise 
        dS_perm[i, ] = ind_S[ , 2] - ind_S[ , 1]
    }
    dS_mean = apply(dS_perm, 2, mean)
    dS_qt = apply(dS_perm, 2, function(x) 
                  quantile(x, c((1 - CI)/2, (1 + CI)/2)))
    out = data.frame(mean = dS_mean, lowerCI = dS_qt[1, ], upperCI = dS_qt[2, ])
    return(out)
}

null_N = function(dat_sp, dat_plot, groups, nperm = 1000, CI = 0.95, ScaleBy = NA) {
    # This function generates null confidence intervals for the deltaSN curve
    # assuming that the two treatments do not differe systematically in N.
    S = ncol(dat_sp)
    dat_sp = dat_sp[dat_plot$group %in% groups, ]
    dat_plot = dat_plot[dat_plot$group %in% groups, ]
    Nplot = min(table(dat_plot$group))
    grp_sads = list()
    for (i in groups) {
        plots = dat_plot$plot[dat_plot$group == i]
        grp_sad = colSums(dat_sp[row.names(dat_sp) %in% plots, ])
        grp_sads[[i]] = rep(1:S, grp_sad)
    }
    # Abundance within each plot
    plot_abus = apply(dat_sp, 1, sum)  
    Nind = min(sapply(grp_sads, length)) 
    deltaSN_perm = matrix(NA, nperm, Nind)
    for (i in 1:nperm) {
        plot_abu_perm = sample(plot_abus)
        sp_draws = sapply(1:nrow(dat_sp), function(x) 
                          sample(rep(1:S, dat_sp[x,]), 
                                 size=plot_abu_perm[x], replace=T))
        dat_sp_perm = t(sapply(1:nrow(dat_sp), function(x) 
                               table(c(1:S, sp_draws[[x]])) - 1 ))
        deltaSN_perm[i, ] = get_deltaSN(dat_sp_perm, dat_plot, groups, ScaleBy,
                                        Nind)
    }
    quant_mean = apply(deltaSN_perm, 2, mean, na.rm = T)
    quant_lower = apply(deltaSN_perm, 2, function(x) 
                        quantile(x, (1 - CI)/2, na.rm = T))
    quant_higher = apply(deltaSN_perm, 2, function(x) 
                         quantile(x, (1 + CI)/2, na.rm = T))
    return(data.frame(mean = quant_mean, lowerCI = quant_lower, upperCI = quant_higher))
}

null_agg = function(dat_sp, dat_plot, groups, nperm = 1000, CI = 0.95) {
    # This function generates null confidence intervals for the deltaSagg curve
    # assuming that plots within each treatment have no spatial structure.
    nplots = table(dat_plot$group)
    grp_index = sapply(groups, function(x) which(dat_plot$group == x))
    deltaSagg_perm = matrix(NA, nperm, min(nplots))
    for (i in 1:nperm) {
        # Shuffle within treatments
        index_shuffle = sapply(grp_index, sample)
        new_order = unlist(index_shuffle)[order(unlist(grp_index))]
        plot_shuffle_loc = dat_plot
        plot_shuffle_loc[unlist(grp_index), c('x', 'y')] = plot_shuffle_loc[new_order, c('x','y')]
        deltaSagg_perm[i, ] = get_deltaSagg(dat_sp, plot_shuffle_loc, groups) 
    }
    quant_mean = apply(deltaSagg_perm, 2, mean, na.rm = T)
    quant_lower = apply(deltaSagg_perm, 2, function(x) 
                        quantile(x, (1 - CI)/2, na.rm = T))
    quant_higher = apply(deltaSagg_perm, 2, function(x) 
                         quantile(x, (1 + CI)/2, na.rm = T))
    return(data.frame(mean = quant_mean, lowerCI = quant_lower, upperCI = quant_higher))
}

table_effect_on_S = function(dat_sp, dat_plot, groups, ScaleBy = NA) {
    # Returns a data frame with the effects of SAD, N, and aggregation on diversity
    # across scales
    # not b/c of spatial ties the values will change every time 
    # this is calculated therefore best pratice may be
    # tst = replicate(20, table_effect_on_S(dat_sp, dat_plot, groups, ScaleBy), simplify=FALSE)
    # plyr::aaply(plyr::laply(tst, as.matrix), c(2, 3), mean)
    nplots = table(dat_plot$group)
    explicit_sample = sapply(groups, function(x) 
                             rarefy_sample_explicit(dat_sp, dat_plot, x, 1:min(nplots)))
    overall = as.numeric(na.omit(explicit_sample[ , 2] - explicit_sample[ , 1]))
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, groups)
    deltaSN = get_deltaSN(dat_sp, dat_plot, groups, ScaleBy) # why is this call diff
    deltaSagg = get_deltaSagg(dat_sp, dat_plot, groups)
    # Rarefy to desired abundances
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    max_level = floor(log10(avg_dens * min(nplots)))
    out = as.data.frame(matrix(NA, 4, max_level))
    row.names(out) = c("overall", "SAD", "N", "aggregation")
    names(out) = as.character(10^(1:max_level))
    for (row in c(2, 3)) {
        deltaS = unlist(list(overall, deltaSsad, deltaSN, deltaSagg)[row])
        out[row, ] = sapply(10^(1:max_level), function(x) 
                            ifelse(length(deltaS) >= x, deltaS[x], NA))
    }
    for (row in c(1, 4)) {
        deltaS = unlist(list(overall, deltaSsad, deltaSN, deltaSagg)[row])
        out_row = pchip(xi=(0:length(deltaS)) * avg_dens, yi=c(0, deltaS), 
                        x=10^(1:min(max_level, floor(log10(length(deltaS) * avg_dens)))))
        out[row, 1:length(out_row)] = out_row
    }
    out = cbind(out, c(overall[length(overall)], deltaSsad[length(deltaSsad)], 
                deltaSN[length(deltaSN)], deltaSagg[length(deltaSagg)]))
    names(out)[max_level + 1] = length(deltaSsad)
    return(out)
}

pairwise_t = function(dat_sp, dat_plot, groups, lower_N = NA) {
    dat_plot_grps = dat_plot[dat_plot$group %in% groups, ]
    dat_sp = dat_sp[match(dat_plot_grps$plot, row.names(dat_sp)), ]
    S_list = rowSums(dat_sp > 0)
    N_list = rowSums(dat_sp)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, ]/N_list[x])^2)))
    if (is.na(lower_N)) {
        rarefied_S_list = apply(dat_sp, 1, function(x) 
                                rarefaction(x, 'indiv', effort = 1:min(N_list)))
    } else {
        # Remove plots with abundance below lower_N in the analysis of rarefied S
        rarefied_S_list = apply(dat_sp, 1, function(x) 
                                if (sum(x) < lower_N)
                                    rep(NA, lower_N)
                                else 
                                    rarefaction(x, 'indiv', effort = 1:lower_N))
        if (any(is.na(rarefied_S_list))) 
            print("Warning: some plots are removed in rarefaction.")
    }
    out = as.data.frame(matrix(NA, 5, 4))
    stats_list = list(rarefied_S_list, N_list, PIE_list, S_list)
    for (i in 1:length(stats_list)) {
        stat = unlist(stats_list[i])
        stat_1 = stat[dat_plot$group == groups[1]]
        stat_2 = stat[dat_plot$group == groups[2]]
        stat_1 = stat_1[!is.na(stat_1)]
        stat_2 = stat_2[!is.na(stat_2)]
        out[ , i] = c(mean(stat_1), sd(stat_1), 
                      mean(stat_2), sd(stat_2), 
                      t.test(stat_1, stat_2)$p.val)
    }
    names(out) = c("S_rarefied", "N", "PIE", "S_raw")
    row.names(out) = c(paste(groups[1], "(mean)", sep = ""), 
                       paste(groups[1], "(sd)", sep = ""), 
                       paste(groups[2], "(mean)", sep = ""), 
                       paste(groups[2], "(sd)", sep = ""), "p_value")
    # Boxplots
    par(mfrow = c(2, 2))  # This is not ideal but I cannot get layout to work in Rstudio
    plot_names = c(paste("Rarified S at N=", 
                         ifelse(is.na(lower_N), min(N_list), lower_N), sep = ""),
                   "N", "PIE", "Raw S")
    plot_names = sapply(1:4, function(x) 
                        paste(plot_names[x], " (p=", round(out[5, x], 6), ")", sep = ""))
    for (i in 1:length(stats_list)) {
        stat = unlist(stats_list[i])
        stat_1 = stat[dat_plot$group == groups[1]]
        stat_2 = stat[dat_plot$group == groups[2]]
        stat_1 = stat_1[!is.na(stat_1)]
        stat_2 = stat_2[!is.na(stat_2)]
        boxplot(stat_1, stat_2, names = c(groups[1], groups[2]), main = plot_names[i])
    }
    return(out)
}

## Functions for plotting
plotEffectS = function(dat_sp, dat_plot, groups, Nperm = 1000,
                       CI = 0.95, ScaleBy = NA) {
    # get individual densities
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    #avg_A = mean(dat_plot$area)
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, groups)
    deltaSN = get_deltaSN(dat_sp, dat_plot, groups, ScaleBy)
    deltaSagg = get_deltaSagg(dat_sp, dat_plot, groups)
    deltaS_list = list(deltaSsad, deltaSN, deltaSagg)
    
    sad_CI = null_sad(dat_sp, dat_plot, groups, Nperm, CI)
    N_CI = null_N(dat_sp, dat_plot, groups, Nperm, CI, ScaleBy)
    agg_CI = null_agg(dat_sp, dat_plot, groups, Nperm, CI)
    null_list = list(sad_CI, N_CI, agg_CI)
    
    par(mfrow = c(1, 3), mar = c(5, 4, 6.5, 2) + 0.1)
    main_list = c("SAD", "N", "Aggregation")
    for (i in 1:3) {
        deltaS = unlist(deltaS_list[i])
        null_CI = unlist(null_list[i], recursive = F)
        if (i %in% c(1, 2)) {
            plot(1:length(null_CI$lowerCI), null_CI$lowerCI, type = "l",
                 col = "grey84", ylab = expression(Delta ~ "S"), 
                 xlab = "Number of Individuals", 
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(1, length(deltaS)), cex.lab = 1.5, cex.axis = 1.5)
            polygon(c(1:length(null_CI$lowerCI), length(null_CI$lowerCI):1), 
                    c(null_CI$lowerCI, rev(null_CI$upperCI)), col = "grey84", border = NA)
            lines(1:length(null_CI$mean), null_CI$mean, type = "l", lwd = 2)
            par(new = T)
            plot((1:length(deltaS))/avg_dens * 1, deltaS, col = "red", lwd = 2, 
                 type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(1, length(deltaS))/avg_dens * 1)
            axis(3, cex.axis = 1.5)
            mtext("Number of Plots", side = 3, line = 2.5)
            title(main = main_list[i], line = 4.5, cex.main = 1.8)
        } else {
            plot(1 * 1:length(null_CI$lowerCI), null_CI$lowerCI, type = "l", 
                 col = "grey84", ylab = expression(Delta ~ "S"), xlab = "Number of Plots",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(0, 1 * length(deltaS)), cex.lab = 1.5, cex.axis = 1.5)
            polygon(c(1 * 1:length(null_CI$lowerCI), 1 * length(null_CI$lowerCI):1), 
                    c(null_CI$lowerCI, rev(null_CI$upperCI)), col = "grey84", border = NA)
            lines(1 * 1:length(null_CI$mean), null_CI$mean, type = "l", lwd = 2)
            par(new = T)
            plot((1:length(deltaS)) * avg_dens, deltaS, col = "red", lwd = 2,
                 type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)), 
                 xlim = c(0, avg_dens * length(deltaS)))
            axis(3, cex.axis = 1.5)
            mtext("Number of Individuals", side = 3, line = 2.5)
            title(main = main_list[i], line = 4.5, cex.main = 1.8)
        }
    }
}

plotSADs = function(dat_sp, dat_plot, col = NA) {
    # TO DO: add check to ensure that col is the same length as treatments
    require(scales)
    par(mfrow = c(1, 1))
    grps = unique(dat_plot$group)
    if (is.na(col)) 
        col = rainbow(length(grps))
    plot(1, type = "n", xlab = "% abundance (log scale)", ylab = "% species", 
         xlim = c(0.01, 1), ylim = c(0, 1), log = "x")
    for (i in 1:length(grps)) {
        col_grp = col[i]
        plots_grp = dat_plot[dat_plot$group == grps[i], 1]
        dat_grp = dat_sp[match(plots_grp, row.names(dat_sp)), ]
        for (j in 1:nrow(dat_grp)) {
            sad_row = as.numeric(sort(dat_grp[j, dat_grp[j, ] != 0]))
            s_cul = 1:length(sad_row)/length(sad_row)
            n_cul = sapply(1:length(sad_row), function(x) sum(sad_row[1:x]) / sum(sad_row))
            lines(n_cul, s_cul, col = alpha(col_grp, 0.5), lwd = 1, type = "l")
        }
    }
    legend("bottomright", grps, col = col, lwd = 2)
}

plotSNpie = function(dat_sp, dat_plot, col = NA) {
    # TO DO: add check to ensure that col is the same length as treatments
    require(rgl)
    grps = unique(dat_plot$group)
    if (is.na(col)) 
        col = rainbow(length(grps))
    S_list = rowSums(dat_sp > 0)
    N_list = rowSums(dat_sp)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, ]/N_list[x])^2)))
    grp_list = as.character(dat_plot$group[match(dat_plot$plot, row.names(dat_sp))])
    col_list = sapply(grp_list, function(x) col[which(grps == x)])
    plot3d(S_list, N_list, PIE_list, "S", "N", "PIE", col = col_list, size = 8)
} 

# Auxillary function: difference between the ind-based rarefaction and the sample-based rarefaction
# Output: a list of two components: $plot_sample_size: levels of n at which the deduction is done
#   $deltaS: a data frame, difference between the two curves at plot_sample_size for each plot
effect_of_N = function(comm_group, env_var_keep, ref_dens, min_plot_group){
  out = list()
  # lumped SAD for each group
  group_keep = unique(env_var_keep)
  keep_group_sad = aggregate(comm_group, by = list(env_var_keep), FUN = sum)
  row.names(keep_group_sad) = keep_group_sad[, 1]
  keep_group_sad = keep_group_sad[, -1]
  
  minN_group = min(apply(keep_group_sad, 1, sum))
  plot_sample_size = round(ref_dens * seq(min_plot_group))
  out$plot_sample_size = plot_sample_size[plot_sample_size <= minN_group]
  
  samp_rare = sapply(group_keep, function(x) rarefaction.sample(comm_group[which(env_var_keep == x), ])[1:length(out$plot_sample_size), 2])
  samp_rare = as.data.frame(samp_rare)
  names(samp_rare) = as.character(group_keep)
  n = c(1, ref_dens * seq(length(out$plot_sample_size)))
  out$effect_N = as.data.frame(matrix(NA, nrow = length(group_keep), ncol = length(out$plot_sample_size)))
  for (i in 1:ncol(samp_rare)){
    col = samp_rare[, i]
    col_name = names(samp_rare)[i]
    S = c(1, col)
    interp_sample = pchip(n, S, out$plot_sample_size)
    row_keep_group_sad = keep_group_sad[which(row.names(keep_group_sad) == col_name), ]
    rare_indiv = as.numeric(rarefy(row_keep_group_sad, n))
    out$effect_N[i, ] = rare_indiv[-1] - interp_sample
  }
  row.names(out$effect_N) = names(samp_rare)
  out
}

get_delta_stats = function(comm, type, env_var, test = c('indiv', 'sampl', 'spat'),
                           log.scale = FALSE, inds = NULL, ref = 'mean', corr = 'spearman',
                           nperm = 1000){
  # Inputs:
  # comm - a 'comm' type object, with attributes ...
  # type - if the environmental variable is 'categorical' (anova-type) or 'continuous' (regression-type)
  # env_var - string, name of the environmental variable
  # test - tests to be included. A single value in c('indiv', 'sampl', 'spat'), or a combination of them as a list. 
  #   The default is test = c('indiv', 'sampl', 'spat') (all three tests).
  # inds - argument for individual-based rarefaction. If given a single number, will be taken as the number of points for
  #   for rarefaction. log.scale will then be employed (see below). If given a list of numbers, will be taken as the levels of N 
  #   for rarefaction. Default (NULL) will result in rarefaction at all possible N.
  # log.scale - If number of points for rarefaction is given, log.scale = TRUE leads to values evenly spaced on log scale,
  #   log.scale = FAlSE (default) leads to values evenly spaced on arithmetic scale.
  # ref - reference density used in converting number of plots to number of individuals. Can take one of three values:
  #   'mean', 'max', 'min'. If 'mean', the average plot-level abundance across all plots are used. If 'min' or 'max', 
  #   the minimum/maximum plot-level abundance is used.
  # corr - which correlation is used for S/delta S vs environmental variable. Can be 'spearman' (default, rank correlation)
  #   or 'pearson' (less recommended because of potential nonlinearity)
  # nperm - number of iterations for null models
  # Output:
  # mobr - a 'mobr' type object, with attributes...
  
  #ToDO: check type %in% c('categorical', 'continuous')
  #TODO: other checks
  # TODO: groups as row names in all outputs?
  mobr = list()  # This is the object with the outputs 
  mobr$type = type
  # Assume: env_var is a string with variable name
  # group_sad is the overall SAD for all plots within a group (level of env_var) lumped together
  plot_count = count(comm$envi, env_var) # Number of plots within each group (env level)
  group_keep = plot_count[which(plot_count$freq >= min_plot), 1] # Groups that will be included in steps 2 & 3
  
  group_sad = aggregate(comm$comm, by = list(comm$envi[, env_var]), FUN = sum)
  group_level = group_sad[ , 1]
  group_sad = group_sad[ , -1]
  minN_group = min(apply(group_sad, 1, sum))
  
  rows_keep_group = which(comm$envi[, env_var] %in% group_keep)
  comm_group = comm$comm[rows_keep_group, ]
  env_var_keep = comm$envi[rows_keep_group, env_var]
  keep_group_sad = aggregate(comm_group, by = list(env_var_keep), FUN = sum)
  row.names(keep_group_sad) = keep_group_sad[, 1]
  keep_group_sad = keep_group_sad[, -1]
  if (ref == 'mean'){
    ref_dens = sum(comm$comm_group) / nrow(comm$comm_group)
  }
  else if (ref == 'max'){
    ref_dens = max(apply(comm$comm_group, 1, sum))
  }
  else {
    ref_dens = min(apply(comm$comm_group, 1, sum))
  }
  
  # TODO: warning if minN_group is too low (e.g., < 20)?

  if (type == 'continuous'){
    if ('indiv' %in% test){
      # 1. Ind-based rarefaction (effect of SAD) vs env_var vs N
      if (comm$analyses$indiv == F){
        print('Error: individual-based rarefaction not allowed by data.')
      }
      else {
        # Assess rarefied S - env relationship at 10 log-even levels of N
        if (length(inds) > 1){
          mobr$ind_sample_size = inds
        }
        else if (is.null(inds)){
          mobr$ind_sample_size = seq(minN_group)
        }
        else {
          if (log.scale = T){
            mobr$ind_sample_size = floor(exp(seq(10) * log(minN_group) / 10))
          }
          else {
            mobr$ind_sample_size = floor(seq(10) * minN_group / 10)
          }
        }
        mobr$ind_rare = as.data.frame(rarefy(group_sad, mobr$ind_sample_size)) # rarefied S for each group 
        mobr$ind_r = apply(mobr$ind_rare, 2, 
                            function(x){cor(x, as.numeric(group_level), method = corr)}) 
        # 1'. Null test for ind-based rarefaction
        sp_extent = unlist(apply(group_sad, 1, function(x) rep(1:ncol(group_sad), x)))
        env_extent = rep(group_level, time = apply(group_sad, 1, sum))
        null_ind_r_mat = matrix(NA, nperm, length(mobr$ind_sample_size))
        for (i in 1:nperm){
          env_shuffle = sample(env_extent)
          sad_shuffle = sapply(group_level, function(x) 
            as.integer(table(factor(sp_extent[env_shuffle == x], levels=1:ncol(group_sad)))))
          perm_ind_rare = as.data.frame(rarefy(sad_shuffle, mobr$ind_sample_size, MARGIN = 2))
          null_ind_r_mat[i, ] = apply(perm_ind_rare, 2, function(x){cor(x, as.numeric(group_level), method = corr)})
        }
        mobr$ind_r_null = apply(null_ind_r_mat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975))) # 95% CI
      }
    # 2. Effect of density vs env_var vs N
    if ('sampl' %in% test){
      if (comm$test$sampl == F){
        print('Error: sample-based analysis not allowed by data.')
      }
      else{
        plot_sample_size = round(ref_dens * seq(min(plot_count$freq)))
        mobr$plot_sample_size = plot_sample_size[plot_sample_size <= minN_group] 
        mobr$samp_rare = sapply(group_keep, 
                                function(x) rarefaction.sample(comm$comm[which(comm$envi[, env_var] == x), ])[1:length(mobr$plot_sample_size), 2])
        mobr$samp_rare = as.data.frame(mobr$samp_rare)
        n = c(1, ref_dens * seq(length(mobr$plot_sample_size))) # Levels of n that match sample-based rarefied S (integer unnecessary)
        names(mobr$samp_rare) = as.character(group_keep)
        mobr$effect_N = as.data.frame(matrix(NA, nrow = length(group_keep), ncol = length(mobr$plot_sample_size)))
        for (i in 1:ncol(mobr$samp_rare)){
          col = mobr$samp_rare[, i]
          col_name = names(mobr$samp_rare)[i]
          S = c(1, col)
          interp_sample = pchip(n, S, mobr$plot_sample_size)
          row_keep_group_sad = keep_group_sad[which(row.names(keep_group_sad) == col_name), ]
          rare_indiv = as.numeric(rarefy(row_keep_group_sad, n))
          mobr$effect_N[i, ] = rare_indiv[-1] - interp_sample
        }
        mobr$sample_r = apply(mobr$effect_N, 2, 
                           function(x){cor(x, as.numeric(group_keep), method = corr)})
        # 2'. Null test for effect of N
        
        
        null_N = function(dat_sp, dat_plot, groups, nperm = 1000, CI = 0.95, ScaleBy = NA) {
          # This function generates null confidence intervals for the deltaSN curve
          # assuming that the two treatments do not differe systematically in N.
          S = ncol(dat_sp)
          dat_sp = dat_sp[dat_plot$group %in% groups, ]
          dat_plot = dat_plot[dat_plot$group %in% groups, ]
          Nplot = min(table(dat_plot$group))
          grp_sads = list()
          for (i in groups) {
            plots = dat_plot$plot[dat_plot$group == i]
            grp_sad = colSums(dat_sp[row.names(dat_sp) %in% plots, ])
            grp_sads[[i]] = rep(1:S, grp_sad)
          }
          # Abundance within each plot
          plot_abus = apply(dat_sp, 1, sum)  
          Nind = min(sapply(grp_sads, length)) 
          deltaSN_perm = matrix(NA, nperm, Nind)
          for (i in 1:nperm) {
            plot_abu_perm = sample(plot_abus)
            sp_draws = sapply(1:nrow(dat_sp), function(x) 
              sample(rep(1:S, dat_sp[x,]), 
                     size=plot_abu_perm[x], replace=T))
            dat_sp_perm = t(sapply(1:nrow(dat_sp), function(x) 
              table(c(1:S, sp_draws[[x]])) - 1 ))
            deltaSN_perm[i, ] = get_deltaSN(dat_sp_perm, dat_plot, groups, ScaleBy,
                                            Nind)
          }
          quant_mean = apply(deltaSN_perm, 2, mean, na.rm = T)
          quant_lower = apply(deltaSN_perm, 2, function(x) 
            quantile(x, (1 - CI)/2, na.rm = T))
          quant_higher = apply(deltaSN_perm, 2, function(x) 
            quantile(x, (1 + CI)/2, na.rm = T))
          return(data.frame(mean = quant_mean, lowerCI = quant_lower, upperCI = quant_higher))
        }
        
      }
    }
    
    # 3. Effect of aggregation vs env_var vs N
    
  }
  
  class(mobr) = 'mobr'
  return(mobr)
}
