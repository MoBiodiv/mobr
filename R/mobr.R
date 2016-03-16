## A rough sketch of what we need in the module Self reminder: need two pieces of
## inputs: 1. a data frame of plot characteristics, with columns: plot,
## group, x, y 2. a data frame where rownames are plot ids, and
## subsequent columns are species abundances

## Functions
require(rareNMtests)
require(pracma)
## Functions to obtain the three kinds of curves: sample-based explicit,
## sample-based implicit, and individual-based rescaled to sample
get_avg_dens = function(dat_sp, dat_plot, ScaleBy) {
    # Auxillary function to obtain average density within plot for rescaling Ask user
    # to specify which group is used as 'standard' for rescaling.
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

rarefy_sample_explicit = function(dat_sp, dat_plot, group) {
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
    return(explicit_S)
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

get_deltaSsad = function(dat_sp, dat_plot, groups) {
    # the first element of the argument groups must be the reference group
    rescaled_ind = sapply(groups, function(x) 
                          rarefy_individual(dat_sp, dat_plot, x))
    #Nind = min(tapply(rowSums(dat_sp),  dat_plot$group, sum))
    Nind = min(length(unlist(rescaled_ind[1])), length(unlist(rescaled_ind[2])))
    deltaSsad = as.numeric(unlist(rescaled_ind[2])[1:Nind] - 
                           unlist(rescaled_ind[1])[1:Nind])
    return(deltaSsad)
}

get_deltaSN = function(dat_sp, dat_plot, groups, ScaleBy = NA, Nind) {
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
                             rarefy_sample_explicit(dat_sp, dat_plot, x)[1:min_nplots])
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
        ind_S = rarefaction(t(sad_shuffle), effort=1:Nind)
        # this step is explicitly pairwise 
        dS_perm[i, ] = ind_S[2, ] - ind_S[1 , ]
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
        sp_draws = sapply(1:nplots, function(x) 
                          sample(rep(1:S, dat_sp[x,]), 
                                 size=plot_abu_perm[x], replace=T))
        dat_sp_perm = t(sapply(1:nplots, function(x) 
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

table_effect_on_S = function(dat_sp, dat_plot, treatment1, treatment2, ScaleBy = NA) {
    # Returns a data frame with the effects of SAD, N, and aggregation on diversity
    # across scales
    nplots = c(nrow(dat_plot[dat_plot[ , 2] == treatment1, ]), 
               nrow(dat_plot[dat_plot[ , 2] == treatment2, ]))
    explicit_sample = sapply(c(treatment1, treatment2), function(x) 
                             rarefy_sample_explicit(dat_sp, dat_plot, x)[1:min(nplots)])
    overall = as.numeric(na.omit(explicit_sample[ , 2] - explicit_sample[ , 1]))
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, treatment1, treatment2)
    deltaSN = get_deltaSN(dat_sp, dat_plot, treatment1, treatment2, ScaleBy)
    deltaSagg = get_deltaSagg(dat_sp, dat_plot, treatment1, treatment2)
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
        out_row = pchip((0:length(deltaS)) * avg_dens, c(0, deltaS), 
                        10^(1:min(max_level, floor(log10(length(deltaS) * avg_dens)))))
        out[row, 1:length(out_row)] = out_row
    }
    out$maxN = c(overall[length(overall)], deltaSsad[length(deltaSsad)], 
                 deltaSN[length(deltaSN)], deltaSagg[length(deltaSagg)])
    return(out)
}

pairwise_t = function(dat_sp, dat_plot, treatment1, treatment2, lower_N = NA) {
    dat_plot_grps = dat_plot[dat_plot[ , 2] %in% c(treatment1, treatment2), ]
    dat_sp = dat_sp[match(dat_plot_grps[ , 1], dat_sp[ , 1]), ]
    S_list = sapply(1:nrow(dat_sp), function(x) 
                    length(which(dat_sp[x, 2:ncol(dat_sp)] != 0)))
    N_list = apply(dat_sp[ , 2:ncol(dat_sp)], 1, sum)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, 2:ncol(dat_sp)]/N_list[x])^2)))
    if (is.na(lower_N)) {
        rarefied_S_list = sapply(1:nrow(dat_sp), function(x) 
                                 as.numeric(rarefaction.individual(dat_sp[x, 2:ncol(dat_sp)], 
                                                                   effort = 1:min(N_list))[2]))
    } else {
        # Remove plots with abundance below lower_N in the analysis of rarefied S
        rarefied_S_list = sapply(1:nrow(dat_sp), function(x) {
                                 if (sum(dat_sp[x, 2:ncol(dat_sp)]) < lower_N)
                                     NA 
                                 else 
                                     as.numeric(rarefaction.individual(dat_sp[x, 2:ncol(dat_sp)],
                                                                       effort = 1:lower_N)[2])})
        if (any(is.na(rarefied_S_list))) 
            print("Warning: some plots are removed in rarefaction.")
    }
    out = as.data.frame(matrix(NA, 5, 4))
    stats_list = list(rarefied_S_list, N_list, PIE_list, S_list)
    for (i in 1:length(stats_list)) {
        stat = unlist(stats_list[i])
        stat_1 = stat[dat_plot[ , 2] == treatment1]
        stat_2 = stat[dat_plot[ , 2] == treatment2]
        stat_1 = stat_1[!is.na(stat_1)]
        stat_2 = stat_2[!is.na(stat_2)]
        out[ , i] = c(mean(stat_1), sd(stat_1), 
                      mean(stat_2), sd(stat_2), 
                      t.test(stat_1, stat_2)$p.val)
    }
    names(out) = c("S_rarefied", "N", "PIE", "S_raw")
    row.names(out) = c(paste(treatment1, "(mean)", sep = ""), 
                       paste(treatment1, "(sd)", sep = ""), 
                       paste(treatment2, "(mean)", sep = ""), 
                       paste(treatment2, "(sd)", sep = ""), "p_value")
    # Boxplots
    par(mfrow = c(2, 2))  # This is not ideal but I cannot get layout to work in Rstudio
    plot_names = c(paste("Rarified S at N=", 
                         ifelse(is.na(lower_N), min(N_list), lower_N), sep = ""),
                   "N", "PIE", "Raw S")
    plot_names = sapply(1:4, function(x) 
                        paste(plot_names[x], " (p=", round(out[5, x], 6), ")", sep = ""))
    for (i in 1:length(stats_list)) {
        stat = unlist(stats_list[i])
        stat_1 = stat[dat_plot[ , 2] == treatment1]
        stat_2 = stat[dat_plot[ , 2] == treatment2]
        stat_1 = stat_1[!is.na(stat_1)]
        stat_2 = stat_2[!is.na(stat_2)]
        boxplot(stat_1, stat_2, names = c(treatment1, treatment2), main = plot_names[i])
    }
    return(out)
}

lch = function(n, k) lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)

## Functions for plotting
plotEffectS = function(dat_sp, dat_plot, groups, Nperm = 1000,
                       CI = 0.95, ScaleBy = NA) {
    # get individual densities
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    avg_A = mean(dat_plot$area)
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
            plot((1:length(deltaS))/avg_dens * avg_A, deltaS, col = "red", lwd = 2, 
                 type = "l", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(1, length(deltaS))/avg_dens * avg_A)
            axis(3, cex.axis = 1.5)
            mtext("Area", side = 3, line = 2.5)
            title(main = main_list[i], line = 4.5, cex.main = 1.8)
        } else {
            plot(avg_A * 1:length(null_CI$lowerCI), null_CI$lowerCI, type = "l", 
                 col = "grey84", ylab = expression(Delta ~ "S"), xlab = "Area",
                 ylim = c(min(deltaS, null_CI$lowerCI), max(deltaS, null_CI$upperCI)),
                 xlim = c(0, avg_A * length(deltaS)), cex.lab = 1.5, cex.axis = 1.5)
            polygon(c(avg_A * 1:length(null_CI$lowerCI), avg_A * length(null_CI$lowerCI):1), 
                    c(null_CI$lowerCI, rev(null_CI$upperCI)), col = "grey84", border = NA)
            lines(avg_A * 1:length(null_CI$mean), null_CI$mean, type = "l", lwd = 2)
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
    grps = unique(dat_plot[ , 2])
    if (is.na(col)) 
        col = rainbow(length(grps))
    plot(1, type = "n", xlab = "% abundance (log scale)", ylab = "% species", 
         xlim = c(0.01, 1), ylim = c(0, 1), log = "x")
    for (i in 1:length(grps)) {
        col_grp = col[i]
        plots_grp = dat_plot[dat_plot[ , 2] == grps[i], 1]
        dat_grp = dat_sp[match(plots_grp, dat_sp[ , 1]), 2:ncol(dat_sp)]
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
    grps = unique(dat_plot[ , 2])
    if (is.na(col)) 
        col = rainbow(length(grps))
    S_list = sapply(1:nrow(dat_sp), function(x) 
                    length(which(dat_sp[x, 2:ncol(dat_sp)] != 0)))
    N_list = apply(dat_sp[ , 2:ncol(dat_sp)], 1, sum)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * 
                      (1 - sum((dat_sp[x, 2:ncol(dat_sp)] / N_list[x])^2)))
    grp_list = as.character(dat_plot[match(dat_plot[ , 1], dat_sp[ , 1]), 2])
    col_list = sapply(grp_list, function(x) col[which(grps == x)])
    plot3d(S_list, N_list, PIE_list, "S", "N", "PIE", col = col_list, size = 8)
} 
