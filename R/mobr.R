## A rough sketch of what we need in the module Self reminder: need two pieces of
## inputs: 1. a data frame of plot characteristics, with columns: plot ID,
## treatment, x, y, area 2. a data frame where the 1st column is plot ID, and
## subsequent columns are species abundances

## Functions
require(rareNMtests)
require(pracma)
## Functions to obtain the three kinds of curves: sample-based explicit,
## sample-based implicit, and individual-based rescaled to sample
get_avg_dens = function(dat_sp, dat_plot, ScaleBy) {
    # Auxillary function to obtain average density within plot for rescaling Ask user
    # to specify which treatment is used as 'standard' for rescaling.
    if (!is.na(ScaleBy)) {
        dat_scale = dat_sp[match(dat_plot[which(dat_plot[ , 2] == ScaleBy), 1],
                                 dat_sp[ , 1]), 2:ncol(dat_sp)]
        avg_dens = mean(apply(dat_scale, 1, sum))
    } else {
        # If unspecified, uses min density among treatments.
        trmts = unique(dat_plot[ , 2])
        avg_dens = min(sapply(trmts, function(x) 
            mean(apply(dat_sp[match(dat_plot[which(dat_plot[ , 2] == x), 1], dat_sp[ , 1]),
                              2:ncol(dat_sp)], 1, sum))))
    }
    return(avg_dens)
}

rarefy_sample_explicit = function(dat_sp, dat_plot, treatment) {
    plot_trmt = dat_plot[dat_plot[ , 2] == treatment, ]
    sp_trmt = dat_sp[match(plot_trmt[ , 1], dat_sp[ , 1]), ]
    explic_loop = matrix(0, nrow(sp_trmt), nrow(sp_trmt))
    pair_dist = as.matrix(dist(plot_trmt[ , 3:4]))
    for (i in 1:nrow(plot_trmt)) {
        focal_site = plot_trmt[i, 1]
        dist_to_site = pair_dist[i, ]
        # Shuffle plots, so that tied grouping is not biased by original order.
        new_order = sample(1:nrow(plot_trmt))  
        # new_order = 1:nrow(plot_trmt) # Alternative: no shuffles
        plots_new = plot_trmt[new_order, 1]
        dist_new = dist_to_site[new_order]
        plots_new_ordered = plots_new[order(dist_new)]
        # Move focal site to the front
        plots_new_ordered = c(focal_site, 
                              plots_new_ordered[plots_new_ordered != focal_site])  
        sp_ordered = sp_trmt[match(sp_trmt[ , 1], plots_new_ordered), 2:ncol(sp_trmt)]
        # 1 for absence, 0 for presence
        sp_bool = as.data.frame(ifelse(sp_ordered[ , 1:ncol(sp_ordered)] == 0, 1, 0)) 
        rich = cumprod(sp_bool)
        explic_loop[ , i] = as.numeric(ncol(dat_sp) - 1 - rowSums(rich))
    }
    explic_S = apply(explic_loop, 1, mean)
    return(explic_S)
}

rarefy_sample_implicit = function(dat_sp, dat_plot, treatment) {
    plot_trmt = dat_plot[dat_plot[ , 2] == treatment, ]
    sp_trmt = dat_sp[match(plot_trmt[ , 1], dat_sp[ , 1]), ]
    sample_S = rarefaction.sample(sp_trmt[ , 2:ncol(sp_trmt)])[ , 2]
    return(sample_S)
}

rarefy_individual = function(dat_sp, dat_plot, treatment) {
    # Returns the rarefaction results at each sample size
    sp_trmt = dat_sp[match(dat_plot[dat_plot[ , 2] == treatment, 1], dat_sp[ , 1]), ]
    trmt_sad = as.numeric(apply(sp_trmt[ , 2:ncol(sp_trmt)], 2, sum))
    trmt_sad = trmt_sad[trmt_sad != 0]
    ind_S = rarefaction.individual(trmt_sad)[ , 2]
    return(ind_S)
}

get_deltaSsad = function(dat_sp, dat_plot, treatment1, treatment2) {
    rescaled_ind = sapply(c(treatment1, treatment2), function(x) 
                          rarefy_individual(dat_sp, dat_plot, x))
    Nind = min(length(unlist(rescaled_ind[1])), length(unlist(rescaled_ind[2])))
    deltaSsad = as.numeric(unlist(rescaled_ind[2])[1:Nind] - 
                           unlist(rescaled_ind[1])[1:Nind])
    return(deltaSsad)
}

get_deltaSN = function(dat_sp, dat_plot, treatment1, treatment2, ScaleBy = NA) {
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    nplots = c(nrow(dat_plot[dat_plot[ , 2] == treatment1, ]), nrow(dat_plot[dat_plot[ , 
        2] == treatment2, ]))
    implicit_sample = sapply(c(treatment1, treatment2), function(x)
                             rarefy_sample_implicit(dat_sp, dat_plot, x)[1:min(nplots)])
    # Use S_diff(1) = 0 for interpolation
    implicit_sample_ind = pchip(c(1, (1:nrow(implicit_sample)) * avg_dens), 
                                c(0, implicit_sample[ , 2] - implicit_sample[ , 1]),
                                1:floor(nrow(implicit_sample) * avg_dens))
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, treatment1, treatment2)
    # add variable called min index to clean up code
    deltaSN = implicit_sample_ind[1:min(length(implicit_sample_ind), 
                  length(deltaSsad))] - deltaSsad[1:min(length(implicit_sample_ind),
                  length(deltaSsad))]
    return(deltaSN)
}

get_deltaSagg = function(dat_sp, dat_plot, treatment1, treatment2) {
    nplots = c(nrow(dat_plot[dat_plot[ , 2] == treatment1, ]), 
               nrow(dat_plot[dat_plot[ , 2] == treatment2, ]))
    implicit_sample = sapply(c(treatment1, treatment2), function(x) 
        rarefy_sample_implicit(dat_sp, dat_plot, x)[1:min(nplots)])
    explicit_sample = sapply(c(treatment1, treatment2), function(x) 
        rarefy_sample_explicit(dat_sp, dat_plot, x)[1:min(nplots)])
    deltaSagg = as.numeric(na.omit(explicit_sample[ , 2] - explicit_sample[ , 1] - 
        (implicit_sample[ , 2] - implicit_sample[ , 1])))
    return(deltaSagg)
}

null_sad = function(dat_sp, dat_plot, treatment1, treatment2, nperm = 1000,
                    CI = 0.95) {
    # This function generates null confidence intervals for the deltaSsad curve
    # assuming that samples of the two treatments come from the same SAD.
    dat_plot_trmts = dat_plot[dat_plot[ , 2] %in% c(treatment1, treatment2), ]
    # Only keep the plots for the two given treatments
    dat_sp = dat_sp[match(dat_plot_trmts[ , 1], dat_sp[ , 1]), ] 
    sp_extend = unlist(sapply(1:nrow(dat_sp), function(x) 
                              rep(2:ncol(dat_sp), dat_sp[x, 2:ncol(dat_sp)])))
    trmt_extend = unlist(sapply(1:nrow(dat_sp), function(x) 
                                rep(dat_plot[dat_plot[ , 1] == dat_sp[x, 1], 2],
                                    sum(dat_sp[x, 2:ncol(dat_sp)]))))
    Nind = min(length(trmt_extend[trmt_extend == treatment1]), 
               length(trmt_extend[trmt_extend == treatment2]))
    deltaSsad_perm = matrix(NA, nperm, Nind)
    for (i in 1:nperm) {
        # Shuffle treatment label of each individual  
        trmt_shuffle = sample(trmt_extend)  
        sad_shuffle = sapply(c(treatment1, treatment2), function(x) 
                             as.numeric(table(factor(
                                 sp_extend[trmt_shuffle == x], levels = 2:ncol(dat_sp)))))
        ind_S = lapply(1:2, function(x) rarefaction.individual(sad_shuffle[ , x])[ , 2])
        deltaSsad_perm[i, ] = as.numeric(unlist(ind_S[2])[1:Nind] - 
                                         unlist(ind_S[1])[1:Nind])
    }
    quant_mean = as.numeric(apply(deltaSsad_perm, 2, mean, na.rm = T))
    quant_lower = as.numeric(apply(deltaSsad_perm, 2, function(x) 
                                   quantile(x, (1 - CI)/2, na.rm = T)))
    quant_higher = as.numeric(apply(deltaSsad_perm, 2, function(x) 
                                    quantile(x, (1 + CI)/2, na.rm = T)))
    return(list(mean = quant_mean, lowerCI = quant_lower, upperCI = quant_higher))
}

null_N = function(dat_sp, dat_plot, treatment1, treatment2, nperm = 1000, CI = 0.95, 
    ScaleBy = NA) {
    # This function generates null confidence intervals for the deltaSN curve
    # assuming that the two treatments do not differe systematically in N.
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    nplots = c(nrow(dat_plot[dat_plot[ , 2] == treatment1, ]), 
               nrow(dat_plot[dat_plot[ , 2] == treatment2, ]))
    dat_plot_trmts = dat_plot[dat_plot[ , 2] %in% c(treatment1, treatment2), ]
    # Only keep the plots for the two given treatments
    dat_sp = dat_sp[match(dat_plot_trmts[ , 1], dat_sp[ , 1]), ]  
    trmt_sads = list()
    for (trmt in c(treatment1, treatment2)) {
        plots = dat_plot[dat_plot[ , 2] == trmt, 1]
        trmt_sad = apply(dat_sp[dat_sp[ , 1] %in% plots, 2:ncol(dat_sp)], 2, sum)
        trmt_sads = c(trmt_sads, list(rep(2:ncol(dat_sp), trmt_sad)))
    }
    # Abundance within each plot
    n_plot = apply(dat_sp[ , 2:ncol(dat_sp)], 1, sum)  
    sad_row = lapply(1:nrow(dat_sp), function(x) 
                     rep(2:ncol(dat_sp), dat_sp[x, 2:ncol(dat_sp)]))
    Nind = min(floor(min(nplots) * avg_dens), 
               length(unlist(trmt_sads[1])), 
               length(unlist(trmt_sads[2])))
    deltaSN_perm = matrix(NA, nperm, Nind)
    for (i in 1:nperm) {
        n_plot_shuffle = sample(n_plot)
        dat_sp_perm = as.data.frame(matrix(0, nrow(dat_sp), ncol(dat_sp)))
        dat_sp_perm[ , 1] = dat_sp[ , 1]
        new_counts = sapply(1:nrow(dat_sp), function(x) as.numeric(table(factor(
                            sample(if (length(unlist(sad_row[x])) > 0) 
                                       unlist(sad_row[x])  # If plot is empty, use treatment-level SAD
                                   else
                                       unlist(trmt_sads[
                                           which(c(treatment1, treatment2) == 
                                                 dat_plot[x, 2])]),
                                   n_plot_shuffle[x], replace = T),
                            levels = 2:ncol(dat_sp)))))
        dat_sp_perm[ , 2:ncol(dat_sp_perm)] = as.data.frame(t(new_counts))
        deltaSN_perm[i, ] = get_deltaSN(dat_sp_perm, dat_plot_trmts,
                                        treatment1, treatment2, ScaleBy)[1:Nind]
    }
    quant_mean = apply(deltaSN_perm, 2, mean, na.rm = T)
    quant_lower = apply(deltaSN_perm, 2, function(x) 
                        quantile(x, (1 - CI)/2, na.rm = T))
    quant_higher = apply(deltaSN_perm, 2, function(x) 
                         quantile(x, (1 + CI)/2, na.rm = T))
    return(list(mean = quant_mean, lowerCI = quant_lower, upperCI = quant_higher))
}

null_agg = function(dat_sp, dat_plot, treatment1, treatment2, nperm = 1000, CI = 0.95) {
    # This function generates null confidence intervals for the deltaSagg curve
    # assuming that plots within each treatment have no spatial structure.
    nplots = c(nrow(dat_plot[dat_plot[ , 2] == treatment1, ]),
               nrow(dat_plot[dat_plot[ , 2] == treatment2, ]))
    trmt_index = list(which(dat_plot[ , 2] == treatment1), 
                      which(dat_plot[ , 2] == treatment2))
    deltaSagg_perm = matrix(NA, nperm, min(nplots))
    for (i in 1:nperm) {
        # Shuffle within treatments
        index_shuffle = list(sample(unlist(trmt_index[1])), sample(unlist(trmt_index[2])))  
        new_order = unlist(index_shuffle)[order(unlist(trmt_index))]
        plot_shuffle_loc = dat_plot
        plot_shuffle_loc[unlist(trmt_index), 3:4] = plot_shuffle_loc[new_order, 3:4]
        deltaSagg_perm[i, ] = get_deltaSagg(dat_sp, plot_shuffle_loc, 
                                            treatment1, treatment2)
    }
    quant_mean = apply(deltaSagg_perm, 2, mean, na.rm = T)
    quant_lower = apply(deltaSagg_perm, 2, function(x) 
                        quantile(x, (1 - CI)/2, na.rm = T))
    quant_higher = apply(deltaSagg_perm, 2, function(x) 
                         quantile(x, (1 + CI)/2, na.rm = T))
    return(list(mean = quant_mean, lowerCI = quant_lower, upperCI = quant_higher))
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
    dat_plot_trmts = dat_plot[dat_plot[ , 2] %in% c(treatment1, treatment2), ]
    dat_sp = dat_sp[match(dat_plot_trmts[ , 1], dat_sp[ , 1]), ]
    S_list = sapply(1:nrow(dat_sp), function(x) 
                    length(which(dat_sp[x, 2:ncol(dat_sp)] != 0)))
    N_list = apply(dat_sp[ , 2:ncol(dat_sp)], 1, sum)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * (1 - sum((dat_sp[x, 2:ncol(dat_sp)]/N_list[x])^2)))
    if (is.na(lower_N)) {
        rarefied_S_list = sapply(1:nrow(dat_sp), function(x) 
                                 as.numeric(rarefaction.individual(dat_sp[x, 2:ncol(dat_sp)], 
                                                                   inds = min(N_list))[2]))
    } else {
        # Remove plots with abundance below lower_N in the analysis of rarefied S
        rarefied_S_list = sapply(1:nrow(dat_sp), function(x) {
                                 if (sum(dat_sp[x, 2:ncol(dat_sp)]) < lower_N)
                                     NA 
                                 else 
                                     as.numeric(rarefaction.individual(dat_sp[x, 2:ncol(dat_sp)],
                                                                       inds = lower_N)[2])})
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

## Functions for plotting
plotEffectS = function(dat_sp, dat_plot, treatment1, treatment2, Nperm = 1000,
                       CI = 0.95, ScaleBy = NA) {
    avg_dens = get_avg_dens(dat_sp, dat_plot, ScaleBy)
    avg_A = mean(dat_plot[ , 5])
    deltaSsad = get_deltaSsad(dat_sp, dat_plot, treatment1, treatment2)
    deltaSN = get_deltaSN(dat_sp, dat_plot, treatment1, treatment2, ScaleBy)
    deltaSagg = get_deltaSagg(dat_sp, dat_plot, treatment1, treatment2)
    deltaS_list = list(deltaSsad, deltaSN, deltaSagg)
    
    sad_CI = null_sad(dat_sp, dat_plot, treatment1, treatment2, Nperm, CI)
    N_CI = null_N(dat_sp, dat_plot, treatment1, treatment2, Nperm, CI, ScaleBy)
    agg_CI = null_agg(dat_sp, dat_plot, treatment1, treatment2, Nperm, CI)
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
    trmts = unique(dat_plot[ , 2])
    if (is.na(col)) 
        col = rainbow(length(trmts))
    plot(1, type = "n", xlab = "% abundance (log scale)", ylab = "% species", 
         xlim = c(0.01, 1), ylim = c(0, 1), log = "x")
    for (i in 1:length(trmts)) {
        col_trmt = col[i]
        plots_trmt = dat_plot[dat_plot[ , 2] == trmts[i], 1]
        dat_trmt = dat_sp[match(plots_trmt, dat_sp[ , 1]), 2:ncol(dat_sp)]
        for (j in 1:nrow(dat_trmt)) {
            sad_row = as.numeric(sort(dat_trmt[j, dat_trmt[j, ] != 0]))
            s_cul = 1:length(sad_row)/length(sad_row)
            n_cul = sapply(1:length(sad_row), function(x) sum(sad_row[1:x]) / sum(sad_row))
            lines(n_cul, s_cul, col = alpha(col_trmt, 0.5), lwd = 1, type = "l")
        }
    }
    legend("bottomright", trmts, col = col, lwd = 2)
}

plotSNpie = function(dat_sp, dat_plot, col = NA) {
    # TO DO: add check to ensure that col is the same length as treatments
    require(rgl)
    trmts = unique(dat_plot[ , 2])
    if (is.na(col)) 
        col = rainbow(length(trmts))
    S_list = sapply(1:nrow(dat_sp), function(x) 
                    length(which(dat_sp[x, 2:ncol(dat_sp)] != 0)))
    N_list = apply(dat_sp[ , 2:ncol(dat_sp)], 1, sum)
    PIE_list = sapply(1:nrow(dat_sp), function(x) 
                      N_list[x]/(N_list[x] - 1) * 
                      (1 - sum((dat_sp[x, 2:ncol(dat_sp)] / N_list[x])^2)))
    trmt_list = as.character(dat_plot[match(dat_plot[ , 1], dat_sp[ , 1]), 2])
    col_list = sapply(trmt_list, function(x) col[which(trmts == x)])
    plot3d(S_list, N_list, PIE_list, "S", "N", "PIE", col = col_list, size = 8)
} 
