library(mobr)

data(inv_comm)
data(inv_plot_attr)
inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)

# run analyses
inv_mob_stats = get_mob_stats(inv_mob_in, 'group', n_perm=999)
inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_group='uninvaded',
                           type='discrete', log_scale=TRUE, n_perm=999,
                           overall_p=TRUE)

# create graphics
pdf('./figs/inv_mob_stats.pdf', height = 7*0.5)
plot(inv_mob_stats)
# additional plot of f0 without outline
plot(inv_mob_stats, 'f_0', outline = F)
dev.off()

pdf('./figs/inv_mob_stats_multi.pdf', width = 7*1.25, height = 7*2)
plot(inv_mob_stats, multi_panel = T)
plot(inv_mob_stats, multi_panel = T, outline = F)
dev.off()

pdf('./figs/inv_delta_stats.pdf')
plot(inv_mob_out, 'invaded', 'uninvaded', leg_loc = 'bottomright')
dev.off()

pdf('./figs/inv_mob_stacked.pdf')
overlap_effects(inv_mob_out, 'invaded', leg_loc = NA)
overlap_effects(inv_mob_out, 'invaded', 'stacked', prop=T, leg_loc = NA)
dev.off()
