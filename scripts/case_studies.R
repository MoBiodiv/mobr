library(R.matlab)
source('./R/mobr.R')
source('./R/mobr_boxplots.R')
# 1) invasion data -------------------------
dat_dir = paste('./data/', 'joninvade.mat', sep = '')

dat_matlab = readMat(dat_dir)
dat_plot = as.data.frame(matrix(NA, length(dat_matlab$x), 5))
dat_sp = as.data.frame(matrix(NA, length(dat_matlab$x), nrow(dat_matlab$comary) + 1))

dat_plot[, 1] = 1:nrow(dat_plot)
dat_plot[, 2] = as.vector(ifelse(unlist(dat_matlab$is.invaded) > 0, 'invaded', 'uninvaded'))
dat_plot[, 3] = unlist(dat_matlab$x)
dat_plot[, 4] = unlist(dat_matlab$y)
dat_plot[, 5] = rep(1, nrow(dat_plot))

comm = t(dat_matlab$comary)

# give reasonable names to env. data set:

spat = dat_plot[,3:4]
env = dat_plot[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))

names(comm$env) = 'groups'


tst.inv = get_delta_stats(comm, 'groups', ref_group='uninvaded',
                          type='discrete', log_scale=TRUE,
                          nperm=100)

save(tst.inv, file='./results/tst.inv.Rdata')
#load('./results/tst.inv.Rdata')

stats <- mob_stats(comm, 'groups')

pdf("./figs/tst.inv.pdf")
plot_samples(stats$samples)
plot_groups(stats$groups)
plot_betaPIE(stats)
plotSADs(comm, 'groups')
plot_rarefy(tst.inv)
plot(tst.inv, same_scale=T)
plot(tst.inv, par_args='mfrow=c(1,1)', same_scale=T)
plot_9_panels(tst.inv, 'invaded', 'uninvaded')
plot_9_panels(tst.inv, 'invaded', 'uninvaded',
              par_args='mfrow=c(1,1)')

dev.off()

pdf('./figs/tst.inv.sum.pdf')
plot(tst.inv, same_scale=T)
plot_9_panels(tst.inv, 'invaded', 'uninvaded')
dev.off()

# 2) Morlaix ------------------

dat_dir = paste('./data/', 'morlaix.mat', sep = '')

dat_matlab = readMat(dat_dir)

dat_plot = as.data.frame(matrix(NA, 6, 5))
dat_sp = as.data.frame(matrix(NA, 6, nrow(dat_matlab$ab) + 1))

dat_plot[, 1] = 1:nrow(dat_plot)
dat_plot[, 2] = c(rep('before', 3), rep('after', 3))
dat_plot[, 3] = rep(1, nrow(dat_plot))
dat_plot[, 4] = rep(1, nrow(dat_plot))
dat_plot[, 5] = rep(1, nrow(dat_plot))

comm = t(dat_matlab$ab[, c(1:3, 6:8)])

env = dat_plot[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env))
names(comm$env) = 'groups'

tst.mor = get_delta_stats(comm, 'groups', ref_group='before',
                          type='discrete', log_scale=TRUE,
                          nperm=100)
save(tst.mor, file='./results/tst.mor.Rdata')
#load('./results/tst.mor.Rdata')

stats <- mob_stats(comm, 'groups')

pdf("./figs/tst.mor.pdf")
plot_samples(stats$samples)
plot_groups(stats$groups)
plot_betaPIE(stats)
plotSADs(comm, 'groups')
plot_rarefy(tst.mor)
plot(tst.mor, same_scale=T)
plot(tst.mor, par_args='mfrow=c(1,1)', same_scale=T)
plot_9_panels(tst.mor, 'after', 'before')
plot_9_panels(tst.mor, 'after', 'before',
              par_args='mfrow=c(1,1)')
dev.off()


pdf('./figs/tst.mor.sum.pdf')
plot(tst.mor, same_scale=T)
plot_9_panels(tst.mor, 'after', 'before')
dev.off()

# 3) Jon's fire data: ---------------------------

dat_unburned = read.csv("./data/fire_data_unburned.csv")
head(dat_unburned)
colnames(dat_unburned)
dim(dat_unburned)

dat_burned = read.csv("./data/fire_data_burned.csv")
head(dat_burned)
colnames(dat_burned)
dim(dat_burned)

names(dat_burned)[2] <- "Treatment"
dat_unburned$Plot <- toupper(dat_unburned$Plot)
dat_burned$Plot <- tolower(dat_burned$Plot)

# rbind rows

library(dplyr)

dat.all <- bind_rows(dat_unburned, dat_burned)
str(dat.all)
names(dat.all)

dat.all[is.na(dat.all)] <- 0
dat.all$Treatment <- as.character(dat.all$Treatment)
dat.all$Treatment[25:48] <- "burned"

## get coordinates:
dat_burned_xy= read.csv("./data/Fire_lat_longs.csv")
names(dat_burned_xy)[1] <- "Plot"

# change first letter into capital
dat_burned_xy$Plot <- tolower(as.character(dat_burned_xy$Plot))

dat.all$Plot = tolower(dat.all$Plot)

dat = merge(dat.all, dat_burned_xy, by='Plot', all=T)

comm = dat[ , 3:23]
spat = dat[ , 24:25]
env = dat$Treatment

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,2],
                                      y=spat[,1]))

names(comm$env) = 'groups'

tst.fire = get_delta_stats(comm, 'groups', ref_group='unburned',
                           type='discrete', log_scale=TRUE,
                           nperm=100)
save(tst.fire, file='./results/tst.fire.Rdata')
#load('./results/tst.fire.Rdata')

stats = mob_stats(comm, 'groups')

pdf("./figs/tst.fire.pdf")
plot_samples(stats$samples)
plot_groups(stats$groups)
plot_betaPIE(stats)
plotSADs(comm, 'groups')
plot_rarefy(tst.fire)
plot(tst.fire, same_scale=T)
plot(tst.fire, par_args='mfrow=c(1,1)', same_scale=T)
plot_9_panels(tst.fire, 'burned', 'unburned')
plot_9_panels(tst.fire, 'burned', 'unburned',
              par_args='mfrow=c(1,1)')
dev.off()

pdf('./figs/tst.fire.sum.pdf')
plot(tst.fire, same_scale=T)
plot_9_panels(tst.fire, 'burned', 'unburned')
dev.off()

# 4) coffee -----------------------------------

dat_coffee = read.csv("./data/coffee_comm.csv")
head(dat_coffee)
colnames(dat_coffee)
dim(dat_coffee)

dat_xy = read.csv("./data/coffee_xy.csv")
dat_xy$Treatment <- rep(c("Natural", "Shaded"), each = 3)
head(dat_xy)
colnames(dat_xy)
dim(dat_xy)


comm = dat_coffee[, -1]
spat = dat_xy[,2:3]
env = dat_xy$Treatment

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))

names(comm$env) = 'groups'

tst.coffee = get_delta_stats(comm, 'groups', ref_group='Natural', 
                             type='discrete', log_scale=TRUE,
                             nperm=1000)

save(tst.coffee, file='./results/tst.coffee.Rdata')
#load('./results/tst.coffee.Rdata')

stats = mob_stats(comm, 'groups')
  
pdf("./figs/tst.coffee.pdf")
plot_samples(stats$samples)
plot_groups(stats$groups)
plot_betaPIE(stats)
plotSADs(comm, 'groups')
plot_rarefy(tst.coffee)
plot(tst.coffee, same_scale = T)
plot(tst.coffee, par_args = 'mfrow=c(1,1)', same_scale = T)
plot_9_panels(tst.coffee, 'Shaded', 'Natural') 
plot_9_panels(tst.coffee, 'Shaded', 'Natural', 
              par_args = 'mfrow=c(1,1)')
dev.off()

pdf('./figs/tst.coffee.sum.pdf')
plot(tst.coffee, same_scale = T)
plot_9_panels(tst.coffee, 'Shaded', 'Natural') 
dev.off()


# 5) Cattle_tank -----------------

dat_cattle_high = read.csv("./data/Cattle_tank_high.csv")
head(dat_cattle_high)
colnames(dat_cattle_high)
dim(dat_cattle_high)
dat_cattle_high[is.na(dat_cattle_high)] <- 0

rownames(dat_cattle_high) <- dat_cattle_high[,1]
dat_cattle_high <- dat_cattle_high[,-1]
dat_cattle_high <- t(dat_cattle_high)


dat_cattle_low = read.csv("./data/Cattle_tank_low.csv")
head(dat_cattle_low)
colnames(dat_cattle_low)
dim(dat_cattle_low)
dat_cattle_low[is.na(dat_cattle_low)] <- 0

rownames(dat_cattle_low) <- dat_cattle_low[,1]
dat_cattle_low <- dat_cattle_low[,-1]
dat_cattle_low <- t(dat_cattle_low)

identical(colnames(dat_cattle_low), colnames(dat_cattle_high))

dat_cattle <- rbind(dat_cattle_low, dat_cattle_high)

## coordinates:
dat_cattle_xy = read.csv("./data/Cattle_tank_xy.csv")
dat_cattle_xy$X = tolower(as.character(dat_cattle_xy$X))
rownames(dat_cattle) <- paste0(dat_cattle_xy[,2], dat_cattle_xy[,1])
rownames(dat_cattle_xy) <- rownames(dat_cattle)

comm = dat_cattle
spat = dat_cattle_xy[,3:4]
env = dat_cattle_xy$X

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))
names(comm$env) = 'groups'

tst.cattle = get_delta_stats(comm, 'groups', ref_group='Low',
                             type='discrete', log_scale=TRUE, 
                             nperm=100)
save(tst.cattle, file='./results/tst.cattle.Rdata')
#load('./results/tst.cattle.Rdata')

stats = mob_stats(comm, 'groups')

pdf("./figs/tst.cattle.pdf")
plot_samples(stats$samples)
plot_groups(stats$groups)
plot_betaPIE(stats)
plotSADs(comm, 'groups')
plot_grp_rads(comm, 'groups', log='y')
plot_rarefy(tst.cattle)
plot(tst.cattle, same_scale=T)
plot(tst.cattle, same_scale=T, par_args = 'mfrow=c(1,1)')
plot_9_panels(tst.cattle, 'high', 'Low')
plot_9_panels(tst.cattle, 'high', 'Low',
              par_args = 'mfrow=c(1,1)')
dev.off()

pdf('./figs/tst.cattle.sum.pdf')
plot(tst.cattle, same_scale=T)
plot_9_panels(tst.cattle, 'high', 'Low')
dev.off()

