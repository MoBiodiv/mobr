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

dat_sp = t(dat_matlab$comary)

# give reasonable names to env. data set:

comm = dat_sp
spat = dat_plot[,3:4]
env = dat_plot[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))

names(comm$env) = 'groups'


tst.inv = get_delta_stats(comm, 'groups', ref_group='uninvaded',
                          type='discrete', log_scale=TRUE,
                          nperm=100)

save(tst.inv, file='./results/tst.inv.Rdata')

stats <- mob_stats(comm, 'groups')

pdf("tst.inv.pdf", height = 5, width = 10)
plot_samples(stats$samples)
plot_groups(stats$groups)
plot_betaPIE(stats)
plot_rarefy(tst.inv)

plot(tst.inv)
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

dat_sp[, 1] = 1:nrow(dat_plot)
dat_sp[, 2:ncol(dat_sp)] = t(dat_matlab$ab[, c(1:3, 6:8)])

comm = dat_sp
env = dat_plot[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env))
names(comm$env) = 'groups'

tst.mor = get_delta_stats(comm, 'groups', ref_group='before',
                          type='discrete', log_scale=TRUE,
                          nperm=100)
save(tst.mor, file='./results/tst.mor.Rdata')

pdf("tst.mor.pdf", height = 5, width = 10)
plot.mobr(tst.mor)
boxplot.comm(comm, "groups")
plot_rarefy(tst.mor)

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

dat.all.2 <- dat.all[,-c(23)]
tail(dat.all.2)

dat.all.2[is.na(dat.all.2)] <- 0
dat.all.2$Treatment <- as.character(dat.all.2$Treatment)
dat.all.2$Treatment[25:48] <- "burned"

## get coordinates:
dat_burned_xy= read.csv("./data/Fire_lat_longs.csv")
names(dat_burned_xy)[1] <- "Plot"

# change first letter into capital
dat_burned_xy$Plot <- as.character(dat_burned_xy$Plot)

dat_burned_xy$Plot[1:24] <- tolower(dat_burned_xy$Plot)[1:24]
dat_burned_xy$Plot[25:48] <- toupper(dat_burned_xy$Plot)[25:48]

match(sort(dat.all.2$Plot), sort(dat_burned_xy$Plot))

sort(dat.all.2$Plot)
sort(dat_burned_xy$Plot)

## which coordinates are potentially strange

dat <- inner_join(dat.all.2[,1:2], dat_burned_xy)

comm = dat.all.2[,3:23]
spat = dat_burned_xy[,2:3]
env = dat.all.2[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1],
                                      y=spat[,2]))

names(comm$env) = 'groups'

tst.fire = get_delta_stats(comm, 'groups', ref_group='unburned',
                           type='discrete', log_scale=TRUE,
                           nperm=100)
save(tst.fire, file='./results/tst.fire.Rdata')

pdf("tst.fire.pdf", height = 5, width = 10)
plot.mobr(tst.fire)
boxplot.comm(comm, "groups")
plot_rarefy(tst.fire)
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


comm = dat_coffee[, 2:22]
spat = dat_xy[,2:3]
env = dat_xy[,4]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))

names(comm$env) = 'groups'

tst.coffee = get_delta_stats(comm, 'groups', ref_group='Natural', type='discrete', log_scale=TRUE, inds=NULL, nperm=1000)

save(tst.coffee, file='./results/tst.coffee.Rdata')

pdf("tst.coffee.pdf", height = 5, width = 10)
plot.mobr(tst.coffee)
boxplot.comm(comm, "groups")
plot_rarefy(tst.coffee)
dev.off()

# 5) Cattle_tank -----------------

dat_cattle_high= read.csv("./data/Cattle_tank_high.csv")
head(dat_cattle_high)
colnames(dat_cattle_high)
dim(dat_cattle_high)
dat_cattle_high[is.na(dat_cattle_high)] <- 0

rownames(dat_cattle_high) <- dat_cattle_high[,1]
dat_cattle_high <- dat_cattle_high[,-1]
dat_cattle_high <- t(dat_cattle_high)


dat_cattle_low= read.csv("./data/Cattle_tank_low.csv")
head(dat_cattle_low)
colnames(dat_cattle_low)
dim(dat_cattle_low)
dat_cattle_low[is.na(dat_cattle_low)] <- 0

rownames(dat_cattle_low) <- dat_cattle_low[,1]
dat_cattle_low <- dat_cattle_low[,-1]
dat_cattle_low <- t(dat_cattle_low)

identical(colnames(dat_cattle_low), colnames(dat_cattle_high))

dat_cattle <- rbind(dat_cattle_low, dat_cattle_high)
rownames(dat_cattle) <- paste0(dat_cattle_xy[,2], dat_cattle_xy[,1])

## coordinates:
dat_cattle_xy= read.csv("./data/Cattle_tank_xy.csv")
rownames(dat_cattle_xy) <- rownames(dat_cattle)

comm = dat_cattle
spat = dat_cattle_xy[,3:4]
env = dat_cattle_xy[,1]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))
names(comm$env) = 'groups'

tst.cattle = get_delta_stats(comm, 'groups', ref_group='Low', type='discrete', log_scale=TRUE, inds=NULL, nperm=100)
save(tst.cattle, file='./results/tst.cattle.Rdata')

pdf("tst.cattle.pdf", height = 5, width = 10)
plot.mobr(tst.cattle)
boxplot.comm(comm, "groups")
plot_rarefy(tst.cattle)
dev.off()
