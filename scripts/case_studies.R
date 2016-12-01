library(R.matlab)
library(mobr)

dat = list()
## load case studies and add them to the list dat

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

dat$inv = make_mob_in(comm, data.frame(groups=env, x=spat[,1], y=spat[,2]))

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
dat$mor = make_mob_in(comm, data.frame(groups=env))

# 3) Jon's fire data: ---------------------------

unburned = read.csv("./data/fire_data_unburned.csv")
burned = read.csv("./data/fire_data_burned.csv")

names(burned)[2] <- "Treatment"
unburned$Plot <- toupper(unburned$Plot)
burned$Plot <- tolower(burned$Plot)

# rbind rows

fire <- dplyr::bind_rows(unburned, burned)
fire[is.na(fire)] <- 0
fire$Treatment <- as.character(fire$Treatment)
fire$Treatment[25:48] <- "burned"

## get coordinates:
burned_xy= read.csv("./data/Fire_lat_longs.csv")
names(burned_xy)[1] <- "Plot"

# change first letter into capital
burned_xy$Plot <- tolower(as.character(burned_xy$Plot))

fire$Plot = tolower(fire$Plot)

fire = merge(fire, burned_xy, by='Plot', all=T)

comm = fire[ , 3:23]
spat = fire[ , 24:25]
env = fire$Treatment

## make community matrix:
dat$fire = make_mob_in(comm, data.frame(groups=env, x=spat[ , 2], y=spat[ , 1]))

# 4) coffee -----------------------------------

dat_coffee = read.csv("./data/coffee_comm.csv")
dat_xy = read.csv("./data/coffee_xy.csv")
dat_xy$Treatment <- rep(c("Natural", "Shaded"), each = 3)

comm = dat_coffee[, -1]
spat = dat_xy[,2:3]
env = dat_xy$Treatment

## make community matrix:
dat$coffee = make_mob_in(comm, data.frame(groups=env, x=spat[,1], y=spat[,2]))

# 5) Cattle_tank -----------------

dat_cattle_high = read.csv("./data/Cattle_tank_high.csv")
dat_cattle_high[is.na(dat_cattle_high)] <- 0

rownames(dat_cattle_high) <- dat_cattle_high[,1]
dat_cattle_high <- dat_cattle_high[,-1]
dat_cattle_high <- t(dat_cattle_high)


dat_cattle_low = read.csv("./data/Cattle_tank_low.csv")
dat_cattle_low[is.na(dat_cattle_low)] <- 0

rownames(dat_cattle_low) <- dat_cattle_low[,1]
dat_cattle_low <- dat_cattle_low[,-1]
dat_cattle_low <- t(dat_cattle_low)

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
dat$cattle = make_mob_in(comm, data.frame(groups=env, x=spat[,1], y=spat[,2]))

# 6) Glade ants 2011 ---------------------------
ants = read.csv('./data/glade_ants/source_files/ants 2011 big restored and natural.csv')
plot_info = read.csv('./data/glade_ants/filtered_files/glade_plot_info.csv')
comm = ants[ , 8:42]
comm = ifelse(is.na(as.matrix(comm)), 0, as.matrix(comm))

dat$ants11 = make_mob_in(comm, data.frame(groups = ants$treatment))

# 6) Glade ants 2014 ---------------------------

ants = read.csv('./data/glade_ants/source_files/ants 2014 big rest and nat.csv')
comm = ants[ , 5:45]
comm = ifelse(is.na(as.matrix(comm)), 0, as.matrix(comm))

dat$ants14 = make_mob_in(comm, data.frame(groups = ants$treatment))


## data analysis --------------------------------------------------------------

ref_groups = c('uninvaded', 'before', 'unburned', 'Natural', 'low', 'NAT', 
               'natural')
trt_groups = c('invaded', 'after', 'burned', 'Shaded', 'high', 'BR', 'restored')


tst = vector('list', length(dat))
names(tst) = names(dat)
for(i in seq_along(dat)) {
    cat(paste('Analysis of', names(dat)[i], 'dataset'))
    tst[[i]] =  get_delta_stats(dat[[i]], 'groups', ref_group = ref_groups[i], 
                                type='discrete', log_scale=TRUE, nperm=100)
    stats = mob_stats(dat[[i]], 'groups')
    pdf(paste('./figs/', names(tst)[i], '_results.pdf', sep=''))
        plot_samples(stats$samples)
        plot_groups(stats$groups)
        plot_betaPIE(stats)
        plot_abu(dat[[i]], 'groups', 'rad', pooled=T, log = 'x')
        plot_abu(dat[[i]], 'groups', 'sad', log='x')
        plot(tst[[i]], trt_groups[i], ref_groups[i], same_scale=T) 
    dev.off()
}

save(tst, file='./results/case_study_results.Rdata')
