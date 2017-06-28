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
plot_info = read.csv('./data/glade_ants/filtered_data/ant_plot_attr.csv')
comm = ants[ , 8:42]
comm = ifelse(is.na(as.matrix(comm)), 0, as.matrix(comm))

dat$ants11 = make_mob_in(comm, data.frame(groups = ants$treatment))

# 6) Glade ants 2014 ---------------------------

ants = read.csv('./data/glade_ants/source_files/ants 2014 big rest and nat.csv')
comm = ants[ , 5:45]
comm = ifelse(is.na(as.matrix(comm)), 0, as.matrix(comm))

dat$ants14 = make_mob_in(comm, data.frame(groups = ants$treatment))

# 7) Portal plants 1989-2002 ----------------------
# species data
plant = read.csv('./data/PortalData/Plants/Portal_plant_1981_2015.csv')

# subset to 1989-present
plant = subset(plant, year >= 1989)
table(plant$year)

# drop any species name with unk in it
plant$species = as.character(plant$species)
plant = subset(plant, !grepl('unk', species))
# dropping names with sp is more complex b/c of "Sida_spin"
sort(unique(plant$species))

# site data
utms = read.csv('./data/PortalData/SiteandMethods/Portal_UTMCoords.csv')
utms = subset(utms, type == 'quadrat')
utms = utms[ , c('plot', 'number', 'east', 'north')]
names(utms) = c('plot', 'quadrat', 'east', 'north')
# fill in two missing sites these are guesses!!
utms = rbind(utms, data.frame(plot=c(1,24),
                              quadrat=c(77, 17),
                              east=c(max(utms[utms$plot==1, 'east']),
                                     max(utms[utms$plot==24, 'east'])),
                              north=c(min(utms[utms$plot==1, 'north']),
                                      max(utms[utms$plot==24, 'north']))))


sites = read.csv('./data/PortalData/SiteandMethods/Portal_plots.csv')
#trts = read.csv('./data/PortalData/SiteandMethods/Portal_plot_treatments.csv')
names(sites) = c('year', 'plot', 'treatment')

sites = subset(sites, year >= 1989 & year< 2005)
trt = with(sites, tapply(plot, list(plot, yr, treatment), length))

dim(trt)
dimnames(trt)
trt[1, , ]
trt[9, , ]
apply(trt, c(1, 3), sum, na.rm=T)
# the above shows that each site in this range 1989 to 2004 had the same 
# treatment for the entire 16 year time span
# one can go out to 26 years of sampling but two plots are then dropped

# subset species data to this range
plant = subset(plant, year >= 1989 & year< 2005)
plant = merge(plant, sites)
plant$id = with(plant, paste(plot, quadrat, year, season, sep='_'))
head(plant)
plant$species
plant_comm = tapply(plant$abundance, list(plant$id, plant$species), sum)
plant_comm = ifelse(is.na(plant_comm), 0, plant_comm)
plant_comm[1:5, 1:5]

id = row.names(plant_comm)
plant_attr = as.data.frame(t(sapply(id, function(x) unlist(strsplit(x, '_')))))
names(plant_attr) = c('plot', 'quadrat', 'year', 'season')
plant_attr = merge(plant_attr, sites, all.x=T, all.y=F)
plant_attr = merge(plant_attr, utms, all.x=T, all.y=F)
names(plant_attr) = c('plot', 'quadrat', 'year', 'season', 'groups',
                      'x', 'y')
row.names(plant_attr) = row.names(plant_comm)

true = with(plant_attr, season == 'summer' & year == 2004)
dat$portal_sum = make_mob_in(plant_comm[true, ], plant_attr[true, ])
true = with(plant_attr, season == 'winter' & year == 2004)
dat$portal_win = make_mob_in(plant_comm[true, ], plant_attr[true, ])


## data analysis --------------------------------------------------------------

ref_groups = c('uninvaded', 'before', 'unburned', 'Natural', 'low', 'NAT', 
               'natural', 'control', 'control')
trt_groups = c('invaded', 'after', 'burned', 'Shaded', 'high', 'BR', 'restored',
               'removal', 'removal')


tst = stats = vector('list', length(dat))
names(tst) = names(stats) = names(dat)

for(i in seq_along(dat)) {
    cat(paste('Analysis of', names(dat)[i], 'dataset'))
    
    tst[[i]] =  get_delta_stats(dat[[i]], 'groups', ref_group = ref_groups[i], 
                                type='discrete', log_scale=TRUE, nperm=200)
    stats[[i]] = get_mob_stats(dat[[i]], 'groups', nperm=200)

    pdf(paste('./figs/', names(tst)[i], '_results.pdf', sep=''))
        plot(stats[[i]], multipanel = T, col=c('red', 'blue')) 
        plot_abu(dat[[i]], 'groups', 'rad', pooled=T, log = 'x')
        plot_abu(dat[[i]], 'groups', 'sad', log='x')
        plot(tst[[i]], trt_groups[i], ref_groups[i], same_scale=T)
        par(mfrow=c(1,1))
        stack_effects(tst[[i]], trt_groups[i])
        stack_effects(tst[[i]], trt_groups[i], prop=TRUE)
    dev.off()
}

save(tst, file='./results/case_study_results.Rdata')


## Gentry continous example ------------------------------------
load('./data/gentry.Rdata')
names(gentry)
lapply(gentry, head)
# record filtering ----------
# remove NA count
gentry$counts = gentry$counts[!is.na(gentry$counts$count), ]
# remove if no line recorded
gentry$counts = gentry$counts[!is.na(gentry$counts$line), ]
# remove if line is recorded as 0 or 11
gentry$counts = gentry$counts[gentry$counts$line != 0 & 
                              gentry$counts$line != 11, ]
# kraft et al (2011) removed Poaceae and ferns
# drop psuedo species
nrow(gentry$counts)
sum(gentry$counts$count)
nrow(gentry$stems)
sum(gentry$stems$stem > 0)

# 
summary(gentry$counts$count)
comm_local = with(gentry$counts, 
              tapply(count, list(site_code, line, species_id), sum))
comm_local = ifelse(is.na(comm_local), 0, comm_local)

# remove sites that did not sample all 10 lines
true = apply(apply(comm_local, 1:2, sum) > 0, 1, sum) == 10
comm_local = comm_local[true, , ]

stat_local = apply(comm_local, 1:2, function(x)
                   c(sum(x), sum(x > 0), diversity(x, index='simpson'),
                     rarefaction(x, 'indiv', effort=5),
                     rarefaction(x, 'indiv', effort=10),
                     rarefaction(x, 'indiv', effort=20),
                     rarefaction(x, 'indiv', effort=40),
                     rarefaction(x, 'indiv', effort=80)))
#dimnames(stat_local)[[1]] = c('N_local', 'S_local', 'PIE_local', 'rare_local')

comm_region = apply(comm_local, c(1, 3), sum)

stat_region = apply(comm_region, 1, function(x)
                    c(sum(x), sum(x > 0), diversity(x, index='simpson'),
                      rarefaction(x, 'indiv', effort=43)))
dimnames(stat_region)[[1]] = c('N_region', 'S_region', 'PIE_region', 'rare_region')
stat_region = data.frame(site=colnames(stat_region), 
                         t(stat_region))
rownames(stat_region) = NULL
head(stat_region)

stat_local  = apply(stat_local, 1, function(x)
                    reshape2::melt(x))
stat_local = data.frame(stat_local[[1]], 
                        S_local = stat_local[[2]][ , 3],
                        PIE_local = stat_local[[3]][ , 3], 
                        rare_5 = stat_local[[4]][ , 3],
                        rare_10 = stat_local[[5]][ , 3],
                        rare_20 = stat_local[[6]][ , 3],
                        rare_40 = stat_local[[7]][ , 3],
                        rare_80 = stat_local[[8]][ , 3])

names(stat_local) = c('site', 'line', 'N_local', 'S_local', 'PIE_local',
                      'rare_5', 'rare_10', 'rare_20', 'rare_40', 'rare_80')
head(stat_local)

gentry_div = merge(stat_local, stat_region)
gentry_div = merge(gentry_div, gentry$sites, 
                   by.x='site', by.y='abbreviation',
                   all.x=T, all.y=F)
#gentry_div2 = aggregate(gentry_div, by=list(gentry_div$line), mean)

pdf('./figs/gentry_panels.pdf', height=7*2)
par(mfrow=c(4,2))
#S
plot(S_local ~ abs(lat), data=gentry_div, ylim=c(1, 300), 
     col='pink', ylab='S')
points(S_region ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat), S_local), col='red'))
with(gentry_div, lines(lowess(abs(lat), S_region), col='blue'))
#logS
plot(log10(S_local) ~ abs(lat), data=gentry_div, ylim=c(0, 3), 
     col='pink', ylab='log10 S')
points(log10(S_region) ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat), log10(S_local)), col='red'))
with(gentry_div, lines(lowess(abs(lat), log10(S_region)), col='blue'))
#rare
plot(rare_5 ~ abs(lat), data=gentry_div, ylim=c(1, 300),
     col='pink', ylab='S rarefied')
points(rare_10 ~ abs(lat), data=gentry_div, col='orange')
points(rare_20 ~ abs(lat), data=gentry_div, col='green3')
points(rare_40 ~ abs(lat), data=gentry_div, col='purple')
points(S_region ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_5)]),
                              rare_5[!is.na(rare_5)]), col='pink'))
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_10)]),
                              rare_10[!is.na(rare_10)]), col='orange'))
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_20)]),
                              rare_20[!is.na(rare_20)]), col='green3'))
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_40)]),
                              rare_40[!is.na(rare_40)]), col='purple'))
with(gentry_div, lines(lowess(abs(lat), S_region), col='dodgerblue',
                       lwd=2))
legend('topright', c('n=N', 'n=40','n=20','n=10','n=5'), 
       col=c('dodgerblue', 'purple', 'green3', 'orange', 'pink'),
       pch=1, bty='n')
#lograre
plot(log10(rare_5) ~ abs(lat), data=gentry_div, ylim=c(0, 3), 
     col='pink', ylab='log10 S rarefied')
points(log10(rare_10) ~ abs(lat), data=gentry_div, col='orange')
points(log10(rare_20) ~ abs(lat), data=gentry_div, col='green3')
points(log10(rare_40) ~ abs(lat), data=gentry_div, col='purple')
points(log10(S_region) ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_5)]),
                              log10(rare_5[!is.na(rare_5)])), col='pink'))
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_10)]),
                              log10(rare_10[!is.na(rare_10)])), col='orange'))
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_20)]),
                              log10(rare_20[!is.na(rare_20)])), col='green3'))
with(gentry_div, lines(lowess(abs(lat[!is.na(rare_40)]),
                              log10(rare_40[!is.na(rare_40)])), col='purple'))
with(gentry_div, lines(lowess(abs(lat), log10(S_region)), col='dodgerblue',
                       lwd=2))
legend('topright', c('n=N', 'n=40','n=20','n=10','n=5'), 
       col=c('dodgerblue', 'purple', 'green3', 'orange', 'pink'),
       pch=1, bty='n')
#N
plot(N_local ~ abs(lat), data=gentry_div, ylim=c(1, 1000), 
     col='pink', ylab='N')
points(N_region ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat), N_local), col='red'))
with(gentry_div, lines(lowess(abs(lat), N_region), col='blue'))
#logN
plot(log10(N_local) ~ abs(lat), data=gentry_div, ylim=c(0, 3), 
     col='pink', ylab='log10 N')
points(log10(N_region) ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat), log10(N_local)), col='red'))
with(gentry_div, lines(lowess(abs(lat), log10(N_region)), col='blue'))
#
plot(PIE_local ~ abs(lat), data=gentry_div, ylim=c(0, 1),
     col='pink',  ylab='PIE')
points(PIE_region ~ abs(lat), data=gentry_div, col='dodgerblue')
with(gentry_div, lines(lowess(abs(lat), PIE_local), col='red'))
with(gentry_div, lines(lowess(abs(lat), PIE_region), col='blue'))
dev.off()

       