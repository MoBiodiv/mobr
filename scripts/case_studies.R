library(mobr)
library(betaC)
library(tidyverse)
library(devtools)
load_all()
# copy 4 functions from beta_cartoon.R 
get_curves <- function(mob_in, group_var) {
  # drop rows with no occurrences
  mob_in <- subset(mob_in, subset = rowSums(mob_in$comm) > 0, type = 'logical')
  groups <- mob_in$env[ , group_var]
  group_lvs <- unique(groups)
  uber <- vegan::rarefy(colSums(mob_in$comm), sample = 1:sum(mob_in$comm),
                        se = TRUE)
  out <- data.frame(scale = 'uber', group = 0, sample = 0,
                    N = attributes(uber)$Subsample,
                    S = uber[1, ], se = uber[2, ])
  # compute group level pooled curves
  N_lv <- rep(0, length(group_lvs))
  r_dat_pool <- array(NA, c(2, length(group_lvs), sum(mob_in$comm)))
  for (i in seq_along(group_lvs)) {
    comm_lv <- mob_in$comm[groups == group_lvs[i], ]
    N_lv[i] <- sum(comm_lv)
    r_dat_pool[ , i, 1:N_lv[i]] <- vegan::rarefy(colSums(comm_lv), sample = 1:N_lv[i],
                                   se = TRUE)
    out <- rbind(out, 
                 data.frame(scale = 'pool',
                            group = group_lvs[i],
                            sample = 0, 
                            N = 1:N_lv[i],
                            S = r_dat_pool[1, i, 1:N_lv[i]],
                            se = r_dat_pool[2, i, 1:N_lv[i]]))
    r_dat <- array(NA, c(2, nrow(comm_lv), max(rowSums(comm_lv))))
    # now sample plots within groups
    for (j in 1:nrow(comm_lv)) {
      Nj <- sum(comm_lv[j, ]) 
      r_dat[ , j, 1:Nj] <- vegan::rarefy(comm_lv[j, ], sample = 1:Nj,
                             se = TRUE)
      out <- rbind(out, 
                 data.frame(scale = 'plot',
                            group = group_lvs[i],
                            sample = j,
                            N = 1:Nj,
                            S = r_dat[1, j, 1:Nj],
                            se = r_dat[2, j, 1:Nj]))
    }
    # average across plots within a group
    Nmin <- min(rowSums(comm_lv))
    out <- rbind(out, 
                 data.frame(scale = 'plot',
                            group = group_lvs[i],
                            sample = 'avg', 
                            N = 1:Nmin,
                            S = colMeans(r_dat[1, , 1:Nmin], na.rm = TRUE),
                            se = colMeans(r_dat[2, , 1:Nmin], na.rm = TRUE)))
  }
  # average pooled curves across groups
  Nmin <- min(N_lv)
  out <- rbind(out, 
               data.frame(scale = 'pool',
                          group = 'avg',
                          sample = 0,
                          N = 1:Nmin,
                          S = colMeans(r_dat_pool[1, , 1:Nmin], na.rm=TRUE),
                          se = colMeans(r_dat_pool[2, , 1:Nmin], na.rm =TRUE)))
  # remove any NaN 
  out <- out[!is.na(out$S), ]
  out$group = as.factor(out$group)
  out$sample = as.factor(out$sample)
  return(out)
}
add_CI <- function(dat, col='pink') {
    polygon(c(dat$N, rev(dat$N)), c(dat$S + dat$se, rev(dat$S - dat$se)),
            col = col, border = NA)
}
plot_curves2 = function(dat, log='', lwd=2, cols = c('red', 'dodgerblue'), CI = FALSE,
                        CI_scale = c('gamma', 'study'), ...) {
    scales = unique(dat$scale)
    plot(S ~ N, data = dat, log = log, type = 'n', ...)
    if (CI) 
        add_CI(subset(dat, subset = scale == CI_scale), col = 'grey80')
    for (i in seq_along(scales)) { 
        lines(S ~ N, data = dat, subset = scale == scales[i], col = cols[i],
              lwd = lwd)
    }
        
}
get_turnover <- function(mob_in, method = 'bray', binary = FALSE) {
    groups <- mob_in$env[ , 'yr']
    sites <- mob_in$env$site
    comm <- mob_in$comm
    # compute all possible pairwise turnovers
    D <- as.matrix(vegdist(comm, method = method, binary = binary))
    # reformat into data.frame
    dat <- data.frame()
    for (i in 1:(nrow(D) - 1)) { # rows
        for (j in (i + 1):ncol(D)) { # columns
            dat <- rbind(dat,
                         data.frame(grain = 1, yr1 = groups[i], yr2 = groups[j],
                                    site1 = sites[i], site2 = sites[j], turn = D[i, j]))
        }
    }
    # compute for 2 sample grain on site id (i.e., spatial beta)
    # compare composition between sites long time
    turn <- as.numeric(vegdist(aggregate(comm, list(mob_in$env$site), sum)[ , -1],
                                            method = method, binary = binary))
    dat <- rbind(dat, 
                 data.frame(grain = 2, yr1 = NA, yr2 = NA, site1 = 1, site2 = 2,
                            turn))
    # compute for 2 sample grain on year id (i.e., temporal beta)
    # that is compare composition between years at regional scale
    turn <- as.numeric(vegdist(aggregate(comm, list(groups), sum)[ , -1],
                                            method = method, binary = binary))
    dat <- rbind(dat, 
                 data.frame(grain = 2, yr1 = 1, yr2 = 2, site1 = NA, site2 = NA,
                            turn))


    # reorder such that yr1, yr2, site 1, site 2, w/in yr, between yr, big space, big time
    dat <- dat[c(1, 6, 2, 5, 7, 8), ]
    return(dat)
}

# Konza Patch-burn experiment Grasshopper community data url:
# http://lter.konza.ksu.edu/content/pbg07-grasshopper-species-abundances-patch-burn-grazing-experiment-konza-prairie
# Metadata ---------------------------------------------------------------------
# Effects of patch-burn grazing management on plant and animal diversity and the
# nature and variety of wildlife habitat are being assessed in two replicate
# management units, each consisting of three pastures (watersheds) designated
# C03A/C03B/C03C and C3SA/C3SB/C3SC. In each patch-burn grazing unit, one
# watershed is burned and two that are left unburned in a given year. The
# burning treatments are rotated annually so that each pasture is burned every
# third year. Each patch-burn grazing unit is paired with an annually-burned
# pasture for comparison with traditional grazing systems (C01A and C1SB).
# Methods: Location of Sampling Stations: Grasshopper density is determined on
# upland topographic locations. C3A, C3B, C3C, C1A, C3SA, C3SB, C3SC, C1B. 4
# sits per watershed, 4 transects per site. Frequency of Sampling: Grasshopper
# abundances are sampled once in late summer (August-September), with each site
# sample twice in a season a week apart. Variable Measured: Number of
# individuals (categorized by instar) for individual grasshopper species.
# Methods: Three watershed units (C3A, C3B, C3C) constitute 'patches' that are
# alternately burned in a 3-year rotation within a single, fenced pasture (i.e.,
# patch-burn grazing). Two additional watersheds are serving as controls: a
# grazed, annually/uniformly-burned watershed (C1A) and an ungrazed,
# annually/uniformly-burned watershed (1D). Grasshopper sampling is done by
# standardized sweeping with canvas beating nets 38 cm in diameter. A sample of
# 250 sweeps (ten sets of 25 sweeps each) is taken at each site (4 independent
# sites per watershed) on each occasion.
# DJM notes --------------------------------------------------------------------
# the watershed codes are a mess, here I subset down to a smaller
# set of Watersheds that are clearly assigned to a specific treatment 
# C3C and C1A are neighboring watersheds the first is in the patch burn rotation
# the later is annual burn. 
# Data import and clean --------------------------------------------------------
dat <- read.csv('http://lter.konza.ksu.edu/sites/default/files/data/PBG072.csv')
head(dat)
table(dat$Watershed)
# I assume codes ("0c1a" and "c01a") and ("0c3c" and "c03c") are equivalent designations
dat_sub = subset(dat, Watershed %in% c("0c1a","c01a", "0c3c", "c03c"),
                 drop = TRUE)
table(dat_sub$Watershed)
dat_sub$group <- with(dat_sub, ifelse(Watershed == '0c1a' | Watershed == 'c01a', 
                                      'control', 'treatment'))
table(dat_sub$group)
head(dat_sub)
table(dat_sub$Comments)
# it appears there are no incomplete surveys in this subset based on comments
summary(dat_sub$Total)
# not clear why some of the totals are NA
dat_sub[is.na(dat_sub$Total), ]
# but it appears from the above that they don't indicated bad data per se
# manually calculating transect totals 
calc_tot <- rowSums(dat_sub[ , paste0("S", 1:10)], na.rm = TRUE)
plot(dat_sub$Total, calc_tot)
abline(a = 0, b = 1)
summary(calc_tot - dat_sub$Total)
dat_sub[dat_sub$Total != calc_tot, ]
# so there a few totals that are miscalculated slightly so use manual totals
dat_sub$Total = calc_tot
# check for unidentified species
table(dat_sub$Species)
# they are coded as Spcode == 41
# drop unidenifed species
dat_sub <- subset(dat_sub, Spcode != 41)
table(dat_sub$Spcode)

head(dat_sub)
table(dat_sub$Recyear)
table(dat_sub$Recmonth)
dat_sub$Date <- with(dat_sub, 
                     as.Date(paste(Recyear, Recmonth, Recday, sep = "-"), 
                             format = "%Y-%m-%d"))
table(dat_sub$Date)
# that indicates that really sampling was just once a year and
# check the Repsite variable
table(dat_sub$group, dat_sub$Repsite)
# change all of Repsite to lowercase
dat_sub$Repsite <- tolower(dat_sub$Repsite)
table(dat_sub$group, dat_sub$Repsite)

# add "sp" to Speciescode to make string
dat_sub$Spcode <- paste0('sp', dat_sub$Spcode)

# recode Repsite so that unique letters occur in each group because
# Repsite "a" is not the same site in control and treatment categories so
# it needs a unique letter id
for (i in 1:nrow(dat_sub)) { 
    if (dat_sub$group[i] == 'treatment') {
        if (dat_sub$Repsite[i] == 'a') tmp <- 'e'
        if (dat_sub$Repsite[i] == 'b') tmp <- 'f'
        if (dat_sub$Repsite[i] == 'c') tmp <- 'g'
        if (dat_sub$Repsite[i] == 'd') tmp <- 'h'
        dat_sub$Repsite[i] <- tmp
    }
}
table(dat_sub$Repsite)

# add simple before after variable and drop year 1
dat_sub <- subset(dat_sub, Recyear != 2010)
dat_sub$Recyear <- with(dat_sub, ifelse(Recyear < 2015, 'past', 'modern'))

# drop years except for first and last year to simply analysis
#table(dat_sub$Recyear)
#dat_sub <- subset(dat_sub, Recyear %in% c(2010, 2018))
#dat_sub$Recyear

# simplify further by dropping Repsite and aggregating data into just
# patch-burn vs annual burn
#dat_sub <- subset(dat_sub, Repsite %in% c('a', 'e'))
#table(dat_sub$group)


# data exploration -------------------------------------------------------------
# make site x species matrix 
# define relevant variables
#core_vars <- c('Recyear', 'Repsite', 'group', 'Spcode', 'Total')
# ignore repsite and recyear instead only use time and group as temporal and spatial indicator
core_vars <- c('Recyear', 'group', 'Spcode', 'Total')
comm <- tidyr::pivot_wider(dat_sub[ , core_vars],
                           names_from = Spcode, 
                           values_from = Total,
                           values_fill = 0,
                           values_fn = sum)
comm

# analyze biodiversity

#dat_mob_in <- make_mob_in(comm[ , -c(1:3)], comm[ , 1:3])
dat_mob_in <- make_mob_in(comm[ , -c(1:2)], comm[ , 1:2])
dat_mob_in

par(mfrow = c(3,2))
# group effect
plot_rarefaction(dat_mob_in, group_var = 'group', ref_level = 'control', 
                 method = 'IBR', log = 'xy')
# spatial effect
plot_rarefaction(dat_mob_in, group_var = 'Repsite', 
                 method = 'IBR', log = 'xy')
# temporal effect
#cols <- colorRampPalette(c('blue','red'))(9)
plot_rarefaction(dat_mob_in, group_var = 'Recyear', 
                 method = 'IBR', log = 'xy')#,
#                 col = cols)


turn_yr <- get_curves(dat_mob_in, 'Recyear')
#turn_si <- get_curves(dat_mob_in, 'Repsite')
turn_gp <- get_curves(dat_mob_in, 'group')
head(turn_yr)
tail(turn_yr)
maxS <- max(turn_yr$S)

pdf('./figs/konza_curves.pdf', width = 7*1.25, height = 7)
## now curves
#sites 1-8 for each year
#yrs <- sort(unique(dat_mob_in$env$Recyear))
yrs <- c('past', 'modern')
#par(mfrow=c(3,3))
par(mfrow = c(1,2))
for (i in seq_along(yrs))
  plot_curves2(subset(turn_yr, 
                   (scale == 'plot' | scale == 'pool') & group == yrs[i] & (sample == '0' | sample == 'avg')),
             ylim = c(1, maxS), log = '', main = paste('spatial beta \nyear =', yrs[i]))
#each site across the 9 years
#sites <- levels(as.factor(dat_mob_in$env$Repsite))
#for (i in seq_along(sites))
#  plot_curves2(subset(turn_si, 
#                   (scale == 'plot' | scale == 'pool') & group == sites[i] & (sample == '0' | sample == 'avg')),
#             ylim = c(1, maxS), log = 'xy', main = paste('site =', sites[i]))
#treatments across the 9 years
par(mfrow=c(1,2))
grps <- levels(as.factor(dat_mob_in$env$group))
for (i in seq_along(grps))
  plot_curves2(subset(turn_gp, 
                   (scale == 'plot' | scale == 'pool') & group == grps[i] & (sample == '0' | sample == 'avg')),
             ylim = c(1, maxS), log = '', main = paste('temporal beta \ngroup =', grps[i]))
dev.off()

pdf('./figs/konza_curves_study_scale.pdf')        
par(mfrow = c(1, 1))
plot_curves2(subset(turn_yr,
                    scale == 'uber' | (scale == 'pool' & group == 'avg')),
             ylim = c(1, maxS), lwd = 3, log = '', main = 'Year Differences in Composition',
             CI = TRUE, CI_scale = 'uber')

#plot_curves2(subset(turn_si,
#                    scale == 'uber' | (scale == 'pool' & group == 'avg')),
#             ylim = c(1, maxS), lwd = 3, log = 'xy', main = 'Spatial Effect',
#             CI = TRUE, CI_scale = 'uber')
plot_curves2(subset(turn_gp,
                    scale == 'uber' | (scale == 'pool' & group == 'avg')),
             ylim = c(1, maxS), lwd = 3, log = '', main = 'Spatial differences in Composition',
             CI = TRUE, CI_scale = 'uber')
dev.off()


# pairwise turnover metrics
D_comm <- vegdist(dat_mob_in$comm, method = 'bray')
D_grp <- dist(as.numeric(as.factor(dat_mob_in$env$group)))
D_yr <- dist(as.numeric(as.factor(dat_mob_in$env$Recyear)))
D_si <- dist(as.numeric(as.factor(dat_mob_in$env$Repsite)))
cor(D_grp, D_comm)
cor(D_yr, D_comm)
cor(D_si, D_comm)
# suggests that yes there is strong change through time but not space
# ordination analysis 
library(vegan)
expl <- data.frame(Recyear = as.factor(dat_mob_in$env$Recyear),
                   Repsite = as.factor(dat_mob_in$env$Repsite),
                   group = as.factor(dat_mob_in$env$group))
rda_yr <- rda(dat_mob_in$comm, expl$Recyear)
rda_si <- rda(dat_mob_in$comm, expl$Repsite)
rda_gp <- rda(dat_mob_in$comm, expl$group)
anova(rda_yr)
anova(rda_si)
anova(rda_gp)
rda_yrsi <- rda(dat_mob_in$comm ~ expl$Recyear + expl$Repsite)
anova(rda_yrsi, by = 'margin')
rda_yrgp <- rda(dat_mob_in$comm ~ expl$Recyear + expl$group)
anova(rda_yrgp, by = 'margin')
# suggests that yes there is strong change through time but not space

# one day would be good to define averaging function across groups
#for(i in seq_along(tstl)) tstm[1:length(tstl[[i]]), i] <- tstl[[i]]


# for newer function usage setup nested dat structure

dat_nest <- comm  %>% 
            nest(comm = starts_with("sp"),
                 time = Recyear, 
                 space = Repsite)

index = c('S', 'S_PIE', 'S_n')
effort = 5
  
dat_nest <- dat_nest %>%
            mutate(div = map(comm, calc_comm_div, index = index, effort = effort))

# year group
div_yr <- comm[ , -(2:3)] %>%
    nest(comm = starts_with("sp")) %>% 
  mutate(div = map(comm, calc_comm_div, index = index, effort = effort)) %>%
    unnest(. , cols = div)
# sites group
div_si <- comm[ , -c(1, 3)] %>%
    nest(comm = starts_with("sp")) %>% 
  mutate(div = map(comm, calc_comm_div, index = index, effort = effort)) %>%
    unnest(. , cols = div)
# group group
div_gp <- comm[ , -(1:2)] %>%
    nest(comm = starts_with("sp")) %>% 
  mutate(div = map(comm, calc_comm_div, index = index, effort = effort)) %>%
    unnest(. , cols = div)

group_sad <- aggregate(comm[ , -c(1:3)], list(comm$group), sum)
group_sad

study_div <- calc_comm_div(group_sad[ , -1], index = index, effort = effort,
                           scales = c('gamma', 'beta'))
study_div
# merge this with other matrices
div_yr <- rbind(div_yr[ , -2], data.frame(Recyear = 'study', study_div))
div_si <- rbind(div_si[ , -2], data.frame(Repsite = 'study', study_div))
div_gp <- rbind(div_gp[ , -2], data.frame(group = 'study', study_div))
div_gp$group = factor(div_gp$group, levels = c('control', 'treatment', 'study'))


pdf('./figs/konza_metrics.pdf')
ggplot(subset(div_yr, index != 'beta_C'), aes(as.factor(Recyear), value)) + 
  geom_boxplot(aes(fill = scale)) + 
  facet_wrap(. ~ index, nrow = 4, ncol = 2, scales = 'free')
boxplot(value ~ Recyear, data = div_yr, subset = index == 'beta_C',
        ylab = 'beta_C')

ggplot(subset(div_si, index != 'beta_C'), aes(Repsite, value)) + 
  geom_boxplot(aes(fill = scale)) + 
  facet_wrap(. ~ index, nrow = 4, ncol = 2, scales = 'free')
boxplot(value ~ Repsite, data = div_si, subset = index == 'beta_C',
        ylab = 'beta_C')

ggplot(subset(div_gp, index != 'beta_C'), aes(group, value)) + 
  geom_boxplot(aes(fill = scale)) + 
  facet_wrap(. ~ index, nrow = 4, ncol = 2, scales = 'free')
boxplot(value ~ group, data = div_gp, subset = index == 'beta_C',
        ylab = 'beta_C')
dev.off()


pdf('./figs/konza_metrics_barplots.pdf')
ggplot(subset(div_yr, index != 'beta_C'), aes(as.factor(Recyear), value)) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge", aes(fill = scale)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = .2) +
  facet_wrap(. ~ index, nrow = 4, ncol = 2, scales = 'free')
barplot(value ~ Recyear, data = subset(div_yr_sum, index == 'beta_C'))

ggplot(subset(div_si, index != 'beta_C'), aes(Repsite, value)) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge", aes(fill = scale)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = .2) +
  facet_wrap(. ~ index, nrow = 4, ncol = 2, scales = 'free')
barplot(value ~ Repsite, data = subset(div_si, index == 'beta_C'))


ggplot(subset(div_gp, index != 'beta_C'), aes(group, value)) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge", aes(fill = scale)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width = .2) +
  facet_wrap(. ~ index, nrow = 4, ncol = 2, scales = 'free')
barplot(value ~ group, data = subset(div_gp, index == 'beta_C'))
dev.off()