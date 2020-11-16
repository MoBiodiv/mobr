library(mobr)
library(betaC)
devtools::load_all()

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

# data exploration -------------------------------------------------------------
# make site x species matrix 
# define relevant variables
core_vars <- c('Recyear', 'Repsite', 'group', 'Spcode', 'Total')
comm <- tidyr::pivot_wider(dat_sub[ , core_vars],
                           names_from = Spcode, 
                           values_from = Total,
                           values_fill = 0,
                           values_fn = sum)
comm

# analyze biodiversity

dat_mob_in <- make_mob_in(comm[ , -c(1:3)], comm[ , 1:3])
dat_mob_in

plot_rarefaction(dat_mob_in, group_var = 'group', ref_level = 'control', 
                 method = 'IBR')
plot_rarefaction(dat_mob_in, group_var = 'group', ref_level = 'control', 
                 method = 'IBR', pooled = FALSE, log = 'xy')
plot_rarefaction(dat_mob_in, group_var = 'group', ref_level = 'control', 
                 method = 'SBR')

# for newer function usage setup nested dat structure

dat_nest <- comm  %>% 
            nest(comm = starts_with("sp"),
                 time = Recyear, 
                 space = Repsite)

index = c('S', 'S_PIE', 'S_n')
effort = 5
  
dat_nest <- dat_nest %>%
            mutate(div = map(dat$comm, calc_comm_div, index = index, effort = effort))

# need to plot results of div now... 
dat_nest

