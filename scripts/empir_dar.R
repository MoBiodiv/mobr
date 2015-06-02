## Purpose: to compute the empirical SARs

source('./scripts/div_functions.R')

dir.create('./results')

print('Computing empirical diveristy area curves...')

fileNames = dir('./data')
commFiles = grep('comms', fileNames)

## First aggregate community files into a single list
dat = vector('list', length(commFiles))
names(dat) = sub('_comms.csv', '', fileNames[commFiles])
for (i in seq_along(commFiles))
    dat[[i]] = as.matrix(read.csv(paste('./data/', fileNames[commFiles[i]], sep='')))

## filter the community files to the finest spatial grain
Amin = unlist(lapply(dat, function(x) unique(x[ , 1])[1]))
for (i in seq_along(dat))
    dat[[i]] = (dat[[i]][dat[[i]][ , 1] == Amin[i], ])

## load empirical constants
shrtnames = read.table('./data/shrtnames.txt', colClasses='character')
grain_fine = read.table('./data/grain_fine.txt')
bisect_fine = read.table('./data/bisect_fine.txt')
bisect_coarse = read.table('./data/bisect_coarse.txt')

proper_order = match(names(dat), shrtnames)
grain_fine = grain_fine[proper_order]
bisect_fine = bisect_fine[proper_order]
bisect_coarse = bisect_coarse[proper_order]

## convert them to multidimensional arrays

Ns = n_pixels_long(bisect_fine)
Ms = n_pixels_wide(bisect_fine)
psp = vector('list',length(dat))
names(psp) = names(dat)
for(i in seq_along(dat))
    psp[[i]] = mat2psp(dat[[i]][ , -(1:3)], dat[[i]][ , 2:3],
                       Ns[i], Ms[i])

## compute SARs, non-movinging window 
grains = sapply(1:length(psp), function(x) 2^(0:bisect_fine[1, x]))

dar = vector('list',length(psp))
names(dar) = names(psp)
for (i in seq_along(psp)) {
    dar[[i]] = getDAR(psp[[i]], grains[[i]], avg = F)
    ## add area m2 column
    dar[[i]] = data.frame(dar[[i]], area = dar[[i]][ , 1] * grain_fine[1, i])
}

## export results
for (i in seq_along(dar)) 
    write.csv(dar[[i]], file=paste('./results/', names(dar)[i], '_empir_dar.csv', 
                                   sep=''), row.names=FALSE)

## create flat file of dar data
for (i in seq_along(dar)) {
    if (i == 1)
        dat = data.frame(site=names(dar)[i], dar[[i]])
    else
        dat = rbind(dat, data.frame(site=names(dar)[i], dar[[i]]))
}

write.csv(dat, file='./results/empir_dars.csv', row.names=FALSE)

print('Computing empirical diversity area curves, complete!')

## generate figures -----------------------------------
library(ggtern)
#library(LSD)
library(scatterplot3d)
library(rgl)

dat_sub = dat[ , c('richness', 'indiv', 'shannon', 'simpson', 'coverage')]

png('./figs/empir_div_pairs_plot.png')
dat_sca = scale(dat_sub)
pairs(as.matrix(dat_sca))
dev.off()

#heatpairs(as.matrix(dat_sub[, c('richness', 'indiv')]))

pdf('./figs/prelim_div_plots_unaveraged.pdf')
dat_noNA = scale(na.omit(dat_sub))
pca = rda(dat_noNA)
plot(pca, display='sp')


dat_tern = dat_noNA[ , -(2:3)]


plot = ggtern()
plot = plot + 
       geom_point(data=data.frame(dat_tern), aes(richness, coverage, simpson))
plot

dev.off()

## average by grain
dat_avg = aggregate(dat[ , -(1:2)], by=list(dat$grains, dat$site), mean, na.rm=T)
names(dat_avg) = c('grain', 'site', names(tst)[-(1:2)])
head(dat_avg)
?scatterplot3d
scatterplot3d(dat_avg$richness, dat_avg$simpson, dat_avg$coverage)


plot3d(dat_avg$richness, dat_avg$simpson, dat_avg$coverage)

## pull out BCI site
pdf('./figs/bci_pairs_plots.pdf')
bci = subset(dat, site == 'bci')
uni_grains = unique(bci$grains)
for(i in seq_along(uni_grains)) {
    pairs(bci[bci$grains == uni_grains[i], -c(1, 2, 8)])
}
dev.off()

plot(richness ~ area, data = bci)
par(add=FALSE)

