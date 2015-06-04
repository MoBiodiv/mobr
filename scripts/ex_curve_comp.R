library(vegan)
source('./scripts/div_functions.R')

sherman = read.csv('./data/filtered_data/sherman_census3_filtered.csv')
bigoak = read.csv('./data/filtered_data/m94_bigoak_1993_filtered.csv')

sherman = subset(sherman, subset= y <= 150 & 
                                  x >= 140)

assign_id = function(x, nsites) {
    x_tmp = abs(x - max(x))
    x_breaks = seq(min(x_tmp), max(x_tmp), length.out=nsites + 1)
    ids = as.numeric(cut(x_tmp, x_breaks, right=T, include.lowest=T))
    return(ids)
}

sherman$id = assign_id(sherman$x, 4)
bigoak$id = assign_id(bigoak$x, 4)

sherman_sads = tapply(sherman$id, INDEX = list(sherman$id, sherman$spcode), 
              length)
sherman_sads = ifelse(is.na(sherman_sads), 0, sherman_sads)

bigoak_sads = tapply(bigoak$id, INDEX = list(bigoak$id, bigoak$SPEC), 
                      length)
bigoak_sads = ifelse(is.na(bigoak_sads), 0, bigoak_sads)

sherman_acc = apply(sherman_sads, 1, function(x) rarefy(x, 1:sum(x)))
bigoak_acc = apply(bigoak_sads, 1, function(x) rarefy(x, 1:sum(x)))

sherman_agg = vector('list', length(unique(sherman$id)))
for(i in sort(unique(sherman$id))) {
    sherman_agg[[i]] = near_neigh_ind(sherman[sherman$id == i,
                                              c('x', 'y', 'spcode')])
}

bigoak_agg = vector('list', length(unique(bigoak$id)))
for(i in sort(unique(bigoak$id))) {
    bigoak_agg[[i]] = near_neigh_ind(bigoak[bigoak$id == i,
                                            c('X', 'Y', 'SPEC')])
}

sherman_minN = min(unlist(lapply(sherman_acc, length)))
sherman_avg = rep(0, sherman_minN)
for(i in 1:length(sherman_acc)) {
    sherman_avg = sherman_avg + sherman_acc[[i]][1:sherman_minN]
}
sherman_avg = sherman_avg / length(sherman_acc)

bigoak_minN = min(unlist(lapply(bigoak_acc, length)))
bigoak_avg = rep(0, bigoak_minN)
for(i in 1:length(bigoak_acc)) {
    bigoak_avg = bigoak_avg + bigoak_acc[[i]][1:bigoak_minN]
}
bigoak_avg = bigoak_avg / length(bigoak_acc)

#sherman_sub = sherman_avg[1:length(bigoak_avg)]
bigoak_sub = bigoak_avg[1:length(sherman_avg)]

get_dS = function(acc_curv) {
    dS = rep(0, length(acc_curv) - 1)
    for(i in seq_along(acc_curv)) {
        dS[i] = acc_curv[i + 1] - acc_curv[i]
    }
    return(dS)
}

bigoak_ds = get_dS(bigoak_sub)
sherman_ds = get_dS(sherman_avg)



pdf('./figs/bigoak_sherman_div.pdf')
boxplot(cbind(table(sherman$id), table(bigoak$id)),
        names=c('sherman', 'bigoak'))
##
plot(sherman_acc[[1]][1,], type='l', lwd=2, col='blue')
lines(sherman_agg[[1]], col='dodgerblue', lwd=2)
lines(bigoak_acc[[1]][1,], lwd=2, col='red')
lines(bigoak_agg[[1]], lwd=2, col='pink')
##
#plot(sherman_avg, type='l', col='dodgerblue',
#     xlim=range(1,length(bigoak_avg), length(sherman_avg)))
#lines(1:length(bigoak_avg), bigoak_avg, col='red')
#abline(v=length(sherman_avg))
#abline(h=bigoak_avg[length(sherman_avg)])
##
plot(sherman_avg, type='l', col='dodgerblue',
     xlim=range(1,length(bigoak_avg), length(sherman_avg)))
lines(1:length(bigoak_avg), bigoak_avg, col='red')
abline(v=length(sherman_avg))
abline(h=bigoak_avg[length(sherman_avg)])
##
plot(sherman_avg, type='l', xlab='Individuals', ylab='richness', 
     col='dodgerblue', lwd=4, ylim=range(sherman_avg, bigoak_sub))
lines(1:length(bigoak_sub), bigoak_sub, col='red', lwd=4)
##
plot(sherman_avg - bigoak_sub, type='l', lwd=4,
     ylab='Delta rich (wet - dry)', 
     xlab='Individuals')
##
plot(sherman_ds - bigoak_ds, type='l', lwd=4, 
     ylab='Delta dS', xlab='Individuals')
dev.off()


