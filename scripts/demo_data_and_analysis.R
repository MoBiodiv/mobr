#source("./R/mobr.R")

# Input data -------------------------------------------------------------------
# the mob input object will just be a matrix 
nrow= 100
nsites = 10
# the starting data object for analysis should allow for a number of 
# different data types. For example
# information on stems
stem_occ = data.frame(site = rep(1:nsites, length.out=nrow),
                      sp = sample(nrow, replace = T),
                      x = 1:nrow,
                      y = 1)
class(stem_occ) = c('data.frame', 'indiv_occ')
stem_var = data.frame(site = 1:nsites, 
                      var = rep(c('trt', 'cnt'), length.out=nsites))


# information on abundances at sites
quad_abu = tapply(stem_occ$sp, list(stem_occ$site, stem_occ$sp),
                  length)
out = data.frame(site=NA, sp=NA, abu=NA)
for(i in 1:nsites) {
    for(j in 1:ncol(quad_abu)) {
        if (!is.na(quad_abu[i, j])) {
            tmp = data.frame(site=i, sp=j, abu=quad_abu[i,j])
            out = rbind(out, tmp)
        }
    }
}
quad_abu = out[-1, ]
class(quad_abu) = c('data.frame', 'site_occ')
quad_var = data.frame(site = 1:nsites, 
                      x = 1:nsites, 
                      y = 1,
                      var = rep(c('trt', 'cnt'), length.out=nsites))

# alternatively this data is likely to be input as a sp x site matrix
quad_abu2 = tapply(quad_abu$abu, list(quad_abu$site, quad_abu$sp), sum)
quad_abu2 = ifelse(is.na(quad_abu2), 0, quad_abu2)
class(quad_abu2) = c('matrix', 'site_by_sp')

# we should write the functions so they can handle these three basic 
# types of indput data. 

# Analysis --------------------------------------------------------------------

# the first analysis to run is an ANOVA on the trt effect



site_summary(stem_occ)
site_summary(quad_abu)
site_summary(quad_abu2)


