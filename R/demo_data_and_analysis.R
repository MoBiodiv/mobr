# the starting data object for analysis should allow for a number of 
# different data types. For example
# information on stems
stem_occ = data.frame(id = rep(1:nsites, length.out=nrow),
                      sp = sample(nrow, replace = T),
                      x = 1:nrow,
                      y = 1)
stem_vard

# information on abundances at sites
quad_abu = data.frame(id = rep(1:nsites, length.out=nrow),
                      sp = sample(nrow, replace = T),
                      abu = rpois(nrow, 5))
quad_vars = data.frame(id = 1:nsites, 
                       x = 1:nsites, 
                       y = 1,
                       var = rep(c('trt', 'cnt'), length.out=nsites))

#


get_sad(stems$sp)
# this 
# this will be formatted into 


# the mob input object will just be a matrix 
# if data on individuals is available then 
nrow= 100
nsites = 10
