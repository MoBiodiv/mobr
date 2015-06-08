## Functions for emp data
## Here it's assumed that the input data is a data frame, with the following columns:
## ID: ID for each individual. If missing, can populate with seq()
## spcode: species code for each individual.
## sample: sample that the individual belongs to
## treatment: treatment of the sample
## x, y: coordinates of the individual for type 1 data, or coordinates
##     of the sample for type 2 data.

get_sample_stats = function(dat_in){
  ## Returns a data frame where each row is a sample, with the following columns:
  ## sample: sample ID
  ## treatment: treatment associated with the sample
  ## N: total number of individuals
  ## obsS: total number of observed species
  ## PIE: 1 - sum(p_i^2) 
  ## rareS: rarefied S at min(N). This value depends on other samples.
  ## MIHS: difference between obsS and rareS. Thus the min value among the rows will always be zero. 
  library(vegan)
  samples = unique(dat_in$sample)
  dat_out = as.data.frame(matrix(nrow = length(samples), ncol = 7))
  names(dat_out) = c('sample', 'treatment', 'N', 'obsS', 'PIE', 'rareS', 'MIHS')
  dat_out$sample = samples
  dat_out$N = as.numeric(table(dat_in$sample))
  minN = min(dat_out$N)
  for (i in 1:dim(dat_out)[1]){
    dat_in_sample = dat_in[dat_in$sample == samples[i], ]
    dat_out$treatment[i] = dat_in_sample$treatment[1]
    dat_out$N[i] = dim(dat_in_sample)[1]
    dat_out$obsS[i] = length(unique(dat_in_sample$spcode))
    sp_counts = as.numeric(table(dat_in_sample$spcode))
    dat_out$PIE[i] = dat_out$N[i] / (dat_out$N[i] - 1) * (1 - sum((sp_counts / sum(sp_counts))^2))
    dat_out$rareS[i] = rarefy(sp_counts, minN)
  }
  dat_out$MIHS = dat_out$obsS - dat_out$rareS
  return(dat_out)
}

initial_tests_S_N = function(dat_sample, )

