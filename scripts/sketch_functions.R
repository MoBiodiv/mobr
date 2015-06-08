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

initial_tests_S_N = function(dat_sample, parametric, plot){
  ## This function takes the output from get_sample_stats() and perform tests on whether samples differ 
  ## by treatments. The test will be ANoVA if treatment is categorical, linear regression if it is 
  ## continuous.
  ## Inputs: 
  ## dat_sample: output from get_sample_stats(), with 6 columns and one sample in each row.
  ## parametric: If TRUE, p-value is read directly from the output of lm(). If FALSE, p-value is 
  ##     computed by shuffling the treatment labels 1,000 times. Note that when the sample size is
  ##     low there may not be 1,000 unique permutations of the treatment labels.
  ## plot: If TRUE, a 2*3 plot is created to illustrate the comparisons of the five variables.\
  ## Output:
  ## A list of the five p-values: p_N, p_obsS, p_PIE, p_rareS, p_MIHS
  p_vec = c()
  f_vec = c()
  for (i in 3:7){
    model = summary(lm(dat_sample[, i]~dat_sample[, 2]))
    f = model$fstat
    f_vec = c(f_vec, f[1])
    p_vec = c(p_vec, pf(f[1], f[2], f[3], lower = F))
  }
  if (parametric == F) { # Non-parametric test
    f_frame = as.data.frame(matrix(nrow = 1000, ncol = 5))
    for (j in 1:1000){
      dat_rand = dat_sample
      dat_rand$treatment = sample(dat_rand$treatment, length(dat_rand$treatment))
      for (i in 3:7){
        model = summary(lm(dat_rand[, i]~dat_rand[, 2]))
        f = model$fstat
        f_frame[j, i-2] = f[1]
      }
    }
    for (i in 1:5){
      p_vec = c(p_vec, length(f_frame[f_frame[, i]>f_vec[i], i]) / 1000)
    }
  }
  if (plot == T){
    # Currently it is assumed that treatment has distinct levels.
    par(mfrow = c(2, 3))
    col_names = c('N', 'Observed S', 'PIE', 'Rarefied S', 'MIH delta-S')
    for (i in 1:5){
      boxplot(dat_sample[, i + 2] ~ dat_sample[, 2], main = col_names[i],
              las = 2)
      mtext(paste('p=', p_vec[i], sep = ''), cex = 0.8)
    }
  }
  return(list(p_N = p_vec[1], p_obs = p_vec[2], p_PIE = p_vec[3], p_rareS = p_vec[4], 
              p_MIHS = p_vec[5]))
}

