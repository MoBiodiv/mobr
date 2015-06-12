source('mobr.R') # Change to local path of mobr.R
## Two types of data structure:
## 1. Quadrat-based data
## For now, assume that this kind of data has the following columns (in this order):
## treatment, plot (plot ID, or pair ID), sample, species ID, abundance
## Samples are nested within plots. If there are no multiple samples (i.e., each plot is 
## surveyed only once), the sample column should be the same as the plot column.

## Import data, manipulate data frame into the above 5 columns
setwd('C:\\Users\\Xiao\\Downloads') # As an example
dat = read.csv('SpeciesCompositionData_ForJonChase.csv')
dat = dat[, c(3, 4, 2, 6, 9)]

##     Step 1. Initial comparison between treatments
dat_stats = get_sample_stats(dat)
##     If the plots are not paired, can do either parametric or non-parametric tests:
dat_p_par = initial_test_S_N(dat_stats, parametric = T, plot = T) # Parametric tests
dat_p_nonpar = inital_test_S_N(dat_stats, parametric = F, plot = T) # Non-parametric tests
##     If the plots are paired, can only do parametric, paired t-tests:
dat_p = initial_test_S_N_pair(dat_stats, plot = T)

##     Step 2. delta-S vs N (Step 4 in Cookbook) (shared between paired and unpaired data)
##     This step is not sensitive to whether plots are paired or not
dat_reform = reform_quad_data_to_sitesp(dat)
plot_deltaS_N(dat_reform)

##     Step 3. Kolmogorovo-Smirnov test on the shape of the SAD (Step 6 in Cookbook)
dat_for_ks = reform_quad_data_for_ks(dat)
##     There are four options: plots are paired or unpaired, and parametric or nonparametric tests are performed.
##     Note that the four options return vectors of different lengths, which may require further interpretation.
p_sad_ks = ks_test_between_treatments(dat_for_ks, T, T) # The second and third inputs can be T or F
plot_scaled_sads(dat_for_ks)

## 2. Individual-based data
## For now, assume that this kind of data has the following columns (in this order):
## treatment, plot (plot ID or pair ID), sample (if only a single sample within plot, this
## line should be the same as label), species ID, x, y

## Import data, manipulate data frame into the above 6 columns

##     Step 1. Initial comparison between treatments
dat_abd = reform_ind_data_to_abd(dat_in)
dat_stats = get_sample_stats(dat_abd)
## Then, do different tests depending on whether plots are paired or not
dat_p_par = initial_test_S_N(dat_stats, parametric = T, plot = T) # Parametric tests for unpaired data
dat_p_nonpar = inital_test_S_N(dat_stats, parametric = F, plot = T) # Non-parametric tests for unpaired data
dat_p = initial_test_S_N_pair(dat_stats, plot = T) # Pair-wise tests for paired data

##     Step 2. delta-S vs N (Step 4 in Cookbook)
dat_reform = reform_quad_data_to_sitesp(dat_abd)
plot_deltaS_N(dat_reform)

##     Step 3. Kolmogorovo-Smirnov test on the shape of the SAD (Step 6 in Cookbook)
##     See Step 3 for quadrat data for details.
dat_for_ks = reform_quad_data_for_ks(dat_abd)
p_sad_ks = ks_test_between_treatments(dat_for_ks, T, T) # The second and third inputs can be T or F
plot_scaled_sads(dat_for_ks)
