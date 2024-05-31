#' ---
#' title: "beta diversity demonstration"
#' author: Dan McGlinn
#' output: pdf_document
#' ---
#' 

#+ setup, include=FALSE, eval=FALSE
install.packages(devtools)
library(devtools)
install_github('mobiodiv/mobr', ref = 'dev')

#' load mobr and example data
library(mobr)
data(inv_comm)

#' Calculate whittaker's beta
calc_comm_div(inv_comm[1:2, ], 'S')

#' Calculate beta for ENS of PIE (beta S_PIE)
calc_comm_div(inv_comm[1:2, ], 'S_PIE')

#' Calculate beta for S given a specific coverage (beta C)
calc_comm_div(inv_comm[1:2, ], 'S_C')

#' Calculate beta for rarefied fiechness (S_n) for 20 individuals
calc_comm_div(inv_comm[1:2, ], 'S_n', effort = 20)

#' More than two sites can be used at a time
calc_comm_div(inv_comm[1:10, ], 'S')

#' It is also possible to just calculate beta diversity
#' but it is generally not recommended to examine beta 
#' without reference to alpha and gamma diversity
calc_beta_div(inv_comm[1:10, ] , 'S')