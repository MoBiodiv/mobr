#' @title Invasive plants dataset
#' @name inv_comm
#' @aliases  inv_comm_attr
#' @description 
#' Plant species counts in invaded and unvaded sites. 
#' 
#' @details
#' \code{inv_comm} is a site-by-species matrix with species counts.
#' 
#' \code{inv_plot_attr} is a data frame with corresponding site variables. 
#' The column \code{group} specifies whether a site is "invaded" or "uninvaded". This variable is considered a "treatment" in the mob framework.
#' The columns \code{x} and \code{y} contain the spatial coordinates of the sites.
#' 
#' The data were adapted from Powell et al (2013).
#' 
#' 
#' @references
#' Powell, K. I., Chase, J. M., & Knight, T. M. (2013). Invasive plants have scale-dependent effects on diversity by altering species-area relationships. science, 339(6117), 316-318.
#' 
#' @examples 
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' @keywords data
NULL

