#' @title Invasive plants dataset
#' @name inv_comm
#' @aliases  inv_plot_attr
#' @description Herbaceous plant species counts sites invaded and uninvaded by
#'   Lonicera maackii (Amur honeysuckle) which is an invasive shrub.
#'
#' @details \code{inv_comm} is a site-by-species matrix with individual counts.
#'
#'   \code{inv_plot_attr} is a data frame with corresponding site variables. The
#'   column \code{group} specifies whether a site is "invaded" or "uninvaded".
#'   This variable is considered a "treatment" in the mob framework. The columns
#'   \code{x} and \code{y} contain the spatial coordinates of the sites.
#'
#'   The data were adapted from Powell et al (2013).
#'
#'
#' @references Powell, K. I., Chase, J. M., & Knight, T. M. (2013). Invasive
#'   plants have scale-dependent effects on diversity by altering species-area
#'   relationships. Science, 339: 316-318.
#'
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr)
#' @keywords data invasion invaded
NULL

#' @title Fire data set
#' @name fire_comm
#' @aliases  fire_plot_attr
#' @description Woody plant species counts in burned and unburned forest sites
#'   in the Missouri Ozarks, USA.
#'
#' @details \code{fire_comm} is a site-by-species matrix with individual counts.
#'
#' \code{fire_plot_attr} is a data frame with corresponding site variables. The
#' column \code{group} specifies whether a site is "burned" or "unburned". This
#' variable is considered a "treatment" in the mob framework. The columns
#' \code{x} and \code{y} contain the spatial coordinates of the sites.
#'
#' The data were adapted from Myers et al (2015).
#'
#'
#' @references Myers, J. A., Chase, J. M., Crandall, R. M., & Jiménez, I.
#' (2015). Disturbance alters beta‐diversity but not the relative importance of
#' community assembly mechanisms. Journal of Ecology, 103: 1291-1299.
#' @examples
#' data(fire_comm)
#' data(fire_plot_attr)
#' fire_mob_in = make_mob_in(fire_comm, fire_plot_attr)
#' @keywords data fire burned
NULL

#' @title Cattle tank data set
#' @name tank_comm
#' @aliases  tank_plot_attr
#' @description Species counts of aquatic macro-invertebrates from experimental
#'   freshwater ponds ("cattle tanks") with two different nutrient treatments.
#'
#' @details \code{tank_comm} is a site-by-species matrix with individual counts.
#'
#' \code{tank_plot_attr} is a data frame with corresponding site variables. The
#' column \code{group} specifies whether a pond has received a "high" or "low"
#' nutrient treatment. The columns \code{x} and \code{y} contain the spatial
#' coordinates of the sites.
#'
#' The data were adapted from Chase (2010).
#'
#'
#' @references Chase, J. M. (2010). Stochastic community assembly causes higher
#'   biodiversity in more productive environments. Science. 328:1388-1391.
#' @examples
#' data(tank_comm)
#' data(tank_plot_attr)
#' tank_mob_in = make_mob_in(tank_comm, tank_plot_attr)
#' @keywords data cattle tanks ponds
NULL
