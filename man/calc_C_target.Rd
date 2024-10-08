% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_C.R
\name{calc_C_target}
\alias{calc_C_target}
\title{Calculate the recommended target coverage value for the computation of beta_C}
\usage{
calc_C_target(x, factor = 2)
}
\arguments{
\item{x}{a site by species abundance matrix}

\item{factor}{numeric. A multiplier for how much larger than total community 
abundance to extrapolate to. Defaults to 2.}
}
\value{
numeric value
}
\description{
Returns the estimated gamma-scale coverage that corresponds to the largest
allowable sample size (i.e. the smallest observed sample size at the alpha
scale multiplied by an extrapolation factor). The default (factor = 2) allows
for extrapolation up to 2 times the observed sample size of the smallest
alpha sample. For factor= 1, only interpolation is applied. Factors larger
than 2 are not recommended.
}
\examples{
data(tank_comm)

# What is the largest possible C that I can use to calculate beta_C
calc_C_target(tank_comm)
}
