#' Calculate expected sample coverage C_hat
#'
#' Returns expected sample coverage of a sample `x` for a smaller than observed
#' sample size `m` (Chao & Jost, 2012). This code was copied from INEXT's internal
#' function \code{iNEXT::Chat.Ind} (Hsieh et al 2016).
#' 
#' @param x integer vector (species abundances)
#' @param m integer a number of individuals that is smaller than observed total
#' community abundance. 
#'
#' @return a numeric value that is the expected coverage. 
#' 
#' @references 
#' Chao, A., and L. Jost. 2012. Coverage-based rarefaction and extrapolation:
#'  standardizing samples by completeness rather than size. Ecology 93:2533–2547.
#'  
#' Anne Chao, Nicholas J. Gotelli, T. C. Hsieh, Elizabeth L. Sander, K. H. Ma,
#'  Robert K. Colwell, and Aaron M. Ellison 2014. Rarefaction and extrapolation
#'  with Hill numbers: a framework for sampling and estimation in species
#'  diversity studies.  Ecological Monographs 84:45-67.
#' 
#' T. C. Hsieh, K. H. Ma and Anne Chao. 2024. 
#'  iNEXT: iNterpolation and EXTrapolation for
#'  species diversity. R package version 3.0.1
#'  URL: http://chao.stat.nthu.edu.tw/wordpress/software-download/.
#' 
#' 
#' @export
#'
#' @examples
#' data(inv_comm)
#' # What is the expected coverage at a sample size of 50 at the gamma scale?
#' Chat(colSums(inv_comm), 50)
Chat <- function (x, m)
{
    x <- x[x > 0]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n -
                                                                    1) / n * f1 ^
                         2 / 2 / f2)
    A <- ifelse(f1 > 0, n * f0.hat / (n * f0.hat + f1), 1)
    Sub <- function(m) {
        if (m < n) {
            xx <- x[(n - x) >= m]
            out <- 1 - sum(xx / n * exp(
                lgamma(n - xx + 1) - lgamma(n -
                                                xx - m + 1) - lgamma(n) + lgamma(n - m)
            ))
        }
        if (m == n)
            out <- 1 - f1 / n * A
        if (m > n)
            out <- 1 - f1 / n * A ^ (m - n + 1)
        out
    }
    sapply(m, Sub)
}

#' Number of individuals corresponding to a desired coverage (inverse C_hat)
#'
#' If you wanted to resample a vector to a certain expected sample coverage, how
#' many individuals would you have to draw? This is C_hat solved for the number
#' of individuals. This code is a modification of INEXT's internal function
#' `invChat.Ind` (Hsieh et al 2016).
#' 
#' @param x integer vector (species abundances)
#' @param C coverage value between 0 and 1
#'
#' @return a numeric value which is the number of individuals for a given
#' level of coverage \code{C}.
#' @references 
#' Chao, A., and L. Jost. 2012. Coverage-based rarefaction and extrapolation:
#'  standardizing samples by completeness rather than size. Ecology 93:2533–2547.
#'  
#' Anne Chao, Nicholas J. Gotelli, T. C. Hsieh, Elizabeth L. Sander, K. H. Ma,
#'  Robert K. Colwell, and Aaron M. Ellison 2014. Rarefaction and extrapolation
#'  with Hill numbers: a framework for sampling and estimation in species
#'  diversity studies.  Ecological Monographs 84:45-67.
#' 
#' T. C. Hsieh, K. H. Ma and Anne Chao. 2024. 
#'  iNEXT: iNterpolation and EXTrapolation for
#'  species diversity. R package version 3.0.1
#'  URL: http://chao.stat.nthu.edu.tw/wordpress/software-download/.
#' @seealso \code{\link{calc_S_C}}
#' @export
#' @importFrom stats optimize
#' @examples
#' data(inv_comm)
#' # What sample size corresponds to an expected sample coverage of 55%?
#' invChat(colSums(inv_comm), 0.55)
#'
invChat <- function (x, C)
{
    m <- NULL
    n <- sum(x)
    refC <- Chat(x, n)
    f <- function(m, C)
        abs(Chat(x, m) - C)
    # for interpolation
    if (refC > C) {
        opt <- stats::optimize(f,
                        C = C,
                        lower = 0,
                        upper = sum(x))
        mm <- opt$minimum
    }
    # for extrapolation
    if (refC <= C) {
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        if (f1 > 0 & f2 > 0) {
            A <- (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
        }
        if (f1 > 1 & f2 == 0) {
            A <- (n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2)
        }
        if (f1 == 1 & f2 == 0) {
            A <- 1
        }
        if (f1 == 0 & f2 == 0) {
            A <- 1
        }
        mm <- (log(n / f1) + log(1 - C)) / log(A) - 1
        mm <- n + mm
        
    }
    if (mm > 2 * n)
        warning(
            "The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias."
        )
    return(mm)
}


#' Calculate species richness for a given coverage level. 
#'
#' This function uses coverage-based rarefaction to compute species richness.
#' Specifically, the metric is computed as the
#' 
#' @param x a site by species matrix or a species abundance distribution
#' @param C_target target coverage between 0 and 1 (default is NULL). If not
#' provided then target coverage is computed by \code{\link{calc_C_target}}
#' @param extrapolate logical. Defaults to TRUE in which case richness is 
#' extrapolated to sample sizes larger than observed in the dataset.
#' @param interrupt logical. Should the function throw an error when \code{C_target}
#'  exceeds the maximum recommendable coverage?
#'
#' @returns numeric value which is the species richness at a specific level of 
#' coverage.
#' @references 
#' Chao, A., and L. Jost. 2012. Coverage-based rarefaction and extrapolation:
#'  standardizing samples by completeness rather than size. Ecology 93:2533–2547.
#'  
#' Anne Chao, Nicholas J. Gotelli, T. C. Hsieh, Elizabeth L. Sander, K. H. Ma,
#'  Robert K. Colwell, and Aaron M. Ellison 2014. Rarefaction and extrapolation
#'  with Hill numbers: a framework for sampling and estimation in species
#'  diversity studies.  Ecological Monographs 84:45-67.
#' 
#' T. C. Hsieh, K. H. Ma and Anne Chao. 2024. 
#'  iNEXT: iNterpolation and EXTrapolation for
#'  species diversity. R package version 3.0.1
#'  URL: http://chao.stat.nthu.edu.tw/wordpress/software-download/.
#' 
#' @seealso \code{\link{invChat}}
#' @export
#'
#' @examples
#' data(tank_comm)
#' # What is species richness for a coverage value of 60%?
#' calc_S_C(tank_comm, C_target = 0.6)
calc_S_C <- function(x,
                   C_target = NULL,
                   extrapolate = TRUE,
                   interrupt = TRUE) {
    x <- as.matrix(x)
    if (any(dim(x) == 1))
      sad <- as.numeric(x)
    else
      sad <- colSums(x)
    if (is.null(C_target))
      C_target <- calc_C_target(x)
    N <- round(invChat(sad, C_target))
    C_max = calc_C_target(x, factor = ifelse(extrapolate, 2, 1))
    if (C_target > C_max & interrupt) {
        if (extrapolate) {
          stop(
            paste0(
              "Coverage exceeds the maximum possible value recommendable for extrapolation (i.e. C_target = ",
              round(C_max, 4),
              "). Reduce the value of C_target."
            )
          )
        } else{
          stop(
            paste0(
              "Coverage exceeds the maximum possible value for interpolation (i.e. C_target = ",
              round(C_max, 4),
              "). Use extrapolation or reduce the value of C_target."
            )
          )
        }
    }
    if (N > 1) {
        S_C = rarefaction(x = sad,
                            method = "IBR",
                            effort = N,
                            extrapolate = extrapolate,
                            quiet_mode = TRUE)
    } else {
        S_C = NA
    }
    attr(S_C, "C") = C_target
    attr(S_C, "N") = N
    return(S_C)
}

#' Calculate the recommended target coverage value for the computation of beta_C 
#'
#' Returns the estimated gamma-scale coverage that corresponds to the largest
#' allowable sample size (i.e. the smallest observed sample size at the alpha
#' scale multiplied by an extrapolation factor). The default (factor = 2) allows
#' for extrapolation up to 2 times the observed sample size of the smallest
#' alpha sample. For factor= 1, only interpolation is applied. Factors larger
#' than 2 are not recommended.
#'
#' @param x a site by species abundance matrix
#' @param factor numeric. A multiplier for how much larger than total community 
#' abundance to extrapolate to. Defaults to 2. 
#'
#' @return numeric value
#' @export
#'
#' @examples
#' data(tank_comm)
#'
#' # What is the largest possible C that I can use to calculate beta_C
#' calc_C_target(tank_comm)
calc_C_target <- function(x, factor = 2) {
    x <- as.matrix(x)
    if (any(dim(x) == 1)) {
        n <- factor * sum(x)
        C_target <- Chat(x, n)
    }
    else {
        n <- min(factor * rowSums(x))
        C_target <- Chat(colSums(x), n)
    }
    return(C_target)
}