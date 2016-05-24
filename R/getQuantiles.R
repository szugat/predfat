#' Computation of the alpha/2 and the 1 - alpha/2 quantiles of the hypoexponential distribution for given rates of exponential distributions
#' 
#' @param rates Vector of rates of exponential distributions
#' @param alpha The alpha/2 and the 1 - alpha/2 quantiles of the hypoexponential distribution are computed
#' @return alpha/2 and the 1 - alpha/2 quantiles of the hypoexponential distribution

getQuantiles <- function(rates, alpha) {
  bLower <- qhypoexp2(p = alpha/2, rate = rates, interval = c(0, 10^10))
  bUpper <- qhypoexp2(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))
  return(c(bLower, bUpper))
}