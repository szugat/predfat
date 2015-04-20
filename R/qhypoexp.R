#' Implicit computation of the p-quantile from the hypoexponential distribution
#' 
#' @param p probability for the quantile in [0,1]
#' @param rates vector of rates of hypoexponential distribution
#' @param interval search interval in which uniroot() searches for the quantile
#' @return p-quantile from the hypoexponential distribution
qhypoexp <- function(p, rates, interval){
  require(sdprisk)
  H <- function(b){
    phypoexp(q = b, rate = rates) - p
  }
  uniroot(H, interval = interval)
}