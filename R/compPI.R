#' Computation of prediction intervals
#' 
#' @param xNew values of influential variable 
#' @param x0 observed values of influential variable
#' @param L integer value indicating the number of events to be predicted (used when xNew not given)
#' @param L_max integer value indicating the maximum number of events in Poisson process when xNew is not given
#' @param theta numeric vector of length four with link function's parameters
#' @param link link function for the negative (!) logarithmized mean of the exponential distribution
#' @param gradient gradient of link function w.r.t. parameters
#' @param type if link or gradient are not given, there is collection of implented link functions and gradients, see \code{\link{linkfun}} and \code{\link{gradLambda}}
#' @param alpha confidence niveau in (0,1)
#' @param I information matrix
#' @param withSolve logical value indicating if inverting of information matrix with solve() or ginv() from MASS
#' @return list with prediction interval bounds and corresponding quantiles of hypoexponential
#' @export
compPI <- function(xNew, x0, L, L_max = 35L, theta, lambda, gradient, type, alpha, I, withSolve = TRUE) {
  if(missing(xNew)) {
    stopifnot(L <= L_max)
    stopifnot(L > length(x0))
    xNew <- numeric(L - length(x0))
    for(i in (length(x0)+1):L){
      xNew[i-length(x0)] <- x0[1] * L_max/(L_max - i)
    }
  }
  
  if (missing(lambda)) {
    lambda <- linkfun(type)
  }
  
  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  if (length(theta) == 2) {
    theta <- c(theta, 0, 1)
  }
  lambdaNew <- exp(lambda(x = xNew, theta = theta))
  
  #fitted <- getFit(newdata = xNew, theta = theta, lambda = lambda)
  
  getInterval <- function(rates, newdata) {
    bLower <- qhypoexp(p = alpha/2, rate = rates, interval = c(0, 10^10))
    bUpper <- qhypoexp(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))
    bDot <- function(y){
      -1/dhypoexp(x = y, rate = rates) * gradH(newdata, y, theta, lambda = lambda, gradient = gradient)
    }
    
    if (withSolve) {
      I1 <- solve(I)
    } else {
      library(MASS)
      I1 <- ginv(I)
    }
    
    bdl <- bDot(bLower)
    v1 <- qnorm(1 - alpha/2) * sqrt(t(bdl) %*% I1 %*% bdl)
    
    bdu <- bDot(bUpper)
    v2 <- qnorm(1 - alpha/2) * sqrt(t(bdu) %*% I1 %*% bdu)
    return(list(interval = c(bLower - v1, bUpper + v2), quantiles = c(bLower, bUpper), v = c(v1, v2)))
  }
  lapply(seq_along(lambdaNew), function(i) getInterval(lambdaNew[1:i], xNew[1:i]))
}
