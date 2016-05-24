#' Prediction interval for the hypoexponential distribution based on delta method
#' 
#' @param rates Rates of the exponential distributions whose sum are hypoexponential distributed
#' @param newdata Values for the link function which belong to the vector rates
#' @param alpha The alpha/2 and 1 - alpha/2 quantiles of the hypoexponential distribution are computed
#' @param I information matrix
#' @param gradient function(x,y,theta,lambda,gradient) gradient of the link function
#' @param theta Vector with parameters of the link function
#' @return Prediction interval for jumps at newdata using the delta method and the hypoexponential distribution
predInt <- function(rates, newdata, alpha, I, link, gradient, theta) {
  bLower <- qhypoexp2(p = alpha/2, rate = rates, interval = c(0, 10^10))
  bUpper <- qhypoexp2(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))
  bDot <- function(y){
    -1/dhypoexp(x = y, rate = rates) * predfat:::gradH(x = newdata, y = y, theta = theta, lambda = link, gradient = gradient)
  }
  
  I1 <- solve(I)
  
  bdl <- bDot(bLower)
  v1 <- qnorm(1 - alpha/2) * sqrt(t(bdl) %*% I1 %*% bdl)
  
  bdu <- bDot(bUpper)
  v2 <- qnorm(1 - alpha/2) * sqrt(t(bdu) %*% I1 %*% bdu)
  return(list(interval = c(bLower - v1, bUpper + v2), quantiles = c(bLower, bUpper), v = c(v1, v2)))
}