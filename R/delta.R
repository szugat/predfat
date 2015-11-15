#' Computation of prediction intervals
#' 
#' @export
#compPI <- function(xNew, x0, L, L_max = 35L, theta, lambda, gradient, type, alpha, I, withSolve = TRUE) {
delta <- function(stresses, deltat, truss, start, toPred, link, gradient, type, plot = FALSE, 
                   xlim, alpha = .05, addTrue = TRUE, withSolve = TRUE, ...) {
  if (missing(link)) {
    link <- linkfun(type)
  }
  
  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  l <- length(stresses[[truss]])
  x0 <- stresses[[truss]][(l - (toPred - 1)):l]
  t0 <- deltat[[truss]][(l - (toPred - 1)):l]
  backup <- stresses
  stresses[[truss]] <- stresses[[truss]][-((l - (toPred - 1)):l)]
  deltat[[truss]] <- deltat[[truss]][-((l - (toPred - 1)):l)]
  
  x <- unlist(stresses)
  t <- unlist(deltat)
  
  estimation <- estML(x = x, t = t, start = start, link = link, gradient = gradient)
  theta <- estimation$optimum$par
  
  lambdaNew <- exp(link(x = x0, theta = theta))
  
  #fitted <- getFit(newdata = xNew, theta = theta, lambda = lambda)
  
  getInterval <- function(rates, newdata) {
    bLower <- qhypoexp(p = alpha/2, rate = rates, interval = c(0, 10^10))
    bUpper <- qhypoexp(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))
    bDot <- function(y){
       -1/dhypoexp(x = y, rate = rates) * predfat:::gradH(newdata, y, theta, lambda = link, gradient = gradient, type = type)
    }
    
    if (withSolve) {
      I1 <- solve(estimation$I)
    } else {
      library(MASS)
      I1 <- ginv(estimation$I)
    }
    
    bdl <- bDot(bLower)
    v1 <- qnorm(1 - alpha/2) * sqrt(t(bdl) %*% I1 %*% bdl)
    
    bdu <- bDot(bUpper)
    v2 <- qnorm(1 - alpha/2) * sqrt(t(bdu) %*% I1 %*% bdu)
    return(list(interval = c(bLower - v1, bUpper + v2), quantiles = c(bLower, bUpper), v = c(v1, v2)))
  }
  lapply(seq_along(lambdaNew), function(i) getInterval(lambdaNew[1:i], x0[1:i]))
}
