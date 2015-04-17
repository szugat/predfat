#' Maximum likelihood estimation of parameters in link function for the exponential distribution
#' 
#' @param x values of influental variable for the link function
#' @param t values of dependent variable
#' @param start starting solution for optim()
#' @param ... further paraments passed to optim()
#' @return optimal solution for the parameter in the link function in optimum and estimated 
#'         information matrix in I
estML <- function(x, t, start, ...){
  ## link function
  link <- function(theta, log = FALSE){
    if (length(start) == 4) {
      h <- -theta[1] + theta[2] * x - theta[3] * x^(-theta[4])
    } else if (length(start) == 2) {
      h <- -theta[1] + theta[2] * x
    }
    if(log){
      return(h)
    } else{
      return(exp(h))
    } 
  }
  
  ## objective function
  loglike <- function(theta){
    sum(t * link(theta = theta, log = FALSE) - link(theta = theta, log = TRUE))
  }
  thetahat <- optim(par = start, fn = loglike, ...)
  information <- estI(x, theta = thetahat$par)
  return(list(optimum = thetahat, I = information))
}