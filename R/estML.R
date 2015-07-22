#' Maximum likelihood estimation of parameters in link function for the exponential distribution
#' 
#' @param link a function(theta, log) used in the likelihood which is \eqn{L(\theta) = \sum t * \lambda(x) - log(\lambda(x))}, where log = TRUE refers to term at the end
#' @param x values of influental variable for the link function
#' @param t values of dependent variable
#' @param start starting solution for optim()
#' @param ... further paraments passed to optim()
#' @return optimal solution for the parameter in the link function in optimum and estimated 
#'         information matrix in I
#' @export
estML <- function(link, type = 1, x, t, start, ...){
  ## link function
  if (missing(link)) {
    link <- linkfun(type)
  }
  
  ## objective function
  loglike <- function(theta){
    sum(t * exp(link(theta = theta, x = x)) - link(theta = theta, x = x))
  }
  thetahat <- optim(par = start, fn = loglike, ...)
  information <- predfat:::estI(x, theta = thetahat$par)
  return(list(optimum = thetahat, I = information))
}