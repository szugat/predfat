#' Maximum likelihood estimation of parameters in link function for the exponential distribution
#' 
#' @param x values of influental variable for the link function
#' @param t values of dependent variable
#' @param delta vector indicating of length length(t) indicating if the respective observation in t is completely observed (1) or censored (0)
#' @param link a function(theta, log) used in the likelihood which is \eqn{L(\theta) = \sum t * \lambda(x) - log(\lambda(x))}
#' @param start starting solution for optim()
#' @param ... further paraments passed to optim()
#' @return optimal solution for the parameter in the link function in optimum and estimated 
#'         information matrix in I
#' @export
estML <- function(x, t, delta, link, type = 1, start, ...){
  ## link function
  if (missing(link)) {
    link <- linkfun(type)
  }
  ## if delta is not given assume that no observations are censored:
  if (missing(delta)) {
    delta <- rep(1, length(t))
  }
  ## objective function
  loglike <- function(theta){
    sum(delta * (t * exp(link(theta = theta, x = x)) - link(theta = theta, x = x)) + # density function
       (1 - delta) * exp(link(theta = theta, x = x)) * t) # survival function
  }
  thetahat <- optim(par = start, fn = loglike, ...)
  information <- predfat:::estI(x, theta = thetahat$par)
  return(list(optimum = thetahat, I = information))
}