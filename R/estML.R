#' Maximum likelihood estimation of parameters in link function for the exponential distribution
#' 
#' @param x values of influental variable for the link function
#' @param t values of dependent variable
#' @param delta vector indicating of length length(t) indicating if the respective observation in t is completely observed (1) or censored (0)
#' @param link a function(theta, log) used in the likelihood which is \eqn{L(\theta) = \sum t * \lambda(x) - log(\lambda(x))}
#' @param gradient [\code{function(x, theta, ...)}]\cr
#'        gradient of link function
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @param start starting solution for optim()
#' @param ... further paraments passed to optim()
#' @return optimal solution for the parameter in the link function in optimum and estimated 
#'         information matrix in I
#' @export
estML <- function(x, t, delta, link, gradient, type = 1, start, ...){
  ## link function
  if (missing(link)) {
    link <- linkfun(type)
  }
  ## if delta is not given assume that no observations are censored:
  if (missing(delta)) {
    delta <- rep(1, length(t))
  }

  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  ## objective function
  loglike <- function(theta){
    sum(delta * (t * exp(link(theta = theta, x = x)) - link(theta = theta, x = x)) + # density function
       (1 - delta) * exp(link(theta = theta, x = x)) * t) # survival function
  }
  thetahat <- optim(par = start, fn = loglike, ...)
  information <- estI(x, theta = thetahat$par, lambda = link, gradient = gradient)
  return(list(optimum = thetahat, I = information))
}