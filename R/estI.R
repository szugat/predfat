#' Estimation of information matrix
#' 
#' @param x values of influental variable for the link function
#' @param theta numeric vector of length four with link function's parameters
#' @param lambda [\code{function(theta, x)}]\cr
#'        link function for exponential distribution
#' @param gradient [\code{function(x, theta, ...)}]\cr
#'        gradient of link function
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @return estimated information matrix
estI <- function(x, theta, lambda, gradient, type){
  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  
  if (missing(link)) {
    lambda <- linkfun(type)
  }
  ## gradient of lambda
  helper <- function(y, theta){
    lambdaValue <- exp(lambda(y, theta))
    lambdaDot <- drop(gradient(y, theta))
    mat <- 1/lambdaValue^2 * (lambdaDot %*% t(lambdaDot))
  }
  ##1/length(x) * Reduce("+", lapply(x, helper, theta = theta))
  Reduce("+", lapply(x, helper, theta = theta))
}