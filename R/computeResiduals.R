#' Residuals of fitted median of exponential distribution used for data depth 
#' 
#' @param theta numeric vector of length two or four with link function's parameters
#' @param x values of influental variable for the link function
#' @param y values of dependent variable
#' @param lambda [\code{function(theta, x)}]\cr
#'        link function for exponential distribution
#' @param gradient [\code{function(x, theta, ...)}]\cr
#'        gradient of link function
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @return vector of residuals
computeResiduals <- function(theta, x, y, lambda, type) {
  if (missing(lambda)) {
    lambda <- linkfun(type)
  }
  y - log(2) * 1 / exp(lambda(x = x, theta = theta))
}