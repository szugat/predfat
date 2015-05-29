#' Residuals of fitted median of exponential distribution used for data depth 
#' 
#' @param theta numeric vector of length two or four with link function's parameters
#' @param x values of influental variable for the link function
#' @param y values of dependent variable
#' @return vector of residuals
computeResiduals <- function(theta, x, y) {
  y - log(2) * 1 / lambda(x = x, theta = theta)
}