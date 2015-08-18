#' Computation of fitted values for new data
#' @param newdata values of influental variable
#' @param theta parameter vector for link function
#' @param lambda [\code{function(theta, x)}]\cr
#'        link function for exponential distribution
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @param plot logical value indicating whether the fitted values should be added to an existing plot
#' @param median logical value indicating whether the median (TRUE) oder the expectation should be fitted (FALSE, default)
#' @return vector of fitted values
getFit <- function(newdata, theta, lambda, type, plot = FALSE, median = FALSE){
  if (missing(lambda)) {
    lambda <- linkfun(type)
  }
  if (length(theta) == 2) {
    theta <- c(theta, 0, 1)
  }
  
  fitted <- lambda(x = newdata, theta = theta)/log(10)
  
  if (median) {
    fitted <- log10(log(2)) + fitted
  }
  if(plot)
    points(newdata, fitted, type = "l", lwd = 2)
  return(fitted)
}
