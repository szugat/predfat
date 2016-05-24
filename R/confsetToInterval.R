#' Computes the quantiles of the hypoexponential distribution for all parameter values in a confidence set
#' 
#' @param cs Matrix of parameter values for the link function of the hypoexponential distribution.
#'            Each row of the matrix is one set of parameters
#' @param lambda function(x, theta) link function to link the value of x to the parameter of the hypoexponential distribution
#' @param test.s Vector of x values for the link function
#' @param alpha The alpha/2 and the 1 - alpha/2 quantiles of the hypoexponential distribution are computed
#' @return A list of all intervals for the parameters in the confidence set

confsetToInterval <- function(cs, lambda, test.s, alpha) {
  lambdas <- apply(cs, 1, function(y) exp(lambda(x = test.s, theta = y)))
  if (!is.matrix(lambdas)) {
    lambdas <- t(lambdas)
  }
  quantiles <- apply(lambdas, 2, function(x) getQuantiles(x, alpha))
  lowerBound <- min(quantiles[1, ])
  upperBound <- max(quantiles[2, ])
  return(list(interval = c(lowerBound, upperBound)))
}