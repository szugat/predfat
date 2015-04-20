#' Estimation of information matrix
#' 
#' @param x values of influental variable for the link function
#' @param theta numeric vector of length four with link function's parameters
#' @return estimated information matrix
estI <- function(x, theta){
  ## gradient of lambda
  helper <- function(y, theta){
    lambdaValue <- lambda(y, theta)
    lambdaDot <- drop(gradLambda(y, theta))
    mat <- 1/lambdaValue^2 * (lambdaDot %*% t(lambdaDot))
  }
  ##1/length(x) * Reduce("+", lapply(x, helper, theta = theta))
  Reduce("+", lapply(x, helper, theta = theta))
}