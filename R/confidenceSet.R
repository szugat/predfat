#' Computation of a confidence set for the parameter theta in a two dimensional linear model based on grid search
#' 
#' @param candidates grid of candidates for confidence set
#' @param theta two-dimensional parameter of linear model
#' @param x values of influental variable for the link function 
#' @param y value of dependent variable
#' @param alpha value in (0,1) defining the level of the test
#' @param lambda [\code{function(theta, x)}]\cr
#'        link function for exponential distribution
#' @param gradient [\code{function(x, theta, ...)}]\cr
#'        gradient of link function
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @param ... further arguments passed to \code{generateCandidates}
#' @return confidence set for the two-dimensional parameter
confidenceSet <- function(candidates, theta, method = c("depth", "chisquared", "LR"), x, t = NULL, alpha = .05, lambda, gradient, type, depthType, ...) {
  method <- match.arg(method)
  if (missing(candidates)) {
    candidates <- generateCandidates(theta, ...)
  }
  
  if (missing(lambda)) {
    lambda <- linkfun(type)
  }
  
  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  
  if (method == "depth") {
    require(rexpar)
    residuals <- apply(candidates, 1, computeResiduals, x = x[order(x)], y = t[order(x)],  lambda = lambda)
    if (length(theta) == 2) {
      testCandidates <- apply(residuals, 2, function(x) dS_lin2_test(dS = dS_lin2(resy = x), alpha = alpha, y = t)$phi)
    }
    
    if (length(theta) == 3) {
      testCandidates <- apply(residuals, 2, depth3Test, type = depthType, alpha = alpha)
    }
  }
  
  if (method == "chisquared") {
    I <- predfat:::estI(x = x, theta = theta, lambda = lambda, gradient = gradient)
    diffs <- apply(candidates, 1, function(x) theta - x)
    testCandidates <- apply(diffs, 2, function(x) t(x) %*% I %*% x > qchisq(1 - alpha, df = length(theta)))
  }
  
  if (method == "LR") {
    ML_like <- sum(t * exp(lambda(theta = theta, x = x)) - lambda(theta = theta, x = x))
    cand_like <- apply(candidates, 1, function(y) sum(t * exp(lambda(theta = y, x = x)) - lambda(theta = y, x = x)))
    teststat <- 2 * (cand_like - ML_like)
    testCandidates <- teststat > qchisq(1 - alpha, df = length(theta))
  }
  candidates[!testCandidates, ]
}