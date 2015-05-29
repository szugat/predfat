#' Computation of a confidence set for the parameter theta in a two dimensional linear model based on data depth
#' 
#' @param theta two-dimensional parameter of linear model
#' @param x values of influental variable for the link function 
#' @param y value of dependent variable
#' @param alpha value in (0,1) defining the level of the test
#' @param ... further arguments passed to \code{generateCandidates}
#' @return confidence set for the two-dimensional parameter
confidenceSet <- function(theta, method = c("depth", "chisquared"), x, t = NULL, alpha = .05, ...) {
  method <- match.arg(method)
  candidates <- generateCandidates(theta, ...)
  
  if (method == "depth") {
    require(rexpar)
    residuals <- apply(candidates, 1, computeResiduals, x = x, y = t)
    testCandidates <- apply(residuals, 2, function(x) dS_lin2_test(dS = dS_lin2(resy = x), alpha = alpha, y = t)$phi)
  }
  
  if (method == "chisquared") {
    I <- predfat:::estI(x = x, theta = theta)
    diffs <- apply(candidates, 1, function(x) theta - x)
    testCandidates <- apply(diffs, 2, function(x) t(x) %*% I %*% x > qchisq(1 - alpha, df = length(theta)))
  }
  
  candidates[!testCandidates, ]
}