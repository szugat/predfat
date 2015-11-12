 ## Derivation of H with respect to theta
#' Derivation of quantile from hypoexponential distribution with respect to theta
#' @param x values of influental variable for the link function
#' @param y values of depedent variable
#' @param theta numeric vector of length four with link function's parameters
#' @param lambda [\code{function(theta, x)}]\cr
#'        link function for exponential distribution
#' @param gradient [\code{function(x, theta, ...)}]\cr
#'        gradient of link function
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @return derivation of quantile from hypoexponential distribution with respect to theta
gradH <- function(x, y, theta, lambda, gradient, type){
  if (missing(lambda)) {
    lambda <- linkfun(type)
  }
  
  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  rates <- exp(lambda(x = x, theta = theta))
  derivations <- sapply(x0, gradient, theta = theta, lambda = link)
  
  ai <- sdprisk:::ratetoalpha(lambdaNew)
  ## derivation ai:
  faktoren <- sapply(seq_along(ai), function(k) {
    prod(sapply(seq_along(ai), function(j) {
      if (j !=k ) {
        lambdaNew[k] / (lambdaNew[k] - lambdaNew[j])
      } else {
        1
      }
    }))
  })
  
  ai_dev <- vector(mode = "list", length(ai))
  for(i in seq_along(ai)) {
    ai_dev[[i]] <- sapply(seq_along(ai), function(k) {
      if(k != i) {
        faktoren[k] *(lambdaNew[k] * derivations[, i] - lambdaNew[i] * derivations[, k]) / (lambdaNew[k] - lambdaNew[i])**2
      } else {
        rep(0, length(theta))
      }
    })
  }
  
  deriv_ai <- sapply(ai_dev, rowSums)
  
  lambdaDot <- as.matrix(gradient(x, theta, lambda = lambda))
  
  ## second addend: a_j(theta) * diff(1 - exp(-b * lambda(j, s_0)), theta)
  bterm <- exp(-y * rates)
  gradBterm <- y * apply(lambdaDot, 2, "*", bterm)
  if(length(ai) > 1) {
    second <- ai %*% gradBterm
  } else {
    second <- ai * gradBterm
  }
  
#   ## first addend: diff(a_j(theta), theta) * (1 - exp(-b * lambda(j, s_0)))
#   tmp <- sapply(seq_along(rates),
#                 function(j) {
#                   colSums(((rates[-j] %*% t(lambdaDot[j,])) - (lambdaDot[-j,] * rates[j]))/(rates[-j] * (rates[-j] - rates[j])))
#                 })
#   first <- drop(tmp %*% (1 - bterm))
  return(drop(first + second))
}