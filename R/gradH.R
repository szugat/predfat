## Derivation of H with respect to theta
#' Derivation of quantile from hypoexponential distribution with respect to theta
#' @param x values of influental variable for the link function
#' @param y values of depedent variable
#' @param theta numeric vector of length four with link function's parameters
#' @return derivation of quantile from hypoexponential distribution with respect to theta
gradH <- function(x, y, theta){
  rates <- lambda(x, theta)
  ratetoalpha <- function(rate) {   # function from sdprisk to compute a_j
    stopifnot(is.numeric(rate), all(is.finite(rate)))
    vapply(X = seq_along(rate),
           FUN = function(i) {
             prod(rate[-i] / (rate[-i] - rate[i]))
           },
           FUN.VALUE = numeric(1L))
  }
  aj <- ratetoalpha(rates)
  
  lambdaDot <- as.matrix(gradLambda(x, theta))
  
  ## second addend: a_j(theta) * diff(1 - exp(-b * lambda(j, s_0)), theta)
  bterm <- exp(-y * rates)
  gradBterm <- y * apply(lambdaDot, 2, "*", bterm)
  if(length(aj) > 1) {
    second <- aj %*% gradBterm
  } else {
    second <- aj * gradBterm
  }
  
  ## first addend: diff(a_j(theta), theta) * (1 - exp(-b * lambda(j, s_0)))
  tmp <- sapply(seq_along(rates),
                function(j) {
                  colSums(((rates[-j] %*% t(lambdaDot[j,])) - (lambdaDot[-j,] * rates[j]))/(rates[-j] * (rates[-j] - rates[j])))
                })
  first <- drop(tmp %*% (1 - bterm))
  return(drop(first + second))
}