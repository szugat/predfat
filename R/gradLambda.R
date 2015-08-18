#' Gradient of link function with respect to parameter theta
#' @param x value of influental variable for the link function
#' @param theta numeric vector of length four with link function's parameters
#' @return gradient of link function at x with respect to theta 
gradLambda <- function(type = 4) {
  lambda <- linkfun(type)
  lambdaValue <- exp(lambda(x, theta))
  if (type == 1 || type == 2) {
    lambdaDot <- cbind(-lambdaValue,
                       x * lambdaValue,
                       -x^(-theta[4]) * lambdaValue,
                       theta[3]*x^(-theta[4]) * lambdaValue * log(x))
    if ((length(theta) == 2) || theta[3] == 0) {
      lambdaDot <- lambdaDot[,1:2]
    }
  }
  
  if (type == 4) {
    lambdaDot <- cbind(-lambdaValue,
                       lambdaValue * x,
                       lambdavalue * log(x))
  }
  return(lambdaDot)
}