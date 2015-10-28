#' Collection of gradient functions of link function for exponential distribution
#' @param type integer (1-6) indicating the desired link function, see \code{\link{linkfun}}
#' @return gradient of link function at x with respect to theta 
gradLambda <- function(type = 4) {
  if (type == 1 || type == 2) {
    res <- function(x, theta, lambda) {
      lambdaValue <- exp(lambda(x = x, theta = theta))
      lambdaDot <- cbind(-lambdaValue,
                         x * lambdaValue,
                         -x^(-theta[4]) * lambdaValue,
                         theta[3]*x^(-theta[4]) * lambdaValue * log(x))
      if ((length(theta) == 2) || theta[3] == 0) {
        lambdaDot <- lambdaDot[,1:2]
      }
      return(lambdaDot)
    }
  }
  
  if (type == 4) {
    res <- function(x, theta, lambda) {
      lambdaValue <- exp(lambda(x = x, theta = theta))
      cbind(-lambdaValue,
            lambdaValue * x,
            lambdaValue * log(x))
    }
  }
  
  if (type == 5) {
    res <- function(x, theta, lambda) {
      lambdaValue <- exp(lambda(x = x, theta = theta))
      cbind(-lambdaValue,
            lambdaValue * log(x))
    }
  }
  return(res)
}