#' Collection of gradient functions of H function for delta-method
#' @param type integer (1-6) indicating the desired link function, see \code{\link{linkfun}}
#' @return gradient of H function at b with respect to theta 
gradBterm <- function(type) {
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
    res <- function(b, x, theta, lambda) {
      lambdaValue <- exp(lambda(x = x, theta = theta))
      fak <- b * exp(lambdaValue - b * exp(lambdaValue))
      cbind(- fak,
            x * fak,
            log(x) * fak)
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