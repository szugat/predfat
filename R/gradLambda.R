#' Gradient of link function with respect to parameter theta
#' @param x value of influental variable for the link function
#' @param theta numeric vector of length four with link function's parameters
#' @return gradient of link function at x with respect to theta 
gradLambda <- function(x, theta){
  lambdaValue <- lambda(x, theta)
  lambdaDot <- cbind(-lambdaValue,
                     x * lambdaValue,
                     -x^(-theta[4]) * lambdaValue,
                     theta[3]*x^(-theta[4]) * lambdaValue * log(x))
  if ((length(theta) == 2) || theta[3] == 0) {
    return(lambdaDot[,1:2])
  } else {
    return(lambdaDot)
  }
}