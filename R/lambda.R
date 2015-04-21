#' Link function for continous four parameter model
#' 
#' @param x value of influental variable
#' @param theta vector of length four, parameter of link function
#' 
#' @return value of the link function
lambda <- function(x, theta){
  if ((length(theta) == 2) || theta[3] == 0) {
    theta <- c(theta, 0, 1)
  }
  exp(-theta[1] + theta[2]*x - theta[3]*x^(-theta[4]))
}