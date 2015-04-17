#' Computation of fitted values for new data
#' @param newdata values of influental variable
#' @param theta parameter vector for link function
#' @param plot logical value indicating whether the fitted values should be added to an existing plot
#' @param median logical value indicating whether the median (TRUE) oder the expectation should be fitted (FALSE, default)
#' @return vector of fitted values
getFit <- function(newdata, theta, plot = FALSE, median = FALSE){
  if (length(theta) == 2) {
    theta <- c(theta, 0, 1)
  }
  
  fitted <- (theta[1] - theta[2]*(newdata) + theta[3] * (newdata)^(-theta[4]))/log(10)
  if (median) {
    fitted <- log10(log(2)) + fitted
  }
  if(plot)
    points(newdata, fitted, type = "l", lwd = 2)
  return(fitted)
}
