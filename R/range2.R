#' Computes the range of a vector
#' 
#' @param x vector for which the range should be computed
#' @param ... parameters passed to base function \code{range}
#' @return returns the range of a vector, but in contrast to base R it returns the difference of \code{range}

range2 <- function(x, ...) {
  rg <- range(x, ...)
  rg[2] - rg[1]
}