#' Check if entries of a vector are contained in an interval
#' 
#' @param x Vector with entries which should be checked
#' @param int Interval given by a vector of two length 2, the first entry is the lower bound, the second entry is the upper bound
#' @return Boolean vector indicating if the elements of x are contained in the interval

inInterval <- function(x, int) {
  (x >= int[1]) & (x <= int[2])
}
