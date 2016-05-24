#' Returns last column of a matrix
#' 
#' @param m matrix
#' @return last column of m
lastColumn <- function(m) {
  m[, ncol(m)]
}