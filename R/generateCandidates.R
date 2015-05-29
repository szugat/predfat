#' Grid generation around a given two-dimensional point needed for confidence sets
#' 
#' Candidate generation for the confidence sets around a given a center with specified range and step width
#' @param center center of the grid, a two-dimensional point
#' @param range two-dimensional vector specifiying the width of the region around the center
#' @param steps numeric value indicating the step width of the grid
#' @return grid around \code{center} using \code{expand.grid}
generateCandidates <- function(center, range = c(0.5, 1.5), steps = center * 0.05) {
  x <- seq(center[1] * range[1], center[1] * range[2], by = steps[1])
  y <- seq(center[2] * range[1], center[2] * range[2], by = steps[2])
  expand.grid(x, y)
}