#' Grid generation around a given two-dimensional point needed for confidence sets
#' 
#' Candidate generation for the confidence sets around a given a center with specified range and step width
#' @param center center of the grid, a two-dimensional point
#' @param range two-dimensional vector specifiying the width of the region around the center
#' @param steps numeric value indicating the step width of the grid
#' @return grid around \code{center} using \code{expand.grid}
generateCandidates <- function(center, range = c(0.5, 1.5), steps = center * 0.05) {
  p <- length(center)
  seqs <- vector(mode = "list", p)
  for(i in seq_along(seqs)) {
    seqs[[i]] <- seq(center[i] * range[1], center[i] * range[2], by = steps[i])
  }
  do.call("expand.grid", seqs)
}