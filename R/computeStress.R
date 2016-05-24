#' Computation of load sharing rule
#' 
#' @param initial Initial stress value when no components are damaged
#' @param nJumps Number of failures for which the load sharing should be computed
#' @param L_max Number of components in the system
#' @return Vector of initial stress and the stress values with respect to the load sharing rule
computeStress <- function(initial, nJumps, L_max = 35) {
  if (nJumps == 0) {
    return(initial)
  }
  return(c(initial, initial * L_max/(L_max - (1:nJumps))))
}