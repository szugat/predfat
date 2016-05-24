#' Returns last element of a list
#' 
#' @param l List
#' @return last element of list l
lastList <- function(l) {
  l[[length(l)]]
}