#' A collection of link functions for the estimation of S-N-Curves. 
#' 
#' @param type integer (1-6) indicating the desired link function
#' @return function(theta, x)
linkfun <- function(type) {
  if (type == 1) { # nonlinear
    res <- function(theta, x) {
      -theta[1] + theta[2] * x - theta[3] * x**(-theta[4])
    }
  }
  
  if (type == 2) { # linear
    res <- function(theta, x) {
      -theta[1] + theta[2] * x
    } 
  }
  
  if (type == 3) { # linear with 1/x
    res <- function(theta, x) {
      -theta[1] + theta[2] * x - theta[3] * 1/x
    }
  }
  
  if (type == 4) { # linear with log(x)
    res <- function(theta, x) {
      -theta[1] + theta[2] * x + theta[3] * log(x)
    }
  }
  
  if (type == 5) { # Basquin (1910)
    res <- function(theta, x) {
      -theta[1] + theta[2] * log(x)
    }
  }
  
  if (type == 6) { # Block und Dreier
    res <- function(theta, x) {
      - (log((x - theta[1])/(theta[2] - theta[1])) / log(theta[3]))**(1/theta[4])
    }
  }
  return(res)
}