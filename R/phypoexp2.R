phypoexp2 <- function(q, rate) {
  if (identical(sum(sdprisk:::ratetoalpha(rate)), 1)) {
    require(sdprisk)
    return(phypoexp(q, rate))
  } else {
    convolution2 <- function(q, rate) {
      require(Matrix)
      b <- length(rate)
      D <- matrix(0, ncol = b, nrow = b)
      D[col(D) - row(D) == 1] <- rate[1:(b-1)]
      A <- diag(-rate) + D
      A <- expm(A * q)
      drop(1 - sum(A[1, ]))
    }
    return(vapply(X = q, FUN = convolution2, rate = rate, FUN.VALUE = numeric(1)))
  }
}