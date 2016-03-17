qhypoexp2 <- function (p, rate, interval = c(0, 1e+10)) {
  stopifnot(all(is.finite(rate)))
  invalid.p <- (p < 0 | p > 1)
  if (any(invalid.p)) {
    warning(paste("Some elements of", sQuote("p"), "are outside of the unit interval."))
    is.na(p) <- invalid.p
  }
  H <- function(x, p) {
    phypoexp2(q = x, rate = rate) - p
  }
  auxfun <- function(p) {
    tryCatch(expr = uniroot(f = H, interval = interval, p = p, 
                            extendInt = "upX", tol = .Machine$double.eps^0.5)$root, 
             error = function(err) NA_real_)
  }
  vapply(X = p, FUN = auxfun, FUN.VALUE = numeric(1L))
}