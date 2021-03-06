#' Compute prediction intervals for state dependent Poisson process for beam experiments
#'
#' @param stresses list of values of the influental variable for indepedent fatigue experiments.
#' @param deltat list of waiting times for indepedent fatigue experiments
#' @param truss index of list entry in stresses and deltat for which the prediction should be made
#' @param start starting value for the optimization, length four
#' @param toPred integer value indicating the number of events to be predicted
#' @param link [\code{function(theta, x)}]\cr
#'        link function for exponential distribution
#' @param gradient [\code{function(x, theta, ...)}]\cr
#'        gradient of link function
#' @param type [\code{integer}]\cr
#'        if link function is not given a collection of given link function is available, see \code{\link{linkfun}}
#' @param plot logical value indicating whether the prediction intervals should be plotted or not
#' @param method one of "depth" (default), "chisquared". Method for generating confidence set of parameter theta
#' @export
piBeam <- function(stresses, deltat, truss, start, toPred, link, gradient, type, method = c("depth", "chisquared", "LR", "naive"), plot = FALSE, 
                   xlim, alpha = .05, addTrue = TRUE, ...) {
  if (missing(link)) {
    link <- linkfun(type)
  }
  
  if (missing(gradient)) {
    gradient <- gradLambda(type)
  }
  method <- match.arg(method)
  l <- length(stresses[[truss]])
  x0 <- stresses[[truss]][(l - (toPred - 1)):l]
  t0 <- deltat[[truss]][(l - (toPred - 1)):l]
  backup <- stresses
  stresses[[truss]] <- stresses[[truss]][-((l - (toPred - 1)):l)]
  deltat[[truss]] <- deltat[[truss]][-((l - (toPred - 1)):l)]
  
  x <- unlist(stresses)
  t <- unlist(deltat)
  
  estimation <- estML(x = x, t = t, start = start, link = link, gradient = gradient)
  theta <- estimation$optimum$par

  getQuantiles <- function(rates) {
    bLower <- qhypoexp2(p = alpha/2, rate = rates, interval = c(0, 10^10))
    bUpper <- qhypoexp2(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))
    return(c(bLower, bUpper))
  }
   
  if (method == "naive") {
    rate <- exp(link(x = x0, theta = theta))
    quantiles <- vapply(seq_along(rate), function(x) getQuantiles(rate[1:x]), FUN.VALUE = numeric(2))
    lowerBounds <- quantiles[1, ]
    upperBounds <- quantiles[2, ]
  }  else {
    confSet <- confidenceSet(theta = theta, x = x, t = t, alpha = alpha, method = method, lambda = link, gradient = gradient, ...)
    lambdas <- apply(confSet, 1, function(y) exp(link(x = x0, theta = y)))
    #lambdas <- apply(confSet, 1,  function(y) exp(-y[1] + y[2]*x0))  
    if (is.vector(lambdas))
      lambdas <- t(lambdas)
    
    quantiles <- apply(lambdas, 2, 
                       function(x) vapply(seq_along(x), function(y) getQuantiles(x[1:y]), FUN.VALUE = numeric(2)))
    lower <- seq(1, toPred * 2, by = 2)
    upper <- lower + 1
    lowerBounds <- vapply(lower, function(x) min(quantiles[x, ]), FUN.VALUE = numeric(1))
    upperBounds <- vapply(upper, function(x) max(quantiles[x, ]), FUN.VALUE = numeric(1))
  }
  tObserved <- sum(deltat[[truss]])
  if (length(tObserved) ==  0)
    tObserved <- rep(0, length(lowerBounds))
  unten <- tObserved + lowerBounds
  oben <- tObserved + upperBounds
  unten[is.na(oben)] <- NA
  # when lower boundary gets smaller set to NA, because this does not make sense
  unten[which(diff(unten) < 0)] <- NA
  oben[is.na(unten)] <- NA
  jumps <- c(0, seq_along(cumsum(c(deltat[[truss]], t0))))
  if (plot) {
    if (missing(xlim)) {
      xlim <- c(0, max(oben, na.rm = TRUE))/1000000 
    }
    
    plot(x = c(0, cumsum(deltat[[truss]]), unten)/1000000, y = jumps,
         lty = 2, type = "s", col = "red", xlab = "Cycle loads in millions", ylim = c(0, length(backup[[truss]])), 
         xlim = xlim, lwd = 2, ylab = "Number of broken tension wires")
    title(bquote((1 - .(alpha))^2 ~ plain("prediction intervals for" ~ .(names(backup[truss])))))
    points(x = c(0, cumsum(deltat[[truss]]), oben)/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
           lty = 2, type = "s", col = "red", lwd = 2)
    
    points(x = c(0, cumsum(deltat[[truss]])/1000000), y = c(0, seq_along(cumsum(deltat[[truss]]))),
           type = "s", lwd = 3)
    
    ## add true jumps:
    if (addTrue) {
    points(x = c(0, cumsum(c(deltat[[truss]], t0)))/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
           lty = 2, type = "s", lwd = 2)
    }
  }
  return(list(theta = theta, PI = rbind(unten, oben), jumps = jumps[(l-toPred + 2):(l+1)],
              lastJump = c(cumsum(deltat[[truss]])[l-toPred], jumps[l-toPred + 1])))
}