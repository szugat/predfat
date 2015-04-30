library(predfat)
library(rexpar)
load("~/Documents/sfb/data/data.RData")

x <- unlist(stress)
t <- unlist(t_jm)

starter <- c(3, 0.005, 30, 2)

## estimate parameters for linear model:
estimation <- estML(x = x, t = t, start = starter[1:2])
theta <- estimation$optimum$par

residuals <- t - log(2) * 1 / predfat:::lambda(x = x, theta = theta)

dS <- dS_lin2(resy = residuals)
dS2 <- dS2_lin2(res = residuals)
dS3 <- dS3_lin2(res = residuals)

test_dS <- dS_lin2_test(dS = dS, alpha = 0.05, y = t)
test_dS2 <- dS2_lin2_test(dS2 = dS2, alpha = 0.05, y = t)
test_dS3 <- dS3_lin2_test(dS3 = dS3, alpha = 0.05, y = t)

## region around theta
generateCandidates <- function(center, range = c(0.5, 1.5), steps = center * 0.05) {
  x <- seq(center[1] * range[1], center[1] * range[2], by = steps[1])
  y <- seq(center[2] * range[1], center[2] * range[2], by = steps[2])
  expand.grid(x, y)
}

computeResiduals <- function(theta, x, y) {
  y - log(2) * 1 / predfat:::lambda(x = x, theta = theta)
}

candidates <- generateCandidates(theta)
residuals <- apply(candidates, 1, computeResiduals, x = x, y = t)

testFullDepth <- apply(residuals, 2, function(x) dS_lin2_test(dS = dS_lin2(resy = x), alpha = 0.05, y = t)$phi)

confidenceSet <- candidates[!testFullDepth, ]


## plot confidence set:
plot(theta[1], theta[2], col = "red")
points(confidenceSet)
points(theta[1], theta[2], col = "red")

## prediction intervall for given confidence set
rates <- sapply(30:600, function(x) apply(confidenceSet, 1, function(y) predfat:::lambda(x = x, theta = y)))
blower <- apply(rates, 2, function(x) sapply(x, function(y)
  predfat:::qhypoexp(p = 0.05/2, rate = y, interval = c(0, 10^10))))
bupper <- apply(rates, 2, function(x) sapply(x, function(y)
  predfat:::qhypoexp(p = 1 - 0.05/2, rate = y, interval = c(0, 10^10))))
lower <- apply(blower, 2, min)
upper <- apply(bupper, 2, max)

predAll <- sapply(30:600, compPI, theta = theta, I = estimation$I, alpha = as.numeric(predfat:::qs(1, -2, 0.05)[2]))
interVALLs <- vapply(predAll, function(x) x$interval, FUN.VALUE = numeric(2))

pdf("plots/PI_fulldepth.pdf", width = 9, height = 6)
plotData(x, t)
fitted <- predfat:::getFit(newdata = 30:600, theta = theta, plot = TRUE, median = TRUE)
points(30:600, log10(lower), lwd = 2, lty = 2, type = "l")
points(30:600, log10(upper), lwd = 2, lty = 2, type = "l")
points(30:600, log10(interVALLs[1,]), lwd = 2, lty = 2, type = "l", col = "red")
points(30:600, log10(interVALLs[2,]), lwd = 2, lty = 2, type = "l", col = "red")
legend("topright", lwd = c(2, 2, 2), lty = c(1, 2, 2), legend = c("Median", "95% PI based on full depth", "95% PI based on delta method"),
       col = c("black", "black", "red"))
dev.off()

## compare to other method for prediction intervals

pdf("plots/")
plotData(x, t)
fitted <- predfat:::getFit(newdata = 30:600, theta = theta, plot = TRUE, median = TRUE)
points(30:600, log10(lower), lwd = 2, lty = 2, type = "l")
points(30:600, log10(upper), lwd = 2, lty = 2, type = "l")

## predict more than one jump for one beam
confidenceSet


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

#' Residuals of fitted median of exponential distribution used for data depth 
#' 
#' @param theta numeric vector of length two or four with link function's parameters
#' @param x values of influental variable for the link function
#' @param y values of dependent variable
#' @return vector of residuals
computeResiduals <- function(theta, x, y) {
  y - log(2) * 1 / predfat:::lambda(x = x, theta = theta)
}

#' Computation of a confidence set for the parameter theta in a two dimensional linear model based on data depth
#' 
#' @param theta two-dimensional parameter of linear model
#' @param x values of influental variable for the link function 
#' @param y value of dependent variable
#' @param alpha value in (0,1) defining the level of the test
#' @param ... further arguments passed to \code{generateCandidates}
#' @return confidence set for the two-dimensional parameter
confidenceSet <- function(theta, x, t, alpha = .05, ...) {
  require(rexpar)
  candidates <- generateCandidates(theta, ...)
  residuals <- apply(candidates, 1, computeResiduals, x = x, y = t)
  testFullDepth <- apply(residuals, 2, function(x) dS_lin2_test(dS = dS_lin2(resy = x), alpha = alpha, y = t)$phi)
  candidates[!testFullDepth, ]
}

piDepth <- function(stresses, deltat, truss, start, toPred, plot = FALSE, xlim) {
  stopifnot(L <= L_max)
  stopifnot(L > length(x0))
  l <- length(stresses[[truss]])
  x0 <- stresses[[truss]][(l - (toPred - 1)):l]
  t0 <- deltat[[truss]][(l - (toPred - 1)):l]
  backup <- stresses
  stresses[[truss]] <- stresses[[truss]][-((l - (toPred - 1)):l)]
  deltat[[truss]] <- deltat[[truss]][-((l - (toPred - 1)):l)]
  
  x <- unlist(stresses)
  t <- unlist(deltat)
  
  estimation <- estML(x = x, t = t, start = start)
  theta <- estimation$optimum$par
  
  confSet <- confidenceSet(theta = theta, x = x, t = t, alpha = alpha)
  lambdas <- apply(confSet, 1,  function(theta) exp(-theta[1] + theta[2]*x0))  
  
  getQuantiles <- function(rates) {
    bLower <- sdprisk:::qhypoexp(p = alpha/2, rate = rates, interval = c(0, 10^10))
    bUpper <- sdprisk:::qhypoexp(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))
    return(c(bLower, bUpper))
  }
  
  quantiles <- apply(lambdas, 2, 
                     function(x) vapply(seq_along(x), function(y) getQuantiles(x[1:y]), FUN.VALUE = numeric(2)))
  lower <- seq(1, toPred * 2, by = 2)
  upper <- lower + 1
  lowerBounds <- vapply(lower, function(x) min(quantiles[x, ]), FUN.VALUE = numeric(1))
  upperBounds <- vapply(upper, function(x) max(quantiles[x, ]), FUN.VALUE = numeric(1))
  
  if (plot) {
    if (missing(xlim)) {
      xlim <- c(0, (max(cumsum(deltat[[truss]])) + max(upperBounds))/1000000)
    }
    unten <- cumsum(deltat[[truss]])[length(deltat[[truss]])]+ lowerBounds
    oben <- cumsum(deltat[[truss]])[length(deltat[[truss]])]+ upperBounds
    
    plot(x = c(0, cumsum(deltat[[truss]]), unten)/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
         lty = 2, type = "s", col = "red", xlab = "Cycle loads in millions", ylim = c(0, length(backup[[truss]])), 
         xlim = xlim, lwd = 2, ylab = "Number of broken tension wires")
    title(bquote((1 - .(alpha))^2 ~ plain("prediction intervals for" ~ .(names(backup[truss])))))
    points(x = c(0, cumsum(deltat[[truss]]), oben)/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
           lty = 2, type = "s", col = "red", lwd = 2)
    
    points(x = c(0, cumsum(deltat[[truss]])/1000000), y = c(0, seq_along(cumsum(deltat[[truss]]))),
           type = "s", lwd = 3)
    
    ## add true jumps:
    points(x = c(0, cumsum(c(deltat[[truss]], t0)))/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
           lty = 2, type = "s", lwd = 2)
  }
}

