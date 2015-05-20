 #' Compute prediction intervals for state dependent Poisson process for beam experiments
#'
#' @param stresses list of values of the influental variable for indepedent fatigue experiments.
#' @param deltat list of waiting times for indepedent fatigue experiments
#' @param truss index of list entry in stresses and deltat for which the prediction should be made
#' @param start starting value for the optimization, length four
#' @param toPred integer value indicating the number of events to be predicted
#' @param plot logical value indicating whether the prediction intervals should be plotted or not
#' @param withSolve logical value indicating whether the information matrix should be inverted with solve() (TRUE) or ginv() from MASS package (FALSE)
#' @param plotControl list with arguments for plot: stairs (logical, plot as step function), xlim, legend (logical, generate legend)
#' @export
predBeam <- function(stresses, deltat, truss, start, toPred, plot = FALSE, 
                    plotControl = list(stairs = TRUE, xlim = c(0,0), legend = TRUE), withSolve = TRUE){ 
  # toPred = 1: predict next jump, toPred = 2: predict next two jumps
  # truss in {1, ..., 8}: prediction for which steel truss
  l <- length(stresses[[truss]])
  x0 <- stresses[[truss]][(l - (toPred - 1)):l]
  t0 <- deltat[[truss]][(l - (toPred - 1)):l]
  backup <- stresses
  stresses[[truss]] <- stresses[[truss]][-((l - (toPred - 1)):l)]
  deltat[[truss]] <- deltat[[truss]][-((l - (toPred - 1)):l)]
  
  x <- unlist(stresses)
  t <- unlist(deltat)
  
  estimation <- estML(x = x, t = t, start = start)
  test <- compPI(x0 = x0, L = length(x0) + toPred, L_max = 35L, theta = estimation$optimum$par, 
                 alpha = 0.05, I = estimation$I, withSolve = withSolve)
  intervals <- sapply(test, function(x) x$interval)
  quantiles <- vapply(test, function(x) x$quantiles, FUN.VALUE = numeric(2))
  jumps <- c(0, seq_along(cumsum(c(deltat[[truss]], t0))))
  PI <- intervals
  PI[PI < 0] <- max(PI[1,])
  lower <- cumsum(deltat[[truss]])[length(deltat[[truss]])]+ PI[1,]
  upper <- cumsum(deltat[[truss]])[length(deltat[[truss]])]+ PI[2,]
  if (plot) {
    if (plotControl$stairs == FALSE) {
      plotData(x, t, xlim = c(0, max(unlist(backup))))
      points(stresses[[truss]], log10(deltat[[truss]]), pch = c(rep(19, l-1), 8), col = "red")
      logint <- log10(intervals)
      logint[is.na(logint)] <- 1
      p <- ncol(intervals)
      for (i in 1:p) {
        points(backup[[truss]][l-(i-1)], logint[1,p - (i - 1)], pch = "_", col = "red")
        points(backup[[truss]][l-(i-1)], logint[2,p - (i - 1)], pch = "_", col = "red")
        segments(x0 = backup[[truss]][l-(i-1)], y0 = logint[1,p - (i - 1)], y1 = logint[2,p - (i-1)], lty = 2, col = "red")
      }
      points(x0, log10(t0), pch = 8, col = "red")
      title(bquote((1 - alpha)^2 ~ plain("prediction intervals for" ~ .(names(backup[truss])))))
      if (plotControl$legend) {
        legend("topright", lty = c(2, NA_integer_), lwd = c(1, NA_integer_),
               pch = c(NA_integer_, 8), legend = c("prediction interval", "true observation"),
               col = c("red", "red"))
      }
    } else {
      if (all(plotControl$xlim == c(0,0))) {
        plotControl$xlim <- c(0, (max(cumsum(deltat[[truss]])) + max(intervals))/1000000)
      }
      
      plot(x = c(0, cumsum(deltat[[truss]]), lower)/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
           lty = 2, type = "s", col = "red", xlab = "cycle loads in millions", ylim = c(0, length(backup[[truss]])), 
           xlim = plotControl$xlim,
           ylab = "Number of broken tension wires")
      title(bquote((1 - alpha)^2 ~ plain("prediction intervals for" ~ .(names(backup[truss])))))
      points(x = c(0, cumsum(deltat[[truss]]), upper)/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
             lty = 2, type = "s", col = "red", lwd = 2)
      
      points(x = c(0, cumsum(deltat[[truss]])/1000000), y = c(0, seq_along(cumsum(deltat[[truss]]))),
             type = "s", lwd = 2)
      
      ## add true jumps:
      points(x = c(0, cumsum(c(deltat[[truss]], t0)))/1000000, y = c(0, seq_along(cumsum(c(deltat[[truss]], t0)))),
             lty = 2, type = "s")
      
      if (plotControl$legend) {
        legend("bottomright", lty = c(1,2,2), lwd = c(1,1,1), col = c("black", "black", "red"),
               legend = c("known and used jumps", "true jumps", "prediction interval"))
      }
    }
  }
  return(list(theta = estimation$optimum$par, predictionIntervals = rbind(lower, upper), deltat = deltat[[truss]],
              t0 = t0, quantiles = quantiles, jumps = jumps[(l-toPred + 2):(l+1)], lastJump = c(cumsum(deltat[[truss]])[l-toPred], jumps[l-toPred + 1])))
}