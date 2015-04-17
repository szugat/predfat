#' Plots data for a S-N-curve
#' 
#' @param x stress values for the x-axis
#' @param t values for y-axis
#' @param ... further arguments passed to plot()
#' @return plot of x and log10(t)
plotData <- function(x, t, ...){
  plot(x, log10(t), ylim = c(0, 10), pch = 19, 
       ylab = expression(log[10](t)), xlab = expression(s %.% L[max]/(L[max]-j)), ...)
}