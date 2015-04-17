#' Plot event times
#' 
#' @param jumptimes Event times to be plotted
#' @param log logical value indicating whether to plot log10(jumptimes) or not, default FALSE
#' @param ... further arguments passed to plot()
#' @return Plot of event times
plotJumps <- function(jumptimes, log = FALSE, ...) {
  if (log) {
    plot(x = c(0, log10(cumsum(jumptimes))), y = c(0, seq_along(jumptimes)), ... )
  } else {
    plot(x = c(0, cumsum(jumptimes)), y = c(0, seq_along(jumptimes)), ... )
  }
}