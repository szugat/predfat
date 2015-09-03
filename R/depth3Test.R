depth3Test <- function(residuals, theta, x, y, type = c("full", "1", "2"), alpha = 0.05) {
  if (missing(residuals)) {
    residuals <- predfat:::computeResiduals(theta = theta, x = x, y = y, type = 4)
    residuals <- residuals[order(x)] # sort for ascending stress ranges
  }
  depth <- depth3Par(res = residuals, type = type)
  K <- 3
  N <- length(residuals)
  if (type == "1") {
    Tn <- sqrt(floor(N/(K+1))) * (depth - (1/2)**K) / sqrt((1/2)**K * (1 - (1/2)**K))
  }
  
  if (type == "2") {
    Tn <- sqrt(N-K) * (depth - (1/2)**K) / (sqrt((1/2)**K * (3 - (1/2)**(K-1) * K - 3 * (1/2)**K)))
  }
  return(Tn < qnorm(alpha))
}