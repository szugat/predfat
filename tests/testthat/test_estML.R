context("Maximum likelihood estimation")

load(system.file("data.RData", package = "predfat"))
x <- unlist(stress[-12])
t <- unlist(t_jm[-12])
startSol <- c(42, 0.013,  0.6)

test_that("estML", {
  estim <- estML(x = x, t = t, type = 5, start = startSol[-2])
  theta <- estim$optimum$par
  I <- estim$I
  expect_is(theta, "numeric") # estimated parameter should be a vector
  expect_is(I, "matrix")      # estimated information matrix should be a matrix
  expect_equal(dim(I), rep(length(theta), 2)) # dimension of I should correspond to theta
})