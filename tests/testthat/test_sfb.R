lolcontext("Check if prediction for SFB data works")
load(system.file("data.RData", package = "predfat"))
# starter <- c(25, 0.005, 2)
x <- unlist(stress)
t <- unlist(t_jm)
# estML(x = x, t = t, type = 4, start = starter)

# ## find starting solution for model type 4:
# theta1 <- 1:50
# theta2 <- seq(0.001, 0.1, by = 0.001)
# theta3 <- seq(0.1, 5, by = 0.5)
# grid <- expand.grid(theta1, theta2, theta3)
# best <- Inf
# solution <- rep(NA, 3)
# for(i in 1:nrow(grid)) {
#   cat(paste("Candidate", i, "\n"))
#   estim <- estML(x = x, t = t, type = 4, start = grid[i, ])
#   if(estim$optimum$value < best) {
#     best <- estim$optimum$value
#     solution <- estim$optimum$par
#     start <- grid[i, ]
#   }
# }

# > start
# Var1  Var2 Var3
# 5642   42 0.013  0.6
# > best
# [1] 1697.063
# > solution
# Var1         Var2         Var3 
# 23.353373008  0.005126494  1.791889051 

startSol <- c(42, 0.013,  0.6)
lambda <- predfat:::linkfun(type = 4)
#plotData(x, t, xlim = c(0, 600))
newdata <- 0:600
estimation <- estML(x = x, t = t, type = 4, start = startSol)
fitted <- -lambda(theta = estimation$optimum$par, x = newdata)/log(10)
#points(newdata, fitted, type = "l", lwd = 2)
fittedMed <- (log(log(2)) - lambda(theta = estimation$optimum$par, x = newdata))/log(10)
points(newdata, fittedMed, type = "l", lwd = 2, col = "red")
fittedExp <- log(2) / exp(lambda(theta = estimation$optimum$par, x = x))
#plot(x, t, pch = 19)
#points(x, fittedExp, pch = 19, col = "red")
## sieht gut aus!

## vorhersage testen:
testPred <- piBeam(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 5, type = 4, method = "chisquared")
testPred
## chisquared sieht erstmal gescheit aus, geht das Plotten noch?
testPred <- piBeam(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 5, type = 4, method = "chisquared", plot = TRUE)
## jup, geht

## Residuen für Tiefe:
theta <- estimation$optimum$par
y <- t
table(sign(y - (log(2) / exp(lambda(x = x, theta = theta)))))
## 69 negativ, 51 positiv. Das ist gut

residuals <- y - log(2) / exp(lambda(x = x, theta = theta))
fullDepth <- depth3Par(res = residuals, type = "full")
## 0.1200747
dS1 <- depth3Par(res = residuals, type = "1")
## 0.03333333
dS2 <- depth3Par(res = residuals, type = "2")
## 0.06837607
## für dS1 und dS2 müssen die Residuen sinnvoll geordnet sein nach den Schwingbreiten
ordnung <- order(x)
dS1 <- depth3Par(res = residuals[ordnung], type = "1")
print(dS1)
## 0.03333333
dS2 <- depth3Par(res = residuals[ordnung], type = "2")
print(dS2)
## 0.06837607#
## Es fällt auf, dass die vereinfachten Tiefen deutlich kleiner sind als die volle Tiefe.

## Tiefetest testen:
depth3Test(residuals = residuals[ordnung], type = "1")
# FALSE
depth3Test(residuals = residuals[ordnung], type = "2")
# FALSE
## Lehnen beide nicht ab, gut

## Confidence Set testen:
alpha <- 0.05
gradient <- predfat:::gradLambda(4)
candidates <- predfat:::generateCandidates(theta)
residuals <- apply(candidates, 1, predfat:::computeResiduals, x = x[order(x)], y = t[order(x)],  lambda = lambda)
depthType <- "1"
testres <- residuals[, 1]
testCandidates <- apply(residuals, 2, function(y) depth3Test(residuals = y, type = "1", alpha = alpha))
testCandidates2 <- apply(residuals, 2, function(y) depth3Test(residuals = y, type = "2", alpha = alpha))

## gesamtprozedur mal ausprobieren:
testPred_depth1 <- piBeam(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 5, type = 4, method = "depth", depthType = "1")
testPred_depth2 <- piBeam(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 5, type = 4, method = "depth", depthType = "2")
par(mfrow = c(2, 1))
testPred_depth1 <- piBeam(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 5, type = 4, method = "depth", depthType = "1", plot = TRUE)
testPred_depth2 <- piBeam(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 5, type = 4, method = "depth", depthType = "2", plot = TRUE)
## Informationsmatrix testen:
## kein NA -> gut
test_that("information matrix", {
  expect_false(any(is.na(estimation$I)))
})

theta <- estimation$optimum$par
candidates <- predfat:::generateCandidates(theta)
ML_like <- sum(t * exp(lambda(theta = theta, x = x)) - lambda(theta = theta, x = x))
cand_like <- apply(candidates, 1, function(y) sum(t * exp(lambda(theta = y, x = x)) - lambda(theta = y, x = x)))
teststat <- 2 * (cand_like - ML_like)
alpha <- Re(predfat:::qs(a = 1, b = -2, c = +0.05)[2])
testCandidates <- teststat > qchisq(1 - alpha, df = length(theta))
confset <- candidates[!testCandidates, ]



## linear testen
lambda <- predfat:::linkfun(2)
gradient <- predfat:::gradLambda(2)
theta <- estML(x = x, t = t, type = 2, start = startSol[1:2])$optimum$par
candidates <- predfat:::generateCandidates(theta)
ML_like <- sum(t * exp(lambda(theta = theta, x = x)) - lambda(theta = theta, x = x))
cand_like <- apply(candidates, 1, function(y) sum(t * exp(lambda(theta = y, x = x)) - lambda(theta = y, x = x)))
teststat <- 2 * (cand_like - ML_like)
alpha <- Re(predfat:::qs(a = 1, b = -2, c = +0.05)[2])
testCandidates <- teststat > qchisq(1 - alpha, df = length(theta))
confset <- candidates[!testCandidates, ]

confSet_LR <- predfat:::confidenceSet(theta = theta, x = x, t = t, alpha = alpha, method = "LR", lambda = lambda, gradient = gradient)
confSet_chi<- predfat:::confidenceSet(theta = theta, x = x, t = t, alpha = alpha, method = "chisquared", lambda = lambda, gradient = gradient)

## Test Basquin
theta <- estML(x = x, t = t, type = 5, start = startSol[c(1,3)])$optimum$par
newdata <- 0:600
lambda <- predfat:::linkfun(5)
fitted <- -lambda(theta = theta, x = newdata)/log(10)
plotData(x, t, xlim = c(0, 600))
points(newdata, fitted, type = "l", lwd = 2)

## Test Delta:
testDelta <- delta(stresses = stress, deltat = t_jm, truss = 9, start = startSol, toPred = 20, type = 4, withSolve = TRUE, alpha = alpha)

load(system.file("data.RData", package = "predfat"))
stresses <- stress
deltat <- t_jm
truss <- 9
startSol <- c(42, 0.013,  0.6)
start <- startSol
toPred <- 1
type <- 4
withSolve <- TRUE
alpha <- Re(predfat:::qs(a = 1, b = -2, c = +0.05)[2])

link <- predfat:::linkfun(4)
gradient <- predfat:::gradLambda(4)

l <- length(stresses[[truss]])
x0 <- stresses[[truss]][(l - (toPred - 1)):l]
t0 <- deltat[[truss]][(l - (toPred - 1)):l]
backup <- stresses
stresses[[truss]] <- stresses[[truss]][-((l - (toPred - 1)):l)]
deltat[[truss]] <- deltat[[truss]][-((l - (toPred - 1)):l)]

x <- unlist(stresses)
t <- unlist(deltat)

estimation <- estML(x = x, t = t, start = startSol, type = 4)
theta <- estimation$optimum$par

lambdaNew <- exp(link(x = x0, theta = theta))

derivations <- sapply(x0, gradient, theta = theta, lambda = link)
ai <- sdprisk:::ratetoalpha(lambdaNew)

faktoren <- sapply(seq_along(ai), function(k) {
  prod(sapply(seq_along(ai), function(j) {
    if (j !=k ) {
      lambdaNew[k] / (lambdaNew[k] - lambdaNew[j])
    } else {
      1
    }
  }))
})

ai_dev <- vector(mode = "list", length(ai))
for(i in seq_along(ai)) {
  ai_dev[[i]] <- sapply(seq_along(ai), function(k) {
      if(k != i) {
        faktoren[k] *(lambdaNew[k] * derivations[, i] - lambdaNew[i] * derivations[, k]) / (lambdaNew[k] - lambdaNew[i])**2
      } else {
        rep(0, length(theta))
      }
    })
}

deriv_ai <- sapply(ai_dev, rowSums)

## Probe für Ableitung nach theta1 für a_1:

(- derivations[1, 2] * lambdaNew[1] + lambdaNew[2] * derivations[1, 1]) / (lambdaNew[2] - lambdaNew[1])**2 *
  lambdaNew[2] / (lambdaNew[2] - lambdaNew[1]) * lambdaNew[2] / (lambdaNew[2] - lambdaNew[3]) * lambdaNew[2] / (lambdaNew[2] - lambdaNew[4]) * lambdaNew[2] / (lambdaNew[2] - lambdaNew[5]) +
  (- derivations[1, 3] * lambdaNew[1] + lambdaNew[3] * derivations[1, 1]) / (lambdaNew[3] - lambdaNew[1])**2 *
  lambdaNew[3] / (lambdaNew[3] - lambdaNew[1]) * lambdaNew[3] / (lambdaNew[3] - lambdaNew[2]) * lambdaNew[3] / (lambdaNew[3] - lambdaNew[4]) * lambdaNew[3] / (lambdaNew[3] - lambdaNew[5]) + 
  (- derivations[1, 4] * lambdaNew[1] + lambdaNew[4] * derivations[1, 1]) / (lambdaNew[4] - lambdaNew[1])**2 *
  lambdaNew[4] / (lambdaNew[4] - lambdaNew[1]) * lambdaNew[4] / (lambdaNew[4] - lambdaNew[2]) * lambdaNew[4] / (lambdaNew[4] - lambdaNew[3]) * lambdaNew[4] / (lambdaNew[4] - lambdaNew[5]) + 
  (- derivations[1, 5] * lambdaNew[1] + lambdaNew[5] * derivations[1, 1]) / (lambdaNew[5] - lambdaNew[1])**2 *
  lambdaNew[5] / (lambdaNew[5] - lambdaNew[1]) * lambdaNew[5] / (lambdaNew[5] - lambdaNew[2]) * lambdaNew[5] / (lambdaNew[5] - lambdaNew[3]) * lambdaNew[5] / (lambdaNew[5] - lambdaNew[4])
  
## kommt 0 raus, gut

## Probe für Ableitung nach theta2 für a_1:

(- derivations[2, 2] * lambdaNew[1] + lambdaNew[2] * derivations[2, 1]) / (lambdaNew[2] - lambdaNew[1])**2 *
  lambdaNew[2] / (lambdaNew[2] - lambdaNew[1]) * lambdaNew[2] / (lambdaNew[2] - lambdaNew[3]) * lambdaNew[2] / (lambdaNew[2] - lambdaNew[4]) * lambdaNew[2] / (lambdaNew[2] - lambdaNew[5]) +
  (- derivations[2, 3] * lambdaNew[1] + lambdaNew[3] * derivations[2, 1]) / (lambdaNew[3] - lambdaNew[1])**2 *
  lambdaNew[3] / (lambdaNew[3] - lambdaNew[1]) * lambdaNew[3] / (lambdaNew[3] - lambdaNew[2]) * lambdaNew[3] / (lambdaNew[3] - lambdaNew[4]) * lambdaNew[3] / (lambdaNew[3] - lambdaNew[5]) + 
  (- derivations[2, 4] * lambdaNew[1] + lambdaNew[4] * derivations[2, 1]) / (lambdaNew[4] - lambdaNew[1])**2 *
  lambdaNew[4] / (lambdaNew[4] - lambdaNew[1]) * lambdaNew[4] / (lambdaNew[4] - lambdaNew[2]) * lambdaNew[4] / (lambdaNew[4] - lambdaNew[3]) * lambdaNew[4] / (lambdaNew[4] - lambdaNew[5]) + 
  (- derivations[2, 5] * lambdaNew[1] + lambdaNew[5] * derivations[2, 1]) / (lambdaNew[5] - lambdaNew[1])**2 *
  lambdaNew[5] / (lambdaNew[5] - lambdaNew[1]) * lambdaNew[5] / (lambdaNew[5] - lambdaNew[2]) * lambdaNew[5] / (lambdaNew[5] - lambdaNew[3]) * lambdaNew[5] / (lambdaNew[5] - lambdaNew[4])

rates <- lambdaNew
bLower <- qhypoexp(p = alpha/2, rate = rates, interval = c(0, 10^10))
bUpper <- qhypoexp(p = 1 - alpha/2, rate = rates, interval = c(0, 10^10))

bterm <- exp(-bLower * rates)
gradBterm <- t(bLower * apply(derivations, 1, "*", bterm))

## erster Summand:
term1 <- rowSums(deriv_ai * rbind(1 - bterm, 1 - bterm, 1 - bterm))
## zweiterSummand:
term2 <- drop(gradBterm %*% ai)