library(plotly)


NV <- c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024)
timeTV <- c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154)
N0 <- 2
fExpression <- expression((k * N0) / (N0 + (k - N0) * exp(- r * timeT)))
zExpression <- expression(N - (k * N0) / (N0 + (k - N0) * exp(- r * timeT)))
sqeExpression <- expression(((k * N0) / (N0 + (k - N0) * exp(- r * timeT)) - N) ** 2)
partialK <- as.expression(D(sqeExpression, "k"))
partialFK <- as.expression(D(fExpression, "k"))
partialSK <- as.expression(D(partialK, "k"))
partialR <- as.expression(D(sqeExpression, "r"))
partialFR <- as.expression(D(fExpression, "r"))
partialSR <- as.expression(D(partialR, "r"))

SquareErrorQ3 <- function(theta, N0, NV, timeTV)
{ # function SquareErrorQ3
  result <- c()
  k <- theta[1]
  r <- theta[2]
  for (i in 1:length(timeTV))
  { # for i in 1:length(timeTV)
    timeT <- timeTV[i]
    N <- NV[i]
    result[i] <- eval(sqeExpression)
  } # for i in 1:length(timeTV)
  return(sum(result))
} # end function SquareErrorQ3

EvalPK <- function(theta, N0, NV, timeTV)
{ # function EvalPK
  result <- c()
  k <- theta[1]
  r <- theta[2]
  for (i in 1:length(timeTV))
  { # for i in 1:length(timeTV)
    timeT <- timeTV[i]
    N <- NV[i]
    result[i] <- eval(partialK)
  } # end for i in 1:length(timeTV)
  return(sum(result))
} # end function EvalPK

GetA <- function(theta, N0, NV, timeTV)
{ # function getA
  n <- length(NV)
  nt <- length(theta)
  k <- theta[1]
  r <- theta[2]
  aM <- matrix(NA, n, nt)
  for (i in 1:n)
  { # for i in 1:n
    timeT <- timeTV[i]
    N <- NV[i]
    aM[i, 1] <- eval(partialFK)
    aM[i, 2] <- eval(partialFR)
  } # end for i in 1:n
  return(aM)
} # end function getA

GetZ <- function(theta, N0, NV, timeTV)
{ # function getZ
  n <- length(NV)
  nt <- length(theta)
  k <- theta[1]
  r <- theta[2]
  zM <- matrix(NA, n, 1)
  for (i in 1:n)
  { # for i in 1:n
    timeT <- timeTV[i]
    N <- NV[i]
    zM[i] <- eval(zExpression)
  } # end for i in 1:n
  return(zM)
} # end function getZ

EvalSPK <- function(theta, N0, NV, timeTV)
{ # function EvalSPK
  result <- c()
  k <- theta[1]
  r <- theta[2]
  for (i in 1:length(timeTV))
  { # for i in 1:length(timeTV)
    timeT <- timeTV[i]
    N <- NV[i]
    result[i] <- eval(partialSK)
  } # end for i in 1:length(timeTV)
  return(sum(result))
} # end function EvalSPK

EvalPR <- function(k, r, N0, NV, timeTV)
{ # function EvalPR
  result <- c()
  k <- k
  r <- r
  for (i in 1:length(timeTV))
  { # for i in 1:length(timeTV)
    timeT <- timeTV[i]
    N <- NV[i]
    result[i] <- eval(partialR)
  } # end for i in 1:length(timeTV)
  return(sum(result))
} # end function EvalPR

EvalSPR <- function(k, r, N0, NV, timeTV)
{ # function EvalSPR
  result <- c()
  k <- k
  r <- r
  for (i in 1:length(timeTV))
  { # for i in 1:length(timeTV)
    timeT <- timeTV[i]
    N <- NV[i]
    result[i] <- eval(partialSR)
  } # end for i in 1:length(timeTV)
  return(sum(result))
} # end function EvalSPR

dummyK <- seq(1000, 1500, length.out = 500)
dummyR <- seq(0, 0.5, length.out = 500)


GaussNewtonQ3 <- function(theta0, N0, NV, timeTV)
{ # function GaussNewtonQ3
  errorBefore <- Inf
  errorAfter <- 1000
  i <- 1
  while (abs(errorBefore - errorAfter) >= 0.0000001 & i <= 200) 
  { # while abs(errorBefore - errorAfter) >= 0.0000001 & i <= 200
    errorBefore <- errorAfter
    zM <- GetZ(theta0, N0, NV, timeTV)
    aM <- GetA(theta0, N0, NV, timeTV)
    update <- solve(t(aM) %*% aM) %*% t(aM) %*% zM
    theta0 <- theta0 + update
    errorAfter <- sum(abs(update))
    i <- i + 1
  } # end while abs(errorBefore - errorAfter) >= 0.0000001 & i <= 200
  print(i)
  return(theta0)
} # end function GaussNewtonQ3


dummyK <- as.matrix(dummyK)
dummyR <- t(as.matrix(dummyR, byrow = TRUE))
dummyS <- matrix(NA, 500, 500)
for (i in 1:length(dummyK))
{ # for i in 1:length(dummyK)
  for (j in 1:length(dummyR))
  { # for j in 1:length(dummyR)
    dummyS[i, j] <- SquareErrorQ3(c(dummyK[i, 1], dummyR[1, j]), N0, NV, timeTV)
  } # end for j in 1:length(dummyR)
} # end for i in 1:length(dummyK)


plot_ly(x=dummyK[, 1],y=dummyR[1, ],z=dummyS, type="surface")
plot_ly(x=dummyK[, 1],y=dummyR[1, ],z=dummyS, type="contour")

