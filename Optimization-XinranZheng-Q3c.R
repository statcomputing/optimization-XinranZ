library(plotly)

NV <- c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024)
timeTV <- c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154)
N0 <- 2

f3cLLExpression <- expression(- log((2 * pi * sd) ** 0.5) - ((log(N) - log((k * N0) / (N0 + (k - N0) * exp(- r * timeT)))) ** 2) / (2 * (sd ** 2)))
q3cpartialSD <- as.expression(D(f3cLLExpression, "sd"))
q3cpartialK <- as.expression(D(f3cLLExpression, "k"))
q3cpartialR <- as.expression(D(f3cLLExpression, "r"))
q3cpartialSSD <- as.expression(D(q3cpartialSD, "sd"))
q3cpartialSK <- as.expression(D(q3cpartialK, "k"))
q3cpartialSR <- as.expression(D(q3cpartialR, "r"))

GetdSD <- function(theta, N0, NV, timeTV)
{ # function GetdSD
  result <- c()
  N0 <- N0
  k <- theta[1]
  r <- theta[2]
  sd <- theta[3]
  for (i in 1:length(NV))
  { # for
    N <- NV[i]
    timeT <- timeTV[i]
    result[i] <- eval(q3cpartialSD)
  } # end for
  return(sum(result))
} # end function GetdSD

GetdK <- function(theta, N0, NV, timeTV)
{ # function GetdK
  result <- c()
  N0 <- N0
  k <- theta[1]
  r <- theta[2]
  sd <- theta[3]
  for (i in 1:length(NV))
  { # for
    N <- NV[i]
    timeT <- timeTV[i]
    result[i] <- eval(q3cpartialK)
  } # end for
  return(sum(result))
} # end function GetdK

GetdR <- function(theta, N0, NV, timeTV)
{ # function GetdR
  result <- c()
  N0 <- N0
  k <- theta[1]
  r <- theta[2]
  sd <- theta[3]
  for (i in 1:length(NV))
  { # for
    N <- NV[i]
    timeT <- timeTV[i]
    result[i] <- eval(q3cpartialR)
  } # end for
  return(sum(result))
} # end function GetdR

GetsdSD <- function(theta, N0, NV, timeTV)
{ # function GetsdSD
  result <- c()
  N0 <- N0
  k <- theta[1]
  r <- theta[2]
  sd <- theta[3]
  for (i in 1:length(NV))
  { # for
    N <- NV[i]
    timeT <- timeTV[i]
    result[i] <- eval(q3cpartialSSD)
  } # end for
  return(sum(result))
} # end function GetsdSD

GetsdK <- function(theta, N0, NV, timeTV)
{ # function GetsdK
  result <- c()
  N0 <- N0
  k <- theta[1]
  r <- theta[2]
  sd <- theta[3]
  for (i in 1:length(NV))
  { # for
    N <- NV[i]
    timeT <- timeTV[i]
    result[i] <- eval(q3cpartialSK)
  } # end for
  return(sum(result))
} # end function GetsdK

GetsdR <- function(theta, N0, NV, timeTV)
{ # function GetsdR
  result <- c()
  N0 <- N0
  k <- theta[1]
  r <- theta[2]
  sd <- theta[3]
  for (i in 1:length(NV))
  { # for
    N <- NV[i]
    timeT <- timeTV[i]
    result[i] <- eval(q3cpartialSR)
  } # end for
  return(sum(result))
} # end function GetsdR

GetUpdate <- function(theta, N0, NV, timeTV)
{ # function GetUpdate
  updateSD <- -1 * GetdSD(theta, N0, NV, timeTV) / GetsdSD(theta, N0, NV, timeTV)
  updateK <- -1 * GetdK(theta, N0, NV, timeTV) / GetsdK(theta, N0, NV, timeTV)
  updateR <- -1 * GetdR(theta, N0, NV, timeTV) / GetsdR(theta, N0, NV, timeTV)
  return(c(updateK, updateR, updateSD))
} # end function GetUpdate

NewtonQ3C <- function(theta0, N0, NV, timeTV)
{ # function NewtonQ3C
  errorBefore <- Inf
  errorAfter <- 1000
  i <- 1
  varK <- c()
  varR <- c()
  varSD <- c()
  while (is.na(abs(errorBefore - errorAfter)) == FALSE & is.nan(abs(errorBefore - errorAfter)) == FALSE & abs(errorBefore - errorAfter) >= 0.0001 & i <= 200) 
  { # while abs(errorBefore - errorAfter) >= 0.0000001 & i <= 200
    errorBefore <- errorAfter
    update <- GetUpdate(theta0, N0, NV, timeTV)
    theta0 <- theta0 + update
    errorAfter <- sum(abs(update))
    varK[i] <- theta0[1]
    varR[i] <- theta0[2]
    varSD[i] <- theta0[3]
    i <- i + 1
  } # end while abs(errorBefore - errorAfter) >= 0.0000001 & i <= 200
  
  return(c(theta0, var(varK), var(varR), var(varSD)))
} # end function NewtonQ3C


