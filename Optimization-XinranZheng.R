library(ggplot2)

xQuestion1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 
                4.24, -2.44, 3.29, 3.71, -2.40, 4.53, -0.07, 
                -1.05, -13.87, -2.53, -1.75)
sampleDataQ2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
                 2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)

llExpression <- expression(log(1 / (pi * (1 + (x - theta) ** 2))))
llExpressionQ2 <- expression(log((1 - cos(x - theta)) / (2 * pi)))
sampleData <- xQuestion1
thetaMomentQ2 <- asin(mean(sampleDataQ2) - pi)

LogLikelihood1 <- function(theta, xEv, llExpression)
{ # function LogLikelihood1
  result <- c()
  for (i in 1:length(xEv)) 
  {
    theta <- theta
    x <- xEv[i]
    result[i] <- eval(llExpression)
  }
  return(sum(result))
} # end function LogLikelihood1

VectorizedFirstD <- function(theta, xV, expressionFD, nameFD)
{ # function VectorizedFirstD
  x <- xV
  theta <- theta
  return(eval(D(expressionFD, nameFD)))
} # end function VectorizedFirstD

VectorizedSecondD <- function(theta, xV, expressionSD, nameSD)
{ # function VectorizedSecondD
  x <- xV
  theta <- theta
  firstD <- as.expression(D(expressionSD, nameSD))
  return(eval(D(firstD, nameSD)))
} # end function VectorizedSecondD

NewtonUpdateLogLikelihood <- function(theta, sampleData, llExpression)
{ # function NewtonUpdateLogLikelihood
  theta <- theta
  gradient <- sapply(X = sampleData, 
                     FUN = VectorizedFirstD, 
                     theta = theta,
                     expressionFD = llExpression,
                     nameFD = "theta")
  hessian <- sapply(X = sampleData, 
                    FUN = VectorizedSecondD, 
                    theta = theta,
                    expressionSD = llExpression, 
                    nameSD = "theta")
  return(-1 * rowSums(t(as.matrix(gradient))) / rowSums(t(as.matrix(hessian))))
} # end function NewtonUpdateLogLikelihood

NewtonRaphson <- function(theta, sampleData, llExpression)
{ # function NewtonRaphson
  error <- Inf
  i <- 1
  while (error >= 0.001 & i <= 500)
  { # while error >= 0.001
    update <- NewtonUpdateLogLikelihood(theta, sampleData, llExpression)
    theta <- theta + update
    error <- sum(abs(update))
    i <- i + 1
    if (error > 500)
    {
      return(NA)
    }
  } # end while error >= 0.001
  return(theta)
} # end function NewtonRaphson

FixedPoint <- function(theta, sampleData, llExpression, a)
{ # function FixedPoint
  error <- Inf
  i <- 1
  while (error >= 0.001 & i <= 500)
  { # while error >= 0.001
    gradient <- sapply(X = sampleData, 
                       FUN = VectorizedFirstD, 
                       theta = theta,
                       expressionFD = llExpression,
                       nameFD = "theta")
    update <- rowSums(t(as.matrix(gradient))) * a
    theta <- theta + update
    error <- sum(abs(update))
    #print(error)
    i <- i + 1
    if (error > 500)
    {
      return(NA)
    }
  } # end while error >= 0.001
  return(theta)
} # end function FixedPoint

FisherScore <- function(theta, sampleData, llExpression)
{ # function FisherScore
  error <- Inf
  i <- 1
  while (error >= 0.001 & i <= 500)
  { # while error >= 0.001
    gradient <- sapply(X = sampleData, 
                       FUN = VectorizedFirstD, 
                       theta = theta,
                       expressionFD = llExpression,
                       nameFD = "theta")
    update <- rowSums(t(as.matrix(gradient))) / (length(sampleData) / 2)
    theta <- theta + update
    error <- sum(abs(update))
    #print(error)
    i <- i + 1
    if (error > 500)
    {
      return(NA)
    }
  } # end while error >= 0.001
  return(theta)
} # end function FisherScore


theta <-  as.matrix(c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38))
dummy <-  seq(-5, 5, length.out = 200)
dummyQ2 <-  seq(-pi, pi, length.out = 200)
logLikelihood <- apply(as.matrix(dummy), 
                       MARGIN = 1, 
                       FUN = LogLikelihood1, 
                       xEv = sampleData, 
                       llExpression = llExpression)
logLikelihoodQ2 <- apply(as.matrix(dummyQ2), 
                         MARGIN = 1, 
                         FUN = LogLikelihood1, 
                         xEv = sampleDataQ2, 
                         llExpression = llExpressionQ2)
outDF <- data.frame(logLikelihood, rownames = dummy)
outDFQ2 <- data.frame(logLikelihoodQ2, rownames = dummyQ2)
thetaStar <- sapply(theta, 
                    FUN = NewtonRaphson, 
                    sampleData = sampleData, 
                    llExpression = llExpression)
optLogLikelihood <- apply(as.matrix(thetaStar), 
                          MARGIN = 1, 
                          FUN = LogLikelihood1, 
                          xEv = sampleData, 
                          llExpression = llExpression)
optDF <- data.frame(optLogLikelihood, rownames = thetaStar)
plotQ1 <- ggplot(outDF, 
                 aes(dummy, logLikelihood)) + 
                 geom_line() + 
                 labs(x = "θ", y = "Log Likelihood Question 1")
plotQ2 <- ggplot(outDFQ2, 
                 aes(dummyQ2, logLikelihoodQ2)) + 
                 geom_line() + 
                 labs(x = "θ", y = "Log Likelihood Quesion 2")
a <- c(1, 0.64, 0.25)
thetaStarFP1 <- sapply(theta,
                       FUN = FixedPoint,
                       sampleData = sampleData,
                       llExpression = llExpression,
                       a = a[1])
thetaStarFP2 <- sapply(theta,
                       FUN = FixedPoint,
                       sampleData = sampleData,
                       llExpression = llExpression,
                       a = a[2])
thetaStarFP3 <- sapply(theta,
                       FUN = FixedPoint,
                       sampleData = sampleData,
                       llExpression = llExpression,
                       a = a[3])
thetaStarFS <- sapply(theta,
                      FUN = FisherScore,
                      sampleData = sampleData,
                      llExpression = llExpression)
thetaStarFSNR <- sapply(thetaStarFS,
                        FUN = NewtonRaphson,
                        sampleData = sampleData,
                        llExpression = llExpression)
thetaStarQ2c <- sapply(thetaMomentQ2,
                       FUN = NewtonRaphson,
                       sampleData = sampleDataQ2,
                       llExpression = llExpressionQ2)
thetaQ2d <- c(-2.7, 2.7)
thetaStarQ2d <- sapply(thetaQ2d,
                       FUN = NewtonRaphson,
                       sampleData = sampleDataQ2,
                       llExpression = llExpressionQ2)
thetaStarQ2e <- sapply(dummyQ2,
                       FUN = NewtonRaphson,
                       sampleData = sampleDataQ2,
                       llExpression = llExpressionQ2)
optLogLikelihood1 <- apply(as.matrix(thetaStarFP1), 
                           MARGIN = 1, 
                           FUN = LogLikelihood1, 
                           xEv = sampleData, 
                           llExpression = llExpression)
optLogLikelihood2 <- apply(as.matrix(thetaStarFP2), 
                           MARGIN = 1, 
                           FUN = LogLikelihood1, 
                           xEv = sampleData, 
                           llExpression = llExpression)
optLogLikelihood3 <- apply(as.matrix(thetaStarFP3), 
                           MARGIN = 1, 
                           FUN = LogLikelihood1, 
                           xEv = sampleData, 
                           llExpression = llExpression)
optLogLikelihood4 <- apply(as.matrix(thetaStarFS), 
                           MARGIN = 1, 
                           FUN = LogLikelihood1, 
                           xEv = sampleData, 
                           llExpression = llExpression)
optLogLikelihood5 <- apply(as.matrix(thetaStarFSNR), 
                           MARGIN = 1, 
                           FUN = LogLikelihood1, 
                           xEv = sampleData, 
                           llExpression = llExpression)
optLogLikelihoodQ2c <- apply(as.matrix(thetaStarQ2c), 
                           MARGIN = 1, 
                           FUN = LogLikelihood1, 
                           xEv = sampleDataQ2, 
                           llExpression = llExpressionQ2)
optLogLikelihoodQ2d <- apply(as.matrix(thetaStarQ2d), 
                             MARGIN = 1, 
                             FUN = LogLikelihood1, 
                             xEv = sampleDataQ2, 
                             llExpression = llExpressionQ2)
optLogLikelihoodQ2e <- apply(as.matrix(thetaStarQ2e), 
                             MARGIN = 1, 
                             FUN = LogLikelihood1, 
                             xEv = sampleDataQ2, 
                             llExpression = llExpressionQ2)
optLogLikelihood1 <- data.frame(optLogLikelihood1, rownames = thetaStarFP1)
optLogLikelihood2 <- data.frame(optLogLikelihood2, rownames = thetaStarFP2)
optLogLikelihood3 <- data.frame(optLogLikelihood3, rownames = thetaStarFP3)
optLogLikelihood4 <- data.frame(optLogLikelihood4, rownames = thetaStarFS)
optLogLikelihood5 <- data.frame(optLogLikelihood5, rownames = thetaStarFSNR)
optLogLikelihoodQ2c <- data.frame(optLogLikelihoodQ2c, rownames = thetaStarQ2c)
optLogLikelihoodQ2d <- data.frame(optLogLikelihoodQ2d, rownames = thetaStarQ2d)
optLogLikelihoodQ2e <- data.frame(optLogLikelihoodQ2e, rownames = thetaStarQ2e)

plot2 <- plotQ1 + geom_point(data = optDF, 
                            mapping = aes(x = thetaStar, y = optDF[[1]]), 
                            color = "red") + 
                            labs(caption = "Newton Raphson") + 
                            theme(plot.caption = element_text(hjust = 0.5))
plot3 <- plotQ1 + geom_point(data = optLogLikelihood1, 
                            mapping = aes(x = thetaStarFP1, y = optLogLikelihood1[[1]]), 
                            color = "blue") + 
                            labs(caption = "Fixed Point, a = 1") + 
                            theme(plot.caption = element_text(hjust = 0.5))
plot4 <- plotQ1 + geom_point(data = optLogLikelihood2, 
                            mapping = aes(x = thetaStarFP2, y = optLogLikelihood2[[1]]), 
                            color = "yellow") + 
                            labs(caption = "Fixed Point, a = 0.64") + 
                            theme(plot.caption = element_text(hjust = 0.5))
plot5 <- plotQ1 + geom_point(data = optLogLikelihood3, 
                            mapping = aes(x = thetaStarFP3, y = optLogLikelihood3[[1]]), 
                            color = "green") + 
                            labs(caption = "Fixed Point, a = 0.25") + 
                            theme(plot.caption = element_text(hjust = 0.5))
plot6 <- plotQ1 + geom_point(data = optLogLikelihood4, 
                            mapping = aes(x = thetaStarFS, y = optLogLikelihood4[[1]]), 
                            color = "brown") + 
                            labs(caption = "Fisher Score") + 
                            theme(plot.caption = element_text(hjust = 0.5))

plot7 <- plotQ1 + geom_point(data = optLogLikelihood5, 
                            mapping = aes(x = thetaStarFSNR, y = optLogLikelihood5[[1]]), 
                            color = "purple") + 
                            labs(caption = "Fisher Score, Refined by Newton Raphson") + 
                            theme(plot.caption = element_text(hjust = 0.5))
plotQ2c <- plotQ2 + geom_point(data = optLogLikelihoodQ2c, 
                               mapping = aes(x = thetaStarQ2c, y = optLogLikelihoodQ2c[[1]]), 
                               color = "red") + 
                               labs(caption = "Newton-Raphson, method-of-moments estimator") + 
                               theme(plot.caption = element_text(hjust = 0.5))
plotQ2d <- plotQ2 + geom_point(data = optLogLikelihoodQ2d, 
                               mapping = aes(x = thetaStarQ2d, y = optLogLikelihoodQ2d[[1]]), 
                               color = "blue") + 
                               labs(caption = "Newton-Raphson, θ = -2.7 and θ = 2.7") + 
                               theme(plot.caption = element_text(hjust = 0.5))


thetaStarQ2e <- round(thetaStarQ2e, 4)
margin <- c(0)
for (i in 1:(length(thetaStarQ2e) - 1))
{
  if (thetaStarQ2e[i] != thetaStarQ2e[i + 1])
  {
    margin <- c(margin, i)
  }
}


intervals <- c()
for (i in 1:length(margin))
{
  intervals <- c(intervals, dummyQ2[margin[i]])
}


plot2
plot3
plot4
plot5
plot6
plot7
plotQ2c
plotQ2d
