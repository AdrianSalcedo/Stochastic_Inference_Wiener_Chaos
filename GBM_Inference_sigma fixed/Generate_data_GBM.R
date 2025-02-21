################################################################################
source('HermiteFunction.R')
source('BridgeChaosFunction.R')
source('Ou-Chaos-functions-bayes-both.r')
source('IntegralItoTrapecio.R')
library(ggplot2)
############### Read the propagator ord. 8 and parameters system ###############
sigmaV <- 0.2
alphaV <- 0.5
delta <- 1/1000
n <- 1000
nb<- 1000
x0<- 1
TF <- delta*n
TiempoC <- seq(0, TF, length.out = n+1)
M <- 10000
par <- matrix(0,M,2)
sigma <- sigmaV
alpha <- alphaV
for (i in 1:M) {
  X <- Milstein_GBM(x0,alphaV,sigmaV,n,delta,TiempoC)
  Yti <- Logrendimiento(X)
  ##Estimador alpha
  Est_sigma <- Estimador_sigma(Yti,delta)
  Est_alpha <- Estimador_alpha(Yti,Est_sigma,delta)
  print(c(Est_alpha,Est_sigma))
  alpha <- Est_alpha
  sigma <- Est_sigma
  par[i,] <-c(alpha,sigma)
  
}
par(mfrow = c(1, 2))
plot(par[,1])
lines(alphaV*rep(1,M), col="red")
plot(par[,2])
lines(sigmaV*rep(1,M), col="red")



write.csv(X,"Data_GBM_alpha0.2_sigma0.3_v2.csv")



###############################################################################

Logrendimiento <- function(X)
{
  n <- length(X)
  x1 <- X[1:(n-1)]
  x2 <- X[2:n]
  Yti <- log(x2/x1)
  return(Yti)
}

Milstein_GBM <- function(x0,alpha,sigma,n,delta,TiempoC)
{
  Sol = matrix(0,1,n+1)
  Wt = matrix(0,1,n+1)
  for (i in 2:(n+1)){
    Wt[i] = Wt[i-1] + sqrt(delta)*rnorm(1,0,1)
  }
  Sol[1] = x0
  for(k in 2:(n+1))
  {
    DeltaW = Wt[k]-Wt[k-1]
    Sol[k] = Sol[k-1] +alpha*Sol[k-1]*delta +sigma*Sol[k-1]*DeltaW +(1/2)*(sigma^2)*Sol[k-1]*(DeltaW^2-delta)
  }
  return(Sol)
}

##Estimador alpha
Estimador_alpha <-function(Yti,sigma,delta)
{
  n <- length(Yti)
  Y <- mean(Yti)
  Estimador1 <- Y/delta + (sigma^2)/2
  return(Estimador1)
}
Estimador_sigma <-function(Yti,delta)
{
  n <- length(Yti)
  Y <- mean(Yti)
  Estimador2 <- sqrt(mean((1/delta)*(Yti-Y)^2))
   return(Estimador2)
}
