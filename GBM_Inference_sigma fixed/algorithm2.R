################################################################################
source('HermiteFunction.R')
source('BridgeChaosFunction.R')
source('Ou-Chaos-functions-bayes-both.r')
source('IntegralItoTrapecio.R')
library('truncnorm')
library(ggplot2)
############### Read the propagator ord. 8 and parameters system ###############
sigma <- 0.3
sigma0 <- 0.5
alpha <- 0.00
delta <- 1/1000
lambda_2 <- 2900
SigmaM <- 0.0000003
n <- 1000
nb<- 1000
# number of bridges
npb <- 10
Mm <- 1
Warming <- 5000
Acepp <- 5000
#X <- t(read.csv(file = 'Data_GBM_alpha0.2_sigma0.3.csv'))
#Xms <- t(as.matrix(read.csv(file = 'Propagator_order8_int1.0_alpha0.2_sigma0.3_MB1000.csv'))[,-1])
################################################################################
Lp <- dim(Xms)[[1]]
Ti<- 1          # Horizonte de tiempo
TF <- delta*n
TiempoC <- seq(0, TF, length.out = n+1)
TL <- length(TiempoC)
#Yms <- mat.or.vec(1,TL)

Xd <- Data_reduction(X,10)            # Reducción de datos (simulando que tenemos un cantidad pequeña de datos reales)
Tred <- Data_reduction(TiempoC,10)    # Reducción de tiempo respecto a los datos
No <- length(Xd)-1                    # número de puntos en el reducido
deltad <- Ti/No                       # tamaño de paso en el proceso reducido
################################################################################
################################################################################
################# Grafica del proceso original y reducido  #####################
par(mfrow = c(1, 2))
plot(X,type = "l")
plot(Xd,type = "l")
################################################################################

######################  Generate the path using WCE ############################
# BM <- rnorm(nb,0,1)
# Xis <- c(fXis_1(nb,TiempoC,BM),fXis_2_j(c(1,2),TiempoC,BM),fXis_2_j(c(1,3),TiempoC,BM),fXis_2_j(c(2,3),TiempoC,BM),fXis_2(nb,TiempoC,BM),fXis_3_12(c(1,2),TiempoC,BM),fXis_3_21(c(1,2),TiempoC,BM),fXis_3(nb,TiempoC,BM),fXis_4(nb,TiempoC,BM),fXis_5(nb,TiempoC,BM),fXis_6(nb,TiempoC,BM),fXis_7(nb,TiempoC,BM),fXis_8(nb,TiempoC,BM))
# X <- apply(Xms[2:Lp,]*Xis,2,sum)+Xms[1,]
# #Xd <- Data_reduction(X,10)
# No <- length(Xd)-1
# deltad <- TF/No
# par(mfrow = c(1, 2))
# plot(X,type = "l")
# plot(Xd,type = "l")

################################################################################
alpha_prior(alpha,SigmaM)
sigma_prior(lambda_2)
parameters = NULL
probs = NULL
parameters = Gibb_OU_chaos_pars(Xd,npb,deltad,alpha,sigma,SigmaM,lambda_2,Warming,Acepp,Nrow,Nbm)
mean_pars = colMeans(parameters)
mean_alpha = mean_pars[1]
mean_sigma = 1/mean_pars[2]
print(c(mean_alpha,mean_sigma))







