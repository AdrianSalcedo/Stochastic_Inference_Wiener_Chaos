#estimacion bayesiana del OU caso 1
source("Ou-Chaos-functions-bayes-both.r")
source("EM_Ou_Estimador_Stefano.r")
source('Ou_Chaos_functions.R')
library('truncnorm')
#parámetros del proceso
alpha_real <- 0.05 # parametro de inicialización
beta_real <- 0.001   # parametro de inicialización
Ti<- 1        # Horizonte de tiempo
TF <- Ti        #
n<- 1000        # número de puntos en el proceso original
delta <- Ti/n   # tamaño de paso del proceso original
TiempoC<-seq(0, Ti, length.out = n+1)
TL<-length(TiempoC)
Cond = TRUE
x0 <- 0         #condicion inicial de SDE
#parametro de la apriori es exponencial, con mean 1/lambda
lambda_1 <- 19.8
lambda_2 <- 0.001
#numero de puntos en el puente
nb <- 10
#tamaño de muestra para alpha
M <-5000
#periodo de calentamiento
K<- 5000
#número de puentes a utilizar en cada iteración
#Nbm <- 1000
Nrow<- 100000

X <- read.csv("PathX_A0.05_B0.001_T1_n1000.csv")
X<-as.numeric(X[,2])
Est_alpha<- Est_stefano_alpha(X,delta)
Est_beta <- Est_stefano_beta(X,Est_alpha,delta)
print(c(Est_alpha,Est_beta))

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
################################################################################
alpha_prior(lambda_1)
beta_prior(lambda_2)
parameters = NULL
probs = NULL
parameters = Gibb_OU_chaos_pars(Xd,nb,deltad,lambda_1,lambda_2,K,M,Nrow,Nbm)
mean_pars = colMeans(parameters)
mean_alpha = mean_pars[1]
mean_beta = sqrt(1/mean_pars[2])
print(c(mean_alpha,mean_beta))
sprintf("alpha =  %f",mean_alpha)
sprintf("beta = %f",mean_beta)
quantile(parameters[,1],prob = c(0.025,0.975))
quantile(sqrt(1/parameters[,2]),prob = c(0.025,0.975))

write.csv(parameters,"Chaos_Data_alpha0.05_beta0.001_T1_lambda1_19.8lambda2_0.001.csv")

layout(matrix(c(1,1), nrow = 1, ncol = 2, byrow = TRUE))
plot(parameters[,1],main="Weight Normal",xlab ="", ylab=expression(alpha))
abline(h=mean_alpha, col="red")
abline(h=alpha_real, col="blue")
hist(parameters[,1],prob = TRUE,main = expression(paste("Densities of ", alpha)), 
     xlab=expression(alpha),col= "white")
lines(density(parameters[,1]), col= "red")
##
plot(sqrt(1/parameters[,2]),main="Gamma",xlab ="", ylab=expression(beta))
abline(h=mean_beta, col="red")
abline(h=beta_real, col="blue")
hist(sqrt(1/parameters[,2]),prob = TRUE,main = expression(paste("Densities of ", beta)), 
     xlab=expression(beta),col= "white")
lines(density(sqrt(1/parameters[,2])), col= "red")




plot(density(1/parameters[,2]),main = expression(paste("Densities of ", beta)), 
     xlab=expression(beta))
legend("topright",legend = c("Truncated Normal"),lty = c(1, 1), col = c(1, 2),
       cex = 0.8,lwd = 2)
layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
