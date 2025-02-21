#estimacion bayesiana del OU caso 1
source("Ou-Chaos-functions-bayes.r")
source("EM_Ou_Estimador_Stefano.r")
source('Ou_Chaos_functions.R')
library('truncnorm')
#parámetros del proceso
alpha<-1.0
beta<-0.8
Ti<- 10
TF <- Ti
n<- 1000
delta <- Ti/n
TiempoC<-seq(0, Ti, length.out = n+1)
TL<-length(TiempoC)
#edo inicial
x0 <- 0
#parametro de la apriori es exponencial, con mean 1/lambda
lambda <- 1.8
#numero de puntos en el puente
nb <- 10
#tamaño de muestra para alpha
M <- 5000
#periodo de calentamiento
K<- 5000
#número de puentes a utilizar en cada iteración
#Nbm <- 1000
Nrow<- 100000

X <- read.csv("Data_Alp1_Bet_0.8_n1000_T10.csv")
X<-as.numeric(X[,2])

#X<- OU_EULER(x0,n,Ti,alpha,beta)
Est_alpha<- Est_stefano_alpha(X,delta)
Est_beta <- Est_stefano_beta(X,Est_alpha,delta)
print(c(Est_alpha,Est_beta))

Xd <- Data_reduction(X,10)
No <- length(Xd)-1
deltad <- Ti/No
par(mfrow = c(1, 2))
plot(X,type = "l")
plot(Xd,type = "l")
alpha_prior(lambda)
alphas=Gibb_OU_chaos(Xd,nb,deltad,beta,lambda,K,M,Nrow,Nbm)
mean_alphas = rowMeans(alphas)
mean_alphas
#write.csv(t(alphas),"Data_Alphas_1_beta1.0_T10_lambda0.01.csv")
quantile(alphas[1,],prob = c(0.025,0.975))
quantile(alphas[2,],prob = c(0.025,0.975))

layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
plot(alphas[1,],main="Truncated Normal",xlab ="", ylab=expression(alpha))
abline(h=mean_alphas[1], col="red")
abline(h=alpha, col="blue")
#plot(alphas[2,],main="Normal",xlab ="Sample size", ylab=expression(alpha))
#abline(h=alpha, col="blue")
#abline(h=mean_alphas[2], col="red")
hist((alphas[1,]),prob = TRUE,main = expression(paste("Densities of ", alpha)), 
     xlab=expression(alpha),col= "white")
lines(density(alphas[1,]), col= "red")
#legend("topright",legend = c("Truncated", "Normal"),lty = c(1, 1), col = c(1, 2),
#       cex = 0.8,lwd = 2)
layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
