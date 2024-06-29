#estimacion bayesiana del OU caso 1
source("Ou-functions-bayes.r")
source("EM_Ou_Estimador_Stefano.r")
library('truncnorm')
#parámetros del proceso
alpha<-1.0
beta<-0.8
t0 <- 0
Ti<- 1
n<- 1000
delta <- (Ti-t0)/n
#edo inicial
x0 <- 0
#parametro de a priori es exponencial, con mean 1/lambda
lambda <- 1.05
#numero de puntos en el puente
nb <- 10
#tamaño de muestra para alpha
M <- 2000
#periodo de calentamiento
K<- 1000
#número de puentes a utilizar en cada iteración
Nrow<- 1

X <- read.csv("Data_alpha1beta0.8.csv")
X<-as.numeric(X[,2])

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
alphas = Gibb_OU_exact(Xd,nb,deltad,beta,lambda,K,M,Nrow)
mean_alphas = rowMeans(alphas[1:2,])
mean_alphas
quantile(alphas[1,],prob = c(0.025,0.975))
quantile(alphas[2,],prob = c(0.025,0.975))

layout(matrix(c(1,3,2,3), nrow = 2, ncol = 2, byrow = TRUE))
plot(alphas[1,],main="Truncated Normal",xlab ="", ylab=expression(alpha))
abline(h=mean_alphas[1], col="red")
abline(h=alpha, col="blue")
plot(alphas[2,],main="Normal",xlab ="Sample size", ylab=expression(alpha))
abline(h=alpha, col="blue")
abline(h=mean_alphas[2], col="red")
plot(density(alphas[1,]),main = expression(paste("Densities of ", alpha)), 
     xlab=expression(alpha)) + lines(density(alphas[2,]), col="red")
legend("topright",legend = c("Truncated", "Normal"),lty = c(1, 1), col = c(1, 2),
       cex = 0.8,lwd = 2)
