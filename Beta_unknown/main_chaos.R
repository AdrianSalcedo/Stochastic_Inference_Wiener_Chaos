#estimacion bayesiana del OU caso 1
source("Ou-Chaos-functions-bayes.r")
source("EM_Ou_Estimador_Stefano.r")
source('Ou_Chaos_functions.R')
library('truncnorm')
#parámetros del proceso
alpha <- 0.5  # parametro fijo
beta0 <- 1.1    # parametro de inicialización
beta <- 0.95   # parametro real
Ti<- 10          # Horizonte de tiempo
TF <- Ti        #
n<- 1000        # número de puntos en el proceso original
delta <- Ti/n   # tamaño de paso del proceso original
TiempoC<-seq(0, Ti, length.out = n+1)
TL<-length(TiempoC)
Cond = TRUE
x0 <- 0         #condicion inicial de SDE
#parametro de la apriori es exponencial, con mean 1/lambda
lambda <- 805
#numero de puntos en el puente
nb <- 10
#tamaño de muestra para alpha
M <- 5000
#periodo de calentamiento
K<- 5000
#número de puentes a utilizar en cada iteración
#Nbm <- 1000
Nrow<- 100000

X <- read.csv("PathX_A0.5_B1.0_T10_n1000.csv")
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
beta_prior(lambda)
betas = NULL
probs = NULL
betas = Gibb_OU_chaos_beta(Xd,nb,deltad,alpha,lambda,K,M,Nrow,Nbm)
mean_betas = mean(1/betas)
mean_betas
quantile(1/betas,prob = c(0.025,0.975))

write.csv(1/betas,"Chaos_Data_alpha0.5_beta1.0_T10_lambda805.csv")

layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))
plot(1/betas,main="Weight Normal",xlab ="", ylab=expression(beta))
abline(h=mean_betas, col="red")
abline(h=1, col="blue")

plot(density(1/betas),main = expression(paste("Densities of ", beta)), 
     xlab=expression(beta))
legend("topright",legend = c("Truncated Normal"),lty = c(1, 1), col = c(1, 2),
       cex = 0.8,lwd = 2)
layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
