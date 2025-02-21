#estimacion bayesiana del OU caso 1
source("Ou-Chaos-functions-bayes.r")
source("EM_Ou_Estimador_Stefano.r")
source('Ou_Chaos_functions.R')
library('truncnorm')
#parámetros del proceso
alpha <- 1.0  # parametro fijo
beta_real <- 0.4 # parametro real
Ti<- 1         # Horizonte de tiempo
TF <- Ti        #
n<- 1000        # número de puntos en el proceso original
delta <- Ti/n   # tamaño de paso del proceso original
TiempoC<-seq(0, Ti, length.out = n+1)
TL<-length(TiempoC)
Cond = TRUE
x0 <- 0         #condicion inicial de SDE
#parametro de la apriori es exponencial, con mean 1/lambda
lambda <- 0.01
#numero de puntos en el puente
nb <- 10
#tamaño de muestra para alpha
M <- 5000
#periodo de calentamiento
K<- 5000
#número de puentes a utilizar en cada iteración
#Nbm <- 1000
Nrow<- 100000

X <- read.csv("PathX_A1.0_B0.4_T1_n1000.csv")
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
#par(mfrow = c(1, 2))
#plot(X,type = "l")
#plot(Xd,type = "l")
################################################################################
################################################################################
beta_prior(lambda)
betas = NULL
probs = NULL
betas = Gibb_OU_chaos_beta(Xd,nb,deltad,alpha,lambda,K,M,Nrow,Nbm)
mean_betas = mean(1/betas^(1/2))
mean_betas
quantile(1/betas^(1/2),prob = c(0.025,0.975))

#write.csv(1/betas^(1/2),"Chaos_Data_alpha1.0_beta0.3_T1_lambda0.0001.csv")

layout(matrix(c(1,2), nrow = 1, ncol = 1, byrow = TRUE))
par(mfrow = c(1, 2))
plot(1/betas^(1/2),main="Weight Normal",xlab ="", ylab=expression(beta))
abline(h=mean_betas, col="red")
abline(h=beta_real, col="blue")

plot(density(1/betas^(1/2)),main = expression(paste("Densities of ", beta)), 
     xlab=expression(beta))
legend("topright",lty = c(1, 1), col = c(1, 2),
       cex = 0.8,lwd = 2)
layout(matrix(c(1,1,2,2), nrow = 2, ncol = 2, byrow = TRUE))
