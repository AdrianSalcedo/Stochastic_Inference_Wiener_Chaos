#estimacion bayesiana del OU caso 1
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
#lambda <- 0.75
#numero de puntos en el puente
#nb <- 10
#tamaño de muestra para alpha
#M <- 2000
#periodo de calentamiento
#K<- 1000
#número de puentes a utilizar en cada iteración
#Nrow<- 1

X<- OU_EULER(x0,n,Ti,alpha,beta)
Est_alpha<- Est_stefano_alpha(X,delta)
Est_beta <- Est_stefano_beta(X,Est_alpha,delta)
print(c(Est_alpha,Est_beta))

write.csv(X,"Data_alpha1beta0.8.csv")
