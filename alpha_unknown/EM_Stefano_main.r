source("EM_Ou_Estimador_Stefano.r")
################################################################################
###simula trayectorias de OU usando Euler
#parametros
#horizonte de tiempo
T=1
#discretización
n=10000
#tiempos
tim=seq(0,T,by=T/n)
delta=T/n
#deriva
alpha=1
#difusion 
beta=0.1
#punto inicial
x0=0
#numero de puntos en el puente
nb=10
#tamaño de muestra para alpha
M=1000
#periodo de calentamiento
K=1
#while ((0.98>Est_alpha && Est_alpha> 1.1) || (0.79>Est_beta && Est_beta>0.81)) {
  X<- OU_EULER(x0,n,T,alpha,beta)
  alpha_cont = OU_alpha(X,delta)
  #alpha_cont
  Est_alpha<- Est_stefano_alpha(X,delta)
  #Est_alpha
  Est_beta <- Est_stefano_beta(X,Est_alpha,delta)
  print(c(Est_alpha,Est_beta))
#}
 
Xd <- Data_reduction(X,10)
No <- length(Xd)-1
deltad <- T/No
#plot data path
par(mfrow = c(1, 2))
plot(X,type = "l")
plot(Xd,type = "l")

#Q<-comp_Bridge(Xd,nb,deltad/nb,1.5,1.3)

parameters <-Gibb_OU_Est(Xd,nb,deltad,Est_alpha,Est_beta,K,M)
rowMeans(parameters)

par(mfrow = c(1, 2))
plot(parameters[1,],main="Alpha",xlab ="Sample size", ylab="Alpha")
abline(h=rowMeans(parameters)[1], col="red")
abline(h=alpha, col="blue")
plot(parameters[2,],main="Beta",xlab ="Sample size", ylab="Beta")
abline(h=rowMeans(parameters)[2], col="red")
abline(h=beta, col="blue")

