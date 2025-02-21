####### My main
source('HermiteFunction.R')
source('BridgeChaosFunction.R')
source('IntegralItoTrapecio.R')
source('Vikngos_Geometric.R')
source('Milstein_codes.R')
library(ggplot2)
#####
sigma=0.3
alpha=0.2
delta=1/1000
n=1000
nb=1000
a=1
b=1
# number of bridges
M=500
#Xms=t(as.matrix(read.csv(file = 'C:/Users/DELL/Desktop/Propagador/Propagator_order12_100MB.csv'))[,-1])
Xms=t(as.matrix(read.csv(file = 'D:/DiffusionBridgesWienerExpantion/DiffusionBridgesWienerExpantion/GeometricBrownian/Orden8/Propagator_order8_int1.0_alpha0.2_sigma0.3_MB1000.csv'))[,-1])
Lp = dim(Xms)[[1]]

TF=delta*n
TiempoC<-seq(0, TF, length.out = n+1)
TL<-length(TiempoC)
M_Y=mat.or.vec(M,TL)
Yms=mat.or.vec(Lp,TL)
Wt <- matrix(0,nrow=M,nb+1)
Chaos_bridges_8 <- Gen_Bridges(Xms,a,b,alpha,sigma,n,delta,nb,M,Wt)
p<-0
while (p<0.05) {
###########
Chau_method<- matrix(0,nrow=M,n+1)
for(j in 1:M){
  Chau_method[j,] <- Milstein_GBM(a,b,alpha,sigma,n,delta,TiempoC)
}
##################
par(mfrow = c(1, 1))
plot(TiempoC,Chau_method[M/2,],type='l',col="blue", main = "Chaos")
lines(TiempoC,Chaos_bridges_8[M/2,],type='l',col="black")
obs=501
p = ks.test(Chaos_bridges_8[,obs],Chau_method[,obs])$p.value
qqplot(Chau_method[,obs],Chaos_bridges_8[,obs],ylab="Lyons and Zheng's approach",xlab="Wiener Chaos approach")
abline(0,1)
print(p)
}


Vikingos_method=M_bridges_MM_GBM(M,a,b,delta,n,alpha,sigma)
pv = ks.test(Vikingos_method[,obs],Chau_method[,obs])$p.value
#par(mfrow = c(1, 1))
plot(TiempoC,Chau_method[M/2,],type='l',col="blue",ylim=c(-1,12), main = "Vikingos")
lines(TiempoC,Vikingos_method[M/2,],type='l',col="black",ylim=c(-1,12))
print(c(p,pv))

ks.test(Chaos_bridges_8[,obs],Vikingos_method[,obs])
