OU_EULER=function(x0,n,T,alpha,sigma)
{
  delta=T/n
  #tray
  X=numeric(n+1)
  X[1]=x0
  for(i in 1:n)
  {
    X[i+1]=X[i]-alpha*X[i]*delta+sigma*rnorm(1,0,sqrt(delta))
  }
  return(X)
}
#EMV aplha
OU_alpha=function(X,delta)
{
  #numero de observaciones
  n=length(X)-1
  #horizonte de tiempo
  num=sum(ITO_int(X,X))
  dem=sum(int(X^2,delta))
  alpha_bar=-num/dem
  return(alpha_bar)
}
##integral de Ito
ITO_int=function(path,X)
{
  n=length(X)-1
  int=numeric(n-1)
  for(i in 1:n)
  {
    int[i]=path[i]*(X[i+1]-X[i])
  }
  
  ito=cumsum(int)
  return(ito)
}
##integral de Ito
int=function(path,delta)
{
  n=length(path)-1
  int=numeric(n-1)
  for(i in 1:n)
  {
    int[i]=(path[i]+path[i+1])*(delta/2)
  }
  int=cumsum(int)
  return(int)
}
##Estimador dado por el Libro Stefano para el OU
Est_stefano_alpha <-function(X,delta)
{
  Numerator <- sum(X[1:(n-1)]*X[2:n])
  Estimator1 <- ifelse(Numerator >0, -(1/delta)*log(Numerator/sum(X [1:(n-1)]^2)),NA)
  return(Estimator1)
}
Est_stefano_beta <-function(X,alpha,delta)
{
  Estimador2 <- sqrt(2*alpha/((n-1)*(1-exp(-2*delta*alpha)))*sum((X[2:n]-X[1:(n-1)]*exp(-delta*alpha))^2))
  return(Estimador2)
}
################################################################################
Data_reduction <- function(X,porcent)
{
  N<-length(X)
  if(porcent<60){
    Xr=X[seq(1,N, by=round(100/porcent))]
  }else{
    Xr=X[seq(1,N, by=100/porcent)]
  }
  return(Xr)
}
################################################################################
################################################################################
Gibb_OU_Est=function(X,nbri,delta,alph,bet,K,M)
{
  alpha <- alph#alpha_prior(lambda)
  alphas <- numeric(M)
  alphas_aux <- numeric(M+K)
  betas_aux <- alphas_aux
  alphas_aux[1]<- alph
  #beta <- beta#alpha_prior(lambda)
  betas <- numeric(M)
  betas_aux <- numeric(M+K)
  betas_aux[1]<- bet
  
  delta_bri <- delta/nbri
  
  for(i in 2:(K+M))
  {
    cb <- comp_Bridge(X,nbri,delta_bri,alphas_aux[i-1],betas_aux[i-1])
    alphas_aux[i] <- Est_stefano_alpha(cb,delta_bri)
    betas_aux[i] <- Est_stefano_beta(cb,alphas_aux[i],delta_bri)
    if(i>K)
    { 
      j=i-K
      alphas[j] <- alphas_aux[i]
      betas[j] <- betas_aux[i]
      #print(c(j,alphas[j]))
      print(c(j,alphas[j],betas[j]))
    }
  }
  return(rbind(alphas,betas))
}
################################################################################
# M_bridges_Exact_OU=function(M,a,b,delta,n,thetaOU,sigmaOU)
# {
#   MOU=matrix(0,nrow=M,ncol=(n+1))
#   for(i in 1:M){
#     MOU[i,]=Exact_OU(a,b,delta,n,thetaOU,sigmaOU)
#   }
#   return(MOU)
# }

comp_Bridge=function(X,nbri,delta,alpha,beta)
{
  #cantidad de puentes
  k=length(X)-1
  cb=matrix(0,nrow=1,ncol=k*(nbri)+1) #numeric(k*(nbri)+1)
  for(s in 1:1){
    for(i in 1:k)
    {
      ini=(i-1)*nbri+1
      fin=i*nbri+1
      cb[s,ini:fin]= Exact_OU(X[i],X[i+1],delta,nbri,alpha,beta)
    }
  }
  cbm <-colMeans(cb)
  return(cbm)
}
##############################################
Exact_OU=function(a,b,delta,n,thetaOU,sigmaOU)
{
  Z=numeric(n+1)
  T=delta*n
  X=OU_EULER(a,n,T,thetaOU,sigmaOU)
  
  ts=seq(0,delta*n,by=delta)
  for(i in 1:(n+1))
  {
    Z[i]=X[i]+(b-X[n+1])*(exp(thetaOU*ts[i])-exp(-thetaOU*ts[i]))/
      (exp(thetaOU*ts[n+1])-exp(-thetaOU*ts[n+1]))
  }
  return(Z)
}