Gibb_OU_chaos_pars<-function(Xd,nb,deltad,lambda_1,lambda_2,K,M,Nrow,Nbm)
{
  alpha<-alpha_prior(lambda_1)
  beta <-beta_prior(lambda_2)
  alphas=numeric(M)
  alphas_aux=numeric(M+K)
  alphas_aux[1]=alpha
  
  betas=numeric(M)
  betas_aux=numeric(M+K)
  betas_aux[1]=beta
  z1 = 1/betas_aux[1]^2
  
  pars=matrix(0,M,2)
  pars_aux=matrix(0,M+K,2)
  pars_aux[1]=alpha
  
  for(i in 2:(K+M))
  {
    #print(i)
    print(c(i,alphas_aux[i-1],sqrt(1/betas_aux[i-1])))
    pars_aux[i,] = posterior(Xd,nb,deltad,betas_aux[i-1],alphas_aux[i-1],Nrow,Nbm,z1)
    alphas_aux[i] = pars_aux[i,1]
    betas_aux[i] = pars_aux[i,2]
    z1 = betas_aux[i]
    if(i>K){
      j=i-K
      pars[j,]=pars_aux[i,]
    }
  }
  return(pars)
}
###########  apriori
alpha_prior=function(lambda_1)
{
  alpha<- rexp(1,lambda_1)
  return(alpha)
}
beta_prior<-function(lambda_2)
{
  beta<- rexp(1,lambda_2)
  return(beta)
}
####### distribucion posterior
posterior<-function(Xd,nb,deltad,beta,alpha,Nrow,Nbm,z1)
{
  delta=deltad/nb
  betas2 = numeric(1)
  betas2_aux = numeric(1000)
  alpha2 = numeric(1)
  betas2_aux[1] = z1    
  Kv9_aux = numeric(1000)
  cb = comp_Bridge(Xd,nb,delta,alpha,sqrt(1/z1))
  n = length(Xd)
  Kv = Ks(cb,n,delta,nb,lambda_1,lambda_2,alpha,sqrt(1/z1))
  #print(Kv[3]/(2*Kv[2]))
  for (i in 2:1000){
    Cond = TRUE
    while(Cond){

      z2 = rgamma(1,shape = (n-1)/2,scale=1/(Kv[2]))
      J5 = Kv[6]-Kv[1]*z2-lambda_1
      J5p = Kv[6]-Kv[1]*z1-lambda_1
      J6 = Kv[3]*z2+Kv[4]*sqrt(z2)+Kv[5]
      J6p = Kv[3]*z1+Kv[4]*sqrt(z1)+Kv[5]
      Weigth_beta = ((1/sqrt(J6p))/(1/sqrt(J6)))*exp(-lambda_2*(sqrt(1/z2)-sqrt(1/z1))+(J5^2)/(4*J6)-(J5p^2)/(4*J6p))
      Prt1 = min(1,1)
      ut = runif(1,0,1)
      Cond1 = is.na(ut<=Prt1)
      betas2_aux[i] <- ifelse(ut<=Prt1,z2,z1)
      z1 <- betas2_aux[i]
      Cond = Cond1
      if (i>999 && Cond == FALSE){
        j = i-999
        betas2[j] = betas2_aux[i]
      }
    }
  }
  J5b = Kv[6] - Kv[1]*z1-lambda_1
  J6b = Kv[3]*z1+Kv[4]*sqrt(z1)+Kv[5]
  ##alpha part
  mean_alpha <- J5b/(2*J6b)
  sd_alpha <- sqrt(1/(2*J6b))
  alpha2 = rtruncnorm(1,a = 0,b=Inf, mean = mean_alpha,sd = sd_alpha)
  
  alpha = alpha2
  beta = betas2
  z1 = beta
  pars = rbind(alpha,beta)
  return(pars)
}
############ Creating Diffusion Bridges between data
comp_Bridge<-function(Xd,nb,delta,alpha,beta)
{
  #cantidad de puentes
  k=length(Xd)-1
  Nrow <- 1
  Nbm <- 1000
  cb=matrix(0,nrow=Nrow,ncol=k*(nb)+1)
  for(s in 1:Nrow){
    for(i in 1:k)
    {
      ini=(i-1)*nbri+1
      fin=i*nbri+1
      cb[s,ini:fin] <- BridgeChaos(Xd[i],Xd[i+1],alpha,beta,nb,delta,Nbm,Nrow,ini,fin)
    }
  }
  cbm <-colMeans(cb)
  return(cbm)
}
#Calcular las constantes Ks
Ks=function(cb,n,delta,nb,lambda_1,lambda_2,alpha,beta_ant)
{
  nc=length(cb)
   C1 <-i_i_s(Xd,beta,nb,Tred)
   
   x1 = Xd[1:(n-1)]
   x2 = Xd[2:n]
   K1 = (Xd[n]^2-Xd[1]^2)/2
   K2 = sum((x2-x1)^2)/(2*deltad)
   K3 = int_path((C1^2)/2,delta)[nc-1]
   K4 = int_path((1/2)*(2*cb*C1-((2*C1^2)/beta_ant)),delta)[nc-1]
   K5 = int_path((1/2)*(2*cb^2+(C1/beta_ant)^2-(2*cb*C1/beta_ant)),delta)[nc-1]
   K6 = (TiempoC[1001]-TiempoC[1])/2
   K = c(K1,K2,K3,K4,K5,K6)
  return(K)
}
################################################################################
i_i_s<-function(Xd,beta,nb,Tred)
{
  ## Input: 
  # X: Datos reducidos, beta: parámetro a estimar, delta: tamaño de paso ajustado
  # nbri: número de puntos entre puente, tred: vector de tiempo de datos reducidos.
  # Ouput: Evaluación de la funcion l^i_s
  
  k=length(Xd)-1
  t=matrix(0,nrow=1,ncol=k*(nb)+1)
  t=TiempoC
  i_is <- matrix(0,nrow=1,ncol=k*(nb)+1)
  i_is[1] = 0
  
  for(i in 1:k)
  {
    ini=(i-1)*nb+1
    fin=i*nb+1
    for (j in ini:fin){
      i_is[j] <- (1/(Tred[i+1]-Tred[i]))*((Tred[i+1]-t[j])*Xd[i]+(t[j]-Tred[i])*Xd[i+1])
    }
  }
  return(i_is)
}
#calcula integrales de trayectorias
int_path=function(path,delta)
{
  n=length(path)-1
  int=numeric(n)
  for(i in 1:n)
  {
    #int[i]=path[i+1]*delta
    int[i]=(path[i]+path[i+1])*(delta/2)
  }
  int_acu=cumsum(int)
  return(int_acu)
}
##################
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
##############################################
Exact_OU=function(a,b,delta,n,thetaOU,sigmaOU)
{
  Z=numeric(n+1)
  T=n*delta
  X=OU_EULER(a,n,T,thetaOU,sigmaOU)
  ts=seq(0,delta*n,by=delta)
  for(i in 1:(n+1))
  {
    Z[i]=X[i]+(b-X[n+1])*(exp(thetaOU*ts[i])-exp(-thetaOU*ts[i]))/(exp(thetaOU*ts[n+1])-exp(-thetaOU*ts[n+1]))
  }
  return(Z)
}
####
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


### variación cuadrática OU
OU_VC=function(X,delta)
{
  #numero de observaciones
  n=length(X)-1
  #horizonte de tiempo
  HT=n*delta
  VC=sqrt(sum((X[2:(n+1)]-X[1:n])^2)/HT)
  return(VC)
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





















