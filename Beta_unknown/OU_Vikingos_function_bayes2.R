Gibb_OU_chaos_beta=function(Xd,nb,deltad,alpha,lambda,K,M,Nrow,Nbm)
{
  #beta=beta_prior(lambda)
  
  betas=numeric(M)
  betas_norm=numeric(M)
  
  betas_aux=numeric(M+K)
  betas_aux[1]=beta0
  z1 = 1/betas_aux[1]
  
  for(i in 2:(K+M))
  {
    #print(i)
    print(c(i,betas_aux[i-1],mean(betas_aux[1:i-1])))
    betas_aux[i]= pos_beta(Xd,nb,deltad,betas_aux[i-1],alpha,Nrow,Nbm,z1)
    z1 = betas_aux[i]
    if(i>K){
      j=i-K
      betas[j]=betas_aux[i]
    }
  }
  #betas = rbind(betas)
  return(betas)
}
###########  paso 1
beta_prior=function(lambda)
{
  beta=rexp(1,lambda)
  return(beta)
}
######  Paso 3      #distribucion posterior de alpha
pos_beta=function(Xd,nb,deltad,beta,alpha,Nrow,Nbm,z1)
{
  delta=deltad/nb
  betas2 = numeric(500)
  betas2_aux = numeric(1000)
  betas2_aux[1] = z1    
  cb = comp_Bridge(Xd,nb,delta,alpha,1/z1)
  n = length(X)
  Kv = Ks(cb,n,delta,nb,980,alpha,1/z1)
  #print(c(Kv[2]/(2*Kv[1]),sqrt(1/(2*Kv[1]))))
  for (i in 2:1000){
    Cond = TRUE
    while(Cond){
      int = numeric(n)
      z2 = rtruncnorm(1,a = 0,b=Inf, mean = Kv[2]/(2*Kv[1]),sd = sqrt(1/(2*Kv[1])))
      #print(c(z1,z2,z2/z1))
      Weigth = exp(-(lambda)*(1/z2-1/z1))*(z1/z2)^(n-1)
      Prt = min(1,Weigth)
      ut = runif(1,0,1)
      Cond = is.na(ut<=Prt)
      #print(c(i,Prt,z1,z2))
      betas2_aux[i] <- ifelse(ut<=Prt,z2,z1)
      probs[i] <- Prt
      if (i>500 && Cond == FALSE){
        j = i-500
        betas2[j] = betas2_aux[i]
        #print(c(j,betas2[j],Cond))
      }
    }
  }
  beta = mean(betas2)
  z1 = beta
  return(beta)
}
####
# pos_beta_norm=function(X,nbri,delta,beta,alpha,Nrow,Nbm)
# {
#   cb=comp_Bridge(X,nbri,delta,alpha,beta)
#   n=length(X)
#   Kv=Ks(cb,n,delta,nbri,lambda,alpha,beta)
#   int=numeric(n)
#   z1 = rgamma(1,n,1/Kv[2])
#   z2 = rgamma(1,n,1/Kv[2])
#   Weigth1 = gamma(n)*(lambda/Kv[2]) * exp(Kv[1]*z1^2+Kv[3])
#   Weigth2 = gamma(n)*(lambda/Kv[2]) * exp(Kv[1]*z2^2+Kv[3])
#   Prt = min(1,Weigth2/Weigth1)
#   ut = runif(1,0,1)
#   ut < Prt
#   if(ut < Prt){
#     beta = 1/z2
#   }else{
#     beta = 1/z1
#   }
#   return(beta)
# }
###########  paso 2  #puente completo
comp_Bridge=function(Xd,nb,delta,alpha,beta)
{
  #cantidad de puentes
  k=length(Xd)-1
  Nrow <- 1
  Nbm <- 1000
  cb=matrix(0,nrow=Nrow,ncol=k*(nb)+1)
  #cb=numeric(k*(nbri)+1)
  for(s in 1:Nrow){
    for(i in 1:k)
    {
      ini=(i-1)*nb+1
      fin=i*nb+1
      cb[s,ini:fin] <- Exact_OU(X[i],X[i+1],delta,nb,alpha,beta)
        #Bridge_MM(Xd[i],Xd[i+1],delta,nb,alpha,beta)
    }
  }
  cbm <-colMeans(cb)
  return(cbm)
}
#Calcular las constantes Ks
Ks=function(cb,n,delta,nb,lambda,alpha,beta)
{
  nc=length(cb)
  iis<-i_i_s(Xd,beta,nb,Tred)
  C1 = 2 * cb
  C2 = C1*iis
  C3 = C2/beta0
  C4 = iis^2
  C5 = 2 * C4/beta0
  C6 = C4 / beta0^2
  x1 = cb[1:(nc-1)]
  x2 = cb[2:nc]
  K1 = -((alpha/2)*(cb[1]^2-cb[nc]^2)-(1/2)*(sum((x2-x1)^2)/delta)-(1/2)*alpha^2*int_path(C4,delta)[nc-1])
  K2 = lambda*beta + (1/2)*alpha^2*int_path(C2-C5,delta)[nc-1]
  K3 = (1/2)*(alpha*(TiempoC[1001]-TiempoC[1])-alpha^2*int_path(cb^2,delta)[nc-1]-alpha^2*int_path(C6-C3,delta)[nc-1])
  K4 = lambda*exp((K2^2)/(4*K1)+K3)
  K=c(K1,K2,K3,K4)
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
################################################################################
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


Bridge_MM=function(a,b,delta,n,theta,sigma)
{
  T= delta*n
  X=OU_EULER(a,n,T,theta,sigma)
  bridge=numeric(n+1)
  ban=0
  while(ban==0){
    Y=rev(OU_EULER(b,n,T,theta,sigma))
    if(X[1]<=Y[1])
    {
      
      for(i in 2:(n+1))
      {
        
        if(X[i]>Y[i])
        {
          
          bridge[1:(i-1)]=X[1:(i-1)]
          bridge[i:(n+1)]=Y[i:(n+1)]
          ban=1
          break
        }
      }
    }
    else{
      
      for(i in 2:(n+1))
      {
        if(X[i]<Y[i])
        {
          
          bridge[1:(i-1)]=X[1:(i-1)]
          bridge[i:(n+1)]=Y[i:(n+1)]
          ban=1
          break
        }
      }
      
    }
  }
  
  return(bridge)
  
}

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