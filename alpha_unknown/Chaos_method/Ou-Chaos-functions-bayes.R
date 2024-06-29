Gibb_OU_chaos=function(X,nbri,delta,beta,lambda,K,M,Nrow,Nbm)
{
  alpha=alpha_prior(lambda)
  
  alphas=numeric(M)
  alphas_norm=numeric(M)
  
  alphas_aux=numeric(M+K)
  alphas_aux[1]=alpha
  
  alphas_aux_norm=numeric(M+K)
  alphas_aux_norm[1]=alpha
  for(i in 2:(K+M))
  {
    #print(i)
    print(c(i,alphas_aux[i-1],alphas_aux_norm[i-1],mean(alphas_aux[1:i-1]),mean(alphas_aux_norm[1:i-1])))
    #print(alphas_aux_norm[i-1])
    alphas_aux[i]=pos_alpha(X,nbri,delta,alphas_aux[i-1],beta,Nrow,Nbm)
    alphas_aux_norm[i]=pos_alpha_norm(X,nbri,delta,alphas_aux_norm[i-1],beta,Nrow,Nbm)
    if(i>K){
      j=i-K

      alphas[j]=alphas_aux[i]
      alphas_norm[j]=alphas_aux_norm[i]
    }
    #alphas[i]=alphas_aux[i]
  }
  alphas = rbind(alphas,alphas_norm)
  return(alphas)
}
###########  paso 1
alpha_prior=function(lambda)
{
  alpha=rexp(1,lambda)
  return(alpha)
}
######  Paso 3      #distribucion posterior de alpha
pos_alpha=function(X,nbri,delta,alpha,beta,Nrow,Nbm)
{
  cb=comp_Bridge(Xd,nbri,delta,alpha,beta)
  n=length(X)
  Kv=Ks(cb,n,delta,nbri,lambda,beta,alpha)
  int=numeric(n)
  alpha=rtruncnorm(1,a = 0,b=Inf, mean = Kv[2]/(2*Kv[1]),sd = sqrt(1/Kv[1]))
  return(alpha)
}
####
pos_alpha_norm=function(X,nbri,delta,alpha,beta,Nrow,Nbm)
{
  cb=comp_Bridge(X,nbri,delta,alpha,beta)
  n=length(X)
  Kv=Ks(cb,n,delta,nbri,lambda,beta,alpha)
  int=numeric(n)
  alpha=rnorm(1, mean = Kv[2]/(2*Kv[1]),sd = sqrt(1/Kv[1]))
  return(alpha)
}
###########  paso 2  #puente completo
comp_Bridge=function(X,nbri,delta,alpha,beta)
{
  #cantidad de puentes
  k=length(X)-1
  Nrow <- 10000
  Nbm <- 100000
  cb=matrix(0,nrow=Nrow,ncol=k*(nb)+1)
  #cb=numeric(k*(nbri)+1)
  for(s in 1:Nrow){
    for(i in 1:k)
    {
      ini=(i-1)*nbri+1
      fin=i*nbri+1
      cb[s,ini:fin] <- BridgeChaos(X[i],X[i+1],alpha,beta,nbri,delta/nbri,Nbm,Nrow,ini,fin)
    }
  }
    cbm <-colMeans(cb)
  return(cbm)
}
#Calcular las constantes Ks
Ks=function(cb,n,delta,nbri,lambda,beta,alpha)
{
  del_bri=delta/nbri
  nc=length(cb)
  #Aqui debo llamarla...
  #
  #
  K1=int_path(cb^2,del_bri)[nc-1]
  K2=n*delta+(cb[1]^2-cb[nc]^2)/(2*beta^2)-2*lambda
  # print(n*delta)
  # print((cb[1]^2-cb[nc]^2)/(2*beta^2))
  # print(2*lambda)
  x1=cb[1:(nc-1)]
  x2=cb[2:nc]
  K3=lambda*exp((-1/(2*beta^2))*(sum((x2-x1)^2)/del_bri)-(n-1)*log(beta))
  K4=K3*exp((K2^2)/(8*K1))
  K=c(K1,K2,K3,K4)
  return(K)
}
#######
i_i_s <- function(X,beta,delta,nbri,tred)
{
  ## Input: 
  # X: Datos reducidos, beta: parámetro a estimar, delta: tamaño de paso ajustado
  # nbri: número de puntos entre puente, tred: vector de tiempo de datos reducidos.
  # Ouput: Evaluación de la funcion i^i_s
  
  k=length(Xd)-1
  t=matrix(0,nrow=1,ncol=k*(nb)+1)
  for(i in 1:k)
  {
    ini=(i-1)*nbri+1
    fin=i*nbri+1
    t[ini:fin] = seq(tred[i],tred[i+1],delta)
  }
  
  i_is <- matrix(0,nrow=1,ncol=k*(nb)+1)
  i_is[1] = 0
  
  for(i in 2:k)
  {
    ini=(i-1)*nbri+1
    fin=i*nbri+1
    Const = (0.9 - 1/beta0)
    i_is[ini:fin] <- (1/(tred[i]-tred[i-1]))*((tred[i]-t[ini:fin])* Const *X[i-1]+(t[ini:fin]-tred[i])* Const * X[fin])
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
