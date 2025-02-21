################################################################################
Gibb_OU_chaos_pars<-function(Xd,npb,deltad,alpha,sigma,SigmaM,lambda_2,Warming,Acepp,Nrow,Nbm)
{
  alpha <- alpha_prior(alpha,SigmaM)
  alphas <- numeric(Acepp)
  alphas_aux <- numeric(Acepp+Warming)
  alphas_aux[1] <- alpha
  
  #sigma <- sigma_prior(lambda_2)
  sigmas<- numeric(Acepp)
  sigmas_aux <- numeric(Acepp+Warming)
  sigmas_aux[1]<- sigma0
  z1 <- sigmas_aux[1]
  
  pars<- matrix(0,Acepp,2)
  pars_aux <- matrix(0,Acepp+Warming,2)
  pars_aux[1,1]<- alpha
  pars_aux[1,2]<- sigma0
  for(i in 2:(Warming+Acepp))
  {
    #print(i)
    print(c(i,alphas_aux[i-1],sigmas_aux[i-1],mean(alphas_aux[1:i-1]),mean(sigmas_aux[1:i-1])))
    pars_aux[i,] <- posterior(Xd,npb,deltad,sigmas_aux[i-1],alphas_aux[i-1],Nrow,Nbm,z1)
    alphas_aux[i] <- pars_aux[i,1]
    sigmas_aux[i] <- pars_aux[i,2]
    z1 <- sigmas_aux[i]
    if(i>Warming){
      j <- i-Warming
      pars[j,] <- pars_aux[i,]
    }
  }
  return(pars)
}
#################################  apriori #####################################
alpha_prior=function(alpha,lambda_1)
{
  alpha <- rnorm(1,alpha,SigmaM)
  return(alpha)
}
sigma_prior<-function(lambda_2)
{
  sigma <- rexp(1,2*lambda_2)
  return(sigma)
}
######################### distribucion posterior ###############################
posterior<-function(Xd,npb,deltad,sigma,alpha,Nrow,Nbm,z1)
{
  delta <- deltad/npb
  sigmas2 <- numeric(500)
  sigmas2_aux <- numeric(1000)
  alphas2 <- numeric(500)
  alphas2_aux <- numeric(1000)
  alphas2_aux[1] <- alpha 
  sigmas2_aux[1] <- z1    
  
  Kv_aux1 <- numeric(1000)
  Kv_aux2 <- numeric(1000)
  
  cb <- comp_Bridge(Xd,npb,delta,alpha,z1)
  n <- length(cb)
  Kv <- Ks(cb,n,delta,npb,lambda_2,alpha,z1)
  for (i in 2:1000){
    Cond = TRUE
    while(Cond){
      x1 <- cb[1:n-1]
      x2 <- cb[2:n]
      ti1 <-TiempoC[1:(n-1)]
      ti2 <-TiempoC[2:n]
      c1 <- log(x2/x1)*(1/z1-1/sigma0)*delta
      c2 <- (ti2*log(x1)-ti1*log(x2))*(1/z1-1/sigma0)*delta
      Vb <- (1/z1^2)*(log(n)-log(cb[1]))-(1/2)*int_path(exp(cb[2:n]+c1*TiempoC[2:n] +c2),delta)[n-2]
      Lb <- (1/z1^2)*(TiempoC[n]-TiempoC[1])
        
      mean_alpha <-(Vb+alpha/SigmaM)*(Lb+1/SigmaM)^(-1) 
      sd_alpha <- sqrt((Lb+1/SigmaM)^(-1))
      #print(mean_alpha)
      alphas2_aux[i] <- rtruncnorm(1,a = 0,b=Inf, mean = mean_alpha,sd = sd_alpha)
      
      ##sigma part
      z2 <- rgamma(1,shape = n, scale = 1/lambda_2)
      Kv_aux1 = -(1/2)*Kv[1]*(1/z1^2)-(1/8)*(z1^2)*(TiempoC[n]-TiempoC[1])
      Kv_aux2 = -(1/2)*Kv[1]*(1/z2^2)-(1/8)*(z2^2)*(TiempoC[n]-TiempoC[1])
      Weigth_sigma = exp(Kv_aux2-Kv_aux1)*((z2/z1)^(-(n-1)/2))*((z1/z2)^(n-1))
      Prt <- min(1,Weigth_sigma)
      ut <- runif(1,0,1)
      Cond <- is.na(ut<=Prt)
      sigmas2_aux[i] <- ifelse(ut<=Prt,z2,z1)
      #print(c(alphas2_aux[i],sigmas2_aux[i]))
      if (i>500 && Cond == FALSE){
        j = i-500
        alphas2[j] = alphas2_aux[i]
        sigmas2[j] = sigmas2_aux[i]
      }
    }
  }
  alpha = mean(alphas2)
  sigma = mean(sigmas2)
  z1 = sigma
  pars = rbind(alpha,sigma)
  return(pars)
}
############ Creating Diffusion Bridges between data
comp_Bridge<-function(Xd,npb,delta,alpha,sigma)
{
  #cantidad de puentes
  k=length(Xd)-1
  Nrow <- 1
  Nbm <- 1000
  cb=matrix(0,nrow=Nrow,ncol=k*(npb)+1)
    for(i in 1:k)
    {
      Xdi = Xd[i]
      Xdf =Xd[i+1]
      ini=(i-1)*npb+1
      fin=i*npb+1 
      TiempoC1 <- TiempoC[ini:fin]
      Xms1 <- Xms[,ini:fin]
      TFc <- tail(TiempoC1,1)
      TLc=length(TiempoC1[ini:fin])
      Lpc <- fin
      Yms = matrix(0,8006,TLc)
      M_Y <- matrix(0,1,TLc)
      cb[1,ini:fin] <- Gen_Bridges_GBM(Xd[i],Xd[i+1],alpha,sigma,npb,delta,Nbm,ini,fin)
    }
  cbm <-colMeans(cb)
  return(cbm)
}
#Calcular las constantes Ks
Ks=function(cb,n,delta,nb,lambda_2,alpha,sigma)
{
  nc=length(cb)
  x1 = cb[1:(nc-1)]
  x2 = cb[2:nc]
  K1 <- -2*alpha*log(cb[nc]/cb[1])+sum((log(x1/x2))^2)/delta + (alpha^2)*(TiempoC[1001]-TiempoC[1])
  K2 <- (1/2)*log(cb[1]/cb[nc])-(1/2)*sum(log(cb))+(1/2)*alpha*(TiempoC[1001]-TiempoC[1])
    K = c(K1,K2)
  return(K)
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





















