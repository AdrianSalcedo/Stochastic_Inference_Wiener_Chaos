##########
# a =Xd[i]
# b=Xd[i+1]
# thetaOU= alpha 
# sigmaOU=beta 
# n= nbri
# delta= delta/nbri
# nb =Nbm 
# M= Nrow
BridgeChaos=function(a,b,thetaOU,sigmaOU,n,delta,nb,M,ini,fin)
{
  OUchaos=NULL
  M_Y=mat.or.vec(M,fin)
  TiempC <- TiempoC[ini:fin]
  TF <-tail(TiempC,1)
  TL<-length(TiempoC[ini:fin])
  #### propagator for |m|=0
  X_0 <- a*exp(-thetaOU*TiempC)
  Y_0=X_0
  int1=c((TiempC[TL]-TiempC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempC[-TL]),X_0),0)
  
  for(i in 1:fin)
  {
    Y_0[i]= a+(b-a)*TiempC[i]/TF+int1[i]
  }
  #### propagator for |m|=1
  cosenos=fcosenos(nb,TiempC)
  senos=fsenos(nb,TiempC)
  Y_1=NULL
  for(i in 1:nb){
    
    X_1=sqrt(2)*sigmaOU*(i*pi*exp(-thetaOU*TiempC)-i*pi*cosenos[i,]+thetaOU*fin*senos[i,])/(sqrt(1/TF)*((i*pi)^2+(thetaOU*TF)^2))
    Y_1=rbind(Y_1,c((fin-TiempC[1:(TL-1)])*Integrate_Xms(1/(fin-TiempC[-TL]),X_1),0))
    
  }
  ##########################################
  for (k in 1:M){
    #print(k)
    
    Xis = rnorm(nb,mean=0,sd=1)
    #Xis=fXis(nb,TiempoC)
    
    
    M_Y[k,]=apply(Y_1*Xis,2,sum)+Y_0
    
  }
  return(M_Y)
}
#########generando Brownianos estandar
My_SBM=function(n,TiempoC)
{
  MBM<-mat.or.vec(n,length(TiempoC))
  l=length(TiempoC)
  delta=TiempoC[2]-TiempoC[1]
  for(i in 1:n){
    for (j in 2:l)
    {
      MBM[i,j]<- MBM[i,j-1]+rnorm(1,0,1)
    }
  }
  return(MBM)
}
############ 

#integral con respecto a Xms
Integrate_Xms=function(path,Xms)
{
  npoints=length(path)
  integral=numeric(npoints)
  for( i in 2:npoints)
  {
    integral[i]=(Xms[i]-Xms[i-1])*(path[i-1])
  }
  integral=cumsum(integral)
  
  return(integral)
  
}
### Calcula las Xis
fXis=function(nb,TiempoC)
{
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)
  Ints=BM
  Xis=sqrt(2/TF)*Ints

  return(Xis)
}

#calcula senos
fsenos=function(n,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  senos<-mat.or.vec(n,TL)
  for(i in 1:n){
    senos[i,]=sin(i*pi*TiempoC/TF)

  }
  return(senos)

}

#calcula cosenos
fcosenos=function(n,TiempoC)
{
  TL<-length(TiempoC)
  TF=tail(TiempoC,1)
  cosenos<-mat.or.vec(n,TL)
  for(i in 1:n){
    cosenos[i,]=cos(i*pi*TiempoC/TF)

  }
  return(cosenos)

}

