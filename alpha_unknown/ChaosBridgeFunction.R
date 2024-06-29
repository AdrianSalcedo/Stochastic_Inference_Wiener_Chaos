##########
BridgeChaos=function(a,b,thetaOU,sigmaOU,n,delta,nb,M)
{
  OUchaos=NULL
  
  M_Y=mat.or.vec(M,TL)
  #### propagator for |m|=0
  X_0 <- a*exp(-thetaOU*TiempoC)
  Y_0=X_0
  int1=c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),X_0),0)
  
  for(i in 1:TL)
  {
    Y_0[i]= a+(b-a)*TiempoC[i]/TF+int1[i]
  }
  
  
  
  
  #### propagator for |m|=1
  cosenos=fcosenos(nb,TiempoC)
  senos=fsenos(nb,TiempoC)
  Y_1=NULL
  for(i in 1:nb){
    
    X_1=sqrt(2*TF)*sigmaOU*(i*pi*exp(-thetaOU*TiempoC)-i*pi*cosenos[i,]+thetaOU*TF*senos[i,])/((i*pi)^2+(thetaOU*TF)^2)
    Y_1=rbind(Y_1,c((TiempoC[TL]-TiempoC[1:(TL-1)])*Integrate_Xms(1/(TF-TiempoC[-TL]),X_1),0))
    
  }
  
  
  for (k in 1:M){
    print(k)
    
    
    Xis=fXis(nb,TiempoC)
    
    
    M_Y[k,]=apply(Y_1*Xis[,TL],2,sum)+Y_0
    
  }
  return(M_Y)
}

#generando Brownianos estandar
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

#integral con respecto a Xms
Integrate_Xms=function(path,Xms)
{
  npoints=length(path)
  integral=numeric(npoints)
  for( i in 2:npoints)
  {
    integral[i]=(Xms[i]-Xms[i-1])*(path[i-1]+path[i])/2.0
  }
  integral=cumsum(integral)
  
  
  return(integral)
  
}
### Calcula las Xis
fXis=function(nb,TiempoC)
{
  senos=fsenos(nb,TiempoC)
  cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  BM=My_SBM(nb,TiempoC)

  Ints=BM
  for( j in 1:nb){
    Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  }

  Xis=sqrt(2/TF)*Ints
  
  return(Xis)
}
#integral de ito por trapecio acumulado
Int_Trap_Acum_2=function(senos,cosenos,j,t,W)
{
  npoints=length(senos)
  integral=numeric(npoints)
  
  delta=t[2]-t[1]
  TF=tail(t,1)
  for( i in 2:npoints)
  {
    integral[i]=(senos[i]*W[i]-senos[i-1]*W[i-1])-(delta*j*pi/TF)*(W[i]*cosenos[i]+W[i-1]*cosenos[i-1])/2
  }
  
  integral=cumsum(integral)
  return(integral)
  
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
