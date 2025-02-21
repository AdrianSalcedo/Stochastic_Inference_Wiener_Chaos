########################## Codigo Orden 8 ######################################
Gen_Bridges_GBM <- function(Xdi,Xdf,alpha,sigma,n,delta,nb,ini,fin){
  ini<- ini
  fin <- fin
   TiempoC1 <- TiempoC[ini:fin]
   Xms1 <- Xms[,ini:fin]
   TFc <- tail(TiempoC1,1)
   TLc=length(TiempoC1[ini:fin])
   Lpc <- fin
   Yms = matrix(0,8006,TLc)
   M_Y <- matrix(0,1,TLc)
  Yms <- Bri_Chaos_GBM(Xms1,Xdi,Xdf,alpha,sigma,n,delta,nb,M,ini,fin) ##Genera propagador del puente
  for (k in 1:Mm){
    BM <- rnorm(nb,0,1)
    #print(k)
    Xis <- c(fXis_1(nb,TiempoC1,BM),fXis_2_j(c(1,2),TiempoC1,BM),fXis_2_j(c(1,3),TiempoC1,BM),fXis_2_j(c(2,3),TiempoC1,BM),fXis_2(nb,TiempoC1,BM),fXis_3_12(c(1,2),TiempoC1,BM),fXis_3_21(c(1,2),TiempoC1,BM),fXis_3(nb,TiempoC1,BM),fXis_4(nb,TiempoC1,BM),fXis_5(nb,TiempoC1,BM),fXis_6(nb,TiempoC1,BM),fXis_7(nb,TiempoC1,BM),fXis_8(nb,TiempoC1,BM))
    M_Y[k,] <- apply(Yms[2:8006,]*Xis,2,sum)+Yms[1,]
  }
  Chaos_bridges <- M_Y #borrar solo de pruebas
  return(M_Y)
}
################################################################################
Bri_Chaos_GBM <- function(Xms1,Xdi,Xdf,alpha,sigma,n,delta,nb,M,ini,fin)
{
  Xms1 <- Xms1
  ini<- ini
  fin <- fin
  TiempoC1 <- TiempoC[ini:fin]
  TiempoC10 <- TiempoC[ini]
  Yms = matrix(0,8006,TLc)
  #### propagator for |m|>0
  for(i in 1:8006-1){
    Yms[i+1,] = c((TiempoC1[TLc]-TiempoC1[1:(TLc-1)])*Integrate_Xms(1/(TFc-TiempoC1[-TLc]),Xms1[1+i,]),0)
  }
  #### propagator for |m|=0
  int1 <- c((TiempoC1[TLc]-TiempoC1[1:(TLc-1)])*Integrate_Xms(1/(TFc-TiempoC1[-TLc]),Xms1[1,]),0)
  for(i in 1:TLc)
  {
    Yms[1,i] <- Xdi+(Xdf-Xdi)*(TiempoC1[i]-TiempoC10)/(TFc-TiempoC10)+int1[i]
  }
  return(Yms)
}

#integral con respecto a Xms
Integrate_Xms <- function(path,Xms)
{
  npoints <- length(path)
  integral<- numeric(npoints)
  for( i in 2:npoints)
  {
    integral[i] <- path[i-1]*(Xms[i]-Xms[i-1])#*(path[i-1]+path[i])/2
  }
  integral <- cumsum(integral)
  return(integral)
  
}
