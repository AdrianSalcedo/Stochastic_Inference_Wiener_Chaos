### Calcula las Xis
fXis_1=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=Ints
  
  
  return(Xis)
}

fXis_2=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(2)*H2(Ints)
  
  return(Xis)
}

fXis_3=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(3))*H3(Ints)
  
  return(Xis)
}

fXis_4=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(4))*H4(Ints)
  
  return(Xis)
}

fXis_5=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
     Ints=BM
  #   for( j in 1:nb){
  #     Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  #   }
  
  Xis=sqrt(factorial(5))*H5(Ints)
  
  return(Xis)
}

fXis_6=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(6))*H6(Ints)
  
  return(Xis)
}



fXis_7=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  # 
  Xis=sqrt(factorial(7))*H7(Ints)
  
  return(Xis)
}


fXis_8=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(8))*H8(Ints)
  
  return(Xis)
}


fXis_9=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  # 
  Xis=sqrt(factorial(9))*H9(Ints)
  
  return(Xis)
}


fXis_10=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(10))*H10(Ints)
  
  return(Xis)
}


fXis_11=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(11))*H11(Ints)
  
  return(Xis)
}

fXis_12=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(12))*H12(Ints)
  
  return(Xis)
}
fXis_13=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(13))*H13(Ints)
  
  return(Xis)
}
fXis_14=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(14))*H14(Ints)
  
  return(Xis)
}
fXis_15=function(nb,TiempoC,BM)
{
  # senos=fsenos(nb,TiempoC)
  # cosenos=fcosenos(nb,TiempoC)
  
  TF=tail(TiempoC,1)
  
  
  Ints=BM
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],j,TiempoC,BM[j,])
  # }
  
  Xis=sqrt(factorial(15))*H15(Ints)
  
  return(Xis)
}
####hermite polynomial

H2=function(x)
{
  h=((x^2)-1)/sqrt(2)
  return(h)
}

H3=function(x)
{
  h=(x*((x^2)-3))/sqrt(factorial(3))
  return(h)
}

H4=function(x)
{
  h=(x^4-6*(x^2)+3)/sqrt(factorial(4))
  return(h)
}

H5=function(x)
{
  h=(x*((x^4)-10*(x^2)+15))/sqrt(factorial(5))
  return(h)
}

H6=function(x)
{
  h=((x^6)-15*(x^4)+45*(x^2)-15)/sqrt(factorial(6))
  return(h)
}


H7=function(x)
{
  h=(x*((x^6)-21*(x^4)+105*(x^2)-105))/sqrt(factorial(7))
  return(h)
}

H8=function(x)
{
  h=((x^8)-28*(x^6)+210*(x^4)-420*(x^2)+105)/sqrt(factorial(8))
  return(h)
}

H9=function(x)
{
  h=(x*((x^8)-36*(x^6)+378*(x^4)-1260*(x^2)+945))/sqrt(factorial(9))
  return(h)
}
H10=function(x)
{
  h=((x^10)-45*(x^8)+630*(x^6)-3150*(x^4)+4725*(x^2)-945)/sqrt(factorial(10))
  return(h)
}
H11=function(x)
{
  h=(x*(((x^10)-55*(x^8)+990*(x^6)-6930*(x^4)+1735*(x^2)-10395)))/sqrt(factorial(11))
  return(h)
}

H12=function(x)
{
  h=((x^12)-66*(x^10)+1485*(x^8)-13860*(x^6)+51975*(x^4)-62370*(x^2)+10395)/sqrt(factorial(12))
  return(h)
}
H13=function(x)
{
  h=(x*((x^12)-78*(x^10)+2145*(x^8)-25740*(x^6)+135135*(x^4)-270270*(x^2)+135135))/sqrt(factorial(13))
  return(h)
}
H14=function(x)
{
  h=((x^14)-91*(x^12)+3003*(x^10)-45045*(x^8)+315315*(x^6)-945945*(x^4)+945945*(x^2)-135135)/sqrt(factorial(14))
  return(h)
}
H15=function(x)
{
  h=(x*((x^14)-105*(x^12)+4095*(x^10)-75075*(x^8)+675675*(x^6)-2837835*(x^4)+4729725*(x^2)-2027025))/sqrt(factorial(15))
  return(h)
}




fXis_2_j=function(js,TiempoC,BM)
{
  nb=length(js)
  # senos=fsenos_j(js,TiempoC)
  # cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  
  Ints=BM
  Xis=Ints[js[1]]*Ints[js[2]]
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
  #   Xis=Xis*Ints[j,]
  # }
  
  return(Xis)
}


fXis_3_12=function(js,TiempoC,BM)
{
  nb=length(js)
  # senos=fsenos_j(js,TiempoC)
  # cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  Ints=BM
  
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
  #   
  # }
  
  Xis=Ints[js[1]]*H2(Ints[js[2]])
  
  
  return(Xis)
}


fXis_3_21=function(js,TiempoC,BM)
{
  nb=length(js)
  # senos=fsenos_j(js,TiempoC)
  # cosenos=fcosenos_j(js,TiempoC)
  
  TF=tail(TiempoC,1)
  Ints=BM
  
  # for( j in 1:nb){
  #   Ints[j,]=Int_Trap_Acum_2(senos[j,],cosenos[j,],js[j],TiempoC,BM[j,])
  #   
  # }
  
  Xis=sqrt(2)*H2(Ints[js[1]])*Ints[js[2]]
  
  
  return(Xis)
}