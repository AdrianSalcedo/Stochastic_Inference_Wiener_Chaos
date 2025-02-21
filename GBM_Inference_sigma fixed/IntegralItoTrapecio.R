#integral de ito por trapecio acumulado
Int_Trap_Acum_2=function(senos,cosenos,j,t,W)
{
  npoints=length(senos)
  integral=numeric(npoints)
  
  delta=t[2]-t[1]
  TF=tail(t,1)
  for( i in 2:npoints)
  {
    integral[i]=(senos[i]*W[i]-senos[i-1]*W[i-1])-(delta*j*pi)*(W[i-1]*cosenos[i-1]+W[i-1]*cosenos[i-1])/(2*TF)
  }
  
  integral[]=sum(integral)
  return(integral)
  
}