# Metodo de Aceptacion-Rechazo Versión Uniformes.
# 
#
#



X_acep = NULL
X_rjct = NULL
j=1
for (i in 1:500){
  Condition = TRUE
  while(Condition){
    
    Y = runif(1,0,5)
    Z = max(c(1.2,1,0.9,0.7))
    Condition = Y > Z
    if (Y<=Z) {
      X_acep[i] = Y
    }
    X_rjct[j] = Y
    j=j+1
  }
}

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

i_is<- i_i_s(Xd,0.8,0.01,10,Tiempored)
##Datos redu, beta a ajustar, delta ajustado, numero de puentes, tiempo reducido.


Tiempored <- Data_reduction(TiempoC,10)
