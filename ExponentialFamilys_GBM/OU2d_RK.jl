using Plots; plotly()
using Random, Distributions
using BenchmarkTools
#System conf
n = 1000
nb  = 1000
Total_Bridge = 1000
t0 = 0
tf = 1
delta = (tf-t0)/n
Tiempo = range(t0,tf, step=delta)
TL=length(Tiempo)

#initial conditions
x01 = 1
x02 = 1
M11 =1.5
M12 =1

M21 =1
M22 =1.5
sigma11 = 1
sigma12 = 0
sigma21 = 0
sigma22 = 1
theta = [1,1]
x0 = [x01,x02]
M = [M11 M12;M21 M22]
Sigma = [sigma11 sigma12;sigma21 sigma22]

M_Y1=zeros(Total_Bridge,TL)
M_Y2=zeros(Total_Bridge,TL)
M2_Y1=zeros(Total_Bridge,TL)
M2_Y2=zeros(Total_Bridge,TL)
@btime begin
  X_0 = OU_2dRK(x0,n)
  Y_0= zeros(2,n+1)
  int1 =zeros(2,n+1)
  int1[1,:]=[Integrate2d_Xms(1 ./ (tf .- Tiempo[1:TL-1]),X_0)[1,:].* (Tiempo[TL].-Tiempo[1:(TL-1)]);0]
  int1[2,:]=[Integrate2d_Xms(1 ./ (tf .- Tiempo[1:TL-1]),X_0)[2,:].* (Tiempo[TL].-Tiempo[1:(TL-1)]);0]

  for i in 1:TL
    Y_0[1,i]= x0[1]+(theta[1]-x0[1])*Tiempo[i]/tf+int1[1,i]
    Y_0[2,i]= x0[2]+(theta[2]-x0[2])*Tiempo[i]/tf+int1[2,i]
  end

  #### propagator for |m|=1
  X_1=[]
  X_2=[]
  Y1x1=zeros(nb,n+1)
  Y1x2=zeros(nb,n+1)
  Y1_aux=zeros(2,n+1)
  Y2x1=zeros(nb,n+1)
  Y2x2=zeros(nb,n+1)
  Y2_aux=zeros(2,n+1)

  for i in 1:nb
      X_1 = XmOU1_2dRK(i,n)
      X_2 = XmOU2_2dRK(i,n)
      Y1_aux[1,:] = [Integrate2d_Xms(1 ./ (tf .- Tiempo[1:TL-1]),X_1)[1,:] .* (Tiempo[TL] .- Tiempo[1:(TL-1)]);0]
      Y1_aux[2,:] = [Integrate2d_Xms(1 ./ (tf .- Tiempo[1:TL-1]),X_1)[2,:] .* (Tiempo[TL] .- Tiempo[1:(TL-1)]);0]

      Y1x1[i,:] =Y1_aux[1,:]
      Y1x2[i,:] =Y1_aux[2,:]
      Y2_aux[1,:] = [Integrate2d_Xms(1 ./ (tf .- Tiempo[1:TL-1]),X_2)[1,:] .* (Tiempo[TL] .- Tiempo[1:(TL-1)]);0]
      Y2_aux[2,:] = [Integrate2d_Xms(1 ./ (tf .- Tiempo[1:TL-1]),X_2)[2,:] .* (Tiempo[TL] .- Tiempo[1:(TL-1)]);0]
      Y2x1[i,:] = Y2_aux[1,:]
      Y2x2[i,:] = Y2_aux[2,:]
  end

  for k in 1:Total_Bridge
      println(k)
      X1is=rand(Normal(0,1),nb)
      X2is=rand(Normal(0,1),nb)
      M_Y1[k,:]=sum(Y1x1.*X2is*X1is[1],dims=1)+sum(Y2x1.*X2is*X1is[2],dims=1)+transpose(Y_0[1,:])
      M_Y2[k,:]=sum(Y1x2.*X2is*X1is[1],dims=1)+sum(Y2x2.*X2is*X1is[2],dims=1)+transpose(Y_0[2,:])
      #M_Y[k,]=Y_0
  end
end


#################################################################################################
p1=plot(Tiempo,M_Y1[1,:],linestyle = :l)
p2=plot(Tiempo,M_Y2[1,:],linestyle = :l)
p3=plot(Tiempo,M_Y1[1,:],M_Y2[1,:],linestyle = :l,linewidth = 3)
p4=plot(Tiempo,M_Y1[2,:],M_Y2[2,:],linestyle = :l,linewidth = 3)
plot(p1,p2)
plot(p3,p4)





function OU_2dRK(x0,n)
      #tray
      X=zeros(2,n+1)
      X[:,1]=x0
      
      for i in 1:n
        RK11 = -(M[1,1]*X[1,i]+M[1,2]*X[2,i]) 
        RK12 = -(M[2,1]*X[1,i]+M[2,2]*X[2,i]) 
        RK21 = -(M[1,1]*(X[1,i]+0.5*delta*RK11)+M[1,2]*(X[2,i]+0.5*delta*RK12)) 
        RK22 = -(M[2,1]*(X[1,i]+0.5*delta*RK11)+M[2,2]*(X[2,i]+0.5*delta*RK12))
        RK31 = -(M[1,1]*(X[1,i]+0.5*delta*RK21)+M[1,2]*(X[2,i]+0.5*delta*RK22))
        RK32 = -(M[2,1]*(X[1,i]+0.5*delta*RK21)+M[2,2]*(X[2,i]+0.5*delta*RK22))
        RK41 = -(M[1,1]*(X[1,i]+delta*RK31)+M[1,2]*(X[2,i]+delta*RK32))
        RK42 = -(M[2,1]*(X[1,i]+delta*RK31)+M[2,2]*(X[1,i]+delta*RK32))
        X[1,i+1] = X[1,i] + (1/6)*delta*(RK11+2*RK21+2*RK31+RK41)
        X[2,i+1] = X[2,i] + (1/6)*delta*(RK12+2*RK22+2*RK32+RK42)
    end
    return X
end
##################
function Integrate2d_Xms(path,Xms)
      npoints=length(path)
      integral = zeros(2,npoints)
      integral1=zeros(npoints)
      integral2=zeros(npoints)
      for i in 2:npoints
        integral1[i]=(Xms[1,i]-Xms[1,i-1])*(path[i-1]+path[i])/2.0
        integral2[i]=(Xms[2,i]-Xms[2,i-1])*(path[i-1]+path[i])/2.0
      end
      integral1=cumsum(integral1)
      integral2=cumsum(integral2)
      integral[1,:] = integral1
      integral[2,:] = integral2
      
      return integral
end

function XmOU1_2dRK(j,n)
      #tray
    X= zeros(2,n+1)
    X[:,1]=[0,0]
      
      for i in 1:n
        RK11 = -(M[1,1]*X[1,i]+M[1,2]*X[2,i]) + sigma11*sqrt(2/tf)*sin(1*pi*Tiempo[i]/tf)+sigma12*sqrt(2/tf)*sin(j*pi*Tiempo[i]/tf)
        RK12 = -(M[2,1]*X[1,i]+M[2,2]*X[2,i]) + sigma21*sqrt(2/tf)*sin(1*pi*Tiempo[i]/tf)+sigma22*sqrt(2/tf)*sin(j*pi*Tiempo[i]/tf)
        RK21 = -(M[1,1]*(X[1,i]+0.5*delta*RK11)+M[1,2]*(X[2,i]+0.5*delta*RK12)) + sigma11*sqrt(2/tf)*sin(1*pi*(Tiempo[i]+delta/2)/tf)+sigma12*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK22 = -(M[2,1]*(X[1,i]+0.5*delta*RK11)+M[2,2]*(X[2,i]+0.5*delta*RK12)) + sigma21*sqrt(2/tf)*sin(1*pi*(Tiempo[i]+delta/2)/tf)+sigma22*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK31 = -(M[1,1]*(X[1,i]+0.5*delta*RK21)+M[1,2]*(X[2,i]+0.5*delta*RK22)) + sigma11*sqrt(2/tf)*sin(1*pi*(Tiempo[i]+delta/2)/tf)+sigma12*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK32 = -(M[2,1]*(X[1,i]+0.5*delta*RK21)+M[2,2]*(X[2,i]+0.5*delta*RK22)) + sigma21*sqrt(2/tf)*sin(1*pi*(Tiempo[i]+delta/2)/tf)+sigma22*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK41 = -(M[1,1]*(X[1,i]+delta*RK31)+M[1,2]*(X[2,i]+delta*RK32)) +sigma11*sqrt(2/tf)*sin(1*pi*(Tiempo[i]+delta)/tf)+sigma12*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta)/tf)
        RK42 = -(M[2,1]*(X[1,i]+delta*RK31)+M[2,2]*(X[1,i]+delta*RK32)) +sigma21*sqrt(2/tf)*sin(1*pi*(Tiempo[i]+delta)/tf)+sigma22*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta)/tf)
        X[1,i+1] = X[1,i] + (1/6)*delta*(RK11+2*RK21+2*RK31+RK41)
        X[2,i+1] = X[2,i] + (1/6)*delta*(RK12+2*RK22+2*RK32+RK42)
    end
      return X
end
######
function XmOU2_2dRK(j,n)
      #tray
      X = zeros(2,n+1)
      X[:,1]=[0,0]
      
    for i in 1:n
        RK11 = -(M[1,1]*X[1,i]+M[1,2]*X[2,i]) + sigma11*sqrt(2/tf)*sin(2*pi*Tiempo[i]/tf)+sigma12*sqrt(2/tf)*sin(j*pi*Tiempo[i]/tf)
        RK12 = -(M[2,1]*X[1,i]+M[2,2]*X[2,i]) + sigma21*sqrt(2/tf)*sin(2*pi*Tiempo[i]/tf)+sigma22*sqrt(2/tf)*sin(j*pi*Tiempo[i]/tf)
        
        RK21 = -(M[1,1]*(X[1,i]+0.5*delta*RK11)+M[1,2]*(X[2,i]+0.5*delta*RK12)) + sigma11*sqrt(2/tf)*sin(2*pi*(Tiempo[i]+delta/2)/tf)+sigma12*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK22 = -(M[2,1]*(X[1,i]+0.5*delta*RK11)+M[2,2]*(X[2,i]+0.5*delta*RK12)) + sigma21*sqrt(2/tf)*sin(2*pi*(Tiempo[i]+delta/2)/tf)+sigma22*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK31 = -(M[1,1]*(X[1,i]+0.5*delta*RK21)+M[1,2]*(X[2,i]+0.5*delta*RK22)) + sigma11*sqrt(2/tf)*sin(2*pi*(Tiempo[i]+delta/2)/tf)+sigma12*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK32 = -(M[2,1]*(X[1,i]+0.5*delta*RK21)+M[2,2]*(X[2,i]+0.5*delta*RK22)) + sigma21*sqrt(2/tf)*sin(2*pi*(Tiempo[i]+delta/2)/tf)+sigma22*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta/2)/tf)
        RK41 = -(M[1,1]*(X[1,i]+delta*RK31)+M[1,2]*(X[2,i]+delta*RK32)) +sigma11*sqrt(2/tf)*sin(2*pi*(Tiempo[i]+delta)/tf)+sigma12*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta)/tf)
        RK42 = -(M[2,1]*(X[1,i]+delta*RK31)+M[2,2]*(X[1,i]+delta*RK32)) +sigma21*sqrt(2/tf)*sin(2*pi*(Tiempo[i]+delta)/tf)+sigma22*sqrt(2/tf)*sin(j*pi*(Tiempo[i]+delta)/tf)
        X[1,i+1] = X[1,i] + (1/6)*delta*(RK11+2*RK21+2*RK31+RK41)
        X[2,i+1] = X[2,i] + (1/6)*delta*(RK12+2*RK22+2*RK32+RK42)
    end
    return X
end