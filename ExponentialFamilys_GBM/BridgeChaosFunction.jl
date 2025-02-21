function Gen_Bridges_GBM(Xdi,Xdf,alpha,sigma,n,delta,nb,ini,fin,TFc,TLc,Xms1,TiempoC1)
  Xms1 = Xms1
  TFc = TFc
  TLc = TLc
  TiempoC1 = TiempoC1
  Yms = Bri_Chaos_GBM(Xms1,Xdi,Xdf,alpha,sigma,n,delta,nb,Mm,ini,fin,TFc,TLc,TiempoC1)
  for k in 1:Mm
    BM=rand(Normal(0,1),nb)
    Xis = vcat([fXis_1(nb,TiempoC1,BM),fXis_2_j([1,2],TiempoC1,BM),fXis_2_j([1,3],TiempoC1,BM),fXis_2_j([2,3],TiempoC1,BM),fXis_2(nb,TiempoC1,BM),fXis_3_12([1,2],TiempoC1,BM),fXis_3_21([1,2],TiempoC1,BM),fXis_3(nb,TiempoC1,BM),fXis_4(nb,TiempoC1,BM),fXis_5(nb,TiempoC1,BM),fXis_6(nb,TiempoC1,BM),fXis_7(nb,TiempoC1,BM),fXis_8(nb,TiempoC1,BM)]...)
    M_Y[k,:]= sum(Yms[2:end,:].*Xis,dims=1)+transpose(Yms[1,:])
  end
  return M_Y
end
###
function Bri_Chaos_GBM(Xms1,Xdi,Xdf,alpha,sigma,n,delta,nb,M,ini,fin,TFc,TLc,TiempoC1)
  TiempoC10 = TiempoC[ini]
    #### propagator for |m|>0
    for i in 1:(8006-1)
      Yms[i+1,:] = [Integrate_Xms(1 ./ (TFc .- TiempoC1[1:TLc-1]), Xms1[1+i,:]) .* (TiempoC1[TLc] .- TiempoC1[1:(TLc-1)]);0]
    end
    #### propagator for |m|=0
    int1 = [Integrate_Xms(1 ./ (TFc .- TiempoC1[1:(TLc-1)]), Xms1[1,:]) .* (TiempoC1[TLc] .- TiempoC1[1:(TLc-1)]);0]
    for i in 1:TLc
      Yms[1,i] = Xdi+(Xdf-Xdi)*(TiempoC1[i]-TiempoC10)/(TFc-TiempoC10)+int1[i]
    end
    return Yms
end
#integral con respecto a Xms
function Integrate_Xms(path,Xms)
  npoints = length(path)
  integral= zeros(npoints)
  for i in 2:npoints
    integral[i] = path[i-1]*(Xms[i]-Xms[i-1])
  end
  integral = cumsum(integral)
  return integral
end 