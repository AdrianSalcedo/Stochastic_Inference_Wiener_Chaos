################################################################################
function Gibb_OU_chaos_pars(Xd,npb,deltad,alpha,sigma,SigmaM,lambda_2,Warming,Acepp,Nrow,Nbm)
    alpha = alpha_prior(alpha,SigmaM)
    alphas = zeros(Acepp)
    alphas_aux = zeros(Acepp+Warming)
    alphas_aux[1] = alpha
    
    pars = zeros(Acepp)
    pars_aux = zeros(Acepp+Warming)
    pars_aux[1] = alpha

    for i in 2:(Warming+Acepp)
        #print(i)
        println([i,alphas_aux[i-1]])
        pars_aux[i] = posterior(Xd,npb,deltad,sigma,alphas_aux[i-1],Nrow,Nbm)
        pars2_aux[i] = pars_aux[i]
        alphas_aux[i] = pars_aux[i]
        if i > Warming
            j = i-Warming
            pars[j] = pars_aux[i]
        end
    end
  return pars
end
#################################  apriori #####################################
function alpha_prior(alpha,lambda_1)
    alpha = rand(truncated(Normal(alpha,SigmaM);lower =0))
    return alpha 
end
function sigma_prior(lambda_2)
      sigma = rand(Exponential(lambda_2),1)
      return sigma
end
######################### distribucion posterior ###############################
function posterior(Xd,npb,deltad,sigma,alpha,Nrow,Nbm)
    delta = deltad/npb
    alphas2 = zeros(500)
    alphas2_aux = zeros(1000)
    alphas2_aux[1] = alpha 

    Kv_aux1 = zeros(1000)
    Kv_aux2 = zeros(1000)
    cb = comp_Bridge(Xd,npb,delta,alpha,sigma)
    n = length(cb)
    n1 = length(Xd)
    #Kv = Ks(cb,n1,delta,npb,lambda_2,alpha,sigma)       
    x1 = cb[1:n-1]
    x2 = cb[2:n]
    ti1 = TiempoC[1:(n-1)]
    ti2 = TiempoC[2:n]
    #Aquí z = 1/sigma^2
    c1 = log.(x2./x1)*(1/sigma-1/sigma)*delta
    c2 = (ti2.*log.(x1).-ti1.*log.(x2))*(1/sigma-1/sigma)*delta
    Vb = (1/sigma^2)*(log(cb[n])-log(cb[1]))-(1/2)*int_path(exp.(cb[2:n].+c1.*TiempoC[2:n] .+c2),delta)[n-2]
    Lb = (1/sigma^2)*(TiempoC[n]-TiempoC[1])
    mean_alpha = (Vb+alpha/SigmaM)*((Lb+1/SigmaM)^(-1)) 
    sd_alpha = sqrt((Lb+1/SigmaM)^(-1))
    #print(mean_alpha)
    alpha = rand(truncated(Normal(mean_alpha,sd_alpha);lower=0))
    alpha = alpha
    pars = alpha
  return pars
end
############ Creating Diffusion Bridges between data
function comp_Bridge(Xd,npb,delta,alpha,sigma)
    #cantidad de puentes
    k=length(Xd)-1
    Nrow = 1
    Nbm = 1000
    cb = zeros(Nrow,k*(npb)+1)
    for i in 1:k
        Xdi = Xd[i]
        Xdf = Xd[i+1]
        ini= (i-1)*npb+1
        fin= i*npb+1 
        TiempoC1 = TiempoC[ini:fin]
        Xms1 = Xms[:,ini:fin]
        TFc = last(TiempoC1)
        TLc = size(TiempoC1)[1]
        Lpc = fin
        Yms = zeros(8006,TLc)
        M_Y = zeros(1,TLc)
        cb[1,ini:fin] = Gen_Bridges_GBM(Xd[i],Xd[i+1],alpha,sigma,npb,delta,Nbm,ini,fin,TFc,TLc,Xms1,TiempoC1)
    end
    #cbm = colMeans(cb)
    return cb
end
#Calcular las constantes Ks
function Ks(cb,n,delta,nb,lambda_2,alpha,sigma)
    #alpha = 0.01
    #cb = comp_Bridge(Xd,npb,delta,alpha,1/sqrt(z1))
    nc = length(cb)
    x1 = Xd[1:(n-1)]
    x2 = Xd[2:n]
    K1 = -alpha*log(Xd[n]/Xd[1])+(1/2)*sum((log.(x1./x2)).^2)/delta+(1/2)*(alpha^2)*(TiempoC[1001]-TiempoC[1])
    K2 = (1/2)*log(Xd[1]/Xd[n])-(1/2)*sum(log.(Xd))+(1/2)*alpha*(TiempoC[1001]-TiempoC[1])
    #println([K1, K2])
    #K3 = -(1/2)*(n-1)*log(sigma)-(1/8)*(sigma^2)*(TiempoC[1001]-TiempoC[1])
    #K4 = (lambda_2*sigma^2)/(K1)
    K = [K1,K2]
    return K
end
################################################################################
#calcula integrales de trayectorias
function int_path(path,delta)
    n = length(path)-1
    int = zeros(n)
    for i in 1:n
        int[i] = (path[i]+path[i+1])*(delta/2)
    end
    int_acu = cumsum(int)
    return int_acu
end
########################################################################
function Data_reduction(X,porcent)
    if porcent<60
        Xr = X[1:Int(round(100/porcent)):end]
    else
        Xr = X[1:100/porcent:N;]
    end
    return Xr
end