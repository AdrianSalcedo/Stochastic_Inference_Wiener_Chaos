-################################################################################
function Gibb_OU_chaos_pars(Xd,npb,deltad,alpha0,sigma,SigmaM,lambda_2,Warming,Acepp,Nrow,Nbm)
    alpha = alpha_prior(alpha0,SigmaM)
    alphas = zeros(Acepp)
    alphas_aux = zeros(Acepp+Warming)
    alphas_aux[1] = alpha
    
    sigma0 = sigma_prior(lambda_2)
    sigmas= zeros(Acepp)
    sigmas_aux = zeros(Acepp+Warming)
    sigmas_aux[1] = (1/sigma0[1])^2
    z1 = sigmas_aux[1]

    pars = zeros(Acepp,2)
    pars_aux = zeros(Acepp+Warming,2)

    pars_aux[1,1] = alpha
    pars_aux[1,2] = 1/sigma0[1]^2
    for i in 2:(Warming+Acepp)
        #print(i)
        println([i,alphas_aux[i-1],1/sqrt(sigmas_aux[i-1])])
        pars_aux[i,:] = posterior(Xd,npb,deltad,sigmas_aux[i-1],alphas_aux[i-1],Nrow,Nbm,z1)
        pars2_aux[i,:] = pars_aux[i,:]
        alphas_aux[i] = pars_aux[i,1]
        sigmas_aux[i] = pars_aux[i,2]
        z1 = sigmas_aux[i]
        if i > Warming
            j = i-Warming
            pars[j,:] = pars_aux[i,:]
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
function posterior(Xd,npb,deltad,sigma_ant,alpha,Nrow,Nbm,z1)
    delta = deltad/npb
    sigmas2 = zeros(1)
    sigmas2_aux = zeros(1000)
    alphas2 = zeros(500)
    alphas2_aux = zeros(1000)
    alphas2_aux[1] = alpha 
    sigmas2_aux[1] = z1    
    Kv_aux1 = zeros(1000)
    Kv_aux2 = zeros(1000)
    cb = comp_Bridge(Xd,npb,delta,alpha,1/sqrt(z1))
    n = length(cb)
    n1 = length(Xd)
    Kv = Ks(cb,n,deltad,npb,lambda_2,alpha,1/sqrt(z1))       
    for i in 2:1000
        Cond = true
        while Cond 
                    
            ##sigma part
            z2 = rand(truncated(Gamma((n1+4)/4,1/Kv[1]);lower=0))
            Kv_aux1 = -(TiempoC[n]-TiempoC[1])/(8*z1)-lambda_2*sqrt(1/z1)
            Kv_aux2 = -(TiempoC[n]-TiempoC[1])/(8*z2)-lambda_2*sqrt(1/z2)

            Weigth_sigma = exp(Kv_aux2-Kv_aux1)
            Prt = min(1,Weigth_sigma)
            ut = rand()
            Cond = ismissing.(ut<=Prt)
            if ut<=Prt
                sigmas2_aux[i] = z2
            else
                sigmas2_aux[i] = z1
            end
            z1 = sigmas2_aux[i]
            #print(c(alphas2_aux[i],sigmas2_aux[i]))
            if i>999 && Cond == false
                j = i-999
                #alphas2[j] = alphas2_aux[i]
                sigmas2[j] = sigmas2_aux[i]
            end
        end
    end
    x1 = cb[1:n-1]
    x2 = cb[2:n]
    ti1 = TiempoC[1:(n-1)]
    ti2 = TiempoC[2:n]
    #Aquí z = 1/sigma^2
    c1 = log.(x2./x1)*(sqrt(z1)-1/sigma_ant)*delta
    c2 = (ti2.*log.(x1).-ti1.*log.(x2))*(sqrt(z1)-1/sigma_ant)*delta
    Vb = (z1)*(log(cb[n])-log(cb[1]))-(1/2)*int_path(exp.(cb[2:n].+c1.*TiempoC[2:n] .+c2),delta)[n-2]
    Lb = (z1)*(TiempoC[n]-TiempoC[1])
    mean_alpha = (Vb+alpha/SigmaM)*((Lb+1/SigmaM)^(-1)) 
    sd_alpha = sqrt((Lb+1/SigmaM)^(-1))
    #print(mean_alpha)
    alpha = rand(truncated(Normal(mean_alpha,sd_alpha);lower=0))
    alpha = alpha
    sigma = mean(sigmas2)
    z1 = sigma
    pars = [alpha,sigma]
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
    #n = length(cb)
    x1 = cb[1:(n-1)]
    x2 = cb[2:n]
    K1 = -alpha*log(cb[n]/cb[1])+(1/2)*sum((log.(x1./x2)).^2)/delta+(1/2)*(alpha^2)*(TiempoC[1001]-TiempoC[1])
    K2 = (1/2)*log(cb[1]/cb[n])-(1/2)*sum(log.(Xd))+(1/2)*alpha*(TiempoC[1001]-TiempoC[1])
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