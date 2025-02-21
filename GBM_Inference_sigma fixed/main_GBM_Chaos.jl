include("HermiteFunction.jl")
include("IntegralItoTrapecio.jl")
include("BridgeChaosFunction.jl")
include("GBM-Gibbs-Bayes.jl")
using Plots; plotly()
using Random, Distributions, CSV, DataFrames
using BenchmarkTools
path = "D:\\DiffusionBridgesWienerExpantion\\Stochastic_Inference_Wiener_Chaos\\GBM_Inference\\"
####################################
############### Read the propagator ord. 8 and parameters system ###############
sigma_real = 0.3
alpha = 0
alpha_real = 0.2
delta = 1/1000
lambda_2 = 100
SigmaM = 0.01
n = 1000
nb =1000
# number of bridges
npb = 10
Mm = 1
Warming = 5000
Nrow = 1
Acepp = 5000
Nbm=nb
pars2_aux = zeros(Acepp+Warming,1)
X = CSV.read(path * "DataGBM.csv",DataFrame)
Xms = CSV.read(path*"Propagator_order8_int1.0_alpha0.2_sigma0.3_MB1000.csv",DataFrame)
X = transpose(Matrix(X[:,2:end]))
Xms=transpose(Matrix(Xms[:,2:end]))
################################################################################
Lp = size(Xms)[1]
Ti= 1          # Horizonte de tiempo
TF = delta*n
TiempoC = [0:delta:TF;]
TL = length(TiempoC)
Yms = zeros(8006,11)

Xd = Data_reduction(X,10)            # Reducción de datos (simulando que tenemos un cantidad pequeña de datos reales)
Tred = Data_reduction(TiempoC,10)    # Reducción de tiempo respecto a los datos
No = length(Xd)-1                    # número de puntos en el reducido
deltad = Ti/No
################# Grafica del proceso original y reducido  #####################
plot(TiempoC,X,type = "l",ls=:dot)
plot!(Tred,Xd,type = "l",)
################################################################################
alpha_prior(alpha,SigmaM)
#sigma_prior(lambda_2)
#@btime begin
    parameters = Gibb_OU_chaos_pars(Xd,npb,deltad,alpha,sigma_real,SigmaM,lambda_2,Warming,Acepp,Nrow,Nbm)
#end
mean_pars = mean(parameters,dims=1)
mean_alpha = mean_pars[1]
print(mean_alpha)
###########################################
p1 = plot((mean_alpha)*ones(Acepp))
plot!(alpha_real*ones(Acepp))
scatter!(parameters[:,1],ms=4)
plot(p1)

##################
p2 = plot((mean_sigma)*ones(Acepp),label="mean")
plot!(sigma_real*ones(Acepp),label="real value")
scatter!(1 ./sqrt.(parameters[:,2]),ms=4,label="estimated parameters")
plot(p1, p2, layout = (1,2),legend=false)
CSV.write(path * "Parameters_GBM_init1_alpha0.2_fixedsigma0.3.csv",DataFrame(pars2_aux,:auto))

p3 = plot((mean_alpha)*ones(Warming+Acepp),label = false)
plot!(alpha_real*ones(Warming+Acepp),label = false)
scatter!(pars2_aux[:,1],ms=4,label = false)

X_path = CSV.read(path * "Parameters_GBM_init1_alpha0.2sigma0.3_V5.csv",DataFrame)
X_path = X_path[1:end,:]
mean_pars_path = mean(Matrix(X_path),dims=1)
mean_alpha_path = mean_pars_path[1]
mean_sigma_path = 1 ./sqrt(mean_pars[2])
print([mean_alpha_path,mean_sigma_path])

p4 = plot((mean_alpha_path)*ones(Warming+Acepp),label = false)
plot!(alpha_real*ones(Warming+Acepp),label = false)
scatter!(Matrix(X_path)[:,1],ms=4,label = false)
p5 = plot((mean_sigma_path)*ones(Acepp),label="mean")
plot!(sigma*ones(Acepp),label="real value")
scatter!(1 ./sqrt.(Matrix(X_path)[:,2]),ms=4,label="estimated parameters")
plot(p4, p5, layout = (1,2),legend=false)