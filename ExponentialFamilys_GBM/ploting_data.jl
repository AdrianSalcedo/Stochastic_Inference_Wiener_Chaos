using Plots, CSV, DataFrames
path = "D:\\DiffusionBridgesWienerExpantion\\Stochastic_Inference_Wiener_Chaos\\GBM_Inference\\"
####################################
############### Read the propagator ord. 8 and parameters system ###############
sigma = 0.3
sigma0 = 0.413
alpha = 0
alpha_real = 0.2
delta = 1/1000
lambda_2 = 0.000002 #0.0001
SigmaM = 0.00000027103
n = 1000
nb =1000
# number of bridges
npb = 10
Mm = 1
Warming = 10000
Nrow = 1
Acepp = 2000
Nbm=nb
pars2_aux = zeros(Acepp+Warming,2)

X_path1 = CSV.read(path * "Parameters_GBM_init1_alpha0.2sigma0.3_V5.csv",DataFrame)
X_path1 = X_path1[1:end,:]
mean_pars_path1 = mean(Matrix(X_path1)[10000:end,:],dims=1)
mean_alpha_path1 = mean_pars_path1[1]
mean_sigma_path1 = 1 ./sqrt(mean_pars_path1[2])
println([mean_alpha_path1,mean_sigma_path1])

p4 = plot((mean_alpha_path1)*ones(Warming+Acepp),label = false,color=:darkblue, dpi=300)
plot!(alpha_real*ones(Warming+Acepp),label = false)
scatter!(Matrix(X_path1)[:,1],ms=4,label = false)
p5 = plot((mean_sigma_path1)*ones(Acepp),label="mean")
plot!(sigma*ones(Acepp),label="real value")
scatter!(1 ./sqrt.(Matrix(X_path1)[:,2]),ms=4,label="estimated parameters")

X_path2 = CSV.read(path * "Parameters_GBM_init1_alpha0.2sigma0.3_V6.csv",DataFrame)
X_path2 = X_path2[1:end,:]
mean_pars_path2 = mean(Matrix(X_path2)[10000:end,:],dims=1)
mean_alpha_path2 = mean_pars_path2[1]
mean_sigma_path2 = 1 ./sqrt(mean_pars_path2[2])
println([mean_alpha_path2,mean_sigma_path2])

p6 = plot((mean_alpha_path2)*ones(Warming+Acepp),label = false,color=:darkblue, dpi=300)
plot!(alpha_real*ones(Warming+Acepp),label = false)
scatter!(Matrix(X_path2)[:,1],ms=4,label = false,color=:red)
p7 = plot((mean_sigma_path2)*ones(Acepp),label="mean")
plot!(sigma*ones(Acepp),label="real value")
scatter!(1 ./sqrt.(Matrix(X_path2)[:,2]),ms=4,label="estimated parameters",color=:red)

X_path3 = CSV.read(path * "Parameters_GBM_init1_alpha0.2sigma0.3_V7.csv",DataFrame)
X_path3 = X_path3[1:end,:]
mean_pars_path3 = mean(Matrix(X_path3)[10000:end,:],dims=1)
mean_alpha_path3 = mean_pars_path3[1]
mean_sigma_path3 = 1 ./sqrt(mean_pars_path3[2])
println([mean_alpha_path3,mean_sigma_path3])

p8 = plot((mean_alpha_path3)*ones(Warming+Acepp),label = false,color=:darkblue, dpi=300)
plot!(alpha_real*ones(Warming+Acepp),label = false)
scatter!(Matrix(X_path3)[:,1],ms=4,label = false,color=:blue)
p9 = plot((mean_sigma_path3)*ones(Acepp),label="mean")
plot!(sigma*ones(Acepp),label="real value")
scatter!(1 ./sqrt.(Matrix(X_path3)[:,2]),ms=4,label="estimated parameters",color=:blue)
plot(p6, p7, layout = (1,2),legend=false)

plot(p4, p5, p6, p7, p8, p9, layout = (3,2),legend=false,size = (700,700))