X <- read.csv("Data_Alphas_1_beta1.0_T10_lambda0.01.csv")


ICIqz <- mean(X$alphas)-1.96*(sd(X$alphas))/sqrt(5000)
ICIqz
ICDer <- mean(X$alphas)+1.96*(sd(X$alphas))/sqrt(5000)
ICDer