X_data0.1 <- read.csv("Data_Alphas_1_beta0.1_T1_lambda1.4.csv",header = T)
X_data0.5 <- read.csv("Data_Alphas_1_beta0.5_T1_lambda1.2.csv", header = T)
X_data0.8 <- read.csv("Data_Alphas_1_beta0.8_T1_lambda1.05.csv", header = T)
X_data1.0 <- read.csv("Data_Alphas_1_beta1.0_T1_lambda0.75.csv", header = T)

layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
plot(density(X_data0.1[,2]),main = expression(paste("Densities of ", alpha)), 
     xlab=expression(alpha), col="black",lwd = 2,ylim=c(0,1))
lines(density(X_data0.5[,2]), col="red",lwd = 2) 
lines(density(X_data0.8[,2]), col="blue",lwd = 2)
lines(density(X_data1.0[,2]), col="darkgreen",lwd = 2)
legend("topright",
       legend = c(expression(paste(beta," = 0.1,",lambda, " = 5.799999")), 
                  expression(paste(beta," = 0.5,",lambda, " = 4")), 
                  expression(paste(beta," = 0.8,",lambda, " = 1.8")), 
                  expression(paste(beta," = 1.0,",lambda, " = 0.01"))),lty = c(1, 1), 
       col = c("black","red","blue","darkgreen"),cex = 1.0,lwd = 3)


plot(alphas[1,],main="Truncated Normal",xlab ="", ylab=expression(alpha))
abline(h=mean_alphas[1], col="red")

quantile(X_data1.0[,2],prob = c(0.025,0.975))
mean(X_data1.0[,2])