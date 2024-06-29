library(ggplot2)
library(dplyr)

X_data0.1 <- read.csv("Chaos_Data_alpha1.0_beta0.1_T10_lambda1000.csv",header = T)
X_data0.5 <- read.csv("Chaos_Data_alpha1.0_beta0.5_T10_lambda800.csv", header = T)
X_data0.8 <- read.csv("Chaos_Data_alpha1.0_beta0.8_T10_lambda830.csv", header = T)
X_data1.0 <- read.csv("Chaos_Data_alpha1.0_beta1.0_T10_lambda830.csv", header = T)



layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
plot(density(X_data0.1[,2]),main = "", xlab=expression(paste(beta, "=0.1")), col="black",lwd = 2)
plot(density(X_data0.5[,2]), col="red",lwd = 2,main = "",xlab = expression(paste(beta, "=0.5")))     
plot(density(X_data0.8[,2]), col="blue",lwd = 2,main = "",xlab = expression(paste(beta, "=0.8")))     
plot(density(X_data1.0[,2]), col="darkgreen",lwd = 2,main = "",xlab = expression(paste(beta, "=1.0")))     

layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
plot(density(X_data0.1[,2]),main = "", xlab=expression(paste(beta, "=0.1")), 
     col="black",lwd = 2,xlim=c(0,1.05), ylim=c(0,250))
lines(density(X_data0.5[,2]), col="red",lwd = 2) 
lines(density(X_data0.8[,2]), col="blue",lwd = 2)
lines(density(X_data1.0[,2]), col="darkgreen",lwd = 2)
legend("topright",
       legend = c(expression(paste(beta," = 0.1,",lambda, " = 5.799999")), 
                  expression(paste(beta," = 0.5,",lambda, " = 4")), 
                  expression(paste(beta," = 0.8,",lambda, " = 1.8")), 
                  expression(paste(beta," = 1.0,",lambda, " = 0.01"))),lty = c(1, 1), 
       col = c("black","red","blue","darkgreen"),cex = 1.0,lwd = 3)


