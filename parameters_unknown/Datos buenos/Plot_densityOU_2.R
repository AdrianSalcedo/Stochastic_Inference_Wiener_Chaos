library(ggplot2)
library(dplyr)

X_data1 <- read.csv("Chaos_Data_alpha0.05_beta0.001_T1_lambda1_19.8lambda2_0.001.csv", header = T)
X_data2 <- read.csv("Chaos_Data_alpha0.1_beta0.01_T1_lambda1_9.8lambda2_0.001.csv",header = T)
X_data3 <- read.csv("Chaos_Data_alpha0.5_beta0.5_T1_lambda1_2.0lambda2_1.0.csv", header = T)
X_data4 <- read.csv("Chaos_Data_alpha1.0_beta0.8_T1_lambda1_0.5lambda2_0.00001.csv", header = T)



layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
plot(density(X_data0.1[,2]),main = "", xlab=expression(paste(beta, "=0.1")), col="black",lwd = 2)
plot(density(X_data0.5[,2]), col="red",lwd = 2,main = "",xlab = expression(paste(beta, "=0.5")))     
plot(density(X_data0.8[,2]), col="blue",lwd = 2,main = "",xlab = expression(paste(beta, "=0.8")))     
plot(density(X_data1.0[,2]), col="darkgreen",lwd = 2,main = "",xlab = expression(paste(beta, "=1.0")))     

layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
plot(density(X_data0.1[,2]),main = "", xlab=expression(paste(beta)), 
     col="black",lwd = 2,xlim=c(0,0.4))
lines(density(X_data0.5[,2]), col="red",lwd = 2) 
lines(density(X_data0.8[,2]), col="blue",lwd = 2)
lines(density(X_data1.0[,2]), col="darkgreen",lwd = 2)
legend("topright",
       legend = c(expression(paste(beta," = 0.01,",lambda, " = 2")), 
                  expression(paste(beta," = 0.1,",lambda, " = 2")), 
                  expression(paste(beta," = 0.2,",lambda, " = 5")), 
                  expression(paste(beta," = 0.3,",lambda, " = 0.0001"))),lty = c(1, 1), 
       col = c("black","red","blue","darkgreen"),cex = 1.0,lwd = 3)

layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
plot(X_data1[,2], col="black",lwd = 2, xlab = " ",ylab = expression( alpha))
abline(h=0.05, col="blue",lwd = 2)
abline(h=mean(X_data1[,2]), col="red")
legend(x= "topleft",inset = c(-0.05, -0.3),bty = "n", legend=c(expression(paste(alpha,",", beta, " real")), expression(paste(alpha,",",beta, " mean"))), 
       fill = c("blue","red"),horiz = TRUE, xpd = TRUE,) 
hist((X_data1[,2]),prob = TRUE,main = "",xlab=expression(alpha),col= "white",ylim = c(0,15))
lines(density(X_data1[,2]), col= "red")

plot(sqrt(1/X_data1[,3]), col="black",lwd = 2, xlab = " ",ylab = expression(beta))
abline(h=0.001, col="blue",lwd = 2)
abline(h=mean(sqrt(1/X_data1[,3])), col="red")

hist((sqrt(1/X_data1[,3])),prob = TRUE,main = "",xlab=expression(beta),col= "white")
lines(density(sqrt(1/X_data1[,3])), col= "red")


quantile(X_data4[,2],prob = c(0.025,0.975))
mean(X_data4[,2])

quantile(sqrt(1/X_data4[,3]),prob = c(0.025,0.975))
mean(sqrt(1/X_data4[,3]))
