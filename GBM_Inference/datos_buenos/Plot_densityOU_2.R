library(ggplot2)
library(dplyr)

X_data0.1 <- read.csv("Parameters_GBM_init1_alpha0.2sigma0.3_10KIte.csv",header = T)
X_data0.5 <- read.csv("Parameters_GBM_init1_alpha0.5sigma0.5_10KIte.csv", header = T)
X_data0.8 <- read.csv("Parameters_GBM_init1_alpha1.0sigma0.5_10KIte.csv", header = T)
X_data1.0 <- read.csv("Parameters_GBM_init1_alpha1.0sigma0.8_10KIte.csv", header = T)



layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))
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

layout(matrix(c(1,2), nrow = 2, ncol = 2, byrow = TRUE))

par(mfrow = c(2, 2))

plot(X_data0.1[5000:10000,1], col="black",lwd = 2, xlab = " ",ylab = expression(alpha))
abline(h=mean(X_data0.1[5000:10000,1]), col="red")
abline(h=0.2, col="blue")
legend(x= "topleft",inset = c(-0.05, -0.6),bty = "n", legend=c(expression(paste(alpha, " real")), expression(paste(alpha, " mean"))), 
       fill = c("blue","red"),horiz = TRUE, xpd = TRUE) 

hist((X_data0.1[5000:10000,1]),prob = TRUE,main = " ", 
     xlab=expression(alpha),col= "white")
lines(density(X_data0.1[5000:10000,1]), col= "red")

plot(1/sqrt(X_data0.1[5000:10000,2]), col="black",lwd = 2, xlab = " ",ylab = expression( beta))
abline(h=mean(1/sqrt(X_data0.1[5000:10000,2])), col="red")
abline(h=0.3, col="blue")
legend(x= "topleft",inset = c(-0.05, -0.6),bty = "n", legend=c(expression(paste(beta, " real")), expression(paste(beta, " mean"))), 
       fill = c("blue","red"),horiz = TRUE, xpd = TRUE,) 

hist((1/sqrt(X_data0.1[5000:10000,2])),prob = TRUE,main = " ", 
     xlab=expression(beta),col= "white", ylim = c(0,30))
lines(density(1/sqrt(X_data0.1[5000:10000,2])), col= "red")

quantile(X_data1.0[,2],prob = c(0.025,0.975))
mean(X_data1.0[,2])