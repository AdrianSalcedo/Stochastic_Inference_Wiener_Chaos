library(ggplot2)
library(dplyr)
beta=1
X_data0.1 <- read.csv("Chaos_Data_alpha0.1_beta1.0_T10_lambda805.csv",header = T)
X_data0.5 <- read.csv("Chaos_Data_alpha0.5_beta1.0_T10_lambda805.csv", header = T)
X_data0.8 <- read.csv("Chaos_Data_alpha0.8_beta1.0_T10_lambda840.csv", header = T)
X_data1.0 <- read.csv("Chaos_Data_alpha1.0_beta1.0_T10_lambda815.csv", header = T)

d1 = density(X_data0.1[,2])
d1$y = d1$y/max(d1$y)

d2 = density(X_data0.5[,2])
d2$y = d2$y/max(d2$y)

d3 = density(X_data0.8[,2])
d3$y = d3$y/max(d3$y)

d4 = density(X_data1.0[,2])
d4$y = d4$y/max(d4$y)

layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
plot(d1,main = "", xlab=expression(paste(beta, "=1.0")), 
     col="black",lwd = 2)
lines(d2, col="red",lwd = 2) 
lines(d3, col="blue",lwd = 2)
lines(d4, col="darkgreen",lwd = 2)
legend("topright",
       legend = c(expression(paste(alpha," = 0.1,", lambda, " = 805")), 
                  expression(paste(alpha," = 0.5,", lambda, " = 805")), 
                  expression(paste(alpha," = 0.8,", lambda, " = 840")), 
                  expression(paste(alpha," = 1.0,", lambda, " = 815"))),lty = c(1, 1), 
       col = c("black","red","blue","darkgreen"),cex = 1.0,lwd = 3)


layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
plot(1/X_data0.5[,2],main="Truncated Normal",xlab ="", ylab=expression(beta),ylim=c(0.99,1.01))
abline(h=mean(1/X_data0.5[,2]), col="red")
abline(h=beta, col="blue")
#plot(alphas[2,],main="Normal",xlab ="Sample size", ylab=expression(alpha))
#abline(h=alpha, col="blue")
#abline(h=mean_alphas[2], col="red")
hist(1/X_data0.5[,2],prob = TRUE,main = expression(paste("Densities of ", beta)), 
     xlab=expression(beta),col= "white", ylim = c(0,600))
lines(density(1/X_data0.5[,2]), col= "red")
