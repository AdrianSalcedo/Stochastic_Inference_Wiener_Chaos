library(ggplot2)
library(dplyr)
## This file reads the CSV document to generate density graphs for 
## alpha, beta parameters

X_data0.1 <- read.csv("Chaos_Data_alpha0.1_beta0.1_T1_lambda1_10lambda2_23850.csv",header = T)

d1 = density(X_data0.1[,2])
d1$y = d1$y/max(d1$y)
d2 = density(1/(X_data0.1[,3]))
d2$y = d2$y/max(d2$y)
#layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))

par(mfrow = c(2, 2))
plot(X_data0.1[,2],main="Truncated Normal",xlab ="", ylab=expression(alpha))
abline(h=mean(X_data0.1[,2]), col="red")
abline(h=0.1, col="blue")
hist(X_data0.1[,2],main="",probability = TRUE,xlab=expression(paste(alpha, "=0.1")), col="white",lwd = 2)
lines(density(X_data0.1[,2]), col= "red")

plot(1/X_data0.1[,3],xlab ="", ylab=expression(beta),ylim = c(0.095,0.105))
abline(h=mean(1/X_data0.1[,3]), col="red")
abline(h=0.1, col="blue")
hist(1/X_data0.1[,3],main = "",probability = TRUE, xlab=expression(paste(beta, "=0.1")), col="white",lwd = 2)
lines(density(1/X_data0.1[,3]), col= "red")
