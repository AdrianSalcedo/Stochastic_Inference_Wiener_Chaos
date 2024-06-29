library(ggplot2)
library(dplyr)
## This file reads the CSV document to generate density graphs for 
## alpha, beta parameters

X_data0.1 <- read.csv("Chaos_Data_alpha1.0_beta1.0_T1_lambda1_0.75lambda2_2350.csv",header = T)

d1 = density(X_data0.1[,2])
d1$y = d1$y/max(d1$y)
d2 = density(1/(X_data0.1[,3]))
d2$y = d2$y/max(d2$y)
#layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE))
par(mfrow = c(1, 2))
plot(d1,main = "", xlab=expression(paste(alpha, "=1.0")), col="black",lwd = 2)
plot(d2,main = "", xlab=expression(paste(beta, "=1.0")),ylab = "", col="red",lwd = 2)