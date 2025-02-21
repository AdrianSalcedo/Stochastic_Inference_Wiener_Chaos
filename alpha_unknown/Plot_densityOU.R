X_data0.1 <- read.csv("Data_Alphas_1_beta0.1_T10_lambda5.799999.csv",header = T)
X_data0.5 <- read.csv("Data_Alphas_1_beta0.5_T10_lambda4.csv", header = T)
X_data0.8 <- read.csv("Data_Alphas_1_beta0.8_T10_lambda1.8.csv", header = T)
X_data1.0 <- read.csv("Data_Alphas_1_beta1.0_T10_lambda0.01.csv", header = T)

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

layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
plot(X_data0.8[,2],main="",xlab ="", ylab=expression(alpha))
abline(h=mean(X_data0.8[,2]), col="red")
abline(h=1.0, col="blue")
legend(x= "topleft",inset = c(-0.05, -0.2),bty = "n", legend=c(expression(paste(alpha, " real")), expression(paste(alpha, " mean"))), 
       fill = c("blue","red"),horiz = TRUE, xpd = TRUE) 
#plot(alphas[2,],main="Normal",xlab ="Sample size", ylab=expression(alpha))
#abline(h=alpha, col="blue")
#abline(h=mean_alphas[2], col="red")
hist((X_data0.8[,2]),prob = TRUE,main = " ", 
     xlab=expression(alpha),col= "white")
lines(density(X_data0.8[,2]), col= "red")
#legend("topright",legend = c("Truncated", "Normal"),lty = c(1, 1), col = c(1, 2),
#       cex = 0.8,lwd = 2)