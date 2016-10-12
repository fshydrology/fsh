x <- c(68.8, 67.3, 70.2, 113.2, 79.2, 61.2, 66.4, 65.1, 115, 67.3, 
102.2, 54.4, 69.3, 54.3, 36, 64.2, 83.4, 64.2, 76.4, 159.4, 62.1, 
78.3, 74.3, 41, 101.6, 85.6, 51.4, 70.3, 81.3, 85.3, 58.4, 66.3, 
91.3, 72.8, 100, 78.4, 61.8, 83.4, 93.4, 99, 133, 101, 109, 88, 
99.6, 74, 94, 99.2, 101.6, 76.6, 84.8, 114.4, 95.8, 65.4, 114.8)

library(ismev)
library(VGAM)

gumfit <- gum.fit(x) # fit Gumbel distribution
parML <- gumfit$mle
vcov <- gumfit$cov # covariance matrix

# gradient of quantile function at F=0.99
h <- c(1,-log(-log(0.99))) 

# standard quantile error at F=0.99
SF <- sqrt(t(h)%*%vcov%*%h) #

#bounds of the confidence interval for quantile F=0.99 (T=100)
LB <- qgumbel(0.99,parML[1],parML[2])+qnorm(0.025)*SF
UB <- qgumbel(0.99,parML[1],parML[2])+qnorm(0.975)*SF