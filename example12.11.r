library(ismev)

pavia <- c(24.2, 31.3, 32.5, 33.5, 20.2, 38.2, 36.7, 35.2, 35.2, 25.3, 
92.3, 30, 25.2, 50.4, 35.7, 40.5, 10.3, 40.2, 8.1, 10.2, 14.2, 
15.3, 40.2, 20.4, 20.2, 32.8, 43.2, 29.8, 42.8, 45, 34.2, 32.8, 
46.3, 31.9, 34.2, 24.3, 24.3, 24.3, 71.4, 37.4, 31.4, 24.3, 43.8, 
58.2, 34.6, 40.2, 20.8, 69, 44, 27.2, 37.2, 36.7, 49, 38.9, 59.6, 
63.3, 41.2, 46.6, 84.2, 29.5, 70.2, 43.7, 36.2, 29.8, 60.2, 28, 
31.4, 38.4, 29.4, 34, 47, 57, 36.5, 84.2, 45, 95.5, 48.5, 38, 
38.6, 26, 27, 58, 27.8, 37.5, 35.2, 27.5, 28.5, 52, 56.8, 80, 
29, 55.2, 48.4, 33.2, 27.4, 27.4, 18.2, 34.2)
# time covariate in a single clumn matrix
t <- matrix(1:length(pavia),ncol=1)

# fitting GEV and Gumbel models 
GEV0 <- gev.fit(pavia)
GEV1 <- gev.fit(pavia,ydat=t,mul=1)
GUM0 <- gum.fit(pavia)
GUM1 <- gum.fit(pavia,ydat=t,mul=1)
GUM2 <- gum.fit(pavia,ydat=t,mul=1,sigl=1,siglink=exp)

# maximum log-likelihood of each model
-GEV0$nllh
-GEV1$nllh
-GUM0$nllh
-GUM1$nllh
-GUM2$nllh

# likelihood ratio test statistic
2*(GUM1$nllh-GUM2$nllh)

# standard errors of GUM2 parameters
GUM2$se

# Rejection region of test
qnorm(0.975)*GUM2$se[4]
