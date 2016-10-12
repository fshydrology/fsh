library(copula) # fitting the copula
library(VGAM) # GEV distribution functions

# Data Input
vol <- c(12.42, 13.77, 14.96, 17.66, 14.56, 12.28, 14.92, 21.12, 10.72, 
15.79, 7.75, 62.61, 24.47, 12.43, 13.61, 6.43, 9.16, 10.28, 5.6, 
14.99, 14.16, 4.89, 3.93, 9.59, 9.24, 10.74, 11.98, 9.48, 19.88, 
21.21, 13.83, 13.23, 13.13, 9.47, 12.7, 13.34, 11.03, 15.17, 
6.74, 18.16, 18.9, 24.33, 19.26, 15.39, 10.54, 16.15, 6.36, 8.22, 
5.08, 16.01, 10.87, 6.54, 22.6, 14.14, 9.06, 12.8, 8.81, 8.35, 
11.75, 18.76, 27.03, 26.08, 40.51, 11.45, 18.24, 11.82, 14.69, 
26.09, 14.58, 7.94, 6.02) #10^6 m^3

peak <- c(79.29, 87.5, 74.76, 158.57, 55.22, 70.51, 52.95, 205.3, 65.41, 
103.07, 33.98, 903.31, 131.67, 38.51, 100.52, 52.39, 97.41, 46.44, 
31.71, 62.86, 64.28, 14.02, 15.86, 28.32, 27.18, 47.01, 51.82, 
33.41, 90.9, 218.04, 80.42, 60.88, 68.53, 56.35, 84.67, 69.94, 
65.98, 38.23, 32.28, 118.65, 104.77, 237.3, 116.67, 95.71, 41.06, 
88.63, 29.45, 45.59, 28.06, 104.21, 68.53, 29.45, 203.6, 97.41, 
38.23, 72.49, 33.13, 45.59, 154.04, 79.29, 288.83, 184.06, 297.33, 
73.06, 96.56, 73.06, 74.19, 297.33, 27.55, 42.76, 95.43) #m3/s

# Distribution parameters
par.vol <- c(10.5777,4.8896,-0.1799)
par.peak <- c(4.2703,0.8553)

# Plot of Figure 5.24
plot(vol,peak)

U <- as.vector(pgev(vol,par.vol[1],par.vol[2],-par.vol[3]))
V <- plnorm(peak,meanlog=par.peak[1],sdlog=par.peak[2])

# Plot of Figure 5.25
plot(U,V)

fcop=fitCopula(copula=gumbelCopula(),cbind(U,V))
 
# Figure 5.26
layout(t(1:2))		  
persp(gumbelCopula(fcop@estimate),dCopula,xlim=c(0,1),
	ylim=c(0,1), xlab='u', ylab='\nv',zlab='\nJoint density')
	title('(a)')
persp(gumbelCopula(fcop@estimate),pCopula,xlim=c(0,1),ylim=c(0,1),
	xlab='u', ylab='\nv',zlab='\nC(u,v)')
	title('(b)')

mymvdist <- mvdc(copula=gumbelCopula(fcop@estimate), 
	margins=c('gev','lnorm'),
	paramMargins=list(list(location=par.vol[1],scale=par.vol[2],
	shape=(-par.vol[3])),list(meanlog=par.peak[1],sdlog=par.peak[2])))

# Figure 5.27
layout(t(1:2))	  
persp(mymvdist,dMvdc,xlim=c(0,40),ylim=c(0,300)
	, xlab='x', ylab='y',zlab='\n\n Joint PDF')
	title('(a)')

persp(mymvdist,pMvdc,xlim=c(0,40),ylim=c(0,300)
	, xlab='x', ylab='y',zlab='\n Joint CDF')
	title('(b)')

