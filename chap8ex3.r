library(lmomco)
library(VGAM)

x=c(444.57, 79.29, 70, 87.5, 74.76, 158.57, 55.22, 70.51, 52.95, 205.3, 
65.41, 103.07, 33.98,903.31,  131.67, 38.51, 100.52, 52.39, 97.41, 
46.44, 31.71, 62.86, 64.28, 14.02, 15.86, 28.32, 27.18, 47.01, 
51.82, 33.41, 90.9, 218.04, 80.42, 60.88, 68.53, 56.35, 84.67, 
69.94, 65.98, 38.23, 32.28, 118.65, 104.77, 237.3, 116.67, 95.71, 
41.06, 88.63, 29.45, 45.59, 28.06, 104.21, 68.53, 29.45, 203.6, 
97.41, 38.23, 72.49, 33.13, 45.59, 154.04, 79.29, 288.83, 184.06, 
297.33, 73.06, 96.56, 73.06, 74.19, 297.33, 27.55, 42.76, 95.43
)

# Parameter estimation
gevMOM <- function (mom) {
	m <- mom$moments[1]	# mean
	s <- mom$sd			# standard deviation
	g <- mom$skew		# skewness coefficient

	if ( g == 1.1396 ) {# subroutine to solve Equation (5.74)
	} else {
		if ( g > 1.1396 ) {
			domain <- c(-0.333,-0.00001)
		} else {
			domain <- c(0.00001,2) 
		}
		Kappa <- uniroot(function (x)( g - sign(x)*
		(-gamma(1+3*x)+3*gamma(1+x)*gamma(1+2*x)-2*gamma(1+x)^3) 
		/(gamma(1+2*x)-gamma(1+x)^2)^(3/2) ), 
		interval=domain, tol=1e-10)$root
	}

	Alpha <- s*abs(Kappa)/( sqrt(gamma(1+2*Kappa)-gamma(1+Kappa)^2) )
	Beta <- m-(Alpha/Kappa)*(1-gamma(1+Kappa))
	res <- c(Beta,Alpha,Kappa)
	return(res)	
}
parMOM <- gevMOM(pmoms(x)) #MOM

parLMOM <- as.numeric(pargev(lmoms(x))$para) #LMOM

max.ll <- optim(par=parMOM,
	fn = function(pars) -sum(dgev(x,pars[1],pars[2],-pars[3],log=TRUE))
	,hessian=TRUE)

parML <- max.ll$par  #ML

# PPCC test statistic
cunnanePP <- ppoints(x,a=0.4)
cor(qgev(cunnanePP,parLMOM[1],parLMOM[2],-parLMOM[3]),sort(x))

# AIC and BIC
aic <- max.ll$value*2+6
bic <- 2*max.ll$value+log(length(x))*3

# Quantile estimates 
T <- c(1.0001,1.2,1.5,2,5,10,20,50,100,200,300,400,500,1000)
F <- 1-1/T

quantile <- qgev(F,parML[1],parML[2],-parML[3])

# quantile standard errors
Beta <- parML[1]
Alpha <- parML[2]
Kappa <- parML[3]
r <- rep(0,length(F))
h <- attr(numericDeriv(quote(qgev(F,Beta,Alpha,-Kappa)),
	theta=c("Beta","Alpha","Kappa")),'gradient')
vcov <- solve(max.ll$hessian)
for (k in 1:length(F)) r[k]=t(h[k,])%*%vcov%*%h[k,]
SF <- sqrt(r)


