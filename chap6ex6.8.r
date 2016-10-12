x <- c(68.8, 67.3, 70.2, 113.2, 79.2, 61.2, 66.4, 65.1, 115, 67.3, 
102.2, 54.4, 69.3, 54.3, 36, 64.2, 83.4, 64.2, 76.4, 159.4, 62.1, 
78.3, 74.3, 41, 101.6, 85.6, 51.4, 70.3, 81.3, 85.3, 58.4, 66.3, 
91.3, 72.8, 100, 78.4, 61.8, 83.4, 93.4, 99, 133, 101, 109, 88, 
99.6, 74, 94, 99.2, 101.6, 76.6, 84.8, 114.4, 95.8, 65.4, 114.8)

library(lmomco) # LMOM estimators
library(ismev) # GEV ML estimators
library(VGAM) # GEV functions

# MOM estimators of the GEV
gevMOM=function (mom) {
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

# Parameter estimation
parMOM <- gevMOM(pmoms(x)) # MOM
parLMOM <- as.numeric(pargev(lmoms(x))$para) # LMOM
parML <- gev.fit(x,quiet=TRUE)$mle # ML

# Exceedance probabilities of quantile x=150
1-pgev(q=150,location=parMOM[1],scale=parMOM[2],shape=-parMOM[3])
1-pgev(q=150,location=parLMOM[1],scale=parLMOM[2],shape=-parLMOM[3])
1-pgev(q=150,location=parML[1],scale=parML[2],shape=parML[3])

# Quantile with F=0.99 (T=100)
qgev(p=0.99,location=parMOM[1],scale=parMOM[2],shape=-parMOM[3])
qgev(p=0.99,location=parLMOM[1],scale=parLMOM[2],shape=-parLMOM[3])
qgev(p=0.99,location=parML[1],scale=parML[2],shape=parML[3])