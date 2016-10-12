data12.5 <-data.frame(
Y=c(3, 2, 1, 3, 3, 0, 2, 3, 4, 1, 4, 2, 2, 5, 3, 1, 2, 0, 1, 0, 
3, 2, 4, 4, 1, 0, 2, 1, 1, 4, 1, 2, 1, 2, 3, 3, 0, 1, 4),
NAO=c(-1.068, 1.598, -0.53, 0.372, -0.126, 1.892, -0.888, -1.218, 
-0.642, -0.54, -0.176, 1.27, -0.196, -2.068, -0.676, -0.36, -0.232, 
1.276, 0.342, 1.18, 0.64, -0.576, 0.338, -0.294, 0.448, 0.308, 
0.53, 1.942, 0.138, -0.536, -0.596, 0.952, -0.056, 1.994, 1.304, 
-0.132, 1.68, 2.044, 1.594))

glm1 <- glm( Y ~ NAO ,data=data12.5,family=poisson(link='log'))
coef(glm1) # regression parameters
vcov(glm1) # covariance matrix

glm0 <- glm( Y ~ 1, data=data12.5,family=poisson(link='log'))

D <- 2*(logLik(glm1)-logLik(glm0)) #likelihood ratio test statistic
qchisq(0.95,1) #critical region

AIC(glm1)
AIC(glm0)


































    if (sum(is.na(Tave)) != 0 & na.rm == FALSE) {
        stop("Data must not contain NAs")
    }
    if (is.ts(Tave) & frequency(Tave) != 12) {
        stop("Data should be a monthly time series (frequency = 12)")
    }
    if (!is.null(ncol(Tave))) {
        if (length(lat) != ncol(Tave)) {
            stop("Longitude of latitudes vector does not coincide with the number of columns in Tave")
        }
    }
    if (!is.ts(Tave)) {
        Tave <- ts(as.matrix(Tave), frequency = 12)
    }
    else {
        Tave <- ts(as.matrix(Tave), frequency = frequency(Tave), 
            start = start(Tave))
    }
    PE <- Tave * NA
    n <- nrow(Tave)
    m <- ncol(Tave)
    c <- cycle(Tave)
    tanLat <- tan(lat/57.2957795)
    tanDelta <- c(-0.37012566, -0.23853358, -0.04679872, 0.16321764, 
        0.32930908, 0.40677729, 0.3747741, 0.239063, 0.04044485, 
        -0.16905776, -0.33306377, -0.40743608)
    tanLatM <- matrix(tanLat, nrow = 12, ncol = length(tanLat), 
        byrow = TRUE) * tanDelta
    tanLatM[tanLatM < {
        -1
    }] <- -1
    tanLatM[tanLatM > 1] <- 1
    omega <- acos(-tanLatM)
    N <- 24 * omega/pi
    days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    K <- N/12 * days/30
    T <- matrix(unlist(by(Tave, c, colMeans, na.rm = na.rm)), 
        ncol = m, byrow = TRUE)
    T[T < 0] <- 0
    J <- colSums((T/5)^1.514, na.rm = na.rm)
    J <- matrix(J, n, m, TRUE)
    J2 <- J * J
    J3 <- J2 * J
    q <- 6.75e-07 * J3 - 7.71e-05 * J2 + 0.01792 * J + 0.49239
    Tave[Tave < 0] <- 0
    PE <- K[c, ] * 16 * (10 * Tave/J)^q
    colnames(PE) <- rep("PET_tho", m)
    return(PE)
}