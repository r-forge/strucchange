pargmaxV <- function(x, xi = 1, phi1 = 1, phi2 = 1)
{
  phi <- xi * (phi2/phi1)^2

  G1 <- function(x, xi = 1, phi = 1)
  {
    x <- abs(x)
    frac <- xi/phi
    rval <- - exp(log(x)/2 - x/8 - log(2*pi)/2) -
              (phi/xi * (phi + 2*xi)/(phi+xi)) * exp((frac * (1 + frac) * x/2) + pnorm(-(0.5 + frac) * sqrt(x), log.p = TRUE)) +
	      exp(log(x/2 - 2 + ((phi + 2 * xi)^2)/((phi + xi)*xi)) + pnorm(-sqrt(x)/2, log.p = TRUE))
    rval
  }

  G2 <- function(x, xi = 1, phi = 1)
  {
    frac <- xi^2/phi
    rval <- 1 + sqrt(frac) * exp(log(x)/2 - (frac*x)/8  - log(2*pi)/2) +
            (xi/phi * (2*phi + xi)/(phi + xi)) * exp(((phi + xi) * x/2) + pnorm(-(phi + xi/2)/sqrt(phi) * sqrt(x), log.p = TRUE)) -
	    exp(log(((2*phi + xi)^2)/((phi+xi)*phi) - 2 + frac*x/2) + pnorm(-sqrt(frac) * sqrt(x)/2 , log.p = TRUE))
    rval
  }

  ifelse(x < 0, G1(x, xi = xi, phi = phi), G2(x, xi = xi, phi = phi))
}

ci.bp <- function(obj, breaks = NULL, alpha = 0.95, het.reg = FALSE, het.err = FALSE)
{
  X <- obj$X
  y <- obj$y
  n <- obj$nobs

  myfun <- function(x, alpha = 0.975, xi = 1, phi1 = 1, phi2 = 1)
    (pargmaxV(x, xi = xi, phi1 = phi1, phi2 = phi2) - alpha)^2

  bp <- breakpoints(obj, breaks = breaks)$breakpoints
  nbp <- length(bp)
  upper <- rep(0, nbp)
  lower <- rep(0, nbp)
  bp <- c(0, bp, n)

  sigma1 <- sigma2 <- sum(lm.fit(X,y)$residuals^2)/n
  Q1 <- Q2 <- crossprod(X)/n
  xi <- 1

  X2 <- X[(bp[1]+1):bp[2],,drop = FALSE]
  y2 <- y[(bp[1]+1):bp[2]]
  fm <- lm.fit(X2, y2)
  beta2 <- fm$coefficients
  if(het.err) sigma2 <- sum(fm$residuals^2)/nrow(X2)
  if(het.reg) Q2 <- crossprod(X2)/nrow(X2)

  for(i in 2:(nbp+1))
  {
    X1 <- X2
    y1 <- y2
    beta1 <- beta2
    sigma1 <- sigma2
    Q1 <- Q2

    X2 <- X[(bp[i]+1):bp[i+1],,drop = FALSE]
    y2 <- y[(bp[i]+1):bp[i+1]]
    fm <- lm.fit(X2, y2)
    beta2 <- fm$coefficients
    delta <- beta2 - beta1
    frac <- as.vector(crossprod(delta, Q1) %*% delta)/sigma1

    if(het.err) sigma2 <- sum(fm$residuals^2)/nrow(X2)
    if(het.reg) {
      Q2 <- crossprod(X2)/nrow(X2)
      xi <- as.vector(crossprod(delta, Q2) %*% delta)/(sigma1*frac)
    }

    upper[i-1] <- optimize(myfun, c(0,n), alpha = (1-(1-alpha)/2), xi = xi, phi1 = sqrt(sigma1), phi2 = sqrt(sigma2))$minimum/frac
    lower[i-1] <- optimize(myfun, c(-n,0), alpha = (1-alpha)/2, xi = xi, phi1 = sqrt(sigma1), phi2 = sqrt(sigma2))$minimum/frac
    ## upper[i-1] <- optimize(myfun, c(0,n), alpha = (1-(1-alpha)/2), xi = xi, phi1 = sqrt(sigma1), phi2 = sqrt(sigma1))$minimum/frac
    ## lower[i-1] <- optimize(myfun, c(-n,0), alpha = (1-alpha)/2, xi = xi, phi1 = sqrt(sigma1), phi2 = sqrt(sigma1))$minimum/frac
    ## Bai & Perron claim: phi2 = sqrt(sigma1), I don't believe this!
  }
  bp <- bp[-c(1,nbp+2)]
  bp <- cbind(bp+floor(lower),bp,bp+ceiling(upper))
  colnames(bp) <- c("lower", "breakpoints", "upper")
  return(bp)
}



Gsimple <- function(x)
{
  G1 <- function(x) {
    1 + 1/sqrt(2*pi) * sqrt(x) * exp(-x/8) - 0.5 * (x + 5) * pnorm(-sqrt(x)/2) + 1.5 * exp(x) * pnorm(-1.5 * sqrt(x))
  }
  ifelse(x < 0, 1 - G1(-x), G1(x))
}

