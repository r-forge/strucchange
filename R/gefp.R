gefp <- function(...,
  fit = glm, scores = estfun, vcov = NULL,
  scale = TRUE, sandwich = TRUE,
  order.by = NULL, parm = NULL, data = list())
{
  fm <- fit(..., data = data)
  psi <- as.matrix(scores(fm))
  n <- nrow(psi)
  k <- ncol(psi)

  if(!is.null(order.by))
  {
    if("formula" %in% class(order.by)) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
    psi <- psi[index, , drop = FALSE]
    z <- z[index]
  } else {
    index <- 1:n
    if(is.ts(psi)) z <- time(psi)
      else z <- index/n
  }

  if("POSIXt" %in% class(z))
    z <- c(z[1] + as.numeric(difftime(z[1], z[2], units = "secs")), z)
  else
    z <- c(z[1] - diff(z[1:2]), z)

  process <- psi/sqrt(n)

  if(is.null(vcov))
    J12 <- root.matrix(crossprod(process))
  else {
    if(sandwich) {
      Q <- chol2inv(chol(summary(fm)$cov.unscaled))
      J12 <- root.matrix((Q %*% vcov(fm) %*% Q)/n)
    } else {
      J12 <- root.matrix(vcov(fm)/n)
    }
  }

  process <- rbind(0, process)
  process <- apply(process, 2, cumsum)

  if(scale) process <- t(chol2inv(chol(J12)) %*% t(process))
    else {
      process <- t(1/diag(J12) * t(process))
      if(length(parm) > 1) stop("limiting process is not a Brownian bridge")
    }

  if(is.ts(psi) & is.null(order.by))
    process <- ts(process, end = end(psi), frequency = frequency(psi))

  if(!is.ts(process))
    process <- ts(process, start = 0, frequency = n-1)

  colnames(process) <- colnames(psi)
  if(!is.null(parm)) process <- process[, parm]

  retval <- list(process = quasits(process, z),
                 nreg = k,
                 nobs = n,
                 call = match.call(),
		 fit = fit,
		 scores = scores,
                 par = NULL,
                 lim.process = "Brownian bridge",
		 type.name = "M-fluctuation test",
                 J12 = J12,
                 time = z)

  class(retval) <- c("gefp", "efp")
  return(retval)
}

plot.gefp <- function(x, alpha = 0.05, parm = NULL,
                     boundary = TRUE, functional = "max", 
		     main = NULL,  ylim = NULL, xlab = "Time",
		     ylab = "empirical fluctuation process",
		     ...)
{
    if(is.null(functional)) fun <- "max"
      else fun <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))
    bound <- boundary(x, alpha = alpha, functional = fun)
    pos <- FALSE
    ave <- FALSE

    if(is.null(main)){
            if(fun == "meanL2")
              main <- paste(x$type.name, "with mean L2 norm")
	    else if(fun == "maxL2")
	      main <- paste(x$type.name, "with max L2 norm")
	    else
	      main <- x$type.name
    }
    
    z <- x$process
    if(!is.null(parm)) z <- z[, parm]

    if(!is.null(functional)) {
      k <- ncol(z)

      switch(functional,
        "max" = {
          if(k > 1) {
            z <- apply(abs(z), 1, max)
            pos <- TRUE
          }
        },
        "range" = { stop("no plot available for range functional") },
        "maxL2" = {
	  if(x$lim.process == "Brownian bridge") {
            z <- rowSums(z^2)
	    pos <- TRUE
	  } else {
	    stop("no test/plot available for mean L2 functional")
	  }
        },

        "meanL2" = {
	  if(x$lim.process == "Brownian bridge") {
            z <- rowSums(z^2)
	    ave <- TRUE
	    pos <- TRUE
	  } else {
	    stop("no test/plot available for mean L2 functional")
	  }
        })
    }
    
    if(is.null(ylim)) {
      ymax <- max(c(z, bound))
      if(pos) ymin <- 0
      else ymin <- min(c(z, -bound))
      ylim <- c(ymin, ymax)
    }

    if(boundary)
        panel <- function(y, ...) {
            lines(y, ...)
            lines(bound, col=2)
            lines(-bound, col=2)
            abline(0,0)
        }
    else
        panel <- function(y, ...) {
            lines(y, ...)
            abline(0,0)
        }

    if(any(attr(z, "class") == "mts"))
        plot(z, main = main, panel = panel, ...)
    else {
        plot(x$time, z, main = main, xlab = xlab, ylab = ylab, ylim = ylim, type = "l", ...)
        if(boundary) {
            lines(x$time, bound, col=2)
            if(!pos) lines(x$time, -bound, col=2)
            if(ave) {
              avez <- ts(rep(mean(z), length(bound)), start = start(bound), frequency = frequency(bound))
              lines(x$time, avez, lty = 2)
            }
        }
        abline(0,0)
    }
}

efpFunctional <- function(functional = list(comp = function(x) max(abs(x)), time = max),
		     boundary = function(x) rep(1, x),
		     computePval = NULL,
		     computeCritval = NULL,
		     lim.process = "Brownian bridge",
		     nobs = 10000, nrep = 50000, nproc = 1:20, h = 0.5)
{		     
  if(is.list(functional)) {
    if(names(functional) == c("comp", "time")) {
      if(identical(functional[[2]], max)) {
        myfun <- function(x) {
	  rval <- apply(x, 1, functional[[1]])
	  functional[[2]](rval/boundary(1:length(rval)/rval))
	}
      } else
        myfun <- function(x) functional[[2]](apply(x, 1, functional[[1]]))
    }
    else if(names(functional) == c("time", "comp"))
      myfun <- function(x) functional[[2]](apply(x, 2, functional[[1]]))
    else  
      stop("`functional' should be a list with elements `comp' and `time'")
  } else {
    myfun <- functional
  }
  
  compStatistic <- function(x)
    myfun(x)
    
  if(is.null(computePval)) {
    if(is.null(nproc)) {
      z <- simulateDist(nobs = nobs, nrep = nrep, nproc = 1,
             h = h, lim.process = lim.process, functional = myfun)
      computePval <- function(x, nproc = 1) 1 - (sum(x <= z)/nrep)^nproc
      if(is.null(computeCritval))
        computeCritval <- function(alpha, nproc = 1) quantile(z, (1-alpha)^(1/nproc))
    } else {
      z <- matrix(rep(0, nrep * length(nproc)), ncol = length(nproc))
      colnames(z) <- as.character(nproc)
      for(i in nproc)
        z[, as.character(i)] <- simulateDist(nobs = nobs, nrep = nrep, nproc = i,
               h = h, lim.process = lim.process, functional = myfun)
      computePval <- function(x, nproc = 1) sum(x > z[, as.character(nproc)])/nrep
      if(is.null(computeCritval))
        computeCritval <- function(alpha, nproc = 1) quantile(z[, as.character(nproc)], 1-alpha)
    }
  }
  if(is.null(computeCritval)) {
    computeCritval <- function(alpha, nproc = 1)
      uniroot(function(y) {computePval(y, nproc = nproc) - alpha}, c(0,100))$root
  }


  plotProcess <- function(x, alpha = 0.05,
    xlab = "Time", ylab = "empirical fluctuation process", main = x$type.name, ...)
  {
    plot(x$time, x$process, type = "l", xlab = xlab, ylab = ylab, main = main, ...)
    abline(h = computeCritval(alpha), col = 2)
    abline(h = computeStatistic(x$process), lty = 2)
  }  

  
  rval <- list(plotProcess = plotProcess,
               computeStatistic = computeStatistic,
	       computePval = computePval,
	       computeCritval = computeCritval)
  class(rval) <- "efpFunctional"
  return(rval)
}

simulateDist <- function(nobs = 1000, nrep = 5000, nproc = 1,
                         lim.process = "Brownian bridge", 
			 h = 0.5, functional = max)
{
  lim.process <- match.arg(lim.process,
    c("Brownian motion", "Brownian motion increments", 
      "Brownian bridge", "Brownian bridge increments"))
  rval <- numeric(nrep)
  
  switch(lim.process,
  
  "Brownian motion" = {
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, cumsum)
      x <- rbind(0, x)/sqrt(nobs)
      rval[i] <- functional(x)
    }
  },
  
  "Brownian motion increments" = {
    nh <- floor(nobs * h)
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, cumsum)
      x <- rbind(0, x)/sqrt(nobs)
      x <- apply(x, 2, function(z) z[-(1:nh)] - z[1:(nobs-nh+1)])
      rval[i] <- functional(x)
    }
  },
  
  "Brownian bridge" = {
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, function(z) cumsum(z - mean(z)))
      x <- rbind(0, x)/sqrt(nobs)
      rval[i] <- functional(x)
    }
  },
  
  "Brownian bridge increments" = {
    nh <- floor(nobs * h)
    for(i in 1:nrep) {
      x <- matrix(rnorm(nproc * nobs), ncol = nproc)
      x <- apply(x, 2, function(z) cumsum(z - mean(z)))
      x <- rbind(0, x)/sqrt(nobs)
      x <- apply(x, 2, function(z) z[-(1:nh)] - z[1:(nobs-nh+1)])
      rval[i] <- functional(x)
    }
  })
  
  return(rval)
}



quasits <- function(x, order.by)
{
  index <- order(order.by)
  order.by <- order.by[index]
  
  if(is.vector(x))
    x <- x[index]
  else if(is.matrix(x))
    x <- x[index, , drop = FALSE]
  else
    stop("`x' has to be a vector a matrix")

  attr(x, "time") <- order.by
  class(x) <- "quasits"
  return(x)
}

time.quasits <- function(x, ...)
{
  attr(x, "time")
}

plot.quasits <- function(x, xlab = "Time",
  ylab = deparse(substitute(x)), type = "l", ...)
{
  plot(time(x), x, xlab = xlab, ylab = ylab, type = type, ...)
}
