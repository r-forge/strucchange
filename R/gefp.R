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
    if(inherits(order.by, "formula")) {
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
      else if(is.quasits(psi)) z <- time(psi)
      else if(is.ts(data)) z <- time(data)
      else if(is.quasits(data)) z <- time(data)
      else z <- index/n
  }

  if(inherits(z, "POSIXt"))
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
                 J12 = J12)

  class(retval) <- c("gefp", "efp")
  return(retval)
}

time.gefp <- function(x, ...)
{
  time(x$process, ...)
}

efpFunctional <- function(functional = list(comp = function(x) max(abs(x)), time = max),
		     boundary = function(x) rep(1, length(x)),
		     computePval = NULL,
		     computeCritval = NULL,
		     lim.process = "Brownian bridge",
		     nobs = 10000, nrep = 50000, nproc = 1:20, h = 0.5,
		     probs = c(0:84/100, 850:1000/1000))
{		     
  probs <- probs[-which(probs == 0)]

  ## compute from functional list the full functional
  ## lambda = myfun
  
  if(is.list(functional)) {
    if(identical(names(functional), c("comp", "time"))) {
      if(identical(functional[[2]], max)) {
        myfun <- function(x) {
	  rval <- apply(as.matrix(x), 1, functional[[1]])
	  functional[[2]](rval/boundary(0:(length(rval)-1)/(length(rval)-1)))
	}
      } else
        myfun <- function(x) functional[[2]](apply(as.matrix(x), 1, functional[[1]]))
    }
    else if(identical(names(functional), c("time", "comp")))
      myfun <- function(x) functional[[2]](apply(as.matrix(x), 2, functional[[1]]))
    else  
      stop("`functional' should be a list with elements `comp' and `time'")
  } else {
    myfun <- functional
  }

  ## setup computeStatistic function
  computeStatistic <- function(x)
    myfun(as.matrix(x)[-1,])

  ## if missing simulate values for
  ## computePval and computeCritval
  
  if(is.null(computePval) & is.null(computeCritval)) {
    if(is.null(nproc)) {
      z <- simulateDist(nobs = nobs, nrep = nrep, nproc = 1,
             h = h, lim.process = lim.process, functional = myfun)
      
      zquant <- c(0, quantile(z, probs = probs))
      rm(z)

      pfun <- approxfun(zquant, 1 - c(0, probs))
      computePval <- function(x, nproc = 1) {
        1 - (1 - ifelse(x > max(zquant), 0, pfun(x)))^nproc
      }
      
      cfun <- approxfun(c(0, probs), zquant)
      computeCritval <- function(alpha, nproc = 1) cfun((1-alpha)^(1/nproc))

    } else {
      z <- matrix(rep(0, nrep * length(nproc)), ncol = length(nproc))
      colnames(z) <- as.character(nproc)
      for(i in nproc)
        z[, as.character(i)] <- simulateDist(nobs = nobs, nrep = nrep, nproc = i,
               h = h, lim.process = lim.process, functional = myfun)

      zquant <- rbind(0, apply(z, 2, function(x) quantile(x, probs = probs)))
      rm(z)
      computePval <- function(x, nproc = 1) {
        if(as.character(nproc) %in% colnames(zquant)) {
          pfun <- approxfun(zquant[, as.character(nproc)], 1 - c(0, probs))
          ifelse(x > max(zquant[, as.character(nproc)]), 0, pfun(x))
	}
	else stop("insufficient simulated values: cannot compute p value")
      }
      computeCritval <- function(alpha, nproc = 1) {
        if(as.character(nproc) %in% colnames(zquant)) {
          cfun <- approxfun(c(0, probs), zquant[, as.character(nproc)])
          cfun(1 - alpha)
        } else stop("insufficient simulated values: cannot compute critical value")
      }
    }
  }

  if(is.null(computeCritval)) {
    computeCritval <- function(alpha, nproc = 1)
      uniroot(function(y) {computePval(y, nproc = nproc) - alpha}, c(0, 1000))$root
  }

  if(is.null(computePval)) {
    computePval <- function(x, nproc = 1)
      uniroot(function(y) {computeCritval(y, nproc = nproc) - x}, c(0, 1))$root
  }


  ## define sensible default plotting method

  if(is.list(functional)) {

    if(identical(names(functional), c("comp", "time"))) {

    ## lambda = lambda_time(lambda_comp(x))
    ## aggregate first over components then over time

      if(identical(functional[[2]], max)) {

    ## special case: lambda = max(lambda_comp(x))
    ## can also use boundary argument: b(t) = critval * boundary(t)
    
        plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
	  xlab = "Time", ylab = NULL, main = x$type.name, ylim = NULL, ...)
	{
          n <- x$nobs
	  bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary(0:n/n)
	  bound <- quasits(bound, time(x))

	  if(aggregate) {
	    proc <- quasits(apply(as.matrix(x$process), 1, functional[[1]]), time(x))
	    
	    if(is.null(ylab)) ylab <- "empirical fluctuation process"
	    if(is.null(ylim)) ylim <- range(c(range(proc), range(bound)))
	    
	    plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
	    abline(0, 0)
	    lines(bound, col = 2)	    
	  } else {
	    panel <- function(x, ...)
	    {
              lines(x, ...)
	      abline(0, 0)
	      if(paste(deparse(functional[[1]]), collapse = "") == "function (x) max(abs(x))") {
	        lines(bound, col = 2)
		lines(-bound, col = 2)
	      }	      
	    }
	    plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
	  }
	}

      } else {

    ## nothing specific known about lambda_time
    ## plot: first aggregate, add critval and statistic

        plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
	  xlab = "Time", ylab = NULL, main = x$type.name, ylim = NULL, ...)
	{
          n <- x$nobs
	  bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary(0:n/n)
	  bound <- quasits(bound, time(x))
          stat <- computeStatistic(x$process)
	  stat <- quasits(rep(stat, length(time(x))), time(x))

	  if(aggregate) {
	    proc <- quasits(apply(as.matrix(x$process), 1, functional[[1]]), time(x))
	    
	    if(is.null(ylab)) ylab <- "empirical fluctuation process"
	    if(is.null(ylim)) ylim <- range(c(range(proc), range(bound), range(stat)))
	    
	    plot(proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
	    abline(0, 0)
	    lines(bound, col = 2)
	    lines(stat, lty = 2)	    
	  } else {
	    panel <- function(x, ...)
	    {
              lines(x, ...)
	      abline(0, 0)
	    }
	    plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
	  }
	}
      }
    }
    
    else if(identical(names(functional), c("time", "comp"))) {

    ## lambda = lambda_comp(lambda_time(x))

      plotProcess <- function(x, alpha = 0.05, aggregate = TRUE,
	  xlab = "Component", ylab = "Statistic", main = x$type.name, ylim = NULL, ...)
      {
        k <- NCOL(x$process)
        bound <- computeCritval(alpha = alpha, nproc = NCOL(x$process)) * boundary(1:k/k)
        stat <- rep(computeStatistic(x$process), k)

        if(aggregate) {
	  proc <- apply(as.matrix(x$process), 2, functional[[1]])
	  
	  xlabels <- colnames(x$process)
	  if(is.null(xlabels)) xlabels <- paste("Series", 1:k)
          if(is.null(ylim)) ylim <- range(c(range(proc), range(bound), range(stat), 0))
	    
	  plot(1:k, proc, xlab = xlab, ylab = ylab, main = main, ylim = ylim, axes = FALSE, type = "h", ...)
	  points(1:k, proc)
	  box()
	  axis(2)
	  axis(1, at = 1:k, labels = xlabels)
	  abline(0, 0)
	  lines(bound, col = 2)
	  if(!identical(functional[[2]], max)) lines(stat, lty = 2)	    
	} else {
	  panel <- function(x, ...)
	  {
            lines(x, ...)
            abline(0, 0)
	  }
          plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
	}      
      }

    }
    
  } else {

    ## lambda = lambda(x)
    ## functional is already the full functional lambda
    ## for plotting: just plot raw process
    plotProcess <- function(x, alpha = 0.05, aggregate = FALSE,
      xlab = "Time", ylab = NULL, main = x$type.name, ...)
    {
      if(aggregate) warning("aggregation not available")
      panel <- function(x, ...) {
        lines(x, ...)
	abline(0, 0)
      }
      plot(x$process, xlab = xlab, ylab = ylab, main = main, panel = panel, ...)
    }  
  }

  
  rval <- list(plotProcess = plotProcess,
               computeStatistic = computeStatistic,
	       computePval = computePval,
	       computeCritval = computeCritval,
	       boundary = boundary)
	       
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

as.quasits <- function(x, ...)
{
  UseMethod("as.quasits")
}

as.quasits.default <- function(x, ...)
{
  rval <- as.vector(x)
  dim(rval) <- dim(x)
  quasits(rval, 1:ifelse(is.null(dim(x)), length(x), dim(x)[1]))
}
  
as.quasits.ts <- function(x, ...)
{
  rval <- as.vector(x)
  dim(rval) <- dim(x)
  quasits(rval, time(x))
}  

as.quasits.irts <- function(x, ...)
{
  quasits(x$value, x$time)
}

plot.quasits <- function(x,
  plot.type = c("multiple", "single"), panel = lines,
  xlab = "Time", ylab = NULL, main = NULL, ylim = NULL,
  oma = c(6, 0, 5, 0), col = 1, lty = 1, nc, ...)
{
  plot.type <- match.arg(plot.type)
  nser <- NCOL(x)
  x.time <- time(x)
  if(is.ts(x.time)) x.time <- as.vector(x.time)

  if(plot.type == "multiple" && nser > 1) {
    if(is.null(main)) main <- deparse(substitute(x))
    if(is.null(ylab)) ylab <- colnames(x)
    if(is.null(ylab)) ylab <- paste("Series", 1:nser)
    ylab <- rep(ylab, length.out = nser)
    col <- rep(col, length.out = nser)
    lty <- rep(lty, length.out = nser)

    panel <- match.fun(panel)
    if(nser > 10) stop("Can't plot more than 10 series")
    if(missing(nc)) nc <- if(nser >  4) 2 else 1
    oldpar <- par("mar", "oma", "mfcol")
    on.exit(par(oldpar))
    par(mar = c(0, 5.1, 0, 2.1), oma = oma)
    nr <- ceiling(nser / nc)
    par(mfcol = c(nr, nc))
    for(i in 1:nser) {
      if(i%%nr==0 || i == nser)
        plot(x.time, x[, i], xlab= "", ylab= ylab[i], type = "n", ...)
      else {      
        plot(x.time, x[, i], axes = FALSE, xlab= "", ylab= ylab[i], type = "n", ...)
        box()
        axis(2, xpd = NA)
      }
      panel(x.time, x[, i], col = col[i], lty = lty[i], ...)
    }
    par(oldpar)
  } else {
    if(is.null(ylab)) ylab <- deparse(substitute(x))
    if(is.null(main)) main <- ""
    if(is.null(ylim)) ylim <- range(x)
    
    col <- rep(col, length.out = nser)
    dummy <- rep(range(x), length.out = length(time(x)))
	    
    plot(x.time, dummy, xlab= xlab, ylab= ylab[1], type = "n", ylim = ylim, ...)
    box()
    y <- as.matrix(x)
    for(i in 1:nser) {
      panel(x.time, y[, i], col = col[i], lty = lty[i], ...)
    }
  }
  title(main)
  return(invisible(x))
}

lines.quasits <- function(x, type = "l", ...)
{
  if(NCOL(x) == 1) lines(time(x), x, type = type, ...)
    else stop("Can't plot multivariate quasi time series object")
}

"[.quasits" <- function(x, i, j, ...)
{
  if(!is.quasits(x))
    stop("method is only for quasits objects")
  x.time <- time(x)
  attr(x, "time") <- NULL
  nclass <- class(x)[-(1:which(class(x) == "quasits"))]
  if(length(nclass) < 1) nclass <- NULL 
  class(x) <- nclass
  if(NCOL(x) < 2) x <- as.matrix(x)
  if(missing(i)) i <- 1:nrow(x)
  if(missing(j)) j <- 1:ncol(x)
  return(quasits(x[i, j, ...], x.time[i]))
}

print.quasits <- function(x, ...)
{
  if(!is.quasits(x))
    stop("method is only for quasits objects")
  x.time <- time(x)
  attr(x, "time") <- NULL
  nclass <- class(x)[-(1:which(class(x) == "quasits"))]
  if(length(nclass) < 1) nclass <- NULL 
  class(x) <- nclass
  cat("Value:\n")
  print(x)
  cat("\nTime:\n")
  print(x.time)
}

is.quasits <- function(object)
  inherits(object, "quasits")
  



plot.gefp <- function(x, alpha = 0.05, functional = efpMax, ...)
{
  ## if(is.character(functional)) {
  ##   functional <- match.arg(functional, c("max", "range", "maxL2", "meanL2"))
  ##   functional <- switch(functional,
  ##     "max" = efpMax,
  ##     "range" = efpRange,
  ##     "maxL2" = efpMaxL2,
  ##     "meanL2" = efpMeanL2)
  ## }
  functional$plotProcess(x, alpha = alpha, ...)
}

sctest.gefp <- function(x, functional = efpMax)
{
  stat <- functional$computeStatistic(x$process)
  names(stat) <- "f(efp)"
  rval <- list(statistic = stat,
               p.value = functional$computePval(stat, NCOL(x$process)),
	       method = x$type.name,
	       data.name = deparse(substitute(x)))
  class(rval) <- "htest"
  return(rval)
}


gbreakpoints <- function(formula, order.by = NULL, h = 0.15, breaks = NULL,
  objective = function(x, y) sum(lm.fit(x,y)$residuals^2), data = list(), ...)
{
  mf <- model.frame(formula, data = data)
  y <- model.response(mf)
  modelterms <- terms(formula, data = data)
  X <- model.matrix(modelterms, data = data)

  n <- nrow(X)
  k <- ncol(X)
  if(is.null(h)) h <- k + 1
  if(h < 1) h <- floor(n*h)
  if(h <= k)
    stop("minimum segment size must be greater than the number of regressors")
  if(is.null(breaks)) breaks <- ceiling(n/h) - 2

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
    X <- X[index, , drop = FALSE]
    y <- y[index]
    z <- z[index]
  } else {
    index <- 1:n
    if(is.ts(y)) z <- time(y)
      else if(is.quasits(y)) z <- time(y)
      else if(is.ts(data)) z <- time(data)
      else if(is.quasits(data)) z <- time(data)
      else z <- index/n
  }

  ## compute ith row of the RSS diagonal matrix, i.e,
  ## the recursive residuals for segments starting at i = 1:(n-h+1)

  RSSi <- function(i)
  {
    rval <- rep(NA, n - i + 1)
    for(j in (k+1):(n-i+1)) {
      rval[j] <- objective(X[i:(i+j-1),,drop = FALSE], y[i:(i+j-1)])
    }
    return(rval)
  }
  RSS.triang <- sapply(1:(n-h+1), RSSi)

  ## function to extract the RSS(i,j) from RSS.triang

  RSS <- function(i,j) RSS.triang[[i]][j - i + 1]

  ## compute optimal previous partner if observation i is the mth break
  ## store results together with RSSs in RSS.table

  ## breaks = 1

  index <- h:(n-h)
  break.RSS <- sapply(index, function(i) RSS(1,i))

  RSS.table <- cbind(index, break.RSS)
  rownames(RSS.table) <- as.character(index)

  ## breaks >= 2

  extend.RSS.table <- function(RSS.table, breaks)
  {
    if((breaks*2) > ncol(RSS.table)) {
      for(m in (ncol(RSS.table)/2 + 1):breaks)
      {
        my.index <- (m*h):(n-h)
        my.RSS.table <- RSS.table[,c((m-1)*2 - 1, (m-1)*2)]
        my.RSS.table <- cbind(my.RSS.table, NA, NA)
        for(i in my.index)
        {
          pot.index <- ((m-1)*h):(i - h)
          break.RSS <- sapply(pot.index, function(j) my.RSS.table[as.character(j), 2] + RSS(j+1,i))
          opt <- which.min(break.RSS)
          my.RSS.table[as.character(i), 3:4] <- c(pot.index[opt], break.RSS[opt])
        }
        RSS.table <- cbind(RSS.table, my.RSS.table[,3:4])
      }
      colnames(RSS.table) <- as.vector(rbind(paste("break", 1:breaks, sep = ""),
                                             paste("RSS", 1:breaks, sep = "")))
    }
    return(RSS.table)
  }

  RSS.table <- extend.RSS.table(RSS.table, breaks)

  ## extract optimal breaks

  extract.breaks <- function(RSS.table, breaks)
  {
    if((breaks*2) > ncol(RSS.table)) stop("compute RSS.table with enough breaks before")
    index <- RSS.table[, 1, drop = TRUE]
    break.RSS <- sapply(index, function(i) RSS.table[as.character(i),breaks*2] + RSS(i + 1, n))
    opt <- index[which.min(break.RSS)]
    if(breaks > 1) {
      for(i in ((breaks:2)*2 - 1))
        opt <- c(RSS.table[as.character(opt[1]),i], opt)
    }
    names(opt) <- NULL
    return(opt)
  }

  opt <- extract.breaks(RSS.table, breaks)

  if(is.ts(data))
      datatsp <- tsp(data)
  else if(is.ts(y))
      datatsp <- tsp(y)
  else
      datatsp <- c(0, 1, n)

  RVAL <- list(breakpoints = opt,
               RSS.table = RSS.table,
	       RSS.triang = RSS.triang,
	       RSS = RSS,
	       extract.breaks = extract.breaks,
	       extend.RSS.table = extend.RSS.table,
	       nobs = n,
	       nreg = k, y = y, X = X,
	       call = match.call(),
	       datatsp = datatsp)
  class(RVAL) <- c("breakpointsfull", "breakpoints")
  RVAL$breakpoints <- breakpoints(RVAL)$breakpoints
  return(RVAL)
}














## deprecated

boundary.gefp <- function(x, ...)
{
  quasits(as.vector(boundary.efp(x, ...)), time(x))
}

##plot.gefp <- function(x, alpha = 0.05, parm = NULL,
plotGEFP <- function(x, alpha = 0.05, parm = NULL,
                     boundary = TRUE, functional = "max", 
		     main = NULL,  ylim = NULL, xlab = "Time", ylab = NULL,
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
      k <- NCOL(z)

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
        mpanel <- function(y, ...) {
            lines(y, ...)
            lines(bound, col=2)
            lines(-bound, col=2)
            abline(0,0)
        }
    else
        mpanel <- function(y, ...) {
            lines(y, ...)
            abline(0,0)
        }
     if(!is.quasits(z)) z <- quasits(z, time(x))
     if(is.null(ylab) & NCOL(z) < 2) ylab <- "empirical fluctuation process"
     plot(z, main = main, xlab = xlab, ylab = ylab, ylim = ylim, panel = mpanel, ...)
}

