breakpoints <- function(obj, ...)
{
  UseMethod("breakpoints")
}

breakpoints.Fstats <- function(obj, ...)
{
  RVAL <- list(breakpoints = obj$breakpoint,
               nobs = obj$nobs,
	       nreg = obj$nreg,
	       call = match.call(),
               datatsp = obj$datatsp)
  class(RVAL) <- "breakpoints"
  return(RVAL)
}

breakpoints.formula <- function(formula, h = 0.15, breaks = NULL,
                                tol = 1e-15, data = list(), ...)
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
  if(is.null(breaks)) breaks <- floor(n/h)

  ## compute ith row of the SSR diagonal matrix, i.e,
  ## the recursive residuals for segments starting at i = 1:(n-h+1)

  SSRi <- function(i)
  {
    ssr <- recresid(X[i:n,,drop = FALSE],y[i:n], tol = tol)
    c(rep(NA, k), cumsum(ssr^2))
  }
  SSR.diag <- sapply(1:(n-h+1), SSRi)

  ## function to extract the SSR(i,j) from SSR.diag

  SSR <- function(i,j) SSR.diag[[i]][j - i + 1]

  ## compute optimal previous partner if observation i is the mth break
  ## store results together with SSRs in SSR.table

  ## breaks = 1

  index <- h:(n-h)
  break.SSR <- sapply(index, function(i) SSR(1,i))

  SSR.table <- cbind(index, break.SSR)
  rownames(SSR.table) <- as.character(index)

  ## breaks >= 2

  extend.SSR.table <- function(SSR.table, breaks)
  {
    if((breaks*2) > ncol(SSR.table)) {
      for(m in (ncol(SSR.table)/2 + 1):breaks)
      {
        my.index <- (m*h):(n-h)
        my.SSR.table <- SSR.table[,c((m-1)*2 - 1, (m-1)*2)]
        my.SSR.table <- cbind(my.SSR.table, NA, NA)
        for(i in my.index)
        {
          pot.index <- ((m-1)*h):(i - h)
          break.SSR <- sapply(pot.index, function(j) my.SSR.table[as.character(j), 2] + SSR(j+1,i))
          opt <- which.min(break.SSR)
          my.SSR.table[as.character(i), 3:4] <- c(pot.index[opt], break.SSR[opt])
        }
        SSR.table <- cbind(SSR.table, my.SSR.table[,3:4])
      }
      colnames(SSR.table) <- as.vector(rbind(paste("break", 1:breaks, sep = ""),
                                             paste("SSR", 1:breaks, sep = "")))
    }
    return(SSR.table)
  }

  SSR.table <- extend.SSR.table(SSR.table, breaks)

  ## extract optimal breaks

  extract.breaks <- function(SSR.table, breaks)
  {
    if((breaks*2) > ncol(SSR.table)) stop("compute SSR.table with enough breaks before")
    index <- SSR.table[, 1, drop = TRUE]
    break.SSR <- sapply(index, function(i) SSR.table[as.character(i),breaks*2] + SSR(i + 1, n))
    opt <- index[which.min(break.SSR)]
    if(breaks > 1) {
      for(i in ((breaks:2)*2 - 1))
        opt <- c(SSR.table[as.character(opt[1]),i], opt)
    }
    names(opt) <- NULL
    return(opt)
  }

  opt <- extract.breaks(SSR.table, breaks)

  if(is.ts(data))
      datatsp <- tsp(data)
  else if(is.ts(y))
      datatsp <- tsp(y)
  else
      datatsp <- c(0, 1, n)

  RVAL <- list(breakpoints = opt,
               SSR.table = SSR.table,
	       SSR.diag = SSR.diag,
	       SSR = SSR,
	       extract.breaks = extract.breaks,
	       extend.SSR.table = extend.SSR.table,
	       nobs = n,
	       nreg = k,
	       call = match.call(),
	       datatsp = datatsp)
  class(RVAL) <- c("breakpointsfull", "breakpoints")
  return(RVAL)
}

breakpoints.breakpointsfull <- function(obj, breaks = 1, ...)
{
  SSR.tab <- obj$extend.SSR.table(obj$SSR.table, breaks)
  breakpoints <- obj$extract.breaks(SSR.tab, breaks)
  RVAL <- list(breakpoints = breakpoints,
               nobs = obj$nobs,
	       nreg = obj$nreg,
	       call = match.call(),
               datatsp = obj$datatsp)
  class(RVAL) <- "breakpoints"
  return(RVAL)
}

print.breakpoints <- function(x, format.times = FALSE, ...)
{
  cat(paste("\t Optimal ", length(x$breakpoints) + 1,
            "-segment partition: \n\n", sep = ""))
  cat("Call:\n")
  print(x$call)
  cat("\nBreakpoints:\n")
  print(x$breakpoints)
  cat("\nBreakdates:\n")
  print(breakdates(x, format.times = format.times))
}

breakdates <- function(obj, format.times = FALSE)
{
  breakdates <- seq(from = obj$datatsp[1], by = 1/obj$datatsp[3], length = obj$nobs)[obj$breakpoints]

  format.time <- function(timevec, freq)
  {
    first <- floor(timevec)
    second <- round((timevec - first)*freq + 1, digits = 0)
    RVAL <- cbind(first, second)
    dummy <- function(x) paste(x[1], "(", x[2], ")", sep = "")
    RVAL <- apply(RVAL, 1, dummy)
    return(RVAL)
  }

  if(format.times) breakdates <- format.time(breakdates, obj$datatsp[3])
  return(breakdates)
}

breakfactor <- function(obj, labels = NULL, ...)
{
  breaks <- obj$breakpoints
  nbreaks <- length(breaks)
  fac <- rep(1:(nbreaks + 1), c(breaks[1], diff(c(breaks, obj$nobs))))
  if(is.null(labels)) labels <- paste("segment", 1:(nbreaks+1), sep = "")
  fac <- factor(fac, labels = labels, ...)
  return(fac)
}

lines.breakpoints <- function(x, lty = 2, ...)
{
  abline(v = breakdates(x), lty = lty, ...)
}

