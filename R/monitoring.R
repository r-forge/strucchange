mefp <- function(obj, ...) UseMethod("mefp")

mefp.formula <-
    function(obj, type = c("ME", "fluctuation"), h=0.25, data=list(),
             alpha=0.05, functional = c("max", "range"),
             period=10, tolerance=.Machine$double.eps^0.5,
             MECritvalTable=monitorMECritvalTable, rescale=FALSE)
{
    type <- match.arg(type)
    functional <- match.arg(functional)
    val <- efp(obj, type=type, h=h, data=data, rescale=rescale)
    val <- mefp(val, alpha=alpha, functional=functional, period=period,
         tolerance=tolerance, MECritvalTable=MECritvalTable,
         rescale=rescale)
    if(length(data) == 0)
        val$data <- NULL
    else
        val$data <- deparse(substitute(data))
    return(val)
}

mefp.efp <-
    function(obj, alpha=0.05, functional = c("max", "range"),
             period=10, tolerance=.Machine$double.eps^0.5,
             MECritvalTable=monitorMECritvalTable, rescale=FALSE)
{
    functional <- match.arg(functional)
    if(! (obj$type %in% c("ME", "fluctuation")))
        stop("efp must be of type `fluctuation' or `ME'")

    if(is.null(as.list(obj$efpcall)$data)){
       data <- NULL
    }
    else{
       data <- as.character(as.list(obj$efpcall)$data)
    }

    ## Bonferroni correction
    elemsiglevel <- alpha / obj$nreg
    histcoef <- obj$coefficients
    histsize <- obj$nobs
    switch(obj$type,
           "ME" = { winsize <- obj$par },
       { winsize <- 1 })
    K <- floor(winsize*obj$nobs)
    sigmahat <- obj$sigma
    Q12s <- obj$Q12 * K/(sigmahat*sqrt(histsize))

    computeEmpProc <- function(newcoef, Q){
        if(is.null(Q)) Q <- Q12s
        else Q <- Q * K/(sigmahat*sqrt(histsize))
        t(Q %*%(newcoef-histcoef))
    }

    logPlus <- function(x) ifelse(x<=exp(1),1,log(x))

    if(obj$type=="ME"){
        dntab <- dimnames(MECritvalTable)
        if(!(winsize %in% dntab[[1]]))
            stop(paste("winsize h =",winsize,"not available, we have:",
                       paste(dntab[[1]], collapse=", ")))
        if(!(period %in% dntab[[2]]))
            stop(paste("period",period,"not available, we have:",
                       paste(dntab[[2]], collapse=", ")))
        critval <- approx(x=as.numeric(dntab[[3]]),
                          y=MECritvalTable[as.character(winsize),
                          as.character(period),,functional],
                          xout=1-elemsiglevel)$y
        if(is.na(critval))
            stop(paste("Necessary significance level per parameter of",
                       elemsiglevel,
                       "\n\toutside of available range",
                       paste(range(1-as.numeric(dntab[[3]])),
                             collapse="-")))

        computeEstims <- function(x, y, k){
            ok <- (k-K+1):k
            retval <- list(coef=NULL, Qr12=NULL)
            retval$coef <- coef(lm.fit(x[ok,,drop=FALSE], y[ok,,drop=FALSE]))
            if(rescale)
                retval$Qr12 <- root.matrix(crossprod(x[ok,,drop=FALSE]))/sqrt(K)
            retval
        }
        border <- function(k){
            critval*sqrt(2*logPlus(k/histsize))
        }
    }
    else{

        mreSize <- function(a){
            -2*(pnorm(a)-a*dnorm(a))
        }
        mreCritval <- function(a){
            abs(2*(pnorm(a)-a*dnorm(a))+elemsiglevel-2)
        }
        critval <- optim(5, mreCritval)$par
        if((mreSize(critval)-elemsiglevel) > tolerance)
            stop("Could not find critical within tolerance")


        computeEstims <- function(x, y, k){
            retval <- list(coef=NULL, Qr12=NULL)
            retval$coef <- coef(lm.fit(x[1:k,,drop=FALSE], y[1:k,,drop=FALSE]))
            if(rescale)
                retval$Qr12 <- root.matrix(crossprod(x[1:k,,drop=FALSE]))/sqrt(k)
            retval
        }
        border <- function(k){
            x <- k/histsize
            sqrt(x*(x-1)*(critval^2 + log(x/(x-1))))
        }
    }

    if(functional=="max"){
        computeStat <- function(empproc){
            max(abs(empproc))
        }
    }
    else if(functional=="range"){
        if(obj$type=="fluctuation")
            stop("Functional `range' not available for recursive estimates")
        else{
            computeStat <- function(empproc){
                max(apply(empproc, 2, function(x) diff(range(x))))
            }
        }
    }



    obj <- list(breakpoint=NA, last=obj$nobs, process=NULL,
                statistic=NULL, histsize=histsize,
                initcall=match.call(), call=match.call(),
                efpcall=obj$call, efpprocess=obj$process,
                computeEstims=computeEstims,
                computeEmpProc=computeEmpProc,
                border=border, computeStat=computeStat,
                functional=functional, alpha=alpha, critval=critval,
                histcoef=histcoef, formula=obj$formula,
                type.name=paste("Monitoring with", obj$type.name),
                data=data)

    class(obj) <- "mefp"
    obj
}

monitor <- function(obj, data=NULL, verbose=TRUE){

    if(!is.na(obj$breakpoint)) return(TRUE)
    if(missing(data)){
        if(is.null(obj$data)){
            data <- list()
        }
        else{
            data <- get(obj$data)
        }
    }

    mf <- model.frame(obj$formula, data=data)
    y <- as.matrix(model.response(mf))
    x <- model.matrix(obj$formula, data = data)

    if(nrow(x) <= obj$last) return(obj)
    if(nrow(x)!=nrow(y))
        stop("response and regressors must have the same number of rows")
    if(ncol(y)!=1)
        stop("multivariate response not implemented yet")
    foundBreak <- FALSE
    for(k in (obj$last+1):nrow(x)){
        newestims <- obj$computeEstims(x,y,k)
        obj$process <- rbind(obj$process,
                             obj$computeEmpProc(newestims$coef, newestims$Qr12))
        stat <- obj$computeStat(obj$process)
        obj$statistic <- c(obj$statistic, stat)
        if(!foundBreak & (stat > obj$border(k))){
            foundBreak <- TRUE
            obj$breakpoint <- k
            if(verbose) cat("Break detected at observation #", k, "\n")
        }
    }
    obj$last <- k
    obj$lastcoef <- newestims$coef
    obj$call <- match.call()
    obj
}

print.mefp <- function(obj){

    cat(obj$type.name, "\n\n")
    cat("Initial call:\n ", deparse(obj$initcall), "\n\n")
    cat("Last call:\n ", deparse(obj$call), "\n\n")
    cat("Significance level   : ", obj$alpha, "\n")
    cat("Critical value       : ", obj$critval, "\n")
    cat("History size         : ", obj$histsize, "\n")
    cat("Last point evaluated : ", obj$last, "\n")
    if(!is.na(obj$breakpoint))
        cat("Structural break at  : ", obj$breakpoint, "\n")
    cat("\nParameter estimate on history :\n");
    print(obj$histcoef)
    if(!is.null(obj$lastcoef)){
        cat("Last parameter estimate :\n");
        print(obj$lastcoef)
    }
}

plot.mefp <- function(obj, boundary=TRUE, functional="max", main=NULL,
                      ylab="empirical fluctuation process", ylim=NULL, ...){

    if(obj$last>obj$histsize){
        proc <- rbind(as.matrix(obj$efpprocess),
                    as.matrix(obj$process))
        proc <- ts(proc, start=start(obj$efpprocess),
                 frequency=frequency(obj$efpprocess))
        bound <- ts(obj$border((obj$histsize+1):obj$last),
                 end = end(proc), frequency=frequency(proc))
        pos <- FALSE
        if(!is.null(functional) && (functional == "max"))
        {
            proc <- ts(apply(abs(proc), 1, 'max'),
                       start = start(proc), frequency = frequency(proc))
            pos <- TRUE
        }
        ymax <- max(c(proc, bound))
        if(pos)
            ymin <- 0
        else
            ymin <- min(c(proc, -bound))
        if(is.null(ylim)) ylim <- c(ymin, ymax)
        if(is.null(main))
            main <- obj$type.name
        if(boundary)
            panel <- function(y, ...)
            {
                lines(y, ...)
                lines(bound, col=2)
                lines(-bound, col=2)
                abline(0,0)
            }
        else
            panel <- function(y, ...)
            {
                lines(y, ...)
                abline(0,0)
            }
        if(any(attr(proc, "class") == "mts"))
            plot(proc, main = main, panel = panel, ...)
        else
        {
            plot(proc, main = main, ylab = ylab, ylim = ylim, ...)
            if(boundary)
            {
                lines(bound, col=2)
                if(!pos) lines(-bound, col=2)
            }
            abline(0,0)
        }
        abline(v=max(time(obj$efpprocess)), lty=2)
    }
    else{
        cat("Nothing monitored yet!\n")
    }
}

boundary.mefp <- function(x)
{
        ts(x$border((x$histsize+1):x$last),
           start = max(time(x$efpprocess))+deltat(x$efpprocess),
           frequency=frequency(x$efpprocess))
}




