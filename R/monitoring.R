monitor <- function(obj, ...) UseMethod("monitor")

monitor.efp <-
    function(obj, alpha=0.05, functional = c("max", "range"),
             period=10, tolerance=.Machine$double.eps^0.5,
             MECritvalTable=monitorMECritvalTable)
{
    functional <- match.arg(functional)
    if(! (obj$type %in% c("ME", "fluctuation")))
        stop("efp must be of type `fluctuation' or `ME'")
    
    ## Bonferroni correction
    elemsiglevel <- alpha / obj$nreg
    histcoef <- obj$coefficients
    histsize <- obj$nobs
    switch(obj$type,
           "ME" = { winsize <- obj$par },
       { winsize <- 1 })
    K <- floor(winsize*obj$nobs)
    Q12s <- obj$Q12 * K/(obj$sigma*sqrt(histsize))

    computeEmpProc <- function(newcoef){
        t(Q12s %*%(newcoef-histcoef))
    }

    logPlus <- function(x){
        klein<-x<=exp(1)
        z <- x*0
        z[klein] <- 1
        z[!klein] <- log(x[!klein])
        z
    }
    
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
        
        computeCoef <- function(x, y, k){
            ok <- (k-K):k
            coef(lm.fit(x[ok,,drop=FALSE], y[ok,,drop=FALSE]))
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
        
            
        computeCoef <- function(x, y, k){
            coef(lm.fit(x[1:k,,drop=FALSE], y[1:k,,drop=FALSE]))
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
                computeCoef=computeCoef,
                computeEmpProc=computeEmpProc,
                border=border, computeStat=computeStat,
                functional=functional, alpha=alpha, critval=critval,
                histcoef=histcoef, formula=obj$formula,
                type.name=paste("Monitoring with", obj$type.name))

    class(obj) <- "mefp"
    obj
}

monitor.mefp <- function(obj, data=NULL, verbose=TRUE){
  
    if(!is.na(obj$breakpoint)) return(TRUE)
    if(missing(data)){
        if(is.null(as.list(obj$efpcall)$data)){
            data <- list()
        }
        else{
            data <- get(as.character(as.list(obj$efpcall)$data))
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
        newcoef <- obj$computeCoef(x,y,k)
        obj$process <- rbind(obj$process,
                             obj$computeEmpProc(newcoef))
        stat <- obj$computeStat(obj$process)
        obj$statistic <- c(obj$statistic, stat)
        if(!foundBreak & (stat > obj$border(k))){
            foundBreak <- TRUE
            obj$breakpoint <- k
            if(verbose) cat("Break detected at observation #", k, "\n")
        }   
    }
    obj$last <- k
    obj$lastcoef <- newcoef
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

plot.mefp <- function(obj, main=NULL, ...){

    if(obj$last>obj$histsize){
        y1 <- rbind(as.matrix(obj$efpprocess),
                    as.matrix(obj$process))
        y1 <- ts(y1, start=start(obj$efpprocess),
                 frequency=frequency(obj$efpprocess))
        y2 <- ts(obj$border((obj$histsize+1):obj$last),
                 end = end(y1), frequency=frequency(y1))
        plot(y1, ty="l", ylim=c(min(y1,-y2), max(y1,y2)), ...)
        lines(y2, col=2)
        lines(-y2, col=2)
        abline(v=max(time(obj$efpprocess)), lty=2)
        abline(h=0)
        if(is.null(main))
            main <- obj$type.name
        title(main)
    }
    else{
        cat("Nothing monitored yet!\n")
    }
}


