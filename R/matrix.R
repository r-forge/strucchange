root.matrix <- function(X)
{
    if((ncol(X)==1)&&(nrow(X)==1)) return(sqrt(X))
    else
    {
        X.eigen <- eigen(X, symmetric=TRUE)
        if(any(X.eigen$values < 0)) stop("matrix is not positive semidefinite")
        sqomega <- sqrt(diag(X.eigen$values))
        V <- X.eigen$vectors
        return(V %*% sqomega %*% t(V))
    }
}

solveCrossprod <- function(X, method = c("qr", "chol", "solve")) {
  switch(match.arg(method),
    "qr" = chol2inv(qr.R(qr(X))),
    "chol" = chol2inv(chol(crossprod(X))),
    "solve" = solve(crossprod(X)))
}


covHC <- function(x, type = c("HC2", "const", "HC", "HC1", "HC3"))
{
  if(is.matrix(x$x))
    X <- x$x
  else {
    mf <- model.frame(x)
    X <- model.matrix(terms(x), mf)    
  }
  res <- residuals(x)

  n <- nrow(X)
  k <- ncol(X)
  Q1 <- summary(x)$cov.unscaled
  sigma2 <- var(res)*(n-1)/(n-k)
  type <- match.arg(type)

  if( type == "const") {
    V <- sigma2 * Q1
  } else {
    if(type == "HC2")
    {
      diaghat <- 1 - diag(X %*% Q1 %*% t(X))
      res <- res/sqrt(diaghat)
    }
    if(type == "HC3")
    {
      diaghat <- 1 - diag(X %*% Q1 %*% t(X))
      res <- res/diaghat
      Xu <- as.vector(t(X) %*% res)
    }
    VX <- res * X
    if(type %in% c("HC", "HC1", "HC2")) { V <- crossprod(crossprod(t(VX), Q1)) }
    if(type == "HC1") {V <- V * (n/(n-k))}
    if(type == "HC3") {V <- Q1 %*% (crossprod(VX) - (outer(Xu,Xu) /n)) %*% Q1 * (n-1)/n}
  }
  return(V)
}

