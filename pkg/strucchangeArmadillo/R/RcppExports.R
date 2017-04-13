#' Computation of recursive residuals in C++
#' @param X design matrix
#' @param y response vector
#' @param start integer (1-based) index of the first observation to compute recursive residuals
#' @param end integer (1-based) index of the last observation to compute recursive residuals
#' @param tol tolerance in the computation of recursive model coefficients
#' @return vector containing the recursive residuals 
#' @seealso \code{\link{recresid}} and \code{\link{recresid.default}}
.sc_cpp_recresid <- function(X, y, start, end, tol, rcond_min) {
  .Call('strucchange_sc_cpp_recresid', X, y, start, end, tol, rcond_min, PACKAGE = 'strucchangeArmadillo')
}

