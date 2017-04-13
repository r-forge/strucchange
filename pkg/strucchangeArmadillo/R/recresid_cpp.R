## Receresid using C++.
recresid_cpp <- function(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), ...)
{
  ## checks and data dimensions
  stopifnot(start > ncol(x) & start <= nrow(x))
  stopifnot(end >= start & end <= nrow(x))
  return(.sc_cpp_recresid(x, y, start, end, tol, sqrt(.Machine$double.eps)))
}

