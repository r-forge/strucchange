\name{recresid_cpp}
\alias{recresid_cpp}
\title{Recursive Residuals Using the C++ Armadillo Library}
\description{
  Function for computing the recursive residuals
  (standardized one step prediction errors) of a linear regression model
  using the C++ Armadillo library.
}
\usage{
recresid_cpp(x, y, start = ncol(x) + 1, end = nrow(x),
  tol = sqrt(.Machine$double.eps)/ncol(x), \dots)
}
\arguments{
  \item{x, y}{specification of the linear regression model
    by a regressor matrix \code{x} and a response variable \code{y}.}
  \item{start, end}{integer. Index of the first and last observation, respectively,
    for which recursive residuals should be computed. By default, the maximal
    range is selected.}
  \item{tol}{numeric. A relative tolerance for precision of recursive
    coefficient estimates, see details.}
  \item{data}{an optional data frame containing the variables in the model. By
     default the variables are taken from the environment which \code{recresid} is
     called from. Specifying \code{data} might also be necessary when applying
     \code{recresid} to a fitted model of class \code{"lm"} if this does not
     contain the regressor matrix and the response.}
  \item{\dots}{\emph{currently not used.}}
}
\details{
  See function \code{\link[strucchange]{recresid}}.
}
\value{
  A vector containing the recursive residuals.
}

\seealso{\code{\link[strucchange]{recresid}}}

\examples{
y <- rnorm(100) + rep(c(0, 2), each = 50)
x <- matrix(1, nrow = 100, ncol = 1)
rr <- recresid_cpp(x, y)
print(rr)
}

\keyword{regression}
