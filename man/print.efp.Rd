\name{print.efp}
\alias{print.efp}
\title{printing empirical fluctuation processes}
\description{print method for objects of class \code{"efp"}}
\usage{
print(x, ...)

\arguments{
  \item{x}{an object of class "\code{efp}".}
  \item{...}{additional parameters to \code{print.ts}}
}

\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{efp}}, \code{\link{plot.efp}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
print(efp(nhtemp ~ 1, type = "OLS-CUSUM"))
}
\keyword{print}
