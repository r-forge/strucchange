\name{print.Fstats}
\alias{print.Fstats}
\title{printing F statistics}
\description{print method for objects of class \code{"Fstats"}}
\usage{
print(x, ...)

\arguments{
  \item{x}{an object of class "\code{Fstats}".}
  \item{...}{additional parameters to \code{print.ts}}
}

\author{Achim Zeileis \email{zeileis@ci.tuwien.ac.at}}

\seealso{\code{\link{Fstats}}, \code{\link{plot.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
print(Fstats(nhtemp ~ 1))
}
\keyword{print}
