\name{plot.mefp}
\alias{plot.mefp}
\title{plot for monitoring empirical fluctuation processes}
\description{Plotting method for objects of class \code{"mefp"})}
\usage{
plot(x, boundary=TRUE, functional="max", main=NULL, ylab="empirical
fluctuation process", ylim=NULL, ...) }

\arguments{
  \item{x}{an object of class \code{"mefp"}.}
  \item{boundary}{logical. If set to \code{FALSE} the boundary will be computed
     but not plotted.}
  \item{functional}{indicates which functional should be applied to the
     process. If set to \code{NULL} a multiple process is plotted.}
  \item{...}{graphical parameters.}
}

\value{\code{\link{mefp}} returns an object of class \code{"mefp"} which 
contains a historic fluctuation process and sets the initialization for
monitoring an empirical fluctuation process in the future. When new data gets
available an \code{mefp} object can be monitored with the function
\code{monitor}, which also returns an object of class \code{"mefp"}. The latter
contains not only the historic but also a new empirical fluctuation process.
Both can be plotted with boundaries for the monitoring period by
\code{plot.mts}.}

\author{Friedrich Leisch}

\seealso{\code{\link{plot.mefp}}, \code{\link{monitor}}}

\examples{
}
\keyword{dplot}
