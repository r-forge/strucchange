\name{mefp.efp}
\alias{mefp.efp}
\title{monitoring empirical fluctuation processes}
\description{Compute an empirical fluctuation process for monitoring.}
\usage{
mefp(obj, alpha=0.05, functional = c("max", "range"), period=10,
     tolerance=.Machine$double.eps^0.5, MECritvalTable=monitorMECritvalTable,
     rescale=FALSE)

}
\arguments{
  \item{obj}{Object of class \code{"efp"}.}
  \item{alpha}{Significance level of the test, i.e., probability of
    type I error.}
  \item{functional}{Compare maximum or range of parameter differences,}
  \item{period}{(only used for moving estimates). Maximum time (relative
    to the history period) that will be monitored. Default is 10 times
    the history period.}
  \item{tolerance}{Tolerance for numeric \code{==} comparisons.}
  \item{MECritvalTable}{(only used for moving estimates) Table of
    critical values, this table is interpolated to get critical values
    for arbitrary \code{alpha}s. By default the pre-computed table
    \code{monitorMECritvalTable} is used.}
  \item{rescale}{logical. If \code{TRUE} the estimates will be standardized by
    the regressor matrix of the corresponding subsample according to Kuan & Chen
    (1994); if \code{FALSE} the whole regressor matrix will be used.
    (only if \code{type} is either \code{"fluctuation"} or \code{"ME"})}

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

\seealso{\code{\link{mefp.efp}}, \code{\link{plot.mefp}}, \code{\link{monitor}}}

\examples{
}
\keyword{dplot}
