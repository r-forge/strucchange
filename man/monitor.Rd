\name{monitor}
\alias{monitor}
\alias{monitor.efp}
\alias{monitor.mefp}
\alias{monitorMECritvalTable}
\alias{plot.mefp}
\alias{print.mefp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Monitoring of Structural Breaks}
\description{
  Online monitoring of structural breaks in a linear regression model. A
  parameter estimate based on a historical sample is compared with
  estimates based on newly arriving data. A sequential test on the
  difference between the two parameter estimates signals structural breaks.
}
\usage{
monitor(obj, ...)
monitor.efp(obj, alpha=0.05, functional = c("max", "range"),
            period=10, tolerance=.Machine$double.eps^0.5,
            MECritvalTable=monitorMECritvalTable)
monitor.mefp(obj, data=NULL, verbose=TRUE)
}
\arguments{
  \item{alpha}{Significance level of the test, i.e., probability of
    type I error.}
  \item{functional}{Compare maximum or range of parameter differences,}
  \item{period}{(only used for moving estimates). Maximum time (relative
    to the history period) that will be monitored. Default is 10 times
    the history period.}
  \item{tolerance}{Tolerance for numeric \code{==} comparisons.}
  \item{verbose}{If \code{TRUE}, signal breaks by text output.}
  \item{MECritvalTable}{(only used for moving estimates) Table of
    critical values, this table is interpolated to get critical values
    for arbitrary \code{alpha}s. By default the pre-computed table
    \code{monitorMECritvalTable} is used.}
}
\details{
  to be written
}
\value{
  An object of class \code{"mefp"}.
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch, Kurt Hornik, and Chung-Ming Kuan. Monitoring
  structural changes with the generalized fluctuation test.
  Econometric Theory, 16:835-854, 2000.
}
\examples{
df1 <- data.frame(y=rnorm(300))
df1[150:300,"y"] <- df1[150:300,"y"]+1

# use the first 50 observations as history period
e1 <- efp(y~1, data=df1[1:50,,drop=FALSE], type="ME", h=0.25)
me1 <- monitor(e1, alpha=0.05)

# and now monitor on all data
me2 <- monitor(me1, data=df1)
plot(me2)
}
\keyword{regression}%-- one or more ...
