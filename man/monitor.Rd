\name{monitor}
\alias{monitor}
\title{Monitoring of structural changes}
\description{
  Online monitoring of structural breaks in a linear regression model. A
  parameter estimate based on a historical sample is compared with
  estimates based on newly arriving data. A sequential test on the
  difference between the two parameter estimates signals structural breaks.
}
\usage{
monitor(obj, data=NULL, verbose=TRUE)
}
\arguments{
  \item{obj}{Object of class \code{mefp}}
  \item{verbose}{If \code{TRUE}, signal breaks by text output.}
}
\details{
  to be written
}
\value{\code{\link{mefp}} returns an object of class \code{mefp} which contains
a historic fluctuation process and sets the initialization for monitoring an
empirical fluctuation process in the future. When new data gets available an
\code{mefp} object can be monitored with the function \code{monitor}, which
also returns an object of class \code{mefp}. The latter contains not only the
historic but also a new empirical fluctuation process. Both can be plotted with
boundaries for the monitoring period by \code{plot.mts}.}

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
me1 <- mefp(e1, alpha=0.05)

# and now monitor on all data
me2 <- monitor(me1, data=df1)
plot(me2)
}
\keyword{regression}
