\name{sctest.Fstats}
\alias{sctest.Fstats}
\title{supF-, aveF- and expF-test}
\description{performs the supF-, aveF- or expF-test}
\usage{
sctest(x, type = c("supF", "aveF", "expF"), asymptotic = FALSE)}
\arguments{
 \item{x}{an object of class "\code{Fstats}".}
 \item{type}{a character string specifying which test will be performed.}
 \item{asymptotic}{logical. Only necessary if \code{x} contains just a single F
     statistic and type is \code{"supF"} or \code{"aveF"}. If then set to
     \code{TRUE} the asymptotic (chi-square) distribution instead of the exact
     (F) distribution will be used to compute the p value.}
}

\details{If \code{x} contains just a single F statistic and type is
\code{"supF"} or \code{"aveF"} the Chow test will be performed.

The original GAUSS code for computing the p values of the supF-, aveF- and
expF-test was written by Bruce Hansen and is available from
\url{http://www.ssc.wisc.edu/~bhansen/}. R port by Achim Zeileis.
}

\value{an object of class "\code{htest}" containing:
\itemize
\item{statistic}{the test statistic}
\item{p.value}{the corresponding p value}
\item{method}{the method used}
\item{data.name}{the data name}}

\references{
Andrews D.W.K. (1993), Tests for parameter instability and structural
change with unknown change point, \emph{Econometrica}, \bold{61}, 821-856.

Hansen B. (1992), Tests for Parameter Instability in Regressions With I(1)
Processes, \emph{Journal of Business & Economic Statistics}, \bold{10}, 321-335.

Hansen B. (1997), Approximate Asymptotic P Values for Structural-Change
Tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67. }

\seealso{\code{\link{Fstats}}, \code{\link{plot.Fstats}}}

\examples{
## Load dataset "nhtemp" with average yearly temperatures in New Haven
data(nhtemp)
## plot the data
plot(nhtemp)

## test the model null hypothesis that the average temperature remains constant
## over the years for potential break points between 1941 (corresponds to
## from = 0.5) and 1962 (corresponds to to = 0.85)
## compute F statistics
fs <- Fstats(nhtemp ~ 1, from = 0.5, to = 0.85)
## plot the F statistics
plot(fs, alpha = 0.01)
## and the corresponding p values
plot(fs, pval = TRUE, alpha = 0.01)
## perform the aveF test
sctest(fs, type = "aveF")
}
\keyword{htest}
