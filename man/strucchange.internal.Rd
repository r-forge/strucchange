\name{pvalue.efp}
\alias{sc.beta.sup}
\alias{sc.beta.ave}
\alias{sc.beta.exp}
\alias{sc.me}
\alias{pvalue.efp}
\alias{pvalue.Fstats}
\alias{monitorMECritvalTable}
\title{p value computation}
\description{functions and tables of critical values and coefficients 
respectively for the computation of p values (for internal use only)}
\usage{
pvalue.efp(x, type, alt, functional = "max", h = NULL, k = NULL)
pvalue.Fstats(x, type = c("supF", "aveF", "expF"), k, lambda)
}
\source{Hansen B. (1997), Approximate Asymptotic P Values for Structural-Change
Tests, \emph{Journal of Business & Economic Statistics}, \bold{15}, 60-67.

Chu C.-S., Hornik K., Kuan C.-M. (1995), The Moving-Estimates Test for
Parameter Stability, \emph{Econometric Theory}, \bold{11}, 669-720.}
\keyword{data}


