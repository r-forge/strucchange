\name{covHC}
\alias{covHC}
\title{Heteroskedasticity-Consistent Covariance Matrix Estimation}
\description{
 Heteroskedasticity-consistent estimation of the covariance matrix of the
 coefficient estimates in a linear regression model.
}
\usage{
covHC(formula, type = c("HC2", "const", "HC", "HC1", "HC3"), data=list())
}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested.}
  \item{type}{a character string specifying the estimation type. For
     details see below.}
 \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which \code{dwtest}
    is called from.}
}
\details{When \code{type = "const"} constant variances are assumed and
  and \code{cov.estim} just gives a diagonal matrix containing the usual
  estimates for the variances of the coefficient estimates.

  All other methods do not assume constant variances and are suitable in case of
  heterskedasticity. For details see the references.
}
\value{A matrix containing the covariance matrix estimate.
}
\references{

MacKinnon J.G., White H. (1985),
Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
\emph{Journal of Econometrics} \bold{29}, 305-325

}

\seealso{\code{\link{lm}}}

\examples{
}
\keyword{htest}
