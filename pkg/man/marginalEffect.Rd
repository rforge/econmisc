\name{marginalEffect}
\Rdversion{1.1}
\alias{coef.marginalEffect}
\alias{marginalEffect}
\alias{marginalEffect.probit}
\alias{print.marginalEffect}
\alias{print.summary.marginalEffect.probit}
\alias{summary.marginalEffect.probit}
\alias{vcov.marginalEffect}
\title{
  marginal effects of non-linear models
}
\description{
  calculate marginal effects of explanatory variables,
  primary for non-linear models.
}
\usage{
marginalEffect(object, ...)
}
\arguments{
  \item{object}{an estimated model, marginal effects of which are to be calculated.
  }
  \item{\dots}{
    additional arguments for other methods.
  }
}
\details{
  Marginal effect measures the impact of change of the
  explanatory variables on the outcome variable.
  In general, this measure differs from the estimated
  coefficients (a notable exception being OLS).

  In most cases, the marginal effect depends on the values of
  explanatory variables and hence is no unique.  The current method
  computes the average of marginal effects of individual observations.
}
\value{
  Object of class 'marginalEffect' (and probably of a more specific
  class as well, such as 'marginalEffect.probit').
  
  There are coef, print, summary, vcov and std methods for extracting
  the effects, standard deviations, and printing model summaries.
}
\author{
  Ott Toomet \email{otoomet@ut.ee}
}
\seealso{
  \code{\link[pkg:sampleSelection]{probit}}
}
\examples{
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
