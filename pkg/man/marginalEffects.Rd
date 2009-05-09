\name{marginalEffects}
\Rdversion{1.1}
\alias{coef.marginalEffects}
\alias{marginalEffect}
\alias{marginalEffects}
\alias{marginalEffects.probit}
\alias{print.marginalEffects}
\alias{print.summary.marginalEffects.probit}
\alias{summary.marginalEffects.probit}
\alias{vcov.marginalEffects}
\title{
  marginal effects of non-linear models
}
\description{
  calculate marginal effects of explanatory variables,
  primary for non-linear models.
}
\usage{
marginalEffects(object, ...)
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

  In most cases, the marginal effects depend on the values of
  explanatory variables and hence are not unique.  The current method
  computes the average of the marginal effects across the individual observations.
}
\value{
  Object of class 'marginalEffects' (and probably of a more specific
  class as well, such as 'marginalEffects.probit').
  
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
library(Ecdat)
data(Participation)
m <- probit(lfp ~ lnnlinc + poly(age, 2) + educ + I(nyc > 0) + I(noc > 0) + foreign,
            data=Participation)
# print the predicted impact on the linear predictor
print(summary(m))
# print the predicted impact on the participation probability
print(marginalEffects(m))
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
