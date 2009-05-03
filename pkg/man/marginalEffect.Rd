\name{marginalEffect}
\Rdversion{1.1}
\alias{marginalEffect}
\alias{marginalEffect.probit}
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
  Ott Toomet \email{otoomet@ut.ee}
}
\seealso{
  \code{\link[pkg:sampleSelection]{probit}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
