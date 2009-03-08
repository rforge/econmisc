\name{interval}
\alias{interval}
\title{Interval Regression}
\description{
  This function estimates interval regression using either common
  intervals to all observations, or observation-specific intervals.
  Normal, logistic, log-log, and Cauchy disturbances are supported.
}
\usage{
interval(formula, data, weights, start, boundaries, ...,
subset, na.action, contrasts = NULL, Hess = FALSE, model = TRUE,
method = c("logistic", "probit", "cloglog", "cauchit"), print.level = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class '"formula"' (or one that can be coerced to
          that class): a symbolic description of the model to be
          fitted.  The left-hand-side variable must be either a factor,
	  decribing the intervals, observations fall to, or a matrix of
	  two columns, describing the interval boundaries by
	  observation.  See details below, see
	  also \code{\link{lm}} for explanation how to write formulas.}
  \item{data}{an optional data frame, list or environment (or object
          coercible by 'as.data.frame' to a data frame) containing the
          variables in the model.  If not found in 'data', the
          variables are taken from 'environment(formula)', typically
          the environment from which 'interval' is called.}
  \item{weights}{an optional vector of weights to be used in the fitting
          process.  Should be 'NULL' or a numeric vector. If non-NULL,
          weighted likelihood is maximized (that
          is, maximising 'sum(w*loglik)').}
  \item{start}{Initial values for the optimization algorithm.  if
    'NULL', initial values are calculated based on interval means.}
  \item{boundaries}{Boundaries for intervals.  See details}
  \item{\dots}{further arguments to \code{\link[maxLik]{maxLik}}}
  \item{subset}{an optional logical vector specifying a subset of observations to be
          used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  The default is set by the 'na.action' setting
          of 'options', and is 'na.fail' if that is unset.  The
          'factory-fresh' default is 'na.omit'.  Another possible value
          is 'NULL', no action.  Value 'na.exclude' can be useful.}
  \item{contrasts}{an optional list. See the 'contrasts.arg' of
          \code{\link{model.matrix.default}}.}
  \item{Hess}{Should Hessian of the model be returned}
  \item{model}{logical for whether the model matrix should be returned.}
  \item{method}{character, distribution of disturbances: either
    "logistic", "probit", "cloglog" or "cauchit".}
  \item{print.level}{output of run-time information: higher level prints
    more.}
}
\details{
  Interval regression is a form of linear regression where one observes
  only intervals (i.e numeric bounds) of the continuous outcome
  variable.  The current implementation assumes known distribution of
  the error term (default is normal) and estimates the model using
  Maximum Likelihood.

  The intervals may be specified in two ways: either common intervals
  for all the observations by using argument 'boundaries', or by
  specifying the response as Nx2 matrix, where columns correspond to the
  lower- and upper bound for the individual observations.  If you
  specify the common boundaries, you may want to give names to the
  components of 'boundaries'.  These names will appear to the estimation
  results.  For observations-specific boundaries, the names are
  generated automatically.
}
\value{
  Object of class "Interval" which inherits from class "maxLik".  The
  latter has several methods, including \code{\link{summary}} and
  \code{\link{coef}}.
  
  Note that the boundaries are passed as fixed parameters to the
  \code{maxLik} estimation routine and hence returned as fixed estimates
  with 0 standard errors.
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\note{A lot of code is borrowed from \code{\link{polr}} from the package MASS
}
\seealso{\code{\link{polr}}}
\examples{
## Example of observation-specific boundaries
## Estimate the willingness to pay for the Kakadu National Park
## Data given in intervals -- 'lower' for lower bound and 'upper' for upper bound.
## Note that dichotomous-coice answers are already coded to 'lower' and 'upper'
library(Ecdat)
data(Kakadu)
## Estimate in log form, change 999 to Inf
lb <- log(Kakadu$lower)
ub <- Kakadu$upper
ub[ub > 998] <- Inf
ub <- log(ub)
y <- cbind(lb, ub)
m <- interval(y ~ sex + log(income) + age + schooling + 
              recparks + jobs + lowrisk + wildlife + future + aboriginal + finben +
              mineparks + moreparks + gov +
              envcon + vparks + tvenv + major, data=Kakadu)
## You may want to compare the results to Werner (1999),
## Journal of Business and Economics Statistics 17(4), pp 479-486
print(summary(m))
##
## Example of common intervals for all the observations
##
library(Ecdat)
data(Bwages)
## calculate an ordinary Mincer-style wage regression
## Note: gross hourly wage rate in EUR
ols <- lm(log(wage) ~ factor(educ) + poly(exper, 2), data=Bwages)
cat("OLS estimate:\n")
print(summary(ols))
## Now we censor the wages to intervals
intervals <- c(0, 5, 10, 15, 25, Inf)
salary <- cut(Bwages$wage, intervals)
int <- interval(salary ~ factor(educ) + poly(exper, 2), data=Bwages, boundaries=log(intervals))
## Note: use logs for the intervals in Euros.  We do not have to
## transform salaris to log form as this does not change the intervals.
## Ignore any warnings
cat("Interval regression:\n")
print(summary(int))
}
\keyword{models}
\keyword{regression}
