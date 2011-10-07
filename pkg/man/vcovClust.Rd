\name{vcovClust}
\alias{vcovClust}
\title{Clustering-robust variance-covariance matrices}
\description{
  Function to calculate clustering-robust variance-covariance matrix.
  Currently works for 1 or 2-way clustering.   
}
\usage{
vcovClust(fm, dfcw=1, strata)
}
\arguments{
  \item{fm}{a fitted model}
  \item{dfcw}{degrees-of-freedom adjustment parameter, 1 (= no adjustment) by default}
  \item{strata}{a vector of length \code{nObs(fm)} or \code{nObs(fm)*2}
    matrix of clustering (strata) indicators
  }
}
\details{
  works for models which have 'vcov' method.
}
\value{
  a numeric matrix, clustering-robust variance-covariance matrix.
}
\references{Cameron A. C., Gelbach J. B. and Miller D. L. "Robust
  Inference with Multi-way Clustering", NBER working paber t0327, 2006}
\author{Ott Toomet \email{otoomet@ut.ee}}
%\seealso{\code{\link[pkg:robust]{covRob}}} % is this relevant?
\examples{
## Explain crime rate in NC by minority percentage (pctmin) and various wages (wxxx).
## For each county there are 7 observations with the same minority percentage
## As minority percentage is county specific, we have to cluster over it
library(Ecdat)
m <- lm(crmrte ~ pctmin + factor(year) +
        polpc + smsa + wcon + wtrd + wmfg + wsta +
        pctymle,
        data=Crime)
vc <- vcovClust(a, strata=Crime$county)
## compare the corresponding standard errors
library(miscTools)
rbind(sqrt(diag(vc)), stdEr(m))
                        # Note standard errors for pctmin and polpc have increased substantially
}
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{cluster}