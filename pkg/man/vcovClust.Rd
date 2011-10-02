\name{vcovClust}
\alias{vcovClust}
\title{Clustering-robust variance-covariance matrices}
\description{
  Function to calculate clustering-robust variance-covariance matrix.
  Currently works for OLS and 1 or 2-way clustering.  Interface not well
  designed. 
}
\usage{
vcovClust(fm, dfcw, strata)
}
\arguments{
  \item{fm}{a fitted model}
  \item{dfcw}{degrees-of-freedom adjustment parameter, 1 (= no adjustment) by default}
  \item{strata}{a vector of length \code{nObs(fm)} or \code{nObs(fm)*2}
    matrix of clustering (strata) indicators
  }
}
\details{
  works for OLS but should be easier to use.  Uses sandwich package.
}
\value{
  a numeric matrix, clustering-robust variance-covariance matrix.
}
\references{Cameron A. C., Gelbach J. B. and Miller D. L. "Robust
  Inference with Multi-way Clustering", NBER working paber t0327, 2006}
\author{Ott Toomet \email{otoomet@ut.ee}}
%\seealso{\code{\link[pkg:robust]{covRob}}} % is this relevant?
\examples{
}
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{cluster}