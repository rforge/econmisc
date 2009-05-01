\name{vcovClust}
\alias{vcovClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clustering-robust variance-covariance matrices}
\description{
  Function to calculate clustering-robust variance-covariance matrix.
  Currently works for OLS and 1 or 2-way clustering.  Interface not well
  designed. 
}
\usage{
vcovClust(fm, dfcw, cluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fm}{a fitted model}
  \item{dfcw}{degrees-of-freedom adjustment parameter, 1 (= no adjustment) by default}
  \item{cluster}{a vector or nObsx2 matrix of clustering indicators}
}
\details{
  works for OLS but should be easier to use.  Uses sandwich package.
}
\value{
  a numeric matrix, clustering-robust variance-covariance matrix.
}
\references{ ~put references to the literature/web site here ~ }
\author{Ott Toomet \email{otoomet@ut.ee}}
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{cluster}