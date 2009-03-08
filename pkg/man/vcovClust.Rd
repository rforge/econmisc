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
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(fm, dfcw, cluster) {
   ## fm       fitted model
   ## dfcw     degrees of freedom correction
   ## cluster  factor: indicator for clusters
   library(sandwich);
   library(lmtest)
   ## if NA-observations are omitted from the model frame, we have to remove the corresponding elements from cluster as well.
   ## we do it by row names of model frame and hope that either 1) cluster has the same row names, or 2) row names are integers,
   ## describing the retained rows of the original data.
   ## Is there a better way to extract the omitted rows?
   dataRows <- row.names(model.frame(fm))
   cluster <- as.matrix(cluster)
                                        # results to Nx1 or Nx2 matrix
   if(!is.null(names(cluster))) {
      cluster <- cluster[dataRows,,drop=FALSE]
   }
   else {
      cluster <- cluster[as.integer(dataRows),,drop=FALSE]
   }
   ## extract the clusters
   cluster1 <- factor(cluster[,1])
   M1 <- length(levels(cluster1))
   N <- length(cluster1)
   K <- fm$rank
   dfc1  <- (M1/(M1-1))*((N-1)/(N-K))  
   u1j  <- apply(estfun(fm),2, function(x) tapply(x, cluster1, sum));
   vcovCL <- dfc1*sandwich(fm, meat=crossprod(u1j)/N)
   if(ncol(cluster) == 2) {
      cluster2 <- factor(cluster[,2])
      cluster12 <- (as.integer(cluster1) - 1)*length(levels(cluster2)) + as.integer(cluster2)
      cluster12 <- factor(cluster12)
      M2 <- length(levels(cluster2))
      M12 <- length(levels(cluster12))
      dfc2  <- (M2/(M2-1))*((N-1)/(N-K))  
      dfc12 <- (M12/(M12-1))*((N-1)/(N-K))  
      u2j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster2,  sum)) 
      u12j  <- apply(estfun(fm), 2, function(x) tapply(x, cluster12, sum)) 
      vc2   <-  dfc2*sandwich(fm, meat=crossprod(u2j)/N )
      vc12  <- dfc12*sandwich(fm, meat=crossprod(u12j)/N)
      vcovCL <- vcovCL + vc2 - vc12
   }
   if(any(diag(vcovCL) < 0)) {
      warning("diagonal of cluster-robust covariance matrix not positive")
   }
   vcovCL*dfcw
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{cluster}