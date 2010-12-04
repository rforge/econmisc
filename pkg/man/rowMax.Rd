\name{rowMax}
\alias{rowMax}
\title{
  Maximum row values of matrix
}
\description{
  This function finds the maximum values of matrix rows
}
\usage{
rowMax(mat)
}
\arguments{
  \item{mat}{
    numeric (including integer) matrix
}
}
\details{
  \code{NA}-s are currently not handled
}
\value{
  numeric vector, the maximum values of each matrix row
}
\author{
  Ott Toomet \email{otoomet@gmail.com}
}

\seealso{
  \code{\link{rowSums}} and \code{\link{colSums}} for matrix sums by
  rows and columns, and \code{\link{rowMeans}} and
  \code{\link{colMeans}} for the corresponding means.
}
\examples{
a <- matrix(1:6, 3, 2)
rowMax(a)
}
\keyword{array}
\keyword{utilities}
