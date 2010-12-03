\name{rcat}
\alias{rcat}
\title{
  Categorical random variables
}
\description{
  This function produces a vector of categorical random variables (in
  range 1 to K).  The probabilities may differ for each draw.
}
\usage{
rcat(n, prob)
}
\arguments{
  \item{n}{
    numeric, how many draws to generate
}
  \item{prob}{
    numeric, probabilities of categories.  The length (or number of
    columns) determines the number of categories which will be numbered
    as 1 to K.  This may be either a vector of length
    K for
    similar probabilities for all the draws, or a n by K matrix, where
    each row specifies the probabilities for the corresponding draw.

    The probabilities are normalized internally to sum to unity by rows.
}
}
\details{
  The function is implemented in C++.
}
\value{
  An integer vector of length n, the random draws from the categories.
}
\author{
  Ott Toomet \email{otoomet@gmail.com}
}
\seealso{
  \code{\link[stats]{rmultinom}} for multinomial random variables.
}
\examples{
prob <- 1:3  # similar probabilities for all the draws
rcat(20, prob)
prob <- matrix(1:18, 6, 3, byrow=TRUE) # probs differ by draws
rcat(6, prob)
}
\keyword{distribution}
