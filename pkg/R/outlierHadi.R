# S-Plus code by Ali S. Hadi, available at
# http://www1.aucegypt.edu/faculty/hadi/Research.html#splus
# ported to R by Arne Henningsen
outlierHadi <- function( data ) {
# -----------------------------------------------------------------
#  Hadi, Ali S. (1994), "A Modification of a Method for the
#  Detection of Outliers in Multivariate Samples," Journal of the
#  Royal Statistical Society (B), 2, 393-396.
# -----------------------------------------------------------------
  X <- as.matrix( data )
  n <- dim(X) [1]
  p <- dim(X) [2]
  h <- trunc((n + p + 1)/2)     
  id <- 1:n
  r <- p
  out <- 0
  cf <- (1 + ((p + 1)/(n - p)) + (2/(n - 1 - (3*p))) )^2
# cf <- (1 + ((p + 1)/(n - p)) + (1/(n - p - h)) )^2
  alpha <- 0.05
  tol <- max(10^-(p+5), 10^-12)
# -----------------------------------------------------------------
# **  Compute Mahalanobis distance
# -----------------------------------------------------------------
  C <- apply(X, 2, mean)
  S <- var(X)
  if (det(S) < tol) stop ()
  D <- mahalanobis(X, C, S)
  mah.out <- 0
  cv <- qchisq(1-(alpha/n), p)
  for (i in 1:n) if (D[i] >= cv) mah.out <- cbind(mah.out, i)
  mah.out <- mah.out[-1]
  mah <- sqrt(D)
  Xbar <- C
  Covariance <- S   #
# ----------------------------------------------------------------
#  **  Step 0
# ----------------------------------------------------------------
#  **  Compute Di(Cm, Sm)
  C <- apply(X, 2, median)
  C <- array(C, dim = c(n, p))
  Y <- X - C
  S <- ((n - 1)^-1)*(t(Y) %*% Y)
  D <- mahalanobis(X, C[1, ], S)
  Z <- sort.list(D)
# ----------------------------------------------------------------
#  **  Compute Di(Cv, Sv)
  repeat {
    Y <- X[Z[1:h], ]
    C <- apply(Y, 2, mean)
    S <- var(Y)
    if (det(S) > tol) {
       D <- mahalanobis(X, C, S)
       Z <- sort.list(D); break }
    else h <- h + 1
    }
# ----------------------------------------------------------------
#  **  Step 1
# ----------------------------------------------------------------
  repeat {
    r <- r + 1
    if ( h < r) break
    Y <- X[Z[1:r],]
    C <- apply(Y, 2, mean)
    S <- var(Y)
    if (det(S) > tol) {
       D <- mahalanobis(X, C, S)
       Z <- sort.list(D) }
    }
# ----------------------------------------------------------------
#  **  Step 3
# ----------------------------------------------------------------
#  **  Compute Di(Cb, Sb)
  repeat {
    Y <- X[Z[1:h],]
    C <- apply(Y, 2, mean)
    S <- var(Y)
    if (det(S) > tol) {
       D <- mahalanobis(X, C, S)
       Z <- sort.list(D)
       if (D[Z[h + 1]] >= (cf*qchisq(1-(alpha/n), p))) {
            out <- Z[(h + 1) : n]
            break }
       else { h <- h + 1
              if (n <= h) break }
       }
    else { h <- h + 1
          if (n <= h) break }
    }
  D <- sqrt(D/cf)
  dst <- cbind(id, mah, D)
  Outliers <- out
  Cb <- C;
  Sb <- S
  Distances <- dst
  result <- list(Xbar = Xbar, Covariance = Covariance, mah.out = mah.out, 
      Outliers = Outliers, Cb = Cb, Sb = Sb, Distances = Distances )
  class( result ) <- "outlierHadi"
  return( result )
}
# ----------------------------------------------------------------
