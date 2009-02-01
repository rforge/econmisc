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
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{method}{ ~~Describe \code{method} here~~ }
  \item{print.level}{ ~~Describe \code{print.level} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
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
function(formula, data, weights, start, boundaries,
                 ..., subset,
                 na.action, contrasts = NULL, Hess = FALSE,
                 model = TRUE,
                 method = c("logistic", "probit", "cloglog", "cauchit"),
                     print.level=0)
{
   ## beta     parameters for the x-s (linear model)
   ## nBeta    number of x-s (including constant)
    logit <- function(p) log(p/(1 - p))
    ## log likelihood
    loglik <- function(theta) {
       beta <- theta[iBeta]
        zeta <- theta[iZeta]
        sd <- theta[iSd]
        eta <- offset
        if (nBeta > 0)
            eta <- eta + drop(x \%*\% beta)
        Pr <- pfun((zeta[boundaryInterval + 1] - eta)/sd) - pfun((zeta[boundaryInterval] - eta)/sd)
       names(Pr) <- NULL
                                        # individual probability to be in interval (zeta[y+1], zeta[y]])
       ll <- wt * log(Pr)
       (ll)
    }
    ## gradient
    gradlik <- function(theta)
    {
        jacobian <- function(theta) { ## dzeta by dtheta matrix
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0 , k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
         }
        beta <- theta[iBeta]
        zeta <- theta[iZeta]
        sigma <- theta[iSd]
        eta <- offset
        if(nBeta > 0)
            eta <- eta + drop(x \%*\% beta)
        normArg1 <- (zeta[boundaryInterval +1] - eta)/sigma
        normArg2 <- (zeta[boundaryInterval] - eta)/sigma
        Pr <- pfun(normArg1) - pfun(normArg2)
        p1 <- dfun(normArg1)
        p2 <- dfun(normArg2)
        g1 <-
                                        # d loglik / d beta
            if(nBeta > 0)
                x * (wt*(p2 - p1)/Pr/sigma)
            else
                matrix(0, 0, nBeta)
        xx <- .polrY1*p1 - .polrY2*p2
        dg.dzeta <- xx * (wt/Pr/sigma)
##         g2 <- - t(xx) \%*\% (wt/Pr/sigma)
##         g2 <- t(g2) \%*\% jacobian(zeta)
##         if(all(Pr > 0))
##             -c(g1, g2)
##         else
##             rep(NA, nBeta+nInterval)
        s1 <- p1*normArg1
        s1[is.infinite(normArg1)] <- 0
                                        # if a boundary is Inf, we get 0*inf type of NaN
        s2 <- p2*normArg2
        s2[is.infinite(normArg2)] <- 0
        dg.dsd <- -(s1 - s2)*(wt/Pr/sigma)
        grad <- cbind(g1, dg.dzeta, dg.dsd)
        grad
     }
    ## ---------- main function -------------
    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    pfun <- switch(method, logistic = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = npcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm,
                   cloglog = dgumbel, cauchit = dcauchy)
    if(is.matrix(eval.parent(m$data)))
        m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$model <- m$boundaries <- m$... <- m$print.level <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    nObs <- nrow(x)
    nBeta <- ncol(x)
    cons <- attr(x, "contrasts") # will get dropped by subsetting
##     if(xint > 0) {
##        x <- x[, -xint, drop=FALSE]
##     }
##     else
##         warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if(!length(wt)) wt <- rep(1, nObs)
    offset <- model.offset(m)
    if(length(offset) <= 1) offset <- rep(0, nObs)
    y <- model.response(m)
    if(is.matrix(y)) {
       ## use the intervals, given for each observation.  We have interval regression and the interval boundaries are fixed
       ordered <- FALSE
       dimnames(y) <- NULL
       lowerBound <- y[,1]
       upperBound <- y[,2]
       ## in case of interval regression, we have to construct a set of intervals and pack them correctly to the
       ## parameter vector
       library(sets)
                                        # we use sets here (sorry, need R >= 2.7)
       intervals <- set()
       for(i in 1:length(lowerBound)) {
          intervals <- intervals | c(lowerBound[i], upperBound[i])
       }
       ## Now make the set to a list to have ordering
       intervals <- as.list(intervals)
       ## Now find which interval each observation falls into
       y <- boundaryInterval <- numeric(length(lowerBound))
                                        # which interval observation falls to 
                                        # y:  in terms of ordered intervals
                                        # boundaryInterval         in terms of ordered boundaries
       for(i in seq(along=intervals)) {
          j <- lowerBound == intervals[[i]][1] & upperBound == intervals[[i]][2]
          boundaryInterval[j] <- 1 + 2*(i - 1)
          y[j] <- i
       }
       boundaries <- unlist(intervals)
                                        # boundaries as a vector (note the joint boundaries are twice
       names(boundaries) <- paste(c("L", "U"), rep(seq(along=intervals), each=2))
       nInterval <- length(boundaries) - 1
    }
    else {
       ## response is given as an ordered factor
       lev <- levels(y)
       if(length(lev) <= 2)
           stop("response must have 3 or more levels")
       if(!is.factor(y))
           stop("response must be a factor or Nx2 matrix of boundaries")
       y <- unclass(y)
       nInterval <- length(lev)
       if(missing(boundaries))
           ordered <- TRUE
       else {
           ordered <- FALSE
           intervals <- vector("list", length(boundaries) - 1)
           for(i in seq(length=length(boundaries) - 1)) {
              intervals[[i]] <- c(boundaries[i], boundaries[i+1])
           }
           boundaryInterval <- 2*y - 1
        }
    }
    Y <- matrix(0, nObs, nInterval + 1)
    .polrY1 <- col(Y) == boundaryInterval + 1
    .polrY2 <- col(Y) == boundaryInterval 
                                        # .polr are markers for which interval the boundaryInterval falls to
    ## starting values
    iBeta <- seq(length=ncol(x))
    iZeta <- nBeta + seq(along=boundaries)
    iSd <- max(iZeta) + 1
    if(missing(start)) {
       if(ordered) {
          ## try logistic/probit regression on 'middle' cut
          q1 <- nInterval \%/\% 2
          y1 <- (y > q1)
          fit <-
              switch(method,
                     "logistic"= glm.fit(x, y1, wt, family = binomial(), offset = offset),
                     "probit" = glm.fit(x, y1, wt, family = binomial("probit"), offset = offset),
                     ## this is deliberate, a better starting point
                     "cloglog" = glm.fit(x, y1, wt, family = binomial("probit"), offset = offset),
                     "cauchit" = glm.fit(x, y1, wt, family = binomial("cauchit"), offset = offset))
          if(!fit$converged)
              warning("attempt to find suitable starting values did not converge")
          coefs <- fit$coefficients
          if(any(is.na(coefs))) {
             warning("design appears to be rank-deficient, so dropping some coefs")
             keep <- names(coefs)[!is.na(coefs)]
             coefs <- coefs[keep]
#          x <- x[, keep[-1], drop = FALSE]
             ## note: we keep the intercept
             nBeta <- ncol(x)
          }
          spacing <- logit((1:nInterval)/(nInterval+1)) # just a guess
          if(method != "logit") spacing <- spacing/1.7
          zetas <- -coefs[2] + spacing - spacing[q1]
          coefs[1] <- 0
          activePar <- c(FALSE, rep(TRUE, length(coef) - 1 + length(zetas) + 1))
                                        # intercept is fixed to 0
          sigma <- 1
       }
       else {
          ## not ordered: estimate OLS on interval middle points
          means <- sapply(intervals, mean)
                                        # we have to put a reasonable value to infinite intervals.
                                        # Pick the average width of the interval and use it as the meanpoint
          widths <- sapply(intervals, function(x) x[2] - x[1])
          meanWidth <- mean(widths[!is.infinite(widths)])
          negInf <- is.infinite(means) & means < 0
          means[negInf] <- sapply(intervals[negInf], function(x) x[2] - meanWidth)
          posInf <- is.infinite(means) & means > 0
          means[posInf] <- sapply(intervals[posInf], function(x) x[1] + meanWidth)
          yMean <- means[y]
          fit <- lm(yMean ~ x - 1)
          coefs <- coef(fit)
          names(coefs) <- gsub("^x", "", names(coefs))
          zetas <- boundaries
          activePar <- c(TRUE, rep(TRUE, length(coefs) - 1), rep(FALSE, length(boundaries)), TRUE)
                                        # intercept is free
          sigma <- sqrt(var(fit$residuals))
       }
       start <- c(coefs, zetas, "sigma"=sigma)
    }
    else
        if(length(start) != nBeta + nInterval)
            stop("'start' is not of the correct length")
    if(print.level > 0) {
       cat("Initial values:\n")
       print(start)
    }
    library(maxLik)
##     compareDerivatives(loglik, gradlik, t0=start)
##     stop()
    res <- maxLik(loglik, gradlik, start=start, method="BHHH", activePar=activePar, iterlim=500, ...)
    return(res)
    beta <- coef(res)[seq_len(nBeta)]
    zeta <- coef(res)[iZeta]
    deviance <- 2 * logLik(res)
    niter <- nIter(res)
    # names(zeta) <- paste(lev[-length(lev)], lev[-1], sep="|")
    if(nBeta > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x \%*\% beta)
    } else {
        eta <- rep(0, nObs)
    }
    cumpr <- matrix(pfun(matrix(zeta, nObs, nInterval - 1, byrow=TRUE) - eta), , nInterval - 1)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                fitted.values = fitted, lev = lev, terms = Terms,
                df.residual = sum(wt) - nBeta - nInterval + 1, edf = nBeta + nInterval - 1, n = sum(wt),
                nobs = sum(wt),
                call = match.call(), method = method,
		convergence = returnMessage(res), niter = niter)
    if(Hess) {
        dn <- c(names(beta), names(zeta))
        H <- hessian(res)
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if(model) fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "polr"
    fit
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
