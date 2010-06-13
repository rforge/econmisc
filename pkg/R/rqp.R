
incidence.csr <- function(s) {
   ## factor-based sparce incidence matrix.
   ## Note: we cannot use model.matrix() for large problems
   if(!is.factor(s))
       s <- as.factor(s)
   library(SparseM)
   new("matrix.csr", ra=rep(1, length(s)),
       ja=as.integer(s),
       ia=as.integer(c(1:length(s), length(s) + 1)),
       dimension=as.integer(c(length(s), length(levels(s))))
       )
}

rqp <- function(formula, subset, strata,
                     taus=(1:3)/4,
                w = taus*(1 - taus)/sum(taus*(1 - taus)),
                     lambda = 1e5,
                     method="fn",
                conf.int="boot", R=100,
                data=sys.frame(sys.parent())
                ) {
   ## Quantile regression for panel data (Koenker, 2004).
   ## conf.int   method for computing confidence intervals ("boot")
   ## R          # of bootstrap replications
   bootStat <- function(iIndividual, index) {
      ## for bootstrapping standard errors
      strata <- iIndividual[index]
      y1 <- y[s %in% strata]
      s1 <- s[s %in% strata]
      X1 <- X[s %in% strata,]
      res <- rq.fit.panel(X1, y1, s1, taus=taus, w=w, lambda=lambda)
                           # $coef: m * p coefficients + n fixed effects
      return(res$coef[iCoef])
   }
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "subset", "strata", "na.action",
                "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   eval(data)
                                        # we need to eval() data here, otherwise the evaluation of the
                                        # model frame will be wrong if called from inside a function
                                        # inside a function (sorry, I can't understand it either :-(
   mf <- eval(mf, envir=parent.frame())
   if (method == "model.frame")
       return(mf)
   else if (method != "fn")
       warning("method = ", method, " is not supported. Using \"fn\"")
   mt <- attr(mf, "terms")
   y <- model.response( mf )
   X <- model.matrix(mt, mf, contrasts)
   s <- as.integer(mf[["(strata)"]])
                           # we force strata to be integer for later
                           # indexing 
   if(is.null(s))
       s <- seq(along=y)
   m <- length(taus)
   p <- ncol(X)
   ##
   res <- rq.fit.panel(X, y, s, taus=taus, w=w, lambda=lambda)
                           # $coef: m * p coefficients + n fixed effects
   iCoef <- matrix(1:(m*p), p, m)
   iFE <- (m*p + 1):length(res$coef)
   coef <- matrix(res$coef[iCoef], p, m)
   row.names(coef) <- colnames(X)
   colnames(coef) <- taus
   if(conf.int=="boot") {
      library(boot)
      b <- boot(data=seq(length=length(unique(s))), bootStat, R=R)
   }
   else
       b <- NULL
   ret <- list(coef=coef,
               effects=res$coef[iFE],
               boot=b
               )
}

rq.fit.panel <- function(X,y,s, taus=(1:3)/4,
                         w = taus*(1 - taus)/sum(taus*(1 - taus)),
                         lambda = 1) {
   ## taus   requested quantiles
   ## w      weight at each quantile
   ## lambda penalization parameter (see Koenker, 2004)
   ## 
   ## See Koenker (2004), Journal of Multivariate Analysis
   ## author: Roger Koenker
   ## 
   ## prototype function for panel data fitting of QR models
   ## 
   ## the matrix X is assumed to contain an intercept
   ## 
   ## the vector s is a strata indicator assumed (so far) to be a one-way
   ## layout
   ## 
   ## NB:
   ## 
   ## 1. The value of the shrinkage parameter lambda is an open research
   ## 	problem in the simplest homogneous settings it should be the
   ## 	ratio of the scale parameters of the fixed effects and the
   ## 	idiocyncratic errors
   ## 
   ## 2. On return the coefficient vector has m*p + n elements where m
   ##	is the number quantiles being estimated, p is the number of
   ##	colums of X, and n is the number of distinct values of s.  The
   ##	first m*p coefficients are the slope estimates, and the last n
   ##	are the "fixed effects"
   ##
   ## 3. Like all shrinkage (regularization) estimators, asymptotic
   ##	inference is somewhat problematic... so the bootstrap is the
   ##	natural first resort.
   ##
   require(SparseM)
   require(quantreg)
   K <- length(w)
   if(K != length(taus))
       stop("length of w and taus must match")
   X <- as.matrix(X)
   p <- ncol(X)
   n <- length(levels(as.factor(s)))
   N <- length(y)
   if(N != length(s)) {
      stop("length of strata s (currently ", length(s),
           ") must match length of outcome y (currently ", length(y),
           ")", sep="")
   }
   if(N != nrow(X)) {
      stop("number of data rows (currently ", nrow(X),
           ") must match length of outcome y (currently ", length(y),
           ")", sep="")
   }
   Z <- incidence.csr(s)
                           # compressed sparse row matrix
   W <- as.matrix.csr(diag(w, nrow=length(w)))
   Fidelity <- cbind(W %x% X,w %x% Z)
                           # 
   Penalty <- cbind(as.matrix.csr(0,n,K*p),
                           # n x Kp matrix of 0-s
                    lambda*as(n,"matrix.diag.csr")
                           # n x n unit matrix
                    )
   D <- rbind(Fidelity,Penalty)
                           # design matrix (Koenker, 2004, p 78)
   y <- c(w %x% y,rep(0,n))
   a <- c((w*(1-taus)) %x% (t(X)%*%rep(1,N)),
          sum(w*(1-taus)) * (t(Z) %*% rep(1,N)) + lambda * rep(1,n))
   rq.fit.sfn(D, y, rhs=a)
}
