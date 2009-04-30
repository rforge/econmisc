
coef.marginalEffect <- function(x, ...)
    x$coefficients

marginalEffect <- function(object, ...)
    UseMethod("marginalEffect")

marginalEffect.probit <- function(object, X="all", ...) {
   ## object  object of class 'probit'
   ## X       explanatory variables
   ##         'all'    mean of individual effects on original explanatory variables
   ## mean of infinitesimal marginal effects
   x <- model.matrix(object)
   NX <- ncol(x)
   beta <- coef(object)
   Sigma <- vcov(object)
   effects <- mean(dnorm(x %*% beta))*beta
   C <- diag(mean(dnorm(x %*% beta)), nrow=NX, ncol=NX) -
       mean((x %*% beta) * dnorm(x %*% beta)^2)*(beta %*% t(beta))
   Sigma1 <- C %*% Sigma %*% t(C)
   me <- list(coefficients=effects, vcov=Sigma1, ML=object)
   class(me) <- c("marginalEffect.probit", "marginalEffect", class(object))
   me
}

print.marginalEffect <- function(x, ...)
    print(summary(x, ...))

print.summary.marginalEffect.probit <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Probit marginal effects/Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$nObs, " observations (", x$N0, " zeros and ", x$N1, " ones) and ",
          x$NActivePar, " free parameters (df = ",
          x$nObs - x$NActivePar, ")\n", sep="")
      cat("Estimates:\n")
      printCoefmat( x$estimate, ... )
   }
   cat("Significance test:\n")
   cat("chi2(", x$LRT$df, ") = ", x$LRT$LRT, " (p=", x$LRT$pchi2, ")\n", sep="")
   cat("--------------------------------------------\n")
}

summary.marginalEffect.probit <- function(object, ...) {
   ## summary for probit -- adds Likelihood Ratio Test to summary.maxLik
   coef <- coef(object)
   stdd <- sqrt(diag(vcov(object)))
   t <- coef/stdd
   p <- 2*pnorm( -abs( t))
   results <- cbind("Estimate"=coef, "Std. error"=stdd, "t value"=t, "Pr(> t)"=p)
   pchi2 <- pchisq(object$ML$LRT$LRT, object$ML$LRT$df, lower.tail=FALSE)
   a <- c(type=object$ML$type,
          iterations=object$ML$iterations,
          code=object$ML$code,
          message=object$ML$message,
          loglik=object$ML$maximum,
          estimate=list(results),
          vcov=list(vcov(object)),
          activePar=object$ML$activePar,
          NActivePar=sum(object$ML$activePar),
          LRT=list(c(object$ML$LRT, pchi2=pchi2)),
          nParam=object$ML$nParam,
          nObs=nObs(object),
          N0=object$ML$param$N0,
          N1=object$ML$param$N1,
          df=object$ML$df)
   class(a) <- c("summary.marginalEffect.probit", class(a))
   a
}

vcov.marginalEffect <- function(x, ...)
    x$vcov
