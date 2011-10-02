
summary.marginalEffects.probit <- function(object, ...) {
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
   class(a) <- c("summary.marginalEffects.probit", class(a))
   a
}
