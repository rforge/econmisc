summary.intreg <- function(object, ...) {
   estimate <- coefTable(coef(object), std(object), object$param$df)
   if(!all(activePar(object)[object$param$index$boundary])) {
      ## we have not estimated the boundaries -> do no print them
      estimate <- estimate[-object$param$index$boundary,,drop=FALSE]
   }
   library(maxLik)
   s <- maxLik:::summary.maxLik(object)
   s$estimate <- estimate
   class(s) <- c("summary.intreg", class(s))
   return(s)
}

print.summary.intreg <- function(x,
                                 digits=max(3, getOption("digits") - 3),
                                 ...) {
   cat("--------------------------------------------\n")
   cat("Interval regression\n")
   cat( "Maximum Likelihood estimation\n" )
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
         cat("Log-Likelihood:", logLik(x), "\n")
      }
   if(!is.null(x$estimate)) {
      cat( x$param$nObs, "observations, " )
      cat( x$param$nParam, "free parameters" )
      cat( " (df = ", x$param$df, ")\n", sep="")
      printCoefmat( x$estimate, signif.legend = TRUE, digits = digits )
   }
   cat("--------------------------------------------\n")
   invisible( x )
}
   
