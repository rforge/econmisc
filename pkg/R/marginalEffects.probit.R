
marginalEffects.probit <- function(object, X="all", ...) {
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
   class(me) <- c("marginalEffects.probit", "marginalEffects", class(object))
   class(me) <- class(me)[-which(class(me) == "probit")]
                           # it is _not_ 'probit' object any more
   me
}
