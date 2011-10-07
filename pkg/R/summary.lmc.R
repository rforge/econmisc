summary.lmc <- function(object, correlation=FALSE,
                       ...) {
   ## Overwritten summary.lm in order to easily print clustering-robust results
   ##
   ## object: a lm-type of object
   if(is.null(object$strata))
       return(stats::summary.lm(object, ...))
   vcov <- vcovClust(object, 1, object$strata)
   ## The rest is borrowed from 'stats::summary.lm'
   z <- object
    p <- z$rank
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/(n - p)
        ans <- z[c("call", "terms")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0L, 4L)
        dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
            "Std. Error", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        return(ans)
    }
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr)) 
        stop("invalid 'lm' object:  no 'terms' nor 'qr' component")
    n <- NROW(Qr$qr)
    rdf <- n - p
    if (is.na(z$df.residual) || rdf != z$df.residual) 
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1L:p
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept")) 
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
##     R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
##     se <- sqrt(diag(R) * resvar)
   se <- sqrt(diag(vcov))
   est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
        rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept")) 
            1L
        else 0L
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- vcov
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
        1)]
    if (correlation) {
        ans$correlation <- vcov/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
#        ans$symbolic.cor <- symbolic.cor
    }
    if (!is.null(z$na.action)) 
        ans$na.action <- z$na.action
    class(ans) <- "summary.lm"
    ans
}
