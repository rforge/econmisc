lmc <- function(formula, strata, data, subset, weights, na.action, ...)
{
   ## Linear model with clustering.  Extracts strata and runs 'lm'
   ##
   ## call 'lm' with similar arguments, except 'strata'
   cl <- match.call()
   is <- match("strata", names(cl), 0)
   cl[[1]] <- as.name("lm")
   cl[[is]] <- NULL
   ##
   res <- eval(cl, parent.frame())
   if(missing(strata)) {
      return(res)
   }
   mf <- match.call(expand.dots=FALSE)
   m <- match(c("formula", "strata", "data", "subset", "weights", "na.action"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   res$strata <- model.extract(mf, "strata")
   class(res) <- c("lmc", class(res))
   res
   }
