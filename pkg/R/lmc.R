lmc <- function(strata, ...)
{
   ## Linear model with clustering.  Extracts strata and runs 'lm'
   ##
   if(missing(strata))
       return(lm(...))
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match("strata", names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(strata, parent.frame())
   res <- lm(...)
   res$strata <- strata
   class(res) <- c("lmc", class(res))
   res
   }
