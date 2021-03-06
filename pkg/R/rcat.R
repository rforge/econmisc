
rcat <- function(n, prob, na.rm=TRUE) {
   ## Categorical random variables.
   ## n      number of variables to generate
   ## prob   probabilities.  Either NxM matrix -- each row corresponds to the probabilities
   ##        of the corresponding draw, or a vector -- equal probabilities for all the draws. 
   ##        The probabilities will be normalized to sum to unity by draws.
   ## na.rm  What to do with NA-s in 'prob'.
   ##        TRUE   remove (treat them as zeros)
   ##        FALSE  signal error
   ## 
   ## Generates a vector of length N of multinomial variables.
   ## The probabilities may differ for each component.
   ## The rows of probabilities should sum to less than one, runif
   ## will be split according to the intervals, generated by
   ## probabilities.
   if(is.null(dim(prob))) {
      prob <- matrix(prob, n, length(prob), byrow=TRUE)
   }
   if(!na.rm & any(is.na(prob)))
       stop("NA in argument 'prob'")
   prob[is.na(prob)] <- 0
   prob <- prob/rowSums(prob)
                           # rows with only NA-s are 0/0 = NaN
   storage.mode(prob) <- "double"
   r <- .Call("do_rcat", prob, na.rm, PACKAGE="econMisc")
   r
}    
