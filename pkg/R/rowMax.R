
rowMax <- function(mat, na.rm=TRUE) {
   ## Max of matrix elements by rows
   mat <- as.matrix(mat)
   if(!is.numeric(mat))
       stop("rowMax only works on _numeric_ matrices")
   r <- .Call("do_rowMax", mat, na.rm, PACKAGE="econMisc")
   r
}    
