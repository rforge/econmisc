vcov.lmc <- function(object,
                       ...) {
   vcov <- vcovClust(object, 1, object$strata)
   vcov
}
