
### std: extract standard errors
std <- function(x, ...)
    UseMethod("std")

std.maxLik <- function(x, ...) {
   sqrt(diag(vcov(x)))
}
