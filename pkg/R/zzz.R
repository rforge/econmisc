.First.lib <- function(libname, pkgname) {
   library.dynam("econMisc", pkgname)
                           # for some reason it does not work with NAMESPACE directive
}
