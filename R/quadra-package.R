#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib quadra, .registration=TRUE
## usethis namespace: end
.onUnload <- function(libpath) {
  library.dynam.unload("quadra", libpath)
}
