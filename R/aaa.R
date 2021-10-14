## usethis namespace: start
#' @importFrom stats approx
#' @importFrom Rcpp sourceCpp
#' @importFrom dipsaus %?<-%
#' @useDynLib raveutils, .registration = TRUE
## usethis namespace: end
NULL

deparse1 <- function(expr, collapse = ' '){
  paste(deparse(expr), collapse = collapse)
}


stopifnot2 <- function(..., msg = 'Condition not satisfied'){
  if(!all(c(...))){
    stop(msg)
  }
}

tempfile2 <- function(
  pattern = "ravetmp-", tmpdir = file.path(tempdir(check = TRUE), "raveutils"), fileext = ""){
  if(!dir.exists(tmpdir)){
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  }
  tempfile(pattern = pattern, tmpdir = tmpdir, fileext = fileext)
}
