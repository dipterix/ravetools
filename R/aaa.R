#' @importFrom dipsaus %?<-%
NULL

gl <- function(..., .envir = parent.frame()){
  glue::glue(..., .envir = .envir)
}

#' @title Wrapper around \code{\link[glue]{glue}}
#' @param ...,.envir passed to \code{\link[glue]{glue}}
#' @export
rave_glue <- gl

#' @title Print colored messages
#' @param ...,.envir passed to \code{\link[glue]{glue}}
#' @param level passed to \code{\link[dipsaus]{cat2}}
#' @param .pal see \code{pal} in \code{\link[dipsaus]{cat2}}
#' @export
catgl <- function(..., .envir = parent.frame(), level = 'DEBUG', .pal){
  if(missing(.pal)){
    dipsaus::cat2(gl(..., .envir = .envir), level = level)
  }else{
    dipsaus::cat2(gl(..., .envir = .envir), level = level, pal = .pal)
  }
}


deparse1 <- function(expr, collapse = ' '){
  paste(deparse(expr), collapse = collapse)
}


stopifnot2 <- function(..., msg = 'Condition not satisfied', immediate. = TRUE){
  if(!all(c(...))){
    rave_fatal(msg, immediate. = immediate.)
  }
}



#' @export
test_farg <- function(fun, arg, dots = TRUE){
  stopifnot2(is.character(arg) || is.numeric(arg), msg = 'test_farg: arg must be either characters or integers')
  fm <- names(formals(fun))
  has_dots <- dots && ('...' %in% fm)
  if(has_dots){
    return(rep(TRUE, length(arg)))
  }
  if(is.character(arg)){
    arg %in% fm
  } else {
    arg <= length(fm)
  }
}