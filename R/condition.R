#' @title Check States of Current RAVE instances or 'shiny' instances
#' @return Logical, for \code{shiny_is_running}, the result
#' indicates whether shiny is running. For \code{rave_is_running},
#' the result shows whether rave has
#' context higher/equaling to \code{'rave_running'}.
#' @name rave-states
NULL

#' @rdname rave-states
#' @export
shiny_is_running <- function(){
  if(requireNamespace('shiny', quietly = TRUE)){
    return(isTRUE(!is.null(shiny::getDefaultReactiveDomain())))
  }
  return(FALSE)
}

#' @rdname rave-states
#' @export
rave_is_running <- function(){
  isTRUE(compare_context(from_rave_context('context'), 'rave_running', strict = FALSE))
}

#' @name rave-cond
#' @title RAVE Internal condition classes
#' @param ... messages of \code{\link{condition}},
#' \code{\link[glue]{glue}} is internally supported. For
#' \code{with_rave_handlers}, passed to
#' \code{\link{withCallingHandlers}}
#' @param class condition class
#' @param call call expression
#' @param immediate. whether to fire condition immediately
#' @param expr expression to evaluate
#' @param debug,info,warn,error,fatal functions to capture signals
#' @return \code{rave_error} returns \code{NULL}, others return
#' generated condition.
NULL

#' @rdname rave-cond
#' @export
rave_condition <- function(..., class = NULL, call = NULL, immediate. = TRUE){
  if(is.null(call)){
    call <- tryCatch({
      eval(quote(match.call()), envir = parent.frame())
    }, error = function(e){
      NULL
    })
  }
  ctx <- rave_context()
  module_id <- ctx$module_id
  if(!length(module_id)){
    module_id <- ''
  } else {
    module_id <- sprintf('[%s]', module_id)
  }
  msg = gl(sprintf('[%s]%s ', ctx$context, module_id), ...)
  cond <- simpleCondition(msg, call = call)
  class(cond) <- c(class, 'raveCondition', class(cond))
  if(immediate.){
    signalCondition(cond)
  }
  invisible(cond)
}

#' @rdname rave-cond
#' @export
rave_warn <- function(..., class = NULL, call = NULL, immediate. = TRUE){
  with_rave_handlers({
    cond <- rave_condition(..., .envir = parent.frame(), class = c(class, 'raveWarning'), call = call, immediate. = immediate.)
  })

  invisible(cond)
}

#' @rdname rave-cond
#' @export
rave_error <- function(..., class = NULL, call = NULL, immediate. = TRUE){
  with_rave_handlers({
    cond <- rave_condition(..., .envir = parent.frame(), class = c(class, 'raveError'), call = call, immediate. = immediate.)
  })
  invisible(cond)
}

#' @rdname rave-cond
#' @export
rave_fatal <- function(..., class = NULL, call = NULL, immediate. = TRUE){
  with_rave_handlers({
    cond <- rave_condition(..., .envir = parent.frame(), class = c(class, 'raveFatal'), call = call, immediate. = immediate.)
  })
  invisible(cond)
}


#' @rdname rave-cond
#' @export
rave_info <- function(..., class = NULL, call = NULL, immediate. = TRUE){
  with_rave_handlers({
    cond <- rave_condition(..., .envir = parent.frame(), class = c(class, 'raveInfo'), call = call, immediate. = immediate.)
  })
  invisible(cond)
}


#' @rdname rave-cond
#' @export
rave_debug <- function(..., class = NULL, call = NULL, immediate. = TRUE){
  with_rave_handlers({
    cond <- rave_condition(..., .envir = parent.frame(), class = c(class, 'raveDebug'), call = call, immediate. = immediate.)
  })
  invisible(cond)
}


default_debug <- function(e){
  dipsaus::cat2(e$message, level = 'DEBUG')
}
default_info <- function(e){
  dipsaus::cat2(e$message, level = 'INFO')
}
default_warn <- function(e){
  dipsaus::cat2(e$message, level = 'WARNING')
}
default_error <- function(e){
  dipsaus::cat2(e$message, level = 'ERROR')
  if(!shiny_is_running()){
    stop(e)
  }
}
default_fatal <- function(e){
  dipsaus::cat2(e$message, level = 'FATAL')
  stop(e)
}

#' @rdname rave-cond
#' @export
with_rave_handlers <- function(expr, debug, info, warn, error, fatal, ...){
  debug %?<-% default_debug
  info %?<-% default_info
  warn %?<-% default_warn
  error %?<-% default_error
  fatal %?<-% default_fatal

  withCallingHandlers(
    expr = expr,
    raveDebug = debug,
    raveInfo = info,
    raveWarning = warn,
    raveError = error,
    raveFatal = fatal,
    ...
  )
}



