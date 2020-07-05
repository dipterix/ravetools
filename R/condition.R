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
#' @param style styles to apply to text; choices are \code{debug}, \code{info},
#' \code{warn}, \code{error}, \code{fatal}, \code{default}, or a function
#' @param expr expression to evaluate
#' @param debug,info,warn,error,fatal functions to capture signals
#' @return \code{rave_error} returns \code{NULL}, others return
#' generated condition.
NULL


# message style presets
STYLE_PRESETS <- list(
  'fast' = cli::make_ansi_style('grey60'),
  'normal' = cli::make_ansi_style('blue'),
  'slow' = cli::make_ansi_style('orange'),
  'superslow' = cli::make_ansi_style('red'),
  'deadlyslow' = cli::make_ansi_style('#763053', bg = TRUE),

  'debug' = cli::make_ansi_style('grey60'),
  'info' = cli::make_ansi_style('#1d9f34'),
  'warn' = cli::make_ansi_style('#ec942c'),
  'error' = cli::make_ansi_style('#f02c2c'),
  'fatal' = cli::make_ansi_style('#763053', bg = TRUE),
  'default' = cli::make_ansi_style('grey20'),

  'grey60' = cli::make_ansi_style('grey60'),
  'grey100' = cli::make_ansi_style('grey100')
)

style_speed <- function(msg, speed, rule){
  stopifnot(length(rule) == 5)
  rule[[5]] <- Inf
  idx <- which(speed < rule)[[1]]

  if(idx == 5){
    STYLE_PRESETS$deadlyslow(STYLE_PRESETS$grey100(msg))
  } else {
    key <- c('fast', 'normal', 'slow', 'superslow', 'deadlyslow')[[idx]]
    STYLE_PRESETS[[key]](msg)
  }
}


#' @rdname rave-cond
#' @export
rave_condition <- local({
  last_fired <- NULL
  function(..., class = NULL, call = NULL, immediate. = TRUE, style = NULL){
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
    now <- Sys.time()
    if(is.null(last_fired)){
      delta <- 0
    } else {
      delta <- dipsaus::time_delta(last_fired, now, units = 'secs')
    }
    last_fired <<- now

    gl_msg <- gl(...)

    if(length(style)){
      if(is.function(style)){
        gl_msg <- style(gl_msg)
      } else if (style %in% c('fatal', 'deadlyslow')) {
        gl_msg <- STYLE_PRESETS[[style]](STYLE_PRESETS$grey100(gl_msg))
      } else if (style %in% names(STYLE_PRESETS)) {
        gl_msg <- STYLE_PRESETS[[style]](gl_msg)
      }
    }

    msg = STYLE_PRESETS$grey60(sprintf(
      '[%s,%s][%s]%s',
      strftime(now, '%H:%M:%S'),
      style_speed((sprintf('+%.3fs', delta)), delta, c(0.1, 1, 5, 10, Inf)),
      ctx$context,
      module_id
    ))


    cond <- simpleCondition(paste(msg, gl_msg), call = call)
    class(cond) <- c(class, 'raveCondition', class(cond))
    if(immediate.){
      signalCondition(cond)
    }
    invisible(cond)
  }
})

#' @rdname rave-cond
#' @export
rave_warn <- function(..., class = NULL, call = NULL, immediate. = TRUE, style = 'warn'){
  with_rave_handlers({
    cond <-
      rave_condition(
        ...,
        .envir = parent.frame(),
        class = c(class, 'raveWarning'),
        call = call,
        immediate. = immediate.,
        style = style
      )
  })

  invisible(cond)
}

#' @rdname rave-cond
#' @export
rave_error <- function(..., class = NULL, call = NULL, immediate. = TRUE, style = 'error'){
  with_rave_handlers({
    cond <-
      rave_condition(
        ...,
        .envir = parent.frame(),
        class = c(class, 'raveError'),
        call = call,
        immediate. = immediate.,
        style = style
      )
  })
  invisible(cond)
}

#' @rdname rave-cond
#' @export
rave_fatal <- function(..., class = NULL, call = NULL, immediate. = TRUE, style = 'fatal'){
  with_rave_handlers({
    cond <-
      rave_condition(
        ...,
        .envir = parent.frame(),
        class = c(class, 'raveFatal'),
        call = call,
        immediate. = immediate.,
        style = style
      )
  })
  invisible(cond)
}


#' @rdname rave-cond
#' @export
rave_info <- function(..., class = NULL, call = NULL, immediate. = TRUE, style = 'info'){
  with_rave_handlers({
    cond <-
      rave_condition(
        ...,
        .envir = parent.frame(),
        class = c(class, 'raveInfo'),
        call = call,
        immediate. = immediate.,
        style = style
      )
  })
  invisible(cond)
}


#' @rdname rave-cond
#' @export
rave_debug <- function(..., class = NULL, call = NULL, immediate. = TRUE, style = 'debug'){
  with_rave_handlers({
    cond <-
      rave_condition(
        ...,
        .envir = parent.frame(),
        class = c(class, 'raveDebug'),
        call = call,
        immediate. = immediate.,
        style = style
      )
  })
  invisible(cond)
}


default_debug <- function(e){
  cat(e$message, '\n')
}
default_info <- function(e){
  cat(e$message, '\n')
}
default_warn <- function(e){
  cat(e$message, '\n')
}
default_error <- function(e){
  cat(e$message, '\n')
  if(!shiny_is_running()){
    stop(e)
  }
}
default_fatal <- function(e){
  cat(e$message, '\n')
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



