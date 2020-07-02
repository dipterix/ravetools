

CTX_STRING = c("default", "rave_module_debug", "rave_compile", "rave_running", "rave_running_local")
CTX_STRING_LEVEL = list(1,2,3,4,5)
names(CTX_STRING_LEVEL) <- CTX_STRING



#' Internal function to compare contexts
#' @param ctx,target context string
#' @param strict whether \code{ctx} equaling to \code{target} is
#' accepted
#' @return Whether \code{ctx} is greater or equal than \code{target}
compare_context <- function(ctx, target, strict = TRUE){
  if(strict){
    CTX_STRING_LEVEL[[ctx]] > CTX_STRING_LEVEL[[target]]
  } else {
    CTX_STRING_LEVEL[[ctx]] >= CTX_STRING_LEVEL[[target]]
  }
}


context_env <- new.env(parent = emptyenv())
context_env$current_context <- 'default'
context_env$package_name <- NULL
context_env$module_id <- NULL
context_env$running_instance <- NULL

set_context <- function(context = CTX_STRING, package, module, instance, force_set = FALSE){
  context = match.arg(context)


  # check if this is under shiny session and rave is running
  if(context_env$current_context == 'rave_running'){
    if(requireNamespace('shiny', quietly = TRUE)){
      if(shiny_is_running()){
        registered_instance <- add_to_session(session = shiny::getDefaultReactiveDomain(),
                                              key = 'rave_instance', val = instance, override = FALSE)

        if(!force_set && !identical(registered_instance, instance)){
          rave_fatal("[ERROR] Trying to set rave_running context with different instance")
          return()
        }
      }
      if(!force_set){
        # we don't set any context in this case
        return()
      }
    }
  }

  switch(
    context,
    'default' = {
      context_env$current_context <- context
      context_env$package_name <- NULL
      context_env$module_id <- NULL
      context_env$running_instance <- NULL
    },
    'rave_module_debug' = {
      force(package)
      force(module)
      context_env$package_name <- package
      context_env$module_id <- module
      context_env$running_instance <- NULL
      context_env$current_context <- context
    },
    'rave_running' = {
      force(package)
      force(module)
      force(instance)
      context_env$package_name <- package
      context_env$module_id <- module
      context_env$running_instance <- instance
      context_env$current_context <- context
      if(!do.call('requireNamespace', list(package = package))){
        rave_warn("Package {package} has not been installed.")
      }
      if(requireNamespace('shiny', quietly = TRUE)){
        if(shiny_is_running()){
          registered_instance <- add_to_session(
            session = shiny::getDefaultReactiveDomain(),
            key = 'rave_instance', val = instance, override = TRUE)
        }
      }
    },
    {
      force(package)
      force(module)
      force(instance)
      context_env$package_name <- package
      context_env$module_id <- module
      context_env$running_instance <- instance
      context_env$current_context <- context
      if(!do.call('requireNamespace', list(package = package))){
        rave_warn("Package {package} has not been installed.")
      }
    }
  )
  invisible()
}

get_context <- function(){
  if(context_env$current_context == 'rave_running'){
    if(requireNamespace('shiny', quietly = TRUE)){
      session_instance <- add_to_session(
        session = shiny::getDefaultReactiveDomain(),
        key = 'rave_instance', val = NULL, override = FALSE)
      if(inherits(session_instance, 'RAVEContainer')){
        return(structure(list(context = 'rave_running',
                              package = session_instance$package,
                              module_id = session_instance$module_id,
                              instance = session_instance),
                         class = 'rave-context'))
      }
    }
  }
  structure(list(context = context_env$current_context,
                 package = context_env$package_name,
                 module_id = context_env$module_id,
                 instance = context_env$running_instance),
            class = 'rave-context')
}




#' Set or Get RAVE context
#' @name rave-context
#' @param context character or missing. If not missing, will set context
#' before returning current context
#' @param package package name, for example, \code{"ravebuiltins"},
#' this is required for \code{rave_module_debug} and above
#' @param module_id,instance RAVE module ID and instance object,
#' required for context \code{rave_running} and \code{rave_running_local}.
#' @param info which context information to obtain
#' @param .force internally used, do not use
#' @return \code{rave_context} returns a \code{rave-context} object,
#' i.e. a list. \code{from_rave_context} returns corresponding
#' context information.
#' @details There are four contexts in RAVE: \code{default},
#' \code{rave_module_debug}, \code{rave_running}, and
#' \code{rave_running_local}. These four contexts are ordered, with
#' \code{default} being the lowest and \code{rave_running_local} the
#' highest. The order only has internal meaning and users don't have to
#' remember it. However, the four contexts correspond to different
#' scenarios.
#'
#' When creating a 'S3' generic function (see
#' \code{\link{rave_context_generics}}), a \code{default}
#' function must be provided to handle the default action of the
#' function.
#'
#' \code{rave_module_debug} is for module developers to use
#' as it minimizes many interactive procedures and dynamically load
#' modules to simplify debugging process. Normally this context is
#' only used in development mode.
#'
#' \code{rave_running} means when start-RAVE command is called, and a
#' shiny app launches within the browser. Many functions might change
#' their behaviors here. For example, data loading process might need
#' to show a modal dialogue instead of asking in the console.
#'
#' \code{rave_running_local} means when module is running without
#' shiny. In this case, and shiny interactive contents need to be
#' switched back to console-based. What differs
#' \code{rave_running_local} and \code{default} is in the prior
#' context, the module has more information about the state and
#' one can always assume the data has been loaded and use them without
#' validating them. Under \code{default} context, nothing can be
#' assumed.
#'
#' @seealso \code{\link{rave_context_generics}}
NULL

#' @rdname rave-context
#' @export
rave_context <- function(context, package = NULL, module_id = NULL, instance = NULL, .force = FALSE){
  if(!missing(context)){
    rave_debug("Setting context to {context}")
    stopifnot2(context %in% CTX_STRING,
               msg = "Invalid context string.")
    set_context(context, package, module_id, instance, force_set = .force)
  }
  return(get_context())
}

#' @rdname rave-context
#' @export
from_rave_context <- function(info = c('context', 'package', 'module_id', 'instance')){
  info <- match.arg(info)
  ctx <- rave_context()
  ctx[[info]]
}


#' @export
with_rave_context <- function(context, expr, package = NULL, module_id = NULL, instance = NULL,
                              quoted = FALSE, env = parent.frame()){
  if(!quoted){
    expr <- substitute(expr)
  }
  if(!missing(context)){
    rave_debug("Setting context to {context}")
    stopifnot2(context %in% CTX_STRING, msg = "Invalid context string.")
  }
  ctx <- get_context()
  on.exit({
    set_context(ctx$context, ctx$package, ctx$module_id, ctx$instance, force_set = TRUE)
  }, add = FALSE)
  set_context(context, package, module_id, instance, force_set = TRUE)
  results <- eval(expr, envir = env)

  set_context(ctx$context, ctx$package, ctx$module_id, ctx$instance, force_set = TRUE)
  on.exit({}, add = FALSE)
  results
}


#' @export
`print.rave-context` <- function(x, ...){
  cat("<RAVE Context> ", x$context, "\n")
  cat("  Package:\t", deparse1(x$package), '\n')
  cat("  Module:\t", deparse1(x$module_id), '\n')
  cat("  Running instance:\n\t")
  base::print(x$instance)
  invisible(x)
}


#' @title Create S3 Generic Functions that Respects RAVE Context
#' @param fun_name generic function name
#' @param args list of arguments to the function. Use \code{\link{alist}}
#' @param env in which environment the function is declared.
#' @return A generic function.
#' @seealso \code{\link{rave-context}}
#'
#' @examples
#'
#'
#' # Make a generic function that reacts under different contexts
#' # The goal is to ask user a comple question. Under debug mode,
#' # Ignore the question and always answer "Yes"
#' ask_for_input <- rave_context_generics(
#'   'ask_for_input', args = alist(msg =, ...=))
#'
#' # write default action. Since under default context, shiny is not running
#' ask_for_input.default <- function(msg, ...){
#'   utils::askYesNo(msg)
#' }
#'
#' ask_for_input.rave_module_debug <- function(msg, ...){
#'   # Do not prompt questions and always answer yes
#'   return(TRUE)
#' }
#'
#' if(interactive()){
#'   rave_context('default')
#'
#'   # will prompt for answer
#'   ask_for_input('Answer yes or no')
#' }
#'
#' # When debugging modules
#' rave_context('rave_module_debug', package = 'base')
#' ask_for_input('Answer yes or no')  # always returns TRUE
#'
#' @export
rave_context_generics <- function(fun_name, args = alist(), env = parent.frame()){

  stopifnot2(is.character(fun_name), msg = 'fun_name must be characters')

  dipsaus::new_function2(args, {
    UseMethod(
      !!fun_name,
      object = structure(
        'context', class = raveutils::from_rave_context('context')
      )
    )
  }, env)

}
