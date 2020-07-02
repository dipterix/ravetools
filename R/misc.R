#' A dummy function that literally does nothing
#' @param ... nothing
#' @export
do_nothing <- function(...){

}

#' Function to clear all elements within environment
#'
#' @param env environment to clean, can be R environment, or
#' \code{\link[dipsaus]{fastmap2}}
#' @param ... ignored
#'
#' @examples
#' \dontrun{
#' env = new.env()
#' env$a = 1
#' print(as.list(env))
#'
#' clear_env(env)
#' print(as.list(env))
#' }
#' @export
clear_env <- function(env, ...){
  if(is.environment(env)){
    if(environmentIsLocked(env)){
      return(invisible())
    }
    nms = names(env)
    nms = nms[!stringr::str_detect(nms, '^\\.__rave')]
    if(isNamespace(env)){
      nms = nms[!nms %in% c(".__NAMESPACE__.", ".__S3MethodsTable__.")]
    }
    rm(list = nms, envir = env)
  } else if(inherits(env, 'fastmap2')){
    .subset2(env, 'remove')(names(env))
  }
  return(invisible())
}


#' @export
add_to_session <- function(
  session, key = 'rave_id',
  val = paste(sample(c(letters, LETTERS, 0:9), 20), collapse = ''),
  override = FALSE
){
  if(!is.null(session)){
    if(override || !exists(key, envir = session$userData)){
      assign(key, val, envir = session$userData)
    }
    return(get(key, envir = session$userData))
  }
  return(NULL)
}



