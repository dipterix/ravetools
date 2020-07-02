#' @export
rave_module_package <- function(){
  if(from_rave_context('context') == 'default'){
    rave_fatal("Cannot call 'rave_module_package' from default context")
  }
  from_rave_context('package')
}


#' @export
rave_module_root_directory <- function(){
  d = NULL
  if( verify_rstudio_version() ){
    d = rstudioapi::getActiveProject()
  }
  pkgname = rave_module_package()

  if(length(d) == 1 && grepl(paste0('/', pkgname, '$'), d)){
    return(d)
  }else{
    # package user
    return(system.file('', package = pkgname))
  }

}


