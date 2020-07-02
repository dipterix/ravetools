SHINY_UPDATE_FUNCTIONS = dipsaus::fastmap2()
SHINY_UPDATE_FUNCTIONS$shiny = list(
  'selectInput' = 'shiny::updateSelectInput'
)
SHINY_UPDATE_FUNCTIONS$dipsaus = list(
  'actionButtonStyled' = 'dipsaus::updateActionButtonStyled',
  'compoundInput2' = 'dipsaus::updateCompoundInput2'
)
SHINY_UPDATE_FUNCTIONS$default = c(
  SHINY_UPDATE_FUNCTIONS$shiny,
  SHINY_UPDATE_FUNCTIONS$dipsaus,
  SHINY_UPDATE_FUNCTIONS$rave,
  SHINY_UPDATE_FUNCTIONS$ravecore,
  SHINY_UPDATE_FUNCTIONS$threeBrain
)

#' @export
guess_shiny_update <- function(call, parse = TRUE){
  # call <- quote(shiny::textInput('asda'))

  funname = call[[1]]
  pkgname = NA

  if(!is.name(funname)){
    # check if first is :: or :::
    if(identical(funname[[1]], quote(`::`)) || identical(funname[[1]], quote(`:::`))){
      pkgname = as.character(funname[[2]])
      # print(pkgname)
      funname = funname[[3]]
    }
  }
  if(!is.name(funname)){
    rave_error('Cannot find shiny update function for function {sQuote(funname)}')
    return(NULL)
  }

  funname <- as.character(funname)
  if(is.na(pkgname)){
    pkgname = 'default'
  }
  update_str = SHINY_UPDATE_FUNCTIONS[[pkgname]][[funname]]

  if(is.null(update_str)){
    update_str = paste0('update', stringr::str_to_upper(stringr::str_sub(funname, end = 1)), stringr::str_sub(funname, start = 2))
    if(pkgname != 'default'){
      fun <- asNamespace(pkgname)[[update_str]]
      if(!is.function(fun)){
        rave_error("Cannot find shiny update function for function {sQuote(funname)}")
      }
      update_str = sprintf('%s:::%s', pkgname, update_str)
    }
  }

  if(parse){
    eval(parse(text = update_str))
  } else{
    update_str
  }
}




SHINY_OUTPUT_FUNCTIONS = dipsaus::fastmap2()
SHINY_OUTPUT_FUNCTIONS$shiny = list(
  'htmlOutput' = 'shiny::renderUI',
  'verbatimTextOutput' = 'shiny::renderPrint',
  'plotOutput' = 'shiny::renderPlot'
)
SHINY_OUTPUT_FUNCTIONS$threeBrain = list(
  'threejsBrainOutput' = 'threeBrain::renderBrain'
)
SHINY_OUTPUT_FUNCTIONS$rave = list(
  'customizedUI' = 'shiny::renderUI'
)
SHINY_OUTPUT_FUNCTIONS$ravecore = SHINY_OUTPUT_FUNCTIONS$rave
SHINY_OUTPUT_FUNCTIONS$DT = list(
  'DTOutput' = 'DT::renderDT',
  'dataTableOutput' = 'DT::renderDataTable'
)
SHINY_OUTPUT_FUNCTIONS$default = c(
  SHINY_OUTPUT_FUNCTIONS$shiny,
  SHINY_OUTPUT_FUNCTIONS$dipsaus,
  SHINY_OUTPUT_FUNCTIONS$rave,
  SHINY_OUTPUT_FUNCTIONS$ravecore,
  SHINY_OUTPUT_FUNCTIONS$threeBrain
)


#' @export
guess_shiny_output <- function(call, parse = TRUE){
  # call <- quote(shiny::verbatimTextOutput('asd'))

  funname = call[[1]]
  pkgname = NA

  if(!is.name(funname)){
    # check if first is :: or :::
    if(identical(funname[[1]], quote(`::`)) || identical(funname[[1]], quote(`:::`))){
      pkgname = as.character(funname[[2]])
      # print(pkgname)
      funname = funname[[3]]
    }
  }
  if(!is.name(funname)){
    rave_error('Cannot find shiny output function for call {sQuote(funname)}')
    return(NULL)
  }

  funname <- as.character(funname)
  if(is.na(pkgname)){
    pkgname = 'default'
  }
  update_str = SHINY_OUTPUT_FUNCTIONS[[pkgname]][[funname]]

  if(is.null(update_str)){
    update_str = paste0('render', stringr::str_to_upper(stringr::str_sub(funname, end = 1)), stringr::str_sub(funname, start = 2))
    update_str = stringr::str_remove(update_str, '[oO]utput$')
    if(pkgname != 'default'){
      fun <- asNamespace(pkgname)[[update_str]]
      if(!is.function(fun)){
        rave_error("Cannot find shiny output function for name {sQuote(funname)}")
      }
      update_str = sprintf('%s:::%s', pkgname, update_str)
    }
  }

  if(parse){
    eval(parse(text = update_str))
  } else{
    update_str
  }
}
