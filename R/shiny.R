#' Register JavaScript Libraries for RAVE-shiny applications
#' @param rmd whether used in R-markdown files
#' @return None
#' @export
register_js <- function(rmd = FALSE){
  if(requireNamespace('shinyalert', quietly = TRUE)){
    shiny::tagList(
      shinyalert::useShinyalert(rmd = rmd),
      shinyjs::useShinyjs(rmd = rmd)
    )
  } else {
    shinyjs::useShinyjs()
  }
}

#' Show Alert in Shiny apps
#' @param title,message title and message of the alert
#' @param type options are \code{"warning"}, \code{"error"},
#' \code{"success"}, and \code{"info"}
#' @param onConfirm R function to call upon confirmation
#' @param closeOnEsc,closeOnClickOutside whether to easy close the modal
#' @param showCancelButton whether to show dismiss button
#' @param confirmButtonText text of confirm button
#' @param html whether not to escape the message
#' @param inputId input ID
#' @param session shiny session
#' @param fancy whether try to use fancy version
#' @export
safe_alert <- function(
  title, message, type = c("warning", "error", "success", "info"),
  onConfirm = NULL, inputId = 'shinyalert',
  closeOnEsc = FALSE, closeOnClickOutside = FALSE, html = FALSE,
  showCancelButton = FALSE, confirmButtonText = 'OK',
  session = shiny::getDefaultReactiveDomain(),
  fancy = TRUE
){
  type <- match.arg(type)
  message <- as.character(message)
  if(shiny_is_running()){
    shiny::withReactiveDomain(session, {
      if(fancy && requireNamespace('shinyalert', quietly = TRUE)){

        if(is.function(onConfirm)){
          callbackR <- function(value){
            on.exit(shinyalert::dismissalert(), add = TRUE, after = TRUE)
            onConfirm(value)
          }
        } else {
          callbackR <- NULL
        }

        shinyalert::shinyalert(
          title = title, text = message, closeOnEsc = closeOnEsc,
          closeOnClickOutside = closeOnClickOutside, html = html, type = type,
          showCancelButton = FALSE, callbackR = callbackR, inputId = inputId,
          confirmButtonText = confirmButtonText, confirmButtonCol = "#AEDEF4"
        )
        if(is.function(callbackR)){
          shinyalert::shinyalert(
            title = 'Running',
            text = 'Executing code, please wait. (This message will be dismissed once finished)',
            closeOnClickOutside = FALSE, showConfirmButton = FALSE, closeOnEsc = TRUE
          )
        }
      } else {
        cancelbtn <- NULL
        if(showCancelButton){
          cancelbtn <- shiny::modalButton('Cancel')
        }
        confirmButtonText
        shiny::showModal(shiny::modalDialog(
          title = title,
          ifelse(isTRUE(html), shiny::HTML(message), message),
          size = 's', easyClose = closeOnClickOutside,
          footer = shiny::tagList(
            dipsaus::actionButtonStyled(session$ns(inputId), confirmButtonText),
            cancelbtn
          )
        ))
        shiny::observeEvent(session$input[[inputId]], {
          dipsaus::updateActionButtonStyled(session, inputId, disabled = TRUE)
          if(is.function(onConfirm)){
            try(onConfirm())
          }
          shiny::removeModal()
        }, once = TRUE, ignoreInit = TRUE, ignoreNULL = TRUE)
      }
    })
  } else {
    fname = list(
      "warning" = 'rave_warn',
      "error" = 'rave_error',
      "success" = 'rave_info',
      "info" = 'rave_debug'
    )[[type]]
    do.call(fname, list(sprintf('%s: %s', title, message)))
  }
  return(invisible(NULL))
}




#' Wrappers to be compatible with Input ID
#' @param inputId character, passed to \code{define_input}
#' @param label label character
#' @param expr expression returning 'html' code
#' @param ... passed to other methods
#' @return \code{plainUI} returns evaluated \code{expr}; \code{plainLabel}
#' returns a label-only input
#' @name customized-ui
NULL

#' @rdname customized-ui
#' @export
plainUI <- function(inputId, expr){
  expr
}

#' @rdname customized-ui
#' @export
plainLabel <- function(inputId, label, ...){
  shiny::div(
    class = "form-group shiny-input-container",
    shiny::tags$label(label, class = "control-label", `for` = inputId),
    shiny::tags$input(id = inputId, type = "text", class = "form-control",
                      style = 'display:none!important;',
                      value = '', placeholder = ''),
    ...
  )
}

SHINY_UPDATE_FUNCTIONS = dipsaus::fastmap2()
SHINY_UPDATE_FUNCTIONS$shiny = list(
  'selectInput' = 'shiny::updateSelectInput'
)
SHINY_UPDATE_FUNCTIONS$dipsaus = list(
  'actionButtonStyled' = 'dipsaus::updateActionButtonStyled',
  'compoundInput2' = 'dipsaus::updateCompoundInput2'
)
SHINY_UPDATE_FUNCTIONS$raveutils = list(
  'plainLabel' = 'shiny::updateTextInput',
  'plainUI' = 'raveutils::do_nothing'
)

SHINY_UPDATE_FUNCTIONS$default = c(
  SHINY_UPDATE_FUNCTIONS$shiny,
  SHINY_UPDATE_FUNCTIONS$dipsaus,
  SHINY_UPDATE_FUNCTIONS$rave,
  SHINY_UPDATE_FUNCTIONS$ravecore,
  SHINY_UPDATE_FUNCTIONS$threeBrain,
  SHINY_UPDATE_FUNCTIONS$raveutils
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
SHINY_OUTPUT_FUNCTIONS$ravecore = list(
  'customizedUI' = 'shiny::renderUI'
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
