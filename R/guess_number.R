# guess a upper-bound non-negative integer

#' Guess the max of non-negative numbers that pass the test
#' @param validator function that returns true if the number passes the test
#' or false if fails
#' @param initial initial value to guess, must be within \code{(min_v, max_v)}
#' @param min_v minimum number that passes the test
#' @param max_v maximum possible number
#' @returns The upper bound of the number that passes the test
#'
#' @examples
#'
#'
#' # 42
#' guess_max_integer(
#'   function(x) {
#'     x <= 42
#'   }
#' )
#'
#' # 42
#' guess_max_integer(
#'   function(x) {
#'     x <= 100
#'   }, max_v = 42
#' )
#'
#' # 42
#' guess_max_integer(
#'   function(x) {
#'     x <= 42
#'   }, min_v = 42
#' )
#'
#' # NA
#' guess_max_integer(
#'   function(x) {
#'     x <= 4
#'   }, min_v = 42
#' )
#'
#'
#'
#' @noRd
guess_max_integer <- function(validator, initial = NA, min_v = 0, max_v = Inf) {

  # DIPSAUS DEBUG START
  # validator <- function(x) { x <= 42 }
  # initial = NA
  # min_v = 0
  # max_v = Inf

  if(!validator(min_v)) { return(NA) }

  if(is.na(initial)) {
    if(is.infinite(max_v)) {
      initial <- max(min_v, 1)
    } else {
      initial <- ceiling((min_v + max_v) / 2.)
    }
  }

  min_v <- min_v
  max_v <- max_v

  if(is.finite( max_v )) {
    if( validator( max_v ) ) {
      return( max_v )
    }
  } else {
    while( validator(initial) ) {
      min_v <- initial
      initial <- initial * 2
    }
    max_v <- initial
    initial <- ceiling((min_v + max_v) / 2.)
  }



  test_number <- function(v) {

    if( min_v + 1 >= max_v ) { return( min_v ) }

    if( validator(v) ) {
      min_v <<- v
    } else {
      max_v <<- v
    }

    v <- ceiling((min_v + max_v) / 2.)

    return(Recall(v))
  }

  test_number(initial)

}
