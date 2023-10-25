Quaternion <- R6::R6Class(
  classname = "Quaternion",
  portable = TRUE,
  cloneable = FALSE,
  lock_objects = FALSE,
  private = list(
    .extern_ptr = NULL
  ),
  active = list(
    is_quaternion = function() { TRUE },
    pointer = function() { private$.extern_ptr }
  ),
  public = list(
    initialize = function() {
      private$.extern_ptr <- Quaternion__new()
    },
    format = function(...) {
      sprintf(
        "<Quaternion>\n%s",
        paste(apply(format(self$to_array()), 1, paste, collapse = " "), collapse = "\n")
      )
    },
    clone2 = function(...) {
      Quaternion$new()$copy(self)
    },
    set = function(x, y, z, w) {
      data <- as.numeric(c(x, y, z, w))[1:4]
      Quaternion__set(private$.extern_ptr, data[[1]], data[[2]], data[[3]], data[[4]])
      self
    },
    to_array = function() {
      Quaternion__to_array(private$.extern_ptr)
    },
    copy = function(m) {
      stopifnot(R6::is.R6(m) && isTRUE(m$is_quaternion))
      Quaternion__copy(private$.extern_ptr, m$pointer)
      self
    }
  )
)

#' Create a \code{Quaternion} instance to store '3D' rotation
#' @description
#' Create instances that mimic the \code{'three.js'} syntax.
#' @param x,y,z,w numeric of length one
#' @returns A \code{Quaternion} instance
#' @seealso \code{\link{new_vector3}}, \code{\link{new_matrix4}}
#' @export
new_quaternion <- function(x = 0, y = 0, z = 0, w = 1) {
  Quaternion$new()$set(x, y, z, w)
}

#' @export
`[.Quaternion` <- function(x, i, ..., drop = TRUE) {
  if(missing(i)) {
    x$to_array()
  } else {
    x$to_array()[i]
  }
}


#' @export
`==.Quaternion` <- function(e1, e2) {
  e2 == e1$to_array()
}

