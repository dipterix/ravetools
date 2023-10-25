Matrix4 <- R6::R6Class(
  classname = "Matrix4",
  portable = TRUE,
  cloneable = FALSE,
  lock_objects = FALSE,
  private = list(
    .extern_ptr = NULL
  ),
  active = list(
    is_matrix4 = function() { TRUE },
    pointer = function() { private$.extern_ptr },
    elements = function() {
      self$to_array()
    }
  ),
  public = list(
    initialize = function() {
      private$.extern_ptr <- Matrix4__new()
    },
    format = function(...) {
      sprintf(
        "<Matrix4>\n%s",
        paste(apply(format(self$to_array()), 1, paste, collapse = " "), collapse = "\n")
      )
    },
    clone2 = function(...) {
      Matrix4$new()$copy(self)
    },
    set = function(n11, ..., byrow = TRUE) {
      if(is.matrix(n11)) {
        data <- n11
      } else {
        data <- matrix(c(n11, ...), nrow = 4, byrow = byrow)
      }
      if(!is.numeric(data)) {
        data <- as.numeric(data)
      }
      if(length(data) != 16) {
        stop("Matrix4$set() requires 16 elements")
      }
      Matrix4__from_array(private$.extern_ptr, data)
      self
    },
    to_array = function() {
      Matrix4__to_array(private$.extern_ptr)
    },
    identity = function() {
      Matrix4__identity(private$.extern_ptr)
      self
    },
    copy = function(m) {
      stopifnot(R6::is.R6(m) && isTRUE(m$is_matrix4))
      Matrix4__copy(private$.extern_ptr, m$pointer)
      self
    },
    copy_position = function(m) {
      stopifnot(R6::is.R6(m) && isTRUE(m$is_matrix4))
      Matrix4__copy_position(private$.extern_ptr, m$pointer)
      self
    },
    extract_basis = function(x_axis = NULL, y_axis = NULL, z_axis = NULL) {
      if(is.null(x_axis)) {
        x_axis <- new_vector3()
      } else {
        stopifnot(R6::is.R6(x_axis) && isTRUE(x_axis$is_vector3))
      }
      if(is.null(y_axis)) {
        y_axis <- new_vector3()
      } else {
        stopifnot(R6::is.R6(y_axis) && isTRUE(y_axis$is_vector3))
      }
      if(is.null(z_axis)) {
        z_axis <- new_vector3()
      } else {
        stopifnot(R6::is.R6(z_axis) && isTRUE(z_axis$is_vector3))
      }
      Matrix4__extract_basis(private$.extern_ptr, x_axis$pointer, y_axis$pointer, z_axis$pointer)
      self
    },
    make_basis = function(x_axis = NULL, y_axis = NULL, z_axis = NULL) {
      if(is.null(x_axis)) {
        x_axis <- new_vector3()
      } else if(is.numeric(x_axis) && length(x_axis) >= 3) {
        x_axis <- new_vector3(x_axis[1], x_axis[2], x_axis[3])
      } else {
        stopifnot(R6::is.R6(x_axis) && isTRUE(x_axis$is_vector3))
      }
      if(is.null(y_axis)) {
        y_axis <- new_vector3()
      } else if(is.numeric(y_axis) && length(y_axis) >= 3) {
        y_axis <- new_vector3(y_axis[1], y_axis[2], y_axis[3])
      } else {
        stopifnot(R6::is.R6(y_axis) && isTRUE(y_axis$is_vector3))
      }
      if(is.null(z_axis)) {
        z_axis <- new_vector3()
      } else if(is.numeric(z_axis) && length(z_axis) >= 3) {
        z_axis <- new_vector3(z_axis[1], z_axis[2], z_axis[3])
      } else {
        stopifnot(R6::is.R6(z_axis) && isTRUE(z_axis$is_vector3))
      }
      Matrix4__make_basis(private$.extern_ptr, x_axis$pointer, y_axis$pointer, z_axis$pointer)
      self
    },
    extract_rotation = function(m) {
      stopifnot(R6::is.R6(m) && isTRUE(m$is_matrix4))
      Matrix4__extract_rotation(private$.extern_ptr, m$pointer)
      self
    },
    look_at = function(eye, target, up) {
      if(is.numeric(eye) && length(eye) >= 3) {
        eye <- new_vector3(eye[1], eye[2], eye[3])
      }
      if(is.numeric(target) && length(target) >= 3) {
        target <- new_vector3(target[1], target[2], target[3])
      }
      if(is.numeric(up) && length(up) >= 3) {
        up <- new_vector3(up[1], up[2], up[3])
      }
      stopifnot(R6::is.R6(eye) && isTRUE(eye$is_vector3))
      stopifnot(R6::is.R6(target) && isTRUE(target$is_vector3))
      stopifnot(R6::is.R6(up) && isTRUE(up$is_vector3))
      Matrix4__look_at(private$.extern_ptr, eye$pointer, target$pointer, up$pointer)
      self
    },
    multiply_matrices = function(a, b) {
      stopifnot(R6::is.R6(a) && isTRUE(a$is_matrix4))
      stopifnot(R6::is.R6(b) && isTRUE(b$is_matrix4))
      Matrix4__multiply_matrices(private$.extern_ptr, a$pointer, b$pointer)
      self
    },
    multiply = function(m) {
      stopifnot(R6::is.R6(m) && isTRUE(m$is_matrix4))
      Matrix4__multiply_matrices(private$.extern_ptr, private$.extern_ptr, m$pointer)
      self
    },
    premultiply = function(m) {
      stopifnot(R6::is.R6(m) && isTRUE(m$is_matrix4))
      Matrix4__multiply_matrices(private$.extern_ptr, m$pointer, private$.extern_ptr)
      self
    },
    multiply_scalar = function(scalar) {
      if(!is.numeric(scalar)) {
        scalar <- as.numeric(scalar[[1]])
      } else {
        scalar <- scalar[[1]]
      }
      Matrix4__multiply_scalar(private$.extern_ptr, scalar)
      self
    },
    determinant = function() {
      Matrix4__determinant(private$.extern_ptr)
    },
    transpose = function() {
      Matrix4__transpose(private$.extern_ptr)
      self
    },
    set_position = function(v, ...) {
      if(R6::is.R6(v) && isTRUE(v$is_vector3)) {
        v <- as.numeric(v[,1])
      } else {
        v <- as.numeric(c(v, ...))[1:3]
      }
      Matrix4__set_position(private$.extern_ptr, v[[1]], v[[2]], v[[3]])
      self
    },
    invert = function() {
      Matrix4__invert(private$.extern_ptr)
      self
    },
    scale = function(v, ...) {
      if(!R6::is.R6(v) || !isTRUE(v$is_vector3)) {
        v <- as.numeric(c(v, ...))[1:3]
        v <- new_vector3(v[[1]], v[[2]], v[[3]])
      }
      Matrix4__scale(private$.extern_ptr, v$pointer)
      self
    },
    get_max_scale_on_axis = function() {
      Matrix4__get_max_scale_on_axis(private$.extern_ptr)
    },
    make_translation = function(x, y = NULL, z = NULL) {
      if(R6::is.R6(x) && isTRUE(x$is_vector3)) {
        x <- as.numeric(x[,1])
      } else {
        x <- as.numeric(c(x, y, z))[1:3]
      }
      Matrix4__make_translation(private$.extern_ptr, x[[1]], x[[2]], x[[3]])
      self
    },
    make_rotation_x = function(theta) {
      theta <- as.double(theta)[[1]]
      Matrix4__make_rotation_x(private$.extern_ptr, theta)
      self
    },
    make_rotation_y = function(theta) {
      theta <- as.double(theta)[[1]]
      Matrix4__make_rotation_y(private$.extern_ptr, theta)
      self
    },
    make_rotation_z = function(theta) {
      theta <- as.double(theta)[[1]]
      Matrix4__make_rotation_z(private$.extern_ptr, theta)
      self
    },
    make_rotation_axis = function(axis, angle) {
      if(!R6::is.R6(axis) || !isTRUE(axis$is_vector3)) {
        axis <- as.numeric(axis)[1:3]
        axis <- new_vector3(axis[[1]], axis[[2]], axis[[3]])
      }
      angle <- as.double(angle)[[1]]
      Matrix4__make_rotation_axis(private$.extern_ptr, axis$pointer, angle)
      self
    },
    make_scale = function(x, y = NULL, z = NULL) {
      if(R6::is.R6(x) && isTRUE(x$is_vector3)) {
        x <- as.numeric(x[,1])
      } else {
        x <- as.numeric(c(x, y, z))[1:3]
      }
      Matrix4__make_scale(private$.extern_ptr, x[[1]], x[[2]], x[[3]])
      self
    },
    make_shear = function(xy, xz, yx, yz, zx, zy) {
      xy <- as.double(xy)[[1]]
      xz <- as.double(xz)[[1]]
      yx <- as.double(yx)[[1]]
      yz <- as.double(yz)[[1]]
      zx <- as.double(zx)[[1]]
      zy <- as.double(zy)[[1]]
      Matrix4__make_shear(private$.extern_ptr, xy, xz, yx, yz, zx, zy)
      self
    },
    make_perpective = function(left, right, top, bottom, near, far) {
      left <- as.double(left)[[1]]
      right <- as.double(right)[[1]]
      top <- as.double(top)[[1]]
      bottom <- as.double(bottom)[[1]]
      near <- as.double(near)[[1]]
      far <- as.double(far)[[1]]
      Matrix4__make_perspective(private$.extern_ptr, left, right, top, bottom, near, far)
      self
    },
    make_orthographic = function(left, right, top, bottom, near, far) {
      left <- as.double(left)[[1]]
      right <- as.double(right)[[1]]
      top <- as.double(top)[[1]]
      bottom <- as.double(bottom)[[1]]
      near <- as.double(near)[[1]]
      far <- as.double(far)[[1]]
      Matrix4__make_orthographic(private$.extern_ptr, left, right, top, bottom, near, far)
      self
    }
  )
)

#' Create a \code{Matrix4} instance for \code{'Affine'} transform
#' @returns A \code{Matrix4} instance
#' @seealso \code{\link{new_vector3}}, \code{\link{new_quaternion}}
#' @export
new_matrix4 <- function() {
  Matrix4$new()
}

#' @export
as.matrix.Matrix4 <- function(x, ...) {
  if(R6::is.R6(x) && isTRUE(x$is_matrix4)) {
    return(x$to_array())
  }
  NextMethod("as.matrix")
}

#' @export
`[.Matrix4` <- function(x, i, ..., drop = TRUE) {
  if(missing(i)) {
    x$to_array()[, ..., drop = drop]
  } else {
    x$to_array()[i, ..., drop = drop]
  }
}

#' @export
dim.Matrix4 <- function(x) {
  if(R6::is.R6(x) && isTRUE(x$is_matrix4)) {
    return(c(4L, 4L))
  }
  NextMethod("dim")
}

#' @export
`==.Matrix4` <- function(e1, e2) {
  e2 == e1$to_array()
}

