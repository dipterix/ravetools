#' @name rgl-call
#' @title Safe ways to call package \code{'rgl'} without requiring \code{'x11'}
#' @description
#' Internally used for example show-cases. Please install package \code{'rgl'}
#' manually to use these functions.
#'
#' @param FUN \code{'rgl'} function name
#' @param ... passed to \code{'rgl'} function
#' @param expr expression within which \code{'rgl'} functions are called
#' @param quoted whether \code{expr} is quoted
#' @param env environment in which \code{expr} is evaluated
#' @param x triangular \code{'mesh3d'} object
#' @param length,lwd,col normal vector length, size, and color
#' @examples
#'
#'
#' # Make sure the example does not run when compiling
#' # or check the package
#' if(FALSE) {
#'
#'   volume <- array(0, dim = c(8,8,8))
#'   volume[4:5, 4:5, 4:5] <- 1
#'   mesh <- mesh_from_volume(volume, verbose = FALSE)
#'
#'   rgl_view({
#'
#'     rgl_call("shade3d", mesh, col = 3)
#'     rgl_plot_normals(mesh)
#'
#'   })
#'
#' }
#'
#'
NULL

check_rgl <- function(strict = NA) {
  rgl_unavailable <- function() {
    msg <- "Package `rgl` is not installed. Please install `rgl` to use this function."
    if(isTRUE(strict)) {
      stop(msg)
    } else if(is.na(strict)) {
      message(msg)
    }
    FALSE
  }
  if(identical(Sys.getenv("RAVETOOLS_RGL_DISABLED"), "1")) { return(rgl_unavailable()) }
  if( getOption("ravetools.rgl.disabled", FALSE) ) { return(rgl_unavailable()) }
  if( system.file(package = "rgl") == "" ) { return(rgl_unavailable()) }
  TRUE
}

#' @rdname rgl-call
#' @export
rgl_call <- function(FUN, ...) {
  check_rgl()

  # Cannot set this option back with on.exit, rgl may throw warning/errors
  # when X11 is unavailable. This option forces rgl to use webgl renderers.
  options(rgl.useNULL = TRUE)
  rgl <- asNamespace("rgl")
  f <- rgl[[FUN]]
  if(!is.function(f)) {
    stop("Function ", FUN, " is not a `rgl` function.")
  }
  rgl[[FUN]](...)
}

#' @rdname rgl-call
#' @export
rgl_view <- function(expr, quoted = FALSE, env = parent.frame()) {
  # Suppress RGL
  if(!quoted) {
    expr <- substitute(expr)
  }
  dev <- rgl_call("open3d")
  on.exit({
    rgl_call("close3d", dev = dev)
  }, add = TRUE, after = TRUE)
  eval(expr, envir = env)
  rgl_call("rglwidget")
}

#' @rdname rgl-call
#' @export
rgl_plot_normals <- function (x, length = 1, lwd = 1, col = 1, ...) {
  if (!"mesh3d" %in% class(x)) {
    stop("please provide object of class mesh3d")
  }
  args <- list(...)
  if ("long" %in% names(args)) {
    length <- args$long
    warning("argument 'long' is deprecated, please use 'length' instead")
  }
  if (is.null(x$normals)) {
    if (!is.null(x$it)) {
      x <- vcg_update_normals(x)
    } else {
      stop("mesh has neither normals nor faces")
    }
  }
  n.mesh <- list()
  lvb <- dim(x$vb)[2]
  vb <- x$vb
  vb.norm <- vb[1:3, , drop = FALSE] + t(length * t(x$normals[1:3, , drop = FALSE]))
  vb <- cbind(vb[1:3, , drop = FALSE], vb.norm)
  vb <- rbind(vb, 1)
  it <- rbind(1:lvb, 1:lvb, (1:lvb) + lvb)
  n.mesh$vb <- vb
  n.mesh$it <- it
  class(n.mesh) <- c("mesh3d", "shape3d")
  rgl_call("wire3d", n.mesh, color = col, lwd = lwd, lit = FALSE)
}
