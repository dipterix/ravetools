## usethis namespace: start
#' @importFrom stats approx
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics image
#' @importFrom R6 R6Class
#' @importFrom R6 is.R6
#' @useDynLib ravetools, .registration = TRUE
## usethis namespace: end
NULL

`%?<-%` <- function(lhs, value) {
  env <- parent.frame()
  lhs <- substitute(lhs)

  isnull <- tryCatch({
    is.null(eval(lhs, envir = env))
  }, error = function(e) {
    return(TRUE)
  })

  if (isnull) {
    eval(as.call(list( quote(`<-`), lhs, value )), envir = env)
  }
}

deparse1 <- function(expr, collapse = " ") {
  paste(deparse(expr), collapse = collapse)
}


stopifnot2 <- function(..., msg = "Condition not satisfied") {
  if (!all(c(...))) {
    stop(msg)
  }
}

rand_string <- function(length = 10) {
  paste(sample(c(letters, LETTERS, 0:9), length, replace = TRUE), collapse = "")
}


tempdir2 <- function(check = TRUE) {
  path <- getOption("ravetools.tempdir",
            default = Sys.getenv("RAVETOOLS_TEMPDIR",
                                 unset = tempdir(check = FALSE)))
  if (check && !dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }
  path
}

tempfile2 <- function(
  pattern = "ravetmp-", tmpdir = file.path(tempdir2(check = TRUE), "ravetools"),
  fileext = "") {
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  }
  if (getOption("ravetools.debug", FALSE)) {
    file.path(tmpdir, sprintf("%s%s%s", pattern, rand_string(), fileext))
  } else {
    tempfile(pattern = pattern, tmpdir = tmpdir, fileext = fileext)
  }

}


.missing_arg <- alist(x = )

.matlab_palette <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
          "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
          "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
          "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
          "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
          "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
          "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
          "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
          "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
          "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
          "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
          "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
          "#AF0000", "#9F0000", "#8F0000", "#800000")

#' 'Matlab' heat-map plot palette
#' @returns vector of 64 colors
#' @export
matlab_palette <- function() {
  .matlab_palette
}

#' @name parallel-options
#' @title Set or get thread options
#' @param n_threads number of threads to set
#' @param stack_size Stack size (in bytes) to use for worker threads. The
#' default used for \code{"auto"} is 2MB on 32-bit systems and 4MB on 64-bit
#' systems.
#' @returns \code{detect_threads} returns an integer of default threads that
#' is determined by the number of \code{CPU} cores; \code{ravetools_threads}
#' returns nothing.
#'
#' @examples
#'
#' detect_threads()
#'
#' ravetools_threads(n_threads = 2)
#'
#' @export
detect_threads <- function() {
  getDefaultNumThreads()
}

#' @rdname parallel-options
#' @export
ravetools_threads <- function(n_threads = "auto", stack_size = "auto") {
  if (identical(n_threads, "auto"))
    n_threads <- -1L
  else if (!is.numeric(n_threads))
    stop("n_threads must be an integer")
  else n_threads <- as.integer(n_threads)
  if (identical(stack_size, "auto"))
    stack_size <- 0L
  else if (!is.numeric(stack_size))
    stop("stack_size must be an integer")
  else stack_size <- as.integer(stack_size)
  if (n_threads == -1L)
    Sys.unsetenv("RAVETOOLS_NUM_THREADS")
  else Sys.setenv(RAVETOOLS_NUM_THREADS = n_threads)
  if (stack_size == 0L)
    Sys.unsetenv("RAVETOOLS_STACK_SIZE")
  else Sys.setenv(RAVETOOLS_STACK_SIZE = stack_size)
  invisible()
}

#' @title Internal function
#' @description
#' Do not call this function directly
#' @param if_interactive,verbose default is \code{TRUE}
#' @returns logical
#' @keywords internal
#' @export
is_not_cran <- function(if_interactive = TRUE, verbose = FALSE) {
  not_cran_flag <- identical(toupper(as.character(Sys.getenv("NOT_CRAN", ""))), "TRUE")
  limit_core_flag <- identical(toupper(Sys.getenv("_R_CHECK_LIMIT_CORES_")), "TRUE")
  interactive_flag <- interactive()
  if (limit_core_flag) {
    if (verbose) {
      message("_R_CHECK_LIMIT_CORES_ is TRUE/true (on CRAN)")
    }
    return(FALSE)
  }
  if (not_cran_flag) {
    if (verbose) {
      message("NOT_CRAN is TRUE/true (not on CRAN)")
    }
    return(TRUE)
  }
  if (interactive_flag) {
    if_interactive <- isTRUE(if_interactive)
    if (verbose) {
      message(sprintf("Session is interactive (%son CRAN)",
                      ifelse(if_interactive, "", "not ")))
    }
    return(if_interactive)
  }
  if (verbose) {
    message("No flag detected, default is on CRAN")
  }
  return(FALSE)
}

#' Get external function from 'RAVE'
#' @description
#' Internal function used for examples relative to 'RAVE' project and should
#' not be used directly.
#' @param name function or variable name
#' @param pkg 'RAVE' package name
#' @param inherit passed to \code{\link{get0}}
#' @param on_missing default value to return of no function is found
#' @returns Function object if found, otherwise \code{on_missing}.
#' @export
internal_rave_function <- function(name, pkg, inherit = TRUE, on_missing = NULL) {
  if (!pkg %in% c("ravecore", "raveio", "ravedash", "ravebuiltins", "rave", "threeBrain",
                 "dipsaus", "filearray", "readNSx", "rpymat", "rpyANTs")) {
    stop("`extern_function`: Package [", pkg, "] is not a RAVE package.")
  }
  if (system.file(package = pkg) == "") { return(on_missing) }
  ns <- asNamespace(pkg)
  get0(x = name, envir = ns, inherits = inherit, ifnotfound = on_missing)
}

#' @export
`print.ravetools-printable` <- function(x, ...) {
  cat(paste(c(format(x, ...), ""), collapse = "\n"))
}

#' @export
`format.ravetools-printable` <- function(x, ...) {
  printable_message <- attr(x, "printable_message")
  if (length(printable_message)) {
    return(printable_message)
  }
  NextMethod("format")
}

#' Left 'Hippocampus' of 'N27-Collin' brain
#'
#' @format
#' A three-mode integer mask array with values of \code{1} ('Hippocampus')
#' and \code{0} (other brain tissues)
#'
"left_hippocampus_mask"


#' Sample stimulation recording
#'
#' @format
#' A list of one-second signal trace and sample rate (30000)
#'
"stimulation_signal"


lapply_async <- function(x, FUN, FUN.args = list(), callback = NULL, ncores = NULL,
                         on_failure = "multisession", ...) {
  if (system.file(package = "raveio") == "") {
    do.call("lapply", c(list(X = x, FUN = FUN), FUN.args))
    ret <- lapply(x, FUN)
  } else {
    raveio <- asNamespace("raveio")
    ret <- raveio$lapply_async(x = x, FUN = FUN, FUN.args = FUN.args, callback = callback, ncores = ncores, on_failure = on_failure, ...)
  }
  ret
}


#' @title Map continuous values to colors
#' @description
#' Linearly maps a numeric vector onto a color ramp, clamping values outside
#' the given range to the range's endpoints.
#'
#' @param values numeric vector of values to map to colors
#' @param clim length-two numeric vector giving the value range to map from;
#' values outside \code{[clim[1], clim[2]]} are clamped to the nearer endpoint
#' before mapping. Default is \code{range(values, na.rm = TRUE)}
#' @param cmap the color ramp to map onto: either a vector of colors (passed
#' to \code{\link[grDevices]{colorRamp}} to build the ramp function) or a
#' function such as the one returned by \code{\link[grDevices]{colorRamp}}
#' that takes a numeric vector with elements in \code{[0, 1]} and returns an
#' \code{n x 3} (or \code{n x 4}, with alpha) matrix of color-channel values
#' in \code{[0, 255]}. Default is \code{grDevices::hcl.colors(11)}
#' @param ... passed to \code{\link[grDevices]{colorRamp}} when \code{cmap} is
#' a vector of colors (for example, \code{alpha = TRUE} to build an
#' alpha-aware ramp)
#'
#' @returns A character vector of \code{'#RRGGBB'} (or \code{'#RRGGBBAA'})
#' color strings, the same length as \code{values}.
#'
#' @examples
#'
#' x <- rnorm(100)
#' col <- color_ramp_continuous(x)
#'
#' plot(x, col = col, pch = 16)
#'
#'
#' # Change color palettes with vector of colors
#' col <- color_ramp_continuous(
#'   x, cmap = c("lightgreen", "white", "pink"))
#' plot(x, col = col, pch = 16)
#'
#' # Using colorRamp
#' col <- color_ramp_continuous(
#'   x, cmap = colorRamp(c("black", "orangered", "orange")))
#' plot(x, col = col, pch = 16)
#'
#' # Using color ramp palette `function(n) { ... }`
#' col <- color_ramp_continuous(
#'   x, cmap = hcl.colors, palette = "Blue-Red 3")
#' plot(x, col = col, pch = 16)
#'
#' # Set range
#' col <- color_ramp_continuous(
#'   x, clim = c(0, 1),
#'   cmap = c("black", "orangered", "orange"))
#' plot(x, col = col, pch = 16)
#'
#'
#'
#' @export
color_ramp_continuous <- function(values, clim = range(values, na.rm = TRUE),
                                  cmap = grDevices::hcl.colors(11), ...) {
  clim[is.finite(clim)]
  if (!length(clim)) {
    stop("Value range is invalid: the value range must be finite values")
  }
  clim <- range(clim, na.rm = TRUE)
  span <- clim[[2]] - clim[[1]]
  if (span <= 0) { span <- 1}
  values_01 <- (values - clim[[1]]) / span

  values_01[values_01 < 0] <- 0
  values_01[values_01 > 1] <- 1

  # cmap is either a vector or a function returned by colorRamp
  if (!is.function(cmap)) {
    cmap <- grDevices::colorRamp(colors = cmap, ...)(values_01)
  } else {
    if (names(formals(cmap))[[1]] == "n") {
      # colorRampPalette
      cmap <- cmap(256, ...)
      return(cmap[round(values_01 * 255) + 1])
    } else {
      cmap <- cmap(values_01, ...)
    }
  }
  # cmap <- cmap(values_01)
  if (ncol(cmap) >= 4) {
    col <- grDevices::rgb(
      red = cmap[, 1],
      green = cmap[, 2],
      blue = cmap[, 3],
      alpha =  cmap[, 4],
      maxColorValue = 255
    )
  } else {
    col <- grDevices::rgb(
      red = cmap[, 1],
      green = cmap[, 2],
      blue = cmap[, 3],
      maxColorValue = 255
    )
  }
  col
}


# color_ramp_
