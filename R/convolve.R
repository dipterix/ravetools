#' @name convolve
#' @title Convolution of \code{1D}, \code{2D}, \code{3D} data via \code{FFT}
#' @description
#' Use the 'Fast-Fourier' transform to compute the convolutions of two data
#' with zero padding. This function is mainly designed for image convolution.
#' For forward and backward convolution/filter, see \code{\link{filtfilt}}.
#'
#' @param x one-dimensional signal vector, two-dimensional image, or
#' three-dimensional volume; numeric or complex
#' @param filter kernel with the same number of dimensions as \code{x}
#' @returns Convolution results with the same length and dimensions as \code{x}.
#' If \code{x} is complex, results will be complex, otherwise results will
#' be real numbers.
#' @details This implementation uses 'Fast-Fourier' transform to perform
#' \code{1D}, \code{2D}, or \code{3D} convolution. Compared to implementations
#' using original mathematical definition of convolution, this approach is
#' much faster, especially for image and volume convolutions.
#'
#' The input \code{x} is zero-padded beyond edges. This is most common in image
#' or volume convolution, but less optimal for periodic one-dimensional signals.
#' Please use other implementations if non-zero padding is needed.
#'
#' The convolution results might be different to the ground truth by a precision
#' error, usually at \code{1e-13} level, depending on the \code{'FFTW3'}
#' library precision and implementation.
#'
#' @examples
#'
#'
#' # ---- 1D convolution ------------------------------------
#' x <- cumsum(rnorm(100))
#' filter <- dnorm(-2:2)
#' # normalize
#' filter <- filter / sum(filter)
#' smoothed <- convolve_signal(x, filter)
#'
#' plot(x, pch = 20)
#' lines(smoothed, col = 'red')
#'
#' # ---- 2D convolution ------------------------------------
#' x <- array(0, c(100, 100))
#' x[
#'   floor(runif(10, min = 1, max = 100)),
#'   floor(runif(10, min = 1, max = 100))
#' ] <- 1
#'
#' # smooth
#' kernel <- outer(dnorm(-2:2), dnorm(-2:2), FUN = "*")
#' kernel <- kernel / sum(kernel)
#'
#' y <- convolve_image(x, kernel)
#'
#' oldpar <- par(mfrow = c(1,2))
#' image(x, asp = 1, axes = FALSE, main = "Origin")
#' image(y, asp = 1, axes = FALSE, main = "Smoothed")
#' par(oldpar)
#'
#'
NULL

#' @rdname convolve
#' @export
convolve_signal <- function(x, filter) {

  len_x <- length(x)
  len_y <- length(filter)

  padded_len <- len_x + len_y + 1L

  junk_len <- ceiling(len_y / 2) - 1

  x_ <- c(as.double(x), rep(0.0, len_y + 1L))
  y_ <- c(as.double(filter), rep(0.0, len_x + 1L))

  a <- fftw_c2c(x_)
  b <- fftw_c2c(y_)

  a <- fftw_c2c(a * b, inverse = TRUE, ret = a) / padded_len
  b <- a[junk_len + seq_len(len_x)]

  if (!is.complex(x)) {
    b <- Re(b)
  }
  b
}

#' @rdname convolve
#' @export
convolve_image <- function(x, filter) {
  # make sure x and filter are matrix
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(filter)) {
    filter <- as.matrix(filter)
  }

  dim_x <- dim(x)
  dim_y <- dim(filter)

  padded_dim <- dim_x + dim_y + 1L

  junk_dim <- ceiling(dim_y / 2) - 1

  x_ <- array(0.0, padded_dim)
  x_[ seq_len(dim_x[[1]]), seq_len(dim_x[[2]]) ] <- x
  y_ <- array(0.0, padded_dim)
  y_[ seq_len(dim_y[[1]]), seq_len(dim_y[[2]]) ] <- filter

  a <- fftw_c2c_2d(x_)
  b <- fftw_c2c_2d(y_)

  a <- fftw_c2c_2d(a * b, inverse = TRUE, ret = a) / prod(padded_dim)
  b <- a[junk_dim[[1]] + seq_len(dim_x[[1]]),
         junk_dim[[2]] + seq_len(dim_x[[2]]),
         drop = FALSE]

  if (!is.complex(x)) {
    b <- Re(b)
  }
  b
}

#' @rdname convolve
#' @export
convolve_volume <- function(x, filter) {
  # make sure x and filter are matrix
  if (!is.array(x)) {
    x <- is.array(x)
  }

  if (!is.array(filter)) {
    filter <- is.array(filter)
  }

  dim_x <- dim(x)
  dim_y <- dim(filter)
  if (length(dim_x) < 3L) {
    dim_x <- c(dim_x, rep(1, 3 - length(dim_x)))
  } else if (length(dim_x) > 3L) {
    if (prod(dim_x) != prod(dim_x[1:3])) {
      stop("`convolve_volume`: x must be an array with 3 dimensions")
    }
    dim_x <- dim_x[1:3]
    dim(x) <- dim_x
  }
  if (length(dim_y) < 3L) {
    dim_y <- c(dim_y, rep(1, 3 - length(dim_y)))
  } else if (length(dim_y) > 3L) {
    if (prod(dim_y) != prod(dim_y[1:3])) {
      stop("`convolve_volume`: filter must be an array with 3 dimensions")
    }
    dim_y <- dim_y[1:3]
    dim(filter) <- dim_y
  }

  padded_dim <- dim_x + dim_y + 1L
  junk_dim <- ceiling(dim_y / 2) - 1L

  x_ <- array(0.0, padded_dim)
  x_[
    seq_len(dim_x[[1]]),
    seq_len(dim_x[[2]]),
    seq_len(dim_x[[3]])
  ] <- x
  y_ <- array(0.0, padded_dim)
  y_[
    seq_len(dim_y[[1]]),
    seq_len(dim_y[[2]]),
    seq_len(dim_y[[3]])
  ] <- filter

  a <- fftw_c2c_3d(x_)
  b <- fftw_c2c_3d(y_)

  a <- fftw_c2c_3d(a * b, inverse = TRUE, ret = a) / prod(padded_dim)
  b <- a[
    junk_dim[[1]] + seq_len(dim_x[[1]]),
    junk_dim[[2]] + seq_len(dim_x[[2]]),
    junk_dim[[3]] + seq_len(dim_x[[3]]),
    drop = FALSE
  ]

  if (!is.complex(x)) {
    b <- Re(b)
  }
  b
}


#' Grow volume mask
#' @param volume volume mask array, must be 3-dimensional array
#' @param x,y,z size of grow along each direction
#' @param threshold threshold after convolution
#' @returns A binary volume mask
#' @examples
#'
#'
#' oldpar <- par(mfrow = c(2,3), mar = c(0.1,0.1,3.1,0.1))
#'
#' mask <- array(0, c(21,21,21))
#' mask[11,11,11] <- 1
#' image(mask[11,,], asp = 1,
#'       main = "Original mask", axes = FALSE)
#' image(grow_volume(mask, 2)[11,,], asp = 1,
#'       main = "Dilated (size=2) mask", axes = FALSE)
#' image(grow_volume(mask, 5)[11,,], asp = 1,
#'       main = "Dilated (size=5) mask", axes = FALSE)
#'
#' mask[11, sample(11,2), sample(11,2)] <- 1
#' image(mask[11,,], asp = 1,
#'       main = "Original mask", axes = FALSE)
#' image(grow_volume(mask, 2)[11,,], asp = 1,
#'       main = "Dilated (size=2) mask", axes = FALSE)
#' image(grow_volume(mask, 5)[11,,], asp = 1,
#'       main = "Dilated (size=5) mask", axes = FALSE)
#'
#' par(oldpar)
#'
#'
#' @export
grow_volume <- function(volume, x, y = x, z = x, threshold = 0.5) {
  filter <- array(1.0, c(x, y, z))
  re <- convolve_volume(volume, filter)
  sel <- re > threshold
  mode(sel) <- "integer"
  sel
}

#' Morphological operations on a volume mask
#' @param volume volume mask array, must be 3-dimensional array; voxels above
#' zero are the foreground and missing values are the background
#' @param operation \code{'erode'} and \code{'dilate'} use a ball that holds
#' the voxels closer than \code{radius + 0.5} voxels, or with three radii the
#' ellipsoid whose semi-axes are \code{radius + 0.5} voxels (the voxels beyond
#' the edge of the volume count as foreground for the erosion and as background
#' for the dilation); \code{'largest_component'} keeps the largest component
#' (all of them if several tie), connected through the faces of the voxels,
#' and leaves its holes; \code{'fill_holes'} fills every background component
#' except the largest one
#' @param radius radius of the ball in voxels, ignored by the other operations;
#' or three radii, one per axis of \code{volume}, so that a distance in
#' millimeters can be honored on voxels that are longer along one axis
#' @returns A binary volume mask: an integer array of the dimension of
#' \code{volume}
#' @examples
#'
#'
#' mask <- array(0L, c(21, 21, 21))
#' mask[6:16, 6:16, 6:16] <- 1L
#' mask[11, 11, 11] <- 0L
#'
#' # a thin bridge to another object
#' mask[17:20, 11, 11] <- 1L
#'
#' # the erosion cuts the bridge, the largest component drops what is left
#' cleaned <- morph_volume(mask, "erode", radius = 1)
#' cleaned <- morph_volume(cleaned, "largest_component")
#' cleaned <- morph_volume(cleaned, "dilate", radius = 1)
#' cleaned <- morph_volume(cleaned, "fill_holes")
#'
#' oldpar <- par(mfrow = c(1, 2), mar = c(0.1, 0.1, 3.1, 0.1))
#' image(mask[, , 11], asp = 1, main = "Original mask", axes = FALSE)
#' image(cleaned[, , 11], asp = 1, main = "Cleaned mask", axes = FALSE)
#' par(oldpar)
#'
#'
#' @export
morph_volume <- function(volume,
                         operation = c("erode", "dilate", "largest_component", "fill_holes"),
                         radius = 1L) {
  operation <- match.arg(operation)
  dims <- dim(volume)
  if (length(dims) != 3L) {
    stop("`morph_volume`: `volume` must be a 3-dimensional array.")
  }
  radius <- suppressWarnings(as.integer(radius))
  if (!length(radius) %in% c(1L, 3L) || anyNA(radius) || any(radius < 0L)) {
    stop("`morph_volume`: `radius` must be a non-negative integer, or three of them (one per axis).")
  }
  radius <- rep_len(radius, 3L)
  foreground <- as.vector(volume > 0)
  foreground[is.na(foreground)] <- FALSE
  re <- switch(
    operation,
    "erode" = n4_mask_morph(foreground, dims, radius, TRUE),
    "dilate" = n4_mask_morph(foreground, dims, radius, FALSE),
    "largest_component" = n4_mask_largest_component(foreground, dims),
    "fill_holes" = n4_mask_cleanup(foreground, dims, 0L, FALSE, 0L)
  )
  array(as.integer(re), dim = dims)
}
