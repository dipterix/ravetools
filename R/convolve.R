#' @export
convolve_signal <- function(x, filter) {

  len_x <- length(x)
  len_y <- length(filter)

  padded_len <- len_x + len_y + 1L

  prepad_len <- ceiling((len_y + 1) / 2) + 1L

  postpad_len <- len_y + 1 - prepad_len

  x_ <- c(rep(0.0, prepad_len), as.double(x), rep(0.0, postpad_len))
  y_ <- c(as.double(filter), rep(0.0, len_x + 1L))

  a <- fftw_r2c(x_)
  b <- fftw_r2c(y_)

  a <- fftw_c2c(a * b, inverse = TRUE, ret = a) / padded_len
  b <- a[prepad_len + postpad_len + seq_len(len_x)]

  Re(b)
}


#' @export
convolve_image <- function(x, filter) {
  # make sure x and filter are matrix
  if(!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if(!is.matrix(filter)) {
    filter <- as.matrix(filter)
  }

  nrow_x <- nrow(x)
  nrow_y <- nrow(filter)
  ncol_x <- ncol(x)
  ncol_y <- ncol(filter)

  padded_row <- nrow_x + nrow_y + 1L
  padded_col <- ncol_x + ncol_y + 1L

  prepad_row <- ceiling((nrow_y + 1) / 2) + 1L
  prepad_col <- ceiling((ncol_y + 1) / 2) + 1L

  postpad_row <- nrow_y + 1 - prepad_row
  postpad_col <- ncol_y + 1 - prepad_col

  x_ <- array(0.0, c(padded_row, padded_col))
  x_[
    prepad_row + seq_len(nrow_x),
    prepad_col + seq_len(ncol_x)
  ] <- x
  y_ <- array(0.0, c(padded_row, padded_col))
  y_[
    seq_len(nrow_y),
    seq_len(ncol_y)
  ] <- filter

  a <- fftw_r2c_2d(x_)
  b <- fftw_r2c_2d(y_)

  a <- fftw_c2c_2d(a * b, inverse = TRUE, ret = a) / (padded_row * padded_col)
  b <- a[prepad_row + postpad_row + seq_len(nrow_x),
         prepad_col + postpad_col + seq_len(ncol_x),
         drop = FALSE]

  Re(b)
}




#' @export
convolve_volume <- function(x, filter) {
  # make sure x and filter are matrix
  if(!is.array(x)) {
    x <- is.array(x)
  }

  if(!is.array(filter)) {
    filter <- is.array(filter)
  }

  dim_x <- dim(x)
  dim_y <- dim(filter)
  if(length(dim_x) < 3L) {
    dim_x <- c(dim_x, rep(1, 3 - length(dim_x)))
  } else if(length(dim_x) > 3L) {
    if(prod(dim_x) != prod(dim_x[1:3])) {
      stop("`convolve_volume`: x must be an array with 3 dimensions")
    }
    dim_x <- dim_x[1:3]
    dim(x) <- dim_x
  }
  if(length(dim_y) < 3L) {
    dim_y <- c(dim_y, rep(1, 3 - length(dim_y)))
  } else if(length(dim_y) > 3L) {
    if(prod(dim_y) != prod(dim_y[1:3])) {
      stop("`convolve_volume`: filter must be an array with 3 dimensions")
    }
    dim_y <- dim_y[1:3]
    dim(filter) <- dim_y
  }

  padded_dim <- dim_x + dim_y + 1L

  prepad_dim <- ceiling((dim_y + 1) / 2) + 1L

  postpad_dim <- dim_y + 1 - prepad_dim

  x_ <- array(0.0, padded_dim)
  x_[
    prepad_dim[[1]] + seq_len(dim_x[[1]]),
    prepad_dim[[2]] + seq_len(dim_x[[2]]),
    prepad_dim[[3]] + seq_len(dim_x[[3]])
  ] <- x
  y_ <- array(0.0, padded_dim)
  y_[
    seq_len(dim_y[[1]]),
    seq_len(dim_y[[2]]),
    seq_len(dim_y[[3]])
  ] <- filter

  a <- fftw_r2c_3d(x_)
  b <- fftw_r2c_3d(y_)

  a <- fftw_c2c_3d(a * b, inverse = TRUE, ret = a) / prod(padded_dim)
  b <- a[
    prepad_dim[[1]] + postpad_dim[[1]] + seq_len(dim_x[[1]]),
    prepad_dim[[2]] + postpad_dim[[2]] + seq_len(dim_x[[2]]),
    prepad_dim[[3]] + postpad_dim[[3]] + seq_len(dim_x[[3]]),
    drop = FALSE
  ]

  Re(b)
}
